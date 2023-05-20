#!flask/bin/python

import sys
import argparse
import re
from flask import Flask, jsonify, abort, make_response, request, url_for, current_app
from datetime import timedelta
from functools import update_wrapper
import time
import MySQLdb
from MySQLdb.constants import FIELD_TYPE
import config
from tree_exporter import get_tree
import json
import itertools
app = Flask(__name__)

KEGG_ALLOWED_FIELDS = set(['percent_identity', 'eval', 'query_percent_alignment', 'subject_percent_alignment'])
PFAM_ALLOWED_FIELDS = set(['eval'])
TIGRFAM_ALLOWED_FIELDS = set(['eval'])

"""DB query related """

def checkDbName(database):
    if database == 'gtdb_bacteria':
        return config.bacterial_database
    elif database == 'gtdb_archaea':
        return config.archaeal_database
    else:
        raise Exception('Unrecognized database %s ' % database)

def getDb(database, username=config.username, password=config.password, host=config.host, port=config.port):
    database = checkDbName(database)
    _my_conv = {FIELD_TYPE.LONG: int, FIELD_TYPE.INT24: int}
    conn = MySQLdb.connect(user=username, passwd=password, port=port,
                           host=host, db=database, conv=_my_conv)
    return conn


def make_error_response(msg='Default error, unknown reasone', code=500):
    return make_response(jsonify({'error': msg}), code)


def bad_request(msg='Bad request'):
    return make_error_response(msg=msg, code=400)


def crossdomain(origin=None, methods=None, headers=None,
                max_age=21600, attach_to_all=True,
                automatic_options=True):
    if methods is not None:
        methods = ', '.join(sorted(x.upper() for x in methods))
    if headers is not None and not isinstance(headers, basestring):
        headers = ', '.join(x.upper() for x in headers)
    if not isinstance(origin, basestring):
        origin = ', '.join(origin)
    if isinstance(max_age, timedelta):
        max_age = max_age.total_seconds()

    def get_methods():
        if methods is not None:
            return methods

        options_resp = current_app.make_default_options_response()
        return options_resp.headers['allow']

    def decorator(f):
        def wrapped_function(*args, **kwargs):
            if automatic_options and request.method == 'OPTIONS':
                resp = current_app.make_default_options_response()
            else:
                resp = make_response(f(*args, **kwargs))
            if not attach_to_all and request.method != 'OPTIONS':
                return resp

            h = resp.headers
            h['Content-Type'] = 'Application/JSON'
            h['Access-Control-Allow-Origin'] = origin
            h['Access-Control-Allow-Methods'] = get_methods()
            h['Access-Control-Max-Age'] = str(max_age)
            if headers is not None:
                h['Access-Control-Allow-Headers'] = headers
            return resp

        f.provide_automatic_options = False
        return update_wrapper(wrapped_function, f)
    return decorator

# /*----------  query treenodes  ----------*/

def getThresholdClause(thresholds):
    '''
    returns a pair (threshold_sql_clause, [threshold value ...])
    threshold_sql_clause is to be included in a where clause,
        note that it will have an AND keyword in front
    threhsold values will be a list to supply to cursor.execute(sql, threshold_values)
    '''
    threshold_clause = ''
    if len(thresholds) > 0:
        threshold_clause = " AND " + (" AND ".join(\
            ["%s %s %%s" % (threshold['fieldname'],
                '<=' if threshold['fieldname']=='eval' else '>=')
                for threshold in thresholds]))
    threshold_values = [th['value'] for th in thresholds]
    return (threshold_clause, threshold_values)

def check_threshold(thresholds, allowed_fields):
    if type(thresholds) != list:
        raise Exception("thresholds must be a list")
    if not thresholds:
        return
    for th in thresholds:
        if 'fieldname' not in th or 'value' not in th:
            raise Exception('`fieldname` or `value` must be in thresholds')
    threshold_fields = set([th['fieldname'] for th in thresholds])
    if not threshold_fields.issubset(allowed_fields):
        bad_fields = threshold_fields - allowed_fields
        raise Exception(str(list(bad_fields)) + " are not allowed in thresholds")
    try:
        for i in xrange(len(thresholds)):
            thresholds[i]['value'] = float(thresholds[i]['value'])
    except Exception as e:
        raise Exception('Value should be a number')



def _get_node_ids_from_tophits(db, top_hit_table, search_colname, search_ids, thresholds):
    """
    search top_hit_table for genomes' `node_id`
    where genomes must contain all search_ids, and each hit
    must have scores satisfying `thresholds`
    """
    search_id_quoted = (",".join(\
        ["'%s'" % (search_id) for search_id in search_ids]))
    search_clause = "%s IN (%s)" % (search_colname, search_id_quoted)
    print(thresholds)
    threshold_clause, threshold_values = getThresholdClause(thresholds)
    sql = """
        SELECT node_id
        FROM gtdb_node gn
        JOIN
            (
            SELECT
                gtdb_id,
                COUNT(DISTINCT {search_colname}) AS num_hit_per_genome
            FROM {top_hit_table}
            WHERE {search_clause} {threshold_clause}
            GROUP BY gtdb_id
            HAVING num_hit_per_genome >= {ids_length}
            ) g
        ON gn.gtdb_id = g.gtdb_id
        """.format(
            ids_length = str(len(search_ids)),
            top_hit_table=top_hit_table,
            search_clause=search_clause,
            search_colname=search_colname,
            threshold_clause=threshold_clause)
    c = db.cursor()
    c.execute(sql, threshold_values)
    rows = c.fetchall()
    if not rows:  # no result
        return []
    node_ids = [r[0] for r in rows]
    return sorted(node_ids)

def getNodeIdFromDomains(db, domains, thresholds=[]):
    return _get_node_ids_from_tophits(db, 'pfam_top_hits', 'pfam_id', domains, thresholds)


def getNodeIdFromKEGG(db, keggs, thresholds=[]):
    return _get_node_ids_from_tophits(db, 'kegg_top_hits', 'kegg_id', keggs, thresholds)


def getNodeIDFromTigrfam(db, domains, thresholds=[]):
    return _get_node_ids_from_tophits(db, 'tigrfam_top_hits', 'tigrfam_id', domains, thresholds)


def getNodeIdFromTaxIds(db, taxids):
    '''
    receive [species_taxid ...]
    return a [node_id ...]
    '''
    # find all child ncbi_taxid with parents taxids
    sql = '''SELECT node_id FROM node_species WHERE species_id IN
            (SELECT get_species(ncbi_taxid) FROM taxonomy WHERE ncbi_taxid IN (%s)
            AND get_species(ncbi_taxid) != 1);''' % ','.join([str(t) for t in taxids])
    c = db.cursor()
    c.execute(sql)
    rows = c.fetchall()
    return [r[0] for r in rows]

def getTaxIdsFromNames(db, names):
    '''
    receive [scientific_names ...]
    return a [taxids ...]
    '''
    # find all taxids using their corresponding scientific names in db
    sql = 'SELECT taxid FROM names WHERE scientific_name IN (%s);' % ','.join(['"%s"' % str(n) for n in names])
    c = db.cursor()
    c.execute(sql)
    rows = c.fetchall()
    return [r[0] for r in rows]

@app.route('/<database>/treeNodes/by/domains', methods=['OPTIONS'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def allowCORSdomainsDomains(database):
    return True

def isPfamAcc(text):
    return bool(re.match('^PF\d{5}$', text))


@app.route('/<database>/treeNodes/by/domains', methods=['POST'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def queryByDomain(database):
    """
    expected input:
    {
        domains: [<string>... ] where <string> is a
            pfam accession e.g. PF00123
    }
    """
    db = getDb(database)
    if not request.json:
        return abort(400)
    if 'domains' not in request.json:
        return bad_request(msg='domains not given')
    domains = request.json.get('domains')
    if not domains or len(domains) == 0:
        return bad_request(msg='empty domains')
    for d in domains:
        if not isPfamAcc(d):
            return bad_request(msg='%s is not a valid pfam accession' % d)
    thresholds = request.json.get('thresholds',[])
    try:
        check_threshold(thresholds, PFAM_ALLOWED_FIELDS)
    except Exception as e:
        return bad_request(msg=str(e))
    hits = getNodeIdFromDomains(db, domains, thresholds=thresholds)
    return json.dumps(hits)


@app.route('/<database>/treeNodes/by/tigrfam', methods=['OPTIONS'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def allowCORSdomainsTigrfam(database):
    return True

def isTigrfamId(text):
    return bool(re.match('^TIGR\d{5}$', text))


@app.route('/<database>/treeNodes/by/tigrfam', methods=['POST'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def queryByTigrfam(database):
    """
    expected input:
    {
        domains: [<string>... ] where <string> is a
            pfam accession e.g. PF00123
    }
    """
    db = getDb(database)
    if not request.json:
        return abort(400)
    if 'domains' not in request.json:
        return bad_request(msg='domains not given')
    domains = request.json.get('domains')
    if not domains or len(domains) == 0:
        return bad_request(msg='empty domains')
    for d in domains:
        if not isTigrfamId(d):
            return bad_request(msg='%s is not a valid tigrfam id' % d)
    thresholds = request.json.get('thresholds',[])
    try:
        check_threshold(thresholds, TIGRFAM_ALLOWED_FIELDS)
    except Exception as e:
        return bad_request(msg=str(e))
    hits = getNodeIDFromTigrfam(db, domains, thresholds=thresholds)
    return json.dumps(hits)



@app.route('/<database>/treeNodes/by/keggs', methods=['OPTIONS'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def allowCORSdomainsKEGG(database):
    return True


def isKeggId(text):
    # '^': start of string
    # 'K\d{5}': kegg ID format
    # '(?:&K\d{5}){0,}': any number of kegg IDs preceded with '&'
    # '?:': prevent storage of capture group
    # '$': end of string
    return bool(re.match('^K\d{5}(?:&K\d{5}){0,}$', text))


@app.route('/<database>/treeNodes/by/keggs', methods=['POST'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def queryByKEGG(database):
    """
    expected input:
    {
        keggs: [<string>... ] where <string> is a
            KEGG id e.g. K00123,
        thresholds:[{
            fieldname: one of ['percent_identity', 'eval', 'query_percent_alignment'],
            value: float
        }]
    }
    """
    db = getDb(database)
    if not request.json:
        return abort(400)
    if 'keggs' not in request.json:
        return bad_request(msg='keggs not given')
    keggs = request.json.get('keggs')
    if len(keggs) == 0:
        return bad_request(msg='empty keggs')
    for d in keggs:
        if not d or not isKeggId(d):
            return bad_request(msg='%s is not a valid kegg id' % d)
    thresholds = request.json.get('thresholds',[])
    try:
        check_threshold(thresholds, KEGG_ALLOWED_FIELDS)
    except Exception as e:
        return bad_request(msg=str(e))
    hits = getNodeIdFromKEGG(db, keggs, thresholds=thresholds)
    return json.dumps(hits)


@app.route('/<database>/treeNodes/by/taxIds', methods=['OPTIONS'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def allowCORSTaxIds(database):
    return True


@app.route('/<database>/treeNodes/by/taxIds', methods=['POST'])
@crossdomain(origin='*')
def queryByTaxIds(database):
    # print("1backend test by steve")
    """
    _db = getDb(database)
    if request.json is None:
        abort(400)
        return
    taxids = request.json
    if len(taxids) == 0:
        return json.dumps([])
    for taxid in taxids:
        if type(taxid) != int:
            return make_error_response(msg='invalid tax ids', code=400)

    return json.dumps(getNodeIdFromTaxIds(_db, taxids))
    """
    _db = getDb(database)
    conn = MySQLdb.connect(user=config.username, passwd=config.password, port=config.port,
                           host=config.host, db="names")
    if request.json is None:
        abort(400)
        return
    taxids = request.json
    if len(taxids) == 0:
        return json.dumps([])

    # handle xml file:
    if type(taxids[0]) == int:
        for taxid in taxids:
            if type(taxid) != int:
                return make_error_response(msg='invalid tax ids in xml file', code=400)
    # handle csv file:

    elif type(taxids[0]) == unicode:
        for taxid in taxids:

            if type(taxid) != unicode:
                return make_error_response(msg='invalid tax ids in csv file', code=400)

        result = getTaxIdsFromNames(conn,taxids)
        result = [int(x) for x in result]
        result = getNodeIdFromTaxIds(_db,result)
        return json.dumps(result)
    return json.dumps(getNodeIdFromTaxIds(_db, taxids))

# /*---------- END of query treenodes   ----------*/


@app.route('/<database>/tree', methods=['OPTIONS'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def allowCORSdomainsTree(database):
    return True


@app.route('/<database>/tree', methods=['GET'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def export_tree(database):
    try:
        db = getDb(database)
    except Exception as e:
        return make_error_response(str(e))
    tree = json.dumps(get_tree(db))
    if tree is None:
        return make_error_response("Tree retrieval failed")
    resp = make_response(tree, 200)
    resp.headers['Cache-Control'] = 'max-age=864000'  # 10 days of cache

    return resp



# /*----------  Auto complete  ----------*/


def checkValidAutocomplete(phrase):
    if not phrase:
        return jsonify([])
    if not phrase or len(phrase) == 0:
        return jsonify([])
    elif re.search('%', phrase):
        return make_error_response(msg='invalid autocomplete phrase, percent sign is not allowed', code=400)
    return None


@app.route('/<database>/pfamDomain/autocomplete', methods=['OPTIONS'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def allowCORSAutocomplete(database):
    return True


@app.route('/<database>/pfamDomain/autocomplete', methods=['GET'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def pfamDomainAutocomplete(database):
    phrase = request.args.get('query')
    err_msg = checkValidAutocomplete(phrase)
    if err_msg:
        return err_msg
    db = getDb(database)
    c = db.cursor(MySQLdb.cursors.DictCursor)
    sql = """SELECT DISTINCT(c.pfamA_acc) AS pfamA_acc, c.description AS description, c.pfamA_id AS pfamA_id
    FROM (SELECT description, pfamA_acc, pfamA_id
	      FROM pfamA
	      WHERE LOWER(description) LIKE %s
          OR LOWER(pfamA_acc) LIKE %s) AS c
    INNER JOIN pfam_top_hits AS h
    ON c.pfamA_acc = h.pfam_id
    LIMIT 10;"""
    searchPhrase = phrase.lower() + '%'  # so that it matches anything starting with `phrase`
    c.execute(sql, ['%' + searchPhrase, searchPhrase])
    rows = c.fetchall()
    return json.dumps(rows)


@app.route('/<database>/taxonomy/autocomplete', methods=['OPTIONS'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def allowTaxonomyCORSAutocomplete(database):
    return True


@app.route('/<database>/taxonomy/autocomplete', methods=['GET'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def taxonomyAutocomplete(database):
    phrase = request.args.get('query')
    err_msg = checkValidAutocomplete(phrase)
    if err_msg:
        return err_msg
    db = getDb(database)
    c = db.cursor(MySQLdb.cursors.DictCursor)
    sql = """SELECT ncbi_taxid as taxId, species FROM taxonomy WHERE LOWER(species)
    LIKE %s AND species NOT LIKE '%%virus%%' AND species NOT LIKE '%%phage%%'
    AND ncbi_taxid IN (SELECT species_id from node_species where species_id IS NOT NULL) LIMIT 10;"""
    searchPhrase = phrase.lower() + '%'  # so that it matches anything starting with `phrase`
    c.execute(sql, [searchPhrase])
    rows = c.fetchall()
    return json.dumps(rows)


@app.route('/<database>/kegg/autocomplete', methods=['OPTIONS'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def allowKeggCORSAutocomplete(database):
    return True


@app.route('/<database>/kegg/autocomplete', methods=['GET'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def keggAutocomplete(database):
    phrase = request.args.get('query')
    err_msg = checkValidAutocomplete(phrase)
    if err_msg:
        return err_msg
    db = getDb(database)
    c = db.cursor(MySQLdb.cursors.DictCursor)
    sql = """SELECT DISTINCT(h.kegg_id) AS keggId, c.definition AS description
    FROM (SELECT kegg_id, definition
	      FROM kegg_definitions
	      WHERE LOWER(definition) LIKE %s
          OR LOWER(kegg_id) LIKE %s) AS c
    INNER JOIN kegg_top_hits AS h
    ON h.kegg_id = c.kegg_id
    LIMIT 10;"""
    searchPhrase = phrase.lower() + '%'  # so that it matches anything starting with `phrase`
    c.execute(sql, ['%' + searchPhrase, searchPhrase])
    rows = c.fetchall()
    return json.dumps(rows)


@app.route('/<database>/tigrfam/autocomplete', methods=['OPTIONS'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def allowTigrfamCORSAutocomplete(database):
    return True


@app.route('/<database>/tigrfam/autocomplete', methods=['GET'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def tigrfamAutocomplete(database):
    phrase = request.args.get('query')
    err_msg = checkValidAutocomplete(phrase)
    if err_msg:
        return err_msg
    db = getDb(database)
    c = db.cursor(MySQLdb.cursors.DictCursor)
    sql = """SELECT DISTINCT(c.tigrfam_id) AS tigrfamId, c.definition AS description
    FROM (SELECT tigrfam_id, definition
	      FROM tigrfam_definitions
	      WHERE LOWER(tigrfam_id) LIKE %s OR
	      LOWER(definition) LIKE %s) AS c
    INNER JOIN tigrfam_top_hits AS h
    ON c.tigrfam_id = h.tigrfam_id
    LIMIT 10;"""
    c.execute(sql, [phrase.lower() + '%', '%'+phrase.lower()+'%'])
    rows = c.fetchall()
    return json.dumps(rows)

# /*----------  END of auto complete  ----------*/


# /*----------  detailed results  ----------*/
def _getPfamScanResults(db, domains, gtdb_ids, size_limit, with_sequence=False, thresholds=[]):
    '''
    return [{
        gtdbId: string,
        geneId: string,
        pfamId: string,
        eval: float,
        bitscore: float,
        sequence: string [optional if with_sequence=True],
    }]
    '''
    domain_clause = ','.join(['\''+d+'\'' for d in domains])
    gtdb_clause = ','.join(['\''+d+'\'' for d in gtdb_ids])
    limit_clause = (' LIMIT ' + str(min(int(size_limit), 999999))) if size_limit else ''
    threshold_clause, threshold_values = getThresholdClause(thresholds)
    sql = ''
    if not with_sequence:
        sql = '''
        SELECT gtdb_id AS gtdbId, gene_id AS geneId,
        pfam_id AS pfamId,eval,bitscore FROM pfam_top_hits
        WHERE
            pfam_id IN ({domain_clause}) AND gtdb_id IN ({gtdb_clause})
            {threshold_clause}
        {limit_clause}
        ;
        '''.format(domain_clause=domain_clause,
            gtdb_clause=gtdb_clause,
            limit_clause=limit_clause,
            threshold_clause=threshold_clause)
    else:
        sql = '''
        SELECT pt.gtdb_id AS gtdbId, pt.gene_id AS geneId,
        pt.pfam_id AS pfamId,pt.eval,pt.bitscore, ps.sequence
        FROM pfam_top_hits pt
        JOIN protein_sequences ps
        ON pt.gene_id = ps.gene_id AND pt.gtdb_id = ps.gtdb_id
        WHERE
            pt.pfam_id IN ({domain_clause}) AND pt.gtdb_id IN ({gtdb_clause})
            {threshold_clause}
        {limit_clause}
        ;
        '''.format(domain_clause=domain_clause,
            gtdb_clause=gtdb_clause,
            limit_clause=limit_clause,
            threshold_clause=threshold_clause)
    c = db.cursor(MySQLdb.cursors.DictCursor)
    c.execute(sql, threshold_values)
    rows = c.fetchall()
    for r in rows:
        r['eval'] = float(r['eval'])
        r['bitscore'] = float(r['bitscore'])
    return rows


def _getKeggResults(db, keggs, gtdb_ids, size_limit=None, with_sequence=False, thresholds=[]):
    '''
    return [{
        gtdbId: string,
        geneId: string,
        keggId: string,
        eval: float,
        percentIdentity: float,
        queryPercentAlignment: float,
        subjectPercentAlignment: float,
        bitscore: float,
        sequence: string [optional if with_sequence=True],
    }]
    '''
    kegg_clause = ','.join(['\''+d+'\'' for d in keggs])
    gtdb_clause = ','.join(['\''+d+'\'' for d in gtdb_ids])
    limit_clause = (' LIMIT ' + str(min(int(size_limit), 999999))) if size_limit else ''
    threshold_clause, threshold_values = getThresholdClause(thresholds)
    sql = ''
    if not with_sequence:
        sql = '''
        SELECT GROUP_CONCAT(kegg_id SEPARATOR '&') AS keggId, gene_id AS geneId,
        gtdb_id AS gtdbId,eval,bitscore,
        kt.percent_identity AS percentIdentity, kt.query_percent_alignment AS queryPercentAlignment,
        kt.subject_percent_alignment AS subjectPercentAlignment,
	nt.gtdb_taxonomy AS taxonomy
        FROM kegg_top_hits kt
	JOIN node_tax nt
	ON nt.gtdb_id = kt.gtdb_id
        WHERE
            kegg_id IN ({kegg_clause}) AND gtdb_id IN ({gtdb_clause})
            {threshold_clause}
        GROUP BY 2,3,4,5,6,7,8
        {limit_clause}
        ;
        '''.format(kegg_clause=kegg_clause,
            gtdb_clause=gtdb_clause,
            limit_clause=limit_clause,
            threshold_clause=threshold_clause)
    else:
        sql = '''
        SELECT GROUP_CONCAT(kt.kegg_id SEPARATOR '&') AS keggId, kt.gene_id AS geneId,
        kt.gtdb_id AS gtdbId,kt.eval,kt.bitscore,
        kt.percent_identity AS percentIdentity, kt.query_percent_alignment AS queryPercentAlignment,
        kt.subject_percent_alignment AS subjectPercentAlignment,
        ps.sequence,
	nt.gtdb_taxonomy AS taxonomy
        FROM kegg_top_hits kt
        JOIN protein_sequences ps
        ON kt.gene_id = ps.gene_id AND kt.gtdb_id = ps.gtdb_id
	JOIN node_tax nt
	ON nt.gtdb_id = kt.gtdb_id
        WHERE
            kt.kegg_id IN ({kegg_clause}) AND kt.gtdb_id IN ({gtdb_clause})
            {threshold_clause}
        GROUP BY 2,3,4,5,6,7,8,9
        {limit_clause}
        ;
        '''.format(kegg_clause=kegg_clause,
        gtdb_clause=gtdb_clause,
        limit_clause=limit_clause,
        threshold_clause=threshold_clause)
    c = db.cursor(MySQLdb.cursors.DictCursor)
    c.execute(sql, threshold_values)
    rows = c.fetchall()

    data_tuple = []
    for r in rows:
	data = {}
	data['keggId'] = r['keggId']
	data['geneId'] = r['geneId']
	data['gtdbId'] = r['gtdbId']
	data['eval'] = float(r['eval'])
	data['bitscore'] = float(r['bitscore'])
	data['percentIdentity'] = float(r['percentIdentity'])
	data['subjectPercentAlignment'] = float(r['subjectPercentAlignment'])
	if with_sequence:
		data['sequence'] = r['sequence']
	
	if r['taxonomy']:
		for tax in r['taxonomy'].split(';'):
			tax_arr = tax.split('__')
			if tax_arr[0] == 'd':
				data['db'] = tax_arr[1]
			elif tax_arr[0] == 'p':
				data['phylum'] = tax_arr[1]
			elif tax_arr[0] == 'c':
				data['class'] = tax_arr[1]
			elif tax_arr[0] == 'o':
				data['order'] = tax_arr[1]
			elif tax_arr[0] == 'f':
				data['family'] = tax_arr[1]
			elif tax_arr[0] == 'g':
				data['genus'] = tax_arr[1]
			else:
				data['species'] = tax_arr[1]
	else:
		data['db'] = 'None'
		data['phylum'] = 'None'
		data['class'] = 'None'
		data['order'] = 'None'
		data['family'] = 'None'
		data['genus'] = 'None'
		data['species'] = 'None'
	data_tuple.append(data)

    data_tuple = tuple(data_tuple)
    return data_tuple


def _getTigrfamScanResults(db, domains, gtdb_ids, size_limit, with_sequence=False, thresholds=[]):
    '''
    return [{
        gtdbId: string,
        geneId: string,
        pfamId: string,
        eval: float,
        bitscore: float,
        sequence: string [optional if with_sequence=True],
    }]
    '''
    domain_clause = ','.join(['\''+d+'\'' for d in domains])
    gtdb_clause = ','.join(['\''+d+'\'' for d in gtdb_ids])
    limit_clause = (' LIMIT ' + str(min(int(size_limit), 999999))) if size_limit else ''
    threshold_clause, threshold_values = getThresholdClause(thresholds)
    sql = ''
    if not with_sequence:
        sql = '''
        SELECT gtdb_id AS gtdbId, gene_id AS geneId,
        tigrfam_id AS tigrfamId,eval,bitscore FROM tigrfam_top_hits
        WHERE
            tigrfam_id IN ({domain_clause}) AND gtdb_id IN ({gtdb_clause})
            {threshold_clause}
        {limit_clause}
        ;
        '''.format(domain_clause=domain_clause,
                gtdb_clause=gtdb_clause,
                limit_clause=limit_clause,
                threshold_clause=threshold_clause)
    else:
        sql = '''
        SELECT top_hits.gtdb_id AS gtdbId, top_hits.gene_id AS geneId,
        top_hits.tigrfam_id AS tigrfamId,top_hits.eval,top_hits.bitscore, ps.sequence
        FROM tigrfam_top_hits top_hits
        JOIN protein_sequences ps
        ON top_hits.gene_id = ps.gene_id AND top_hits.gtdb_id = ps.gtdb_id
        WHERE
            top_hits.tigrfam_id IN ({domain_clause}) AND top_hits.gtdb_id IN ({gtdb_clause})
            {threshold_clause}
        {limit_clause}
        ;
        '''.format(domain_clause=domain_clause,
                gtdb_clause=gtdb_clause,
                limit_clause=limit_clause,
                threshold_clause=threshold_clause)
    c = db.cursor(MySQLdb.cursors.DictCursor)
    c.execute(sql, threshold_values)
    rows = c.fetchall()
    for r in rows:
        r['eval'] = float(r['eval'])
        r['bitscore'] = float(r['bitscore'])
    return rows


@app.route('/<database>/pfamScanResults', methods=['OPTIONS'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def allowCORSdomainsPfamScanResults(database):
    return True


def is_gtdb_id(gtdb_id):
    return bool(re.match('^[A-Z_0-9.]+$', gtdb_id))


@app.route('/<database>/pfamScanResults', methods=['POST'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def getPfamScanResults(database):
    """
    expected input:
    {
        domains: [<string>... ] where <string> is a
            pfam accession e.g. PF00123,
        gtdbIds: ["UBA12301"]
        withSequence: true|false,
        sizeLimit: <int>
    }
    """
    db = getDb(database)
    if not request.json:
        return abort(400)
    if 'domains' not in request.json:
        return bad_request(msg='domains not given')
    domains = request.json.get('domains')
    gtdb_ids = request.json.get('gtdbIds')
    with_sequence = request.json.get('withSequence')
    size_limit = request.json.get('sizeLimit')
    if size_limit < 0:
        return bad_request(msg='size limit should be > 0')
    if len(domains) == 0:
        return bad_request(msg='empty domains')
    if len(gtdb_ids) == 0:
        return bad_request(msg='empty gtdbIds')
    for d in domains:
        if not isPfamAcc(d):
            return bad_request(msg='%s is not a valid pfam accession' % d)
    for g in gtdb_ids:
        if not is_gtdb_id(g):
            return bad_request(msg='%s is not a valid gtdb id' % g)
    thresholds = request.json.get('thresholds', [])
    try:
        check_threshold(thresholds, PFAM_ALLOWED_FIELDS)
    except Exception as e:
        return bad_request(msg=str(e))
    hits = _getPfamScanResults(db, domains, gtdb_ids, size_limit,
        with_sequence=with_sequence, thresholds=thresholds)
    return json.dumps(hits)


@app.route('/<database>/keggResults', methods=['OPTIONS'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def allowCORSdomainsKeggResults(database):
    return True


@app.route('/<database>/keggResults', methods=['POST'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def getKeggResults(database):
    """
    expected input:
    {
        keggs: [<string>... ]
        gtdbIds: ["UBA12301" ...]
        withSequence: true|false
        sizeLimit: <int>
    }
    """
    db = getDb(database)
    if not request.json:
        return abort(400)
    if 'keggs' not in request.json:
        return bad_request(msg='keggs not given')
    keggs = request.json.get('keggs')
    gtdb_ids = request.json.get('gtdbIds')
    with_sequence = request.json.get('withSequence')
    size_limit = request.json.get('sizeLimit')
    thresholds = request.json.get('thresholds', [])
    try:
        check_threshold(thresholds, KEGG_ALLOWED_FIELDS)
    except Exception as e:
        return bad_request(msg=str(e))
    if size_limit and size_limit < 0:
        return bad_request(msg='size limit should be > 0')
    if len(keggs) == 0:
        return bad_request(msg='empty keggs')
    if len(gtdb_ids) == 0:
        return bad_request(msg='empty gtdbIds')
    for k in keggs:
        if not isKeggId(k):
            return bad_request(msg='%s is not a valid kegg id' % k)
    for g in gtdb_ids:
        if not is_gtdb_id(g):
            return bad_request(msg='%s is not a valid gtdb id' % g)
    hits = _getKeggResults(db,
        keggs,
        gtdb_ids,
        size_limit,
        with_sequence=with_sequence,
        thresholds=thresholds)
    return json.dumps(hits)


@app.route('/<database>/tigrfamResults', methods=['OPTIONS'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def allowCORSdomainsTigrfamResults(database):
    return True


@app.route('/<database>/tigrfamResults', methods=['POST'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def getTigrfamResults(database):
    """
    expected input:
    {
        domains: [<string>... ]
        gtdbIds: ["UBA12301" ...]
        withSequence: true|false
        sizeLimit: <int>
    }
    """
    db = getDb(database)
    if not request.json:
        return abort(400)
    if 'domains' not in request.json:
        return bad_request(msg='domains not given')
    domains = request.json.get('domains')
    gtdb_ids = request.json.get('gtdbIds')
    with_sequence = request.json.get('withSequence')
    size_limit = request.json.get('sizeLimit')
    if size_limit and size_limit < 0:
        return bad_request(msg='size limit should be > 0')
    if len(domains) == 0:
        return bad_request(msg='empty domains')
    if len(gtdb_ids) == 0:
        return bad_request(msg='empty gtdbIds')
    for d in domains:
        if not isTigrfamId(d):
            return bad_request(msg='%s is not a valid tigrfam id' % d)
    for g in gtdb_ids:
        if not is_gtdb_id(g):
            return bad_request(msg='%s is not a valid gtdb id' % g)
    thresholds = request.json.get('thresholds', [])
    try:
        check_threshold(thresholds, TIGRFAM_ALLOWED_FIELDS)
    except Exception as e:
        return bad_request(msg=str(e))
    hits = _getTigrfamScanResults(db, domains, gtdb_ids, size_limit,
        with_sequence=with_sequence, thresholds=thresholds)
    return json.dumps(hits)

# /*----------  end of detailed results  ----------*/


@app.route('/<database>/version', methods=['OPTIONS'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def allowCORSdomainsVersion(database):
    return True


@app.route('/<database>/version', methods=['GET'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def getVersion(database):
    db = getDb(database)
    sql = '''SELECT config_param, file_name
            FROM {db_name}.db_config_data_files
            WHERE config_param IN ('json_tree', 'pfamA_sql'); '''\
            .format(db_name=checkDbName(database))
    c = db.cursor(MySQLdb.cursors.DictCursor)
    c.execute(sql)
    rows = c.fetchall()
    version_info = {}
    print rows
    for r in rows:
        config_param = r['config_param']
        file_name = r['file_name']
        if config_param == 'json_tree':
            m = re.search('r(\d+)[^\d]', file_name)
            if not m:
                return make_error_response('GTDB json tree file name does not contain version info')
            version_info['gtdb'] = m.group(1)
        elif config_param == 'pfamA_sql':
            m = re.search('v(\d+)[^\d]', file_name)
            if not m:
                return make_error_response('Pfam SQL file name does not contain version info')
            version_info['pfam'] = m.group(1)
    return json.dumps(version_info)

if __name__ == '__main__':
    f = open('log.txt', 'a')
    timestr = time.strftime("%Y.%m.%d-%H:%M:%S")
    f.write(timestr+'\n')
    f.write(sys.version)
    f.write('running\n')
    f.close()

    app.run(host="0.0.0.0", port=5001, debug=True, threaded=True)
