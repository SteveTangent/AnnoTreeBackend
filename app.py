#!flask/bin/python3

import os
import sys
import argparse
import re
from flask import Flask, jsonify, abort, make_response, request, url_for, current_app
from datetime import timedelta
from functools import update_wrapper
import time

import mysql.connector
from mysql.connector import FieldType

# import MySQLdb
# from MySQLdb.constants import FIELD_TYPE
import config
# from tree_exporter import get_tree
from new_tree_exporter import get_tree
import json
import itertools


app = Flask(__name__)

KEGG_ALLOWED_FIELDS = {'percent_identity', 'eval', 'query_percent_alignment', 'subject_percent_alignment'}
PFAM_ALLOWED_FIELDS = {'eval'}
TIGRFAM_ALLOWED_FIELDS = {'eval'}

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
    _my_conv = {FieldType.LONG: int, FieldType.INT24: int}
    # conn = mysql.connector.connect(user=username, passwd=password, port=port,
    #                        host=host, db=database, conv=_my_conv)

    conn = mysql.connector.connect(user=username, passwd=password,
                           host=host, db=database)
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
    if headers is not None and not isinstance(headers, str):
        headers = ', '.join(x.upper() for x in headers)
    if not isinstance(origin, str):
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
        for i in range(len(thresholds)):
            thresholds[i]['value'] = float(thresholds[i]['value'])
    except Exception as e:
        raise Exception('Value should be a number')



def _get_node_ids_from_tophits(db, top_hit_table, search_colname, search_ids, thresholds):
    """
    search top_hit_table for genomes' `node_id`
    where genomes must contain all search_ids, and each hit
    must have scores satisfying `thresholds`
    """
    # print(search_ids,search_colname)
    search_id_quoted = (",".join(\
        ["'%s'" % (search_id) for search_id in search_ids]))
    search_clause = "%s IN (%s)" % (search_colname, search_id_quoted)


    """
    kegg should be kegg_id, pfam and tigrfam are subfamily of interpro, so it should be parsed separately
    node_list is a list of sets, we search each individual id first, and then use python to find the intersection
    among sets
    """
    node_list = []
    c = db.cursor()

    for each_id in search_ids:
        sql = ""

        if search_colname == 'kegg_id':
            sql = f'''
            SELECT distinct ka.node_id
            FROM kegg_annotation ka
            WHERE ka.kegg_id = "{each_id}";
            '''
        elif search_colname == 'pfam_id':
            sql = f'''
            SELECT distinct ia.node_id
            FROM interpro_annotation ia
            WHERE ia.member_id = "{each_id}";
            '''

        # now should support kegg, interpro and tigrfam
        elif search_colname == 'tigrfam_id':
            cur_database = check_if_valid_id(each_id)

            if cur_database == 'KEGG':
                sql = f'''
                SELECT distinct ka.node_id
                FROM kegg_annotation ka
                WHERE ka.kegg_id = "{each_id}";
                '''
            elif cur_database == 'InterPro':
                sql = f'''
                SELECT distinct ia.node_id
                FROM interpro_annotation ia
                WHERE ia.interpro_id = "{each_id}";
                '''
            else:
                sql = f'''
                SELECT distinct ia.node_id
                FROM interpro_annotation ia
                WHERE ia.member_id = "{each_id}";
                '''


        # c.execute(sql, threshold_values)
        c.execute(sql)
        rows = c.fetchall()
        cur_set = set([x[0] for x in rows])
        node_list.append(cur_set)

    intersection = set.intersection(*node_list)

#     if search_colname == 'kegg_id':
#         sql = f"""
#         SELECT node_id
# FROM node n
# JOIN
#     (
#     SELECT
#         gtdb_id,
#         COUNT(DISTINCT kegg_id) AS num_hit_per_genome
#     FROM kegg_annotation
#     WHERE kegg_id in ({search_id_quoted})
#     GROUP BY gtdb_id
#     HAVING num_hit_per_genome >= {str(len(search_ids))}
#     ) g
# ON n.gtdb_id = g.gtdb_id
# """
#     elif search_colname == 'pfam_id':
#         sql = f'''
# SELECT node_id
# FROM node n
# JOIN(
# SELECT
#     ia.gtdb_id,
#     COUNT(DISTINCT ia.member_id) AS num_hit_per_genome
# FROM
#     interpro_annotation ia
# WHERE ia.member_id in ({search_id_quoted})
# GROUP BY ia.gtdb_id
# HAVING num_hit_per_genome >= {str(len(search_ids))}
#     ) g
# ON n.gtdb_id = g.gtdb_id
#     '''
#     elif search_colname == 'tigrfam_id':
#
#         sql = f'''
# SELECT node_id
# FROM node n
# JOIN(
# SELECT
#     ia.gtdb_id,
#     COUNT(DISTINCT ia.interpro_id) AS num_hit_per_genome
# FROM
#     interpro_annotation ia
# WHERE ia.interpro_id in ({search_id_quoted})
# GROUP BY ia.gtdb_id
# HAVING num_hit_per_genome >= {str(len(search_ids))}
#     ) g
# ON n.gtdb_id = g.gtdb_id
#     '''
#
#     print(sql)

    # c = db.cursor()
    # c.execute(sql)
    # rows = c.fetchall()
    # if not rows:  # no result
    #     return []
    # node_ids = [r[0] for r in rows]
    # print(node_ids)
    node_ids = list(intersection)
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
    # sql = '''SELECT node_id FROM node_species WHERE species_id IN
    #         (SELECT get_species(ncbi_taxid) FROM taxonomy WHERE ncbi_taxid IN (%s)
    #         AND get_species(ncbi_taxid) != 1);''' % ','.join([str(t) for t in taxids])


    search_ids = ','.join([str(t) for t in taxids])
    sql = f'''
    SELECT node_id 
FROM node AS n 
JOIN ncbi_meta AS nm ON nm.accession = n.gtdb_id 
WHERE nm.ncbi_taxid IN ({search_ids});
'''
    # print(sql)
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
    # sql = 'SELECT taxid FROM names WHERE scientific_name IN (%s);' % ','.join(['"%s"' % str(n) for n in names])

    search_names = ','.join(['"%s"' % str(n) for n in names])
    sql = f'''
    SELECT DISTINCT ncbi_taxid AS taxid FROM ncbi_meta WHERE ncbi_organism_name in ({search_names})
    '''
    # print(sql)
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
    '''
    Steve:
     change tigrfam to interpro
    '''
    return bool(re.match('^TIGR\d{5}$', text))
    # return bool(re.match('^IPR\d{6}$', text))


@app.route('/<database>/treeNodes/by/tigrfam', methods=['POST'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def queryByTigrfam(database):
    """
    expected input:
    {
        domains: [<string>... ] where <string> is a
            pfam accession e.g. PF00123
    }

    Steve add:
    now support interpro instead of tigrfam
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
        # if not isTigrfamId(d):
        if not check_if_valid_id(d):
            return bad_request(msg='%s is not a valid tigrfam id' % d)
            # return bad_request(msg='%s is not a valid interpro id' % d)
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
    # _db = getDb(database)
    # nodeId = request.json
    # print(nodeId)
    #
    # return json.dumps(nodeId)
    # conn = mysql.connector.connect(user=config.username, passwd=config.password, port=config.port,
    #                        host=config.host, db="names")
    # conn = MySQLdb.connect(user=config.username, passwd=config.password, port=config.port,
    #                        host=config.host, db="names")

    db = getDb(database)
    if request.json is None:
        abort(400)
        return
    taxids = request.json
    print(taxids)
    if len(taxids) == 0:
        return json.dumps([])

    # handle xml file:
    if type(taxids[0]) == int:
        for taxid in taxids:
            if type(taxid) != int:
                return make_error_response(msg='invalid tax ids in xml file', code=400)
    # handle csv file:

    elif type(taxids[0]) == str:
        for taxid in taxids:

            if type(taxid) != str:
                return make_error_response(msg='invalid tax ids in csv file', code=400)

    #     result = getTaxIdsFromNames(conn,taxids)
    #     result = [int(x) for x in result]
    #     result = getNodeIdFromTaxIds(_db,result)
    #     return json.dumps(result)
    # return json.dumps(getNodeIdFromTaxIds(_db, taxids))

        result = getTaxIdsFromNames(db,taxids)
        result = [int(x) for x in result]
        result = getNodeIdFromTaxIds(db,result)
        return json.dumps(result)
    return json.dumps(getNodeIdFromTaxIds(db, taxids))

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

    c = db.cursor(dictionary=True)
    # c = db.cursor(MySQLdb.cursors.DictCursor)
    # sql = """SELECT DISTINCT(c.pfamA_acc) AS pfamA_acc, c.description AS description, c.pfamA_id AS pfamA_id
    # FROM (SELECT description, pfamA_acc, pfamA_id
	#       FROM pfamA
	#       WHERE LOWER(description) LIKE %s
    #       OR LOWER(pfamA_acc) LIKE %s) AS c
    # INNER JOIN pfam_top_hits AS h
    # ON c.pfamA_acc = h.pfam_id
    # LIMIT 10;"""
    print(phrase)

    if len(phrase) < 2:
        return jsonify([])

    sql = ""
    if len(phrase) >= 2 and phrase[:2].lower() == 'pf':
        sql = f'''
            SELECT DISTINCT(member_id) AS pfamA_acc, id.description AS description
            FROM interpro_def id
            WHERE member_id LIKE '{phrase}%' LIMIT 10
        '''
    else:

        sql = f'''
            SELECT DISTINCT(member_id) AS pfamA_acc, id.description AS description
            FROM interpro_def id
            WHERE member_id LIKE 'PF%' AND id.description LIKE '%{phrase}%' LIMIT 10
        '''

    searchPhrase = phrase.lower() + '%'  # so that it matches anything starting with `phrase`
    # c.execute(sql, ['%' + searchPhrase, searchPhrase])
    # c.execute(sql, [searchPhrase])
    c.execute(sql)
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
    c = db.cursor(dictionary=True)
    # c = db.cursor(MySQLdb.cursors.DictCursor)
    # sql = """SELECT ncbi_taxid as taxId, species FROM taxonomy WHERE LOWER(species)
    # LIKE %s AND species NOT LIKE '%%virus%%' AND species NOT LIKE '%%phage%%'
    # AND ncbi_taxid IN (SELECT species_id from node_species where species_id IS NOT NULL) LIMIT 10;"""

    searchPhrase = '%' + phrase.lower() + '%'  # so that it matches anything starting with `phrase`
    sql = f'''
    SELECT ncbi_taxid AS taxId, ncbi_organism_name AS species
    FROM ncbi_meta WHERE ncbi_organism_name LIKE '%{searchPhrase}%' LIMIT 10;
    '''

    c.execute(sql)
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
    c = db.cursor(dictionary=True)
    # sql = """SELECT DISTINCT(h.kegg_id) AS keggId, c.definition AS description
    # FROM (SELECT kegg_id, definition
	#       FROM kegg_definitions
	#       WHERE LOWER(definition) LIKE %s
    #       OR LOWER(kegg_id) LIKE %s) AS c
    # INNER JOIN kegg_top_hits AS h
    # ON h.kegg_id = c.kegg_id
    # LIMIT 10;"""
#     sql = '''
#     SELECT DISTINCT(kegg_id) AS keggId, kegg_definition AS description
# FROM kegg_annotation
# WHERE kegg_id LIKE %s
# OR kegg_definition LIKE %s LIMIT 10
#     '''
    search_phrase_kegg = phrase.lower() + '%'
    search_phrase_def = '%' + phrase.lower() + '%'
    print(phrase)
    sql = f'''
        SELECT DISTINCT(kegg_id) AS keggId, kegg_definition AS description
    FROM kegg_def
    WHERE kegg_id LIKE '{search_phrase_kegg}' 
    OR kegg_definition LIKE '{search_phrase_def}' LIMIT 10
        '''
    # print(sql)
    # searchPhrase = phrase.lower() + '%'  # so that it matches anything starting with `phrase`
    # c.execute(sql, ['%' + searchPhrase, searchPhrase])
    c.execute(sql)
    rows = c.fetchall()
    # for r in rows:
    #     print(r)
    return json.dumps(rows)


@app.route('/<database>/tigrfam/autocomplete', methods=['OPTIONS'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def allowTigrfamCORSAutocomplete(database):
    return True

# dictionary: prefix,database_name
member_prefix = {
    'PF': 'Pfam',
    'PTHR': 'PANTHER',
    'G3DSA:': 'CATH-Gene3D',
    'cd': 'CDD',
    'MF_': 'HAMAP',
    'TIGR': 'NCBIfam',
    'PIRSF': 'PIRSF',
    'PR': 'PRINTS',
    'PS': 'PROSITE',
    'SFLDF': 'SFLD',
    'SM': 'SMART',
    'SSF': 'SUPERFAMILY',
    'NF': 'NCBIfam',


    'K': 'KEGG',
    'IPR': 'InterPro'

}

def check_if_valid_id(cur_phrase):
    cur_phrase = cur_phrase.lower()
    prefix_list = member_prefix.keys()
    for prefix in prefix_list:
        cur_prefix = prefix.lower()
        # print(cur_prefix,cur_phrase)
        # if cur_phrase.startswith(cur_prefix) and cur_phrase[len(cur_prefix)].isdigit():
        if cur_phrase.startswith(cur_prefix):
            if len(cur_phrase) > len(cur_prefix) and cur_phrase[len(cur_prefix)].isdigit():
                return member_prefix[prefix]
    return ''

@app.route('/<database>/tigrfam/autocomplete', methods=['GET'])
@crossdomain(origin='*', headers=['content-type', 'accept'])
def tigrfamAutocomplete(database):
    phrase = request.args.get('query')
    err_msg = checkValidAutocomplete(phrase)
    if err_msg:
        return err_msg
    db = getDb(database)
    c = db.cursor(dictionary=True)
    # sql = """SELECT DISTINCT(c.tigrfam_id) AS tigrfamId, c.definition AS description

    sql = ""
    # if len(phrase) >= 4 and phrase[:4].lower() == 'tigr':
    #     sql = f'''
    #         SELECT DISTINCT(member_id) AS tigrfamId, id.description AS description
    #         FROM interpro_def id
    #         WHERE member_id LIKE '{phrase}%' LIMIT 10
    #     '''
    # else:
    #     sql = f'''
    #         SELECT DISTINCT(member_id) AS tigrfamId, id.description AS description
    #         FROM interpro_def id
    #         WHERE member_id LIKE 'tigr%' AND id.description LIKE '%{phrase}%' LIMIT 10
    #     '''

    print(phrase)
    cur_database = check_if_valid_id(phrase)

    if len(cur_database) != 0:
        # sql = f'''
        #     SELECT DISTINCT(member_id) AS tigrfamId, id.description AS description
        #     FROM interpro_def id
        #     WHERE member_id LIKE '{phrase}%' LIMIT 15
        # '''

        # since we are combine all databases
        if cur_database == 'KEGG':
            sql = f'''
                        SELECT DISTINCT(kegg_id) AS tigrfamId, kegg_definition AS description
                      FROM kegg_def
                      WHERE kegg_id LIKE '{phrase}%' LIMIT 15
                          '''
        elif cur_database == 'InterPro':
            sql = f'''
                      SELECT DISTINCT(interpro_id) AS tigrfamId, id.description AS description
                      FROM interpro_def id
                      WHERE interpro_id LIKE '{phrase}%' LIMIT 15
                  '''
        else:
            sql = f'''
                       SELECT DISTINCT(member_id) AS tigrfamId, id.description AS description
                       FROM interpro_def id
                       WHERE member_id LIKE '{phrase}%' LIMIT 15
                   '''
    else:
        sql = f'''
            (
                SELECT kegg_id AS tigrfamId, kegg_definition AS description
                FROM kegg_def
                WHERE kegg_definition LIKE '%{phrase}%'
                LIMIT 15
            )
            UNION
            (
                SELECT member_id AS tigrfamId, description AS description
                FROM interpro_def    
                WHERE description LIKE '%{phrase}%'
                LIMIT 15
            )
            LIMIT 15;'''

    # c.execute(sql, [phrase.lower() + '%', '%'+phrase.lower()+'%'])

    c.execute(sql)
    rows = c.fetchall()

    if rows:
        for item in rows:
            cur_database = check_if_valid_id(item['tigrfamId'])
            item['description'] = f"{item['description']} [{cur_database}]"

    # print(rows)
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
        sql = f'''
    SELECT ia.member_id AS pfamId, 
    ia.gtdb_id AS gtdbId,
    ia.gene_id AS geneId,
    t.taxonomy AS tax
FROM 
    interpro_annotation ia 
JOIN
	taxonomy t ON ia.gtdb_id = t.gtdb_id
WHERE 
    ia.member_id IN ({domain_clause}) AND ia.gtdb_id IN ({gtdb_clause})
    {limit_clause}
        '''
    else:

        sql = f'''
    SELECT ia.member_id AS pfamId, 
    ia.gtdb_id AS gtdbId,
    ia.gene_id AS geneId,
    c.seq AS sequence,
    t.taxonomy AS tax
FROM 
    interpro_annotation ia 
JOIN 
    coordinates c ON ia.gene_id = c.gene_id
JOIN
	taxonomy t ON ia.gtdb_id = t.gtdb_id
WHERE 
    ia.member_id IN ({domain_clause}) AND ia.gtdb_id IN ({gtdb_clause})
    {limit_clause}
        '''

    c = db.cursor(dictionary=True)
    c.execute(sql)
    rows = c.fetchall()
    # for r in rows:
    #     r['eval'] = float(r['eval'])
    #     r['bitscore'] = float(r['bitscore'])
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
    # print(db)
    # print(kegg_clause,gtdb_clause)
    # print('limit',limit_clause)
    # print(with_sequence)
    sql = ''
    if not with_sequence:
        sql = '''
                SELECT 
            ka.kegg_id AS keggId, 
            ka.gene_id AS geneId,
            ka.gtdb_id AS gtdbId,
            t.taxonomy AS tax
        FROM 
            kegg_annotation ka
        JOIN 
            taxonomy t ON t.gtdb_id = ka.gtdb_id
        WHERE
            ka.kegg_id IN ({kegg_clause}) 
        AND 
            ka.gtdb_id IN (
               {gtdb_clause}
            )
        {limit_clause};
                '''.format(kegg_clause=kegg_clause,
                           gtdb_clause=gtdb_clause,
                           limit_clause=limit_clause)
    else:
#         sql = '''
#         SELECT
#     GROUP_CONCAT(ka.kegg_id SEPARATOR '&') AS keggId,
#     ka.gene_id AS geneId,
#     ka.gtdb_id AS gtdbId,
#     t.taxonomy AS tax,
#     c.seq AS sequence
# FROM
#     kegg_annotation ka
# JOIN
#     taxonomy t ON t.gtdb_id = ka.gtdb_id
# JOIN
#     coordinates c ON c.gtdb_id = ka.gtdb_id
# WHERE
#     ka.kegg_id IN ({kegg_clause})
# AND
#     ka.gtdb_id IN (
#        {gtdb_clause}
#     )
# # GROUP BY
# #     ka.gene_id,ka.gtdb_id,t.taxonomy,c.seq
# {limit_clause};
#         '''.format(kegg_clause=kegg_clause,
#                    gtdb_clause=gtdb_clause,
#                    limit_clause=limit_clause)
            sql = '''
                SELECT 
            ka.kegg_id AS keggId, 
            ka.gene_id AS geneId,
            ka.gtdb_id AS gtdbId,
            t.taxonomy AS tax,
            c.seq AS sequence
        FROM 
            kegg_annotation ka
        JOIN 
            taxonomy t ON t.gtdb_id = ka.gtdb_id
        JOIN 
            coordinates c ON c.gtdb_id = ka.gtdb_id
        WHERE
            ka.kegg_id IN ({kegg_clause}) 
        AND 
            ka.gtdb_id IN (
               {gtdb_clause}
            )
        {limit_clause};
                '''.format(kegg_clause=kegg_clause,
                           gtdb_clause=gtdb_clause,
                           limit_clause=limit_clause)

    c = db.cursor(dictionary=True)
    # c.execute(sql, threshold_values)
    c.execute(sql)
    rows = c.fetchall()

    data_tuple = []
    for r in rows:
        data = {}
        data['keggId'] = r['keggId']
        data['geneId'] = r['geneId']
        data['gtdbId'] = r['gtdbId']
        # data['eval'] = float(r['eval'])
        # data['bitscore'] = float(r['bitscore'])
        # data['percentIdentity'] = float(r['percentIdentity'])
        # data['subjectPercentAlignment'] = float(r['subjectPercentAlignment'])
        if with_sequence:
            data['sequence'] = r['sequence']

        # if r['taxonomy']:
        if r['tax']:
            # for tax in r['taxonomy'].split(';'):
            for tax in r['tax'].split(';'):
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
    # print(domain_clause)

    # cur_database = check_if_valid_id(phrase)
    if not with_sequence:
        # sql = f'''
        # SELECT ia.member_id AS tigrfamId,
        #     ia.gtdb_id AS gtdbId,
        #     ia.gene_id AS geneId,
        #     t.taxonomy AS tax
        # FROM
        #     interpro_annotation ia
        # JOIN
        #     taxonomy t ON ia.gtdb_id = t.gtdb_id
        # WHERE
        #     ia.member_id IN ({domain_clause}) AND ia.gtdb_id IN ({gtdb_clause})
        #     {limit_clause}
        #         '''
        sql = f'''(SELECT ia.member_id AS SearchId, 
            ia.gtdb_id AS gtdbId,
            ia.gene_id AS geneId,
            t.taxonomy AS tax
        FROM 
            interpro_annotation ia 
        JOIN
            taxonomy t ON ia.gtdb_id = t.gtdb_id
        WHERE 
            ia.member_id IN ({domain_clause}) AND ia.gtdb_id IN ({gtdb_clause})
            {limit_clause})
        UNION
        (SELECT ia.interpro_id AS SearchId, 
            ia.gtdb_id AS gtdbId,
            ia.gene_id AS geneId,
            t.taxonomy AS tax
        FROM 
            interpro_annotation ia 
        JOIN
            taxonomy t ON ia.gtdb_id = t.gtdb_id
        WHERE 
            ia.interpro_id IN ({domain_clause}) AND ia.gtdb_id IN ({gtdb_clause})
            {limit_clause})
        UNION
        (  SELECT 
            ka.kegg_id AS SearchId, 
            ka.gene_id AS geneId,
            ka.gtdb_id AS gtdbId,
            t.taxonomy AS tax
        FROM 
            kegg_annotation ka
        JOIN 
            taxonomy t ON t.gtdb_id = ka.gtdb_id
        WHERE
            ka.kegg_id IN ({domain_clause}) 
        AND 
            ka.gtdb_id IN (
               {gtdb_clause}
            )
        {limit_clause})
        {limit_clause}
        '''
    else:

        # sql = f'''
        # SELECT ia.member_id AS tigrfamId,
        #     ia.gtdb_id AS gtdbId,
        #     ia.gene_id AS geneId,
        #     c.seq AS sequence,
        #     t.taxonomy AS tax
        # FROM
        #     interpro_annotation ia
        # JOIN
        #     coordinates c ON ia.gene_id = c.gene_id
        # JOIN
        #     taxonomy t ON ia.gtdb_id = t.gtdb_id
        # WHERE
        #     ia.member_id IN ({domain_clause}) AND ia.gtdb_id IN ({gtdb_clause})
        #     {limit_clause}
        #         '''
        sql = f'''        (SELECT ia.member_id AS SearchId, 
            ia.gtdb_id AS gtdbId,
            ia.gene_id AS geneId,
            c.seq AS sequence,
            t.taxonomy AS tax
        FROM 
            interpro_annotation ia 
        JOIN 
            coordinates c ON ia.gene_id = c.gene_id
        JOIN
            taxonomy t ON ia.gtdb_id = t.gtdb_id
        WHERE 
            ia.member_id IN ({domain_clause}) AND ia.gtdb_id IN ({gtdb_clause})
            {limit_clause})
        UNION
        (SELECT ia.interpro_id AS SearchId, 
            ia.gtdb_id AS gtdbId,
            ia.gene_id AS geneId,
            c.seq AS sequence,
            t.taxonomy AS tax
        FROM 
            interpro_annotation ia 
        JOIN 
            coordinates c ON ia.gene_id = c.gene_id
        JOIN
            taxonomy t ON ia.gtdb_id = t.gtdb_id
        WHERE 
            ia.interpro_id IN ({domain_clause}) AND ia.gtdb_id IN ({gtdb_clause})
            {limit_clause})
        UNION
        (  SELECT 
            ka.kegg_id AS SearchId, 
            ka.gene_id AS geneId,
            ka.gtdb_id AS gtdbId,
            t.taxonomy AS tax,
            c.seq AS sequence
        FROM 
            kegg_annotation ka
        JOIN 
            taxonomy t ON t.gtdb_id = ka.gtdb_id
        JOIN 
            coordinates c ON c.gtdb_id = ka.gtdb_id
        WHERE
            ka.kegg_id IN ({domain_clause}) 
        AND 
            ka.gtdb_id IN ({gtdb_clause})
            {limit_clause})
        {limit_clause}
        '''

    # print(sql)
    c = db.cursor(dictionary=True)
    c.execute(sql)
    rows = c.fetchall()
    # for r in rows:
    #     r['eval'] = float(r['eval'])
    #     r['bitscore'] = float(r['bitscore'])
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
        # if not isTigrfamId(d):
        if not check_if_valid_id(d):
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
    c = db.cursor(dictionary=True)
    c.execute(sql)
    rows = c.fetchall()
    version_info = {}
    # print(rows)
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
    # current_working_directory = os.getcwd()
    # log_path = current_working_directory +'/log.txt'
    # print(log_path)
    # f = open(log_path, 'a')
    # timestr = time.strftime("%Y.%m.%d-%H:%M:%S")
    # f.write(timestr+'\n')
    # f.write(sys.version)
    # f.write('running\n')
    # f.close()

    app.run(host="0.0.0.0", port=5050, debug=True, threaded=True)
