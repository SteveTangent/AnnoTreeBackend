import MySQLdb
'''
    export gtdb tree from database to a dictionary representation
'''

def get_rows(db):
    '''
    given a database connection, get all rows from node table
    joined with relevant information
    '''
    sql = """
    SELECT n.id, n.level, n.parent_id AS parentId, n.is_leaf AS isLeaf, n.length,
            r.rank, t.gtdb_taxonomy AS gtdbTaxonomy
        FROM node n
        LEFT JOIN node_gtdb_ranks r ON r.node_id = n.id
        LEFT JOIN node_tax t ON t.gtdb_id = n.gtdb_id
        WHERE (r.is_highest = 1 OR r.node_id IS NULL)
    """.format()
    c = db.cursor(MySQLdb.cursors.DictCursor)
    c.execute(sql)
    rows = c.fetchall()
    return rows

def get_taxonomy_levels(taxonomy_str):
    d = {
        'd': 'domain',
        'p': 'phylum',
        'c': 'class',
        'o': 'order',
        'f': 'family',
        'g': 'genus',
        's': 'species'
    }
    rank_and_names = taxonomy_str.split(';')
    ranks = []
    for r in rank_and_names:
        temp_rank = d[r.split('__')[0]]
        temp_rank_name = r.split('__')[1]
        ranks.append({'rank': temp_rank, 'rank_name': temp_rank_name})
    return ranks

def process_rows(rows):
    '''
        process each row to a node
    '''
    for r in rows:
        gtdb_taxonomy = r['gtdbTaxonomy']
        if gtdb_taxonomy:
            r['gtdbTaxonomy'] = get_taxonomy_levels(gtdb_taxonomy)
        r['id'] = int(r['id'])
        r['parentId'] = int(r['parentId']) if r['parentId'] else None
        r['length'] = float(r['length']) if r['length'] is not None else None
        r['children'] = []
        r['isLeaf'] = int(r['isLeaf'])


def count_tree(node):
    count = sum([count_tree(c) for c in node['children']])
    node['num_child'] = count
    if len(node['children']) == 0:
        return 1
    return count


def link_tree(rows):
    '''
        link tree node according to parent id
    '''
    m = {r['id']:r for r in rows}
    root = None
    for r in rows:
        if not r['parentId']:
            root = r
            continue
        parent = m[r['parentId']]
        parent['children'].append(r)
    return root


def get_tree(db):
    '''
        tree json looks like this, each dictionary is a node
        {
            id: <int>
            level: string
            parentId: int or None,
            isLeaf: string or None,
            length: float,
            rank: string,
            gtdbTaxonomy: [{rank: string, rank_name:string} ...],
            num_child: int,
            children: [ node  ...]
        }
    '''
    rows = get_rows(db)
    process_rows(rows)
    root = link_tree(rows)
    count_tree(root)
    return root