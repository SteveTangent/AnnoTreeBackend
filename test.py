import MySQLdb

mydb = mysql.connector.connect(
  host="localhost",
  port=3306,
  user="annotree",
  password="Ann0tr33",
  bacterial_database='gtdb_bacteria_r95',
  archaeal_database='gtdb_archaea_r95'
)

mycursor = mydb.cursor()

sql = '''
        SELECT GROUP_CONCAT(kegg_id SEPARATOR '&') AS keggId, gene_id AS geneId,
        gtdb_id AS gtdbId,eval,bitscore,
        kt.percent_identity AS percentIdentity, kt.query_percent_alignment AS queryPercentAlignment,
        kt.subject_percent_alignment AS subjectPercentAlignment
        FROM kegg_top_hits kt
        GROUP BY 2,3,4,5,6,7,8
        limit 50
        ;
        '''

mycursor.execute(sql)

myresult = mycursor.fetchall()

for x in myresult:
  print(x)
