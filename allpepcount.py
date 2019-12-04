import json
accessions = []
i=0
for line in open('accessions.json'):
	try:
		data = json.loads(line)
	except Exception:
		pass
	try:
		if data['allpeptides'] == True:
			accessions.append(data['accession'])
			i+=1
	except Exception:
		pass
print(i)
# print(accessions) 	