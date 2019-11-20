import json
import os
path = 'Data/Images/'
# path = '/data/ProteomeToolsRaw/Images/'
try:
	os.remove(path+'metadata_filtered.json')
except Exception:
	print('no filtered version exist')

lines_seen = set()
outfile = open(path+'metadata_filtered.json','w')
for line in open(path+'metadata.json','r'):
	data = json.loads(line)
	#Add Filter here#
	# if float(data['m/z']) > 424.95 and float(data['m/z']) < 425:
	if line not in lines_seen:
		outfile.write(line)
		lines_seen.add(line)
outfile.close