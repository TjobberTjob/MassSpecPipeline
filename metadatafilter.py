import json
import os

path = '/data/ProteomeToolsRaw/Images/'
try:
	os.system('rm '+path+'metadata_filtered.json')
except Exception:
	print('no filtered version exist')

lines_seen = set()
outfile = open(path+'metadata_filtered.json','w')
for line in open(path+'metadata.json','r'):
	data = json.loads(line)

	if float(data['m/z']) > 424.5 and float(data['m/z']) < 425: #Add filter here
		if line not in lines_seen:
			outfile.write(line)
			lines_seen.add(line)
outfile.close