import json
import sys
import os

def filter(path, file):
	try:
		os.remove(path+str(filetofilter)+'_filtered.json')
		print('Removing old filtered version')
	except Exception:
		pass

	lines_seen = set()
	outfile = open(path+str(filetofilter)+'_filtered.json','w')
	for line in open(path+str(filetofilter)+'.json','r'):
		data = json.loads(line)
		##### ADD FILTER HERE #####
		try:
			if data['allpeptides'] and data['allpeptides'] != 'Error' and 'raw' in data['filetypes']:# and data['Modification'] != 'No PTMs are included in the dataset':
		###########################
				if line not in lines_seen:
					outfile.write(line)
					lines_seen.add(line)
		except Exception:
			pass
	outfile.close()


if __name__ == '__main__':
	# path = 'Data/metadata/'
	path = '/data/ProteomeToolsRaw/metadata/'
	
	filetofilter = sys.argv[1]

	filter(path = path, file = filetofilter)