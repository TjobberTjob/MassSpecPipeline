import json
import numpy as np
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
		print(data)
		quit()
		##### ADD FILTER HERE #####
		if str(data['size']) == '[166, 66, 4]':
		# if data['allpeptides'] and 'raw' in data['filetypes']:# and data['Modification'] != 'No PTMs are included in the dataset':
		###########################
			if line not in lines_seen:
				outfile.write(line)
				lines_seen.add(line)
	outfile.close()


if __name__ == '__main__':
	#Read datapath from config file
	with open('config.json') as json_file:
   		data = json.load(json_file)
   	
	datapath = data['path']+'metadata/'
	
	filetofilter = sys.argv[1]

	filter(datapath, filetofilter)