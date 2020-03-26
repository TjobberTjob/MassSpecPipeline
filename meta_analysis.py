import json

if __name__ == '__main__':
	
	import sys

	cmd = sys.argv[1]

	if cmd == 'meta':
		for line in open('Data/metadata/accessions.json'):
			print json.loads(line)

