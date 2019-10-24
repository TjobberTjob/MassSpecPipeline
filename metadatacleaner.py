import json

path = sys.argv[1]

lines_seen = set()
outfile = open(path+"metadata_cleaned.json","w")
for line in open(path+"metadata.json","r"):
	if line not in lines_seen:
		outfile.write(line)
		lines_seen.add(line)
outfile.close