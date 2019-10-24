import json
import os

path = '/data/ProteomeToolsRaw/Images/'
try:
	os.system('rm '+path+"metadata_cleaned.json")
except:
	print("no cleaned version exist")

lines_seen = set()
outfile = open(path+"metadata_cleaned.json","w")
for line in open(path+"metadata.json","r"):
	if line not in lines_seen:
		outfile.write(line)
		lines_seen.add(line)
outfile.close