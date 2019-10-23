datapath = "/data/ProteomeToolsRaw/tryptic/"
import os 
import glob
import re
files = [os.path.dirname(p) for p in glob.glob(datapath+"/*/*")]
print(files)
for f in files:
	print(re.match("/", f).end())
