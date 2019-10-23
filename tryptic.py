datapath = "/data/ProteomeToolsRaw/tryptic/"
import os 
import glob
import re
files = [os.path.dirname(p) for p in glob.glob(datapath+"/*/*")]

for f in files:
	print(f[31:])
