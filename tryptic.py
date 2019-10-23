datapath = "/data/ProteomeToolsRaw/tryptic/"
import os 
import glob
import re
files = [os.path.dirname(p) for p in glob.glob(datapath+"/*/*")]

for f in files:
	file = f[31:]+".zip"
	subprocess.run('unzip -j '+f+file+' txt/allPeptides.txt -d '+f,shell = True)