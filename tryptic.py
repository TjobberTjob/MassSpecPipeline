datapath = "/data/ProteomeToolsRaw/tryptic/"
import os 
import glob
import re
import subprocess
files = [os.path.dirname(p) for p in glob.glob(datapath+"/*/*")]

for f in files:
	file = f[31:]+".zip"
	if os.path.exists(f+'/allPeptides.txt'):
		continue
	subprocess.run('unzip -j '+f+"/"+file+' allPeptides.txt -d '+f,shell = True)