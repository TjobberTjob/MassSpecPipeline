datapath = "/data/ProteomeToolsRaw/tryptic/"
import os 
import glob
import re
import pandas as pd
import subprocess
files = [os.path.dirname(p) for p in glob.glob(datapath+"/*/*")]

for f in files:
	file = f[31:]+".zip"
	if os.path.exists(f+'/allPeptides.txt'):
		continue
	subprocess.run('unzip -j '+f+"/"+file+' allPeptides.txt -d '+f,shell = True)

for f in files:
	df = pd.read_csv(f+'/allPeptides.txt', sep = '\t')
	df.head()
	#df2 = df.loc[df['Sequence'] != ' ',]
	#pd.DataFrame.to_csv(df2,f+'/allPeptides.txt')