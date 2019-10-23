datapath = "/data/ProteomeToolsRaw/tryptic"
if __name__ == '__main__':
	import shutil
	import pandas as pd
	import urllib3
	import csv
	import matplotlib.cm as cm
	import matplotlib.pyplot as plt
	import matplotlib as mpl
	from ftplib import FTP
	from datetime import datetime
	from pyteomics import mzml,mzid,mgf
	from pathlib import Path
	import glob
	import json
	import math
	import numpy as np
	import subprocess
	import os 
	from os.path import join

files = [os.path.dirname(p) for p in glob.glob(datapath+"/*/*")]
files = np.unique(files)
for f in files:
	file = f[31:]+".zip"
	if os.path.exists(f+'/*.zip'):
		print(file)



#files = [os.path.dirname(p) for p in glob.glob(datapath+"/*/*")]
#files = np.unique(files)
#for f in files:
#	file = f[31:]+".raw"
#	if os.path.exists(f+file):
#		os.system(f+file+" "+2)
#	subprocess.run('unzip -j '+f+"/"+file+' allPeptides.txt -d '+f,shell = True)
#
#files = np.unique(files)
#
#for f in files:
#	file = f[31:]+".zip"
#	if f+file == '/data/ProteomeToolsRaw/tryptic/TUM_first_pool_97_01_01_2xIT_2xHCD-1h-R2-tryptic':
#		ss = ','
#	else:
#		ss = '\t'
#	df = pd.read_csv(f+'/allPeptides.txt', sep = ss)
#	print('/data/ProteomeToolsRaw/01650b_BG7-TUM_first_pool_97_01_01-2xIT_2xHCD-1h-R2/file.zip')
#	try:
#		shutil.move(f+'/'+file, '/data/ProteomeToolsRaw/'+df.iloc[0,0]+'/file.zip')
#	except:
#		print("moved")
