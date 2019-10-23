datapath = "/data/ProteomeToolsRaw/tryptic/"
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

for f in files:
	file = f[31:]+".zip"
	if os.path.exists(f+'/allPeptides.txt'):
		continue
	subprocess.run('unzip -j '+f+"/"+file+' allPeptides.txt -d '+f,shell = True)

files = np.unique(files)

for f in files:
	file = f[31:]+".zip"
	print(f+'/allPeptides.txt')
	df = pd.read_csv(f+'/allPeptides.txt', sep = '\t')
	print('mv '+f+file+' '+'/data/ProteomeToolsRaw/'+df.iloc[0,0]+'/file.zip')