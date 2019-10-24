datapath = "/data/ProteomeToolsRaw/"
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

files = glob.glob(datapath+"*.zip")
files = np.unique(files)
for f in files:
	subprocess.run('unzip -j '+f+' allPeptides.txt -d '+datapath,shell = True)
	
	#with ZipFile(f, 'r') as zipObj:
	#	zipObj.extract(allPeptides.txt, 'allPeptides.txt')
	df = pd.read_csv(datapath+'allPeptides.txt', sep = "\t")
	#df2 = pd.read_csv(datapath+'allPeptides.txt', sep = ",")
	print(str(len(df)) + df.iloc[0,0])

	#if not os.path.exists(datapath+df.iloc[0,0]+'/file.zip'):
	#	print(os.path.exists(datapath+df.iloc[0,0]+'/file.zip')+ df.iloc[0,0])
	#	#shutil.move(f, datapath+df.iloc[0,0]+'/file.zip')
	os.system('rm '+datapath+'allPeptides.txt')

