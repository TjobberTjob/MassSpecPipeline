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
	df = pd.read_csv(datapath+'allPeptides.txt', sep = "\t")
	name = str(df.iloc[0,0])
	if name == "01625b_GA1-TUM_first_pool_1_01_01-2xIT_2xHCD-1h-R1":
		print('mv '+f+' '+datapath+name+'/file.zip')
	#if not os.path.exists(datapath+name+'/file.zip'):
		#shutil.copyfile(f, datapath+name+'/file.zip')
	os.system('rm '+datapath+'allPeptides.txt')

