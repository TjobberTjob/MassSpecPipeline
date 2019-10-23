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

files = [os.path.dirname(p) for p in glob.glob(datapath+"/*/*")]
files = np.unique(files)
for f in files:
	file = f[23:]+".raw"
	print(f+"/"+file)
	if os.path.exists(f+"/"+file):
		print(f+"/"+file+" "+f+"file.raw")
	else:
		continue
