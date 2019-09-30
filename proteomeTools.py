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

def download(url):
	#Start Process
	Path(datapath+filename+'/'+year).touch()
	
	#Download the files from the urls and inset into right folder
	http = urllib3.PoolManager()
	urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
	start = datetime.now()
	print(datapath+filename+'/'+filename+'.raw')
	print('mv '+datapath+filename+'/'+filename+'.raw '+datapath+filename+'/file.raw')
	quit()
	if os.path.exists(datapath+filename+'/'+filename+'.raw')
		print('mv 'datapath+filename+'/'+filename+'.raw 'datapath+filename+'/file.raw')
		quit()

	#Check if Raw file exists
	if os.path.exists(datapath+filename+'/file.raw'):
		print('Parsed files already exists')
	else:
		print('Starting rawfile download of '+filename)
		with http.request('GET', url+ '.raw', preload_content=False) as r, open(datapath+filename+'/file.raw', 'wb') as out_file:
			shutil.copyfileobj(r, out_file)
		end1 = datetime.now()
		diff1 = end1 - start
		print('Rawfile downloaded \ntime: '+str(diff1))
	#Check if txt file exists
	if os.path.exists(datapath+filename+'/allPeptides.txt'):
		print('Txt file already exists')
	else:
		print('Starting zipfile download of '+filename)
		with http.request('GET', url+'.zip', preload_content=False) as r, open(datapath+filename+'/file.zip', 'wb') as out_file:		
			shutil.copyfileobj(r, out_file)
		subprocess.run('unzip -j Data/'+filename+'/file.zip txt/allPeptides.txt -d Data/'+filename+'/',shell = True)
		os.remove(datapath+filename+'/file.zip')
		#Fix the allpep.txt
		df = pd.read_csv(datapath+filename+'/allPeptides.txt', sep = '\t')
		df2 = df.loc[df['Sequence'] != ' ',]
		pd.DataFrame.to_csv(df2,datapath+filename+'/allPeptides.txt')
		#finish the DL
		end2 = datetime.now()
		diff2 = end2 - start
		print('Files downloaded and handled \ntime: '+str(diff2))
		
def formatFile(filename):
	#Check whether the docker file is implemented or not
	output = subprocess.check_output('docker image ls',shell = True)
	if 'thermorawparser' not in str(output):
		os.chdir('ThermoRawFileParser/')
		subprocess.run('docker build --no-cache -t thermorawparser .', shell= True)
		os.chdir('..')

	if os.path.exists(datapath+filename+'/file.mzML') or os.path.exists(datapath+filename+'/mzML.json'):
		print('Already formatted to mzML')
	else:
		print('Formatting '+filename+' to mzML')
		subprocess.run('docker run -v \"'+os.getcwd()+'/Data:/data_input\" -i -t thermorawparser  mono /home/biodocker/bin/bin/x64/Debug/ThermoRawFileParser.exe -i=/data_input/'+filename+'/file.raw -o=/data_input/'+filename+'/ -f=1 -m=1', shell=True)
		os.remove(datapath+filename+'/file-metadata.txt')
		# os.remove(datapath+filename+'/file.raw')
		
def process_ms1(spectrum):
	#Scan information
	scan_info = spectrum['scanList']
	if scan_info['count'] !=1:
		print('Scan with more than one value - not designed for this, please review data')
		quit()
			
	#Time
	scan_time = scan_info['scan'][0]['scan start time']
	#mass to charge (m/z)
	mz = spectrum['m/z array']
	
	#ion intensity
	intensity = spectrum['intensity array']
	return {'scan_time':scan_time,'intensity':intensity.tolist(),'mz':mz.tolist()}

def internalmzML(filename):
	if os.path.exists(datapath+filename+'/mzML.json'):
		print('mzML data already extracted')
	else:
		print('Extracting data from mzML file')
		data = mzml.MzML(datapath+filename+'/file.mzml')

		stats = {'ms-levels':{1:0}}

		#Extracted data
		extracted = {'ms1':{}}
		i= 0 
		#Extract the necessary data from spectra
		for spectrum in data:
			
			i+=1
			if i%1000 == 0 :
				print(i)
		
			#Scan id
			scan_id = int(spectrum['id'].split('scan=')[1])
			
			#Deal with ms level 1 spectra
			ms1_spectrum = process_ms1(spectrum)
			extracted['ms1'][scan_id] = {'mz':ms1_spectrum['mz'],'intensity':ms1_spectrum['intensity'],'scan_time':ms1_spectrum['scan_time']}

		f = open(datapath+filename+'/mzML.json','w')
		f.write(json.dumps(extracted))
		f.close()
		# os.remove(datapath+filename+'/file.mzml')

def full_image(interval,resolution,filename,show=False):
	if os.path.exists(datapath+filename+'/'+str(resolution['x'])+'x'+str(resolution['y'])+'.png'):
		print('Full image already exists')

	else:

		print('Creating full image')		
		mzml = json.load(open(datapath+filename+'/mzML.json'))
		# #Define the intervals for the given resolution
		x_d = (float(interval['mz']['max']) - float(interval['mz']['min']))/resolution['x']
		y_d = (float(interval['rt']['max']) - float(interval['rt']['min']))/resolution['y']
		#Create the initial array.
		#elements are given by (x,y)
		ms1_array = {}
		#print('Intervals',x_d,y_d)
		#Make sure the intervals can define actual bins
		if x_d == 0 or y_d == 0:
			print ('Nil interval spacing- please review resolution and interval of mz/rt')
			quit()
			
		#Collect inverval statistics.
		stats = { 
			'x' : {},#x-axis, m/z
			'y' : [],#y-axis rt
			'clashed_cells':0,#Number of pixels in which there are more than one value
		}
	
		#Get sorted list of scan ids.
		scan_ids = []
		for scan_id in mzml['ms1']:
			scan_ids.append(int(scan_id))
		intensitys = []
		for intensity in mzml['ms1']:
			intensity = intensity['intensity']
			intensitys.append(int(intensity))

		for scan_id in sorted(scan_ids):
			scan_id = str(scan_id)
			#Get the intervals
			scan_time = float(mzml['ms1'][scan_id]['scan_time'])
			if scan_time < interval['rt']['min'] or scan_time > interval['rt']['max']:
				continue
			stats['y'].append(scan_time)
		#Calculate the y axis.	
			y_n = int((scan_time - interval['rt']['min'])/y_d)
			#print (scan_time,y_n)
			i = 0
			for mz_elem in mzml['ms1'][scan_id]['mz']:
				if mz_elem < interval['mz']['min'] or mz_elem > interval['mz']['max']:
					continue
				stats['x'][mz_elem] = 0
				x_n = int((mz_elem - interval['mz']['min'])/x_d)
				_key = (x_n,y_n)
				#Current strategy for collapsing the intensity values is taking their logs
				intensity_val = math.log(mzml['ms1'][scan_id]['intensity'][i])
				try:
					ms1_array[_key].append(intensity_val)
				except KeyError:
					ms1_array[_key] = [intensity_val]
				i+=1
	
		#Create the final image.
		image = []
		for y_i in range(0,resolution['y']):
			row = []
			for x_i in range(0,resolution['x']):
				_key = (x_i,y_i)
				try:
					#Current strategy for normalizing intensity is mean.
					intensity = np.mean(ms1_array[_key])
					if len(ms1_array[_key])>1:
						stats['clashed_cells']+=1
				except KeyError:
					intensity = 0.0
				row.append(intensity)
			image.append(row)
		image = image[::-1]
		image = np.ma.masked_less(image,0.01)
		#Setup colormap
		colMap = cm.jet
		colMap.set_bad('darkblue')
		#Save or show image
		plt.imshow(image,cmap=colMap,extent = [interval['mz']['min'], interval['mz']['max'], interval['rt']['min'], interval['rt']['max']],aspect = 'auto',vmax = 16,vmin = 6)
		plt.tight_layout()
		plt.xlabel('m/z', fontsize=12)
		plt.ylabel('Retention time - Minutes', fontsize=12)
		plt.axis([interval['mz']['min'], interval['mz']['max'], interval['rt']['min'], interval['rt']['max']])
		plt.colorbar(extend = 'both')
		plt.tight_layout()
		print('Image created')
		if show == True:
			plt.show()
		elif show == False:
			plt.savefig(datapath+filename+'/'+str(resolution['x'])+'x'+str(resolution['y'])+'.png')		

def sub_images(wash_out,resolution,filename):
	print('Creating subimages')

	df = pd.read_csv(datapath+filename+'/allPeptides.txt')
	mzml = json.load(open(datapath+filename+'/mzML.json'))
	
	maxint = 0 #GLOBAL MAXIMA
	for f in mzml['ms1']:
		maxint = max(maxint,max(mzml['ms1'][str(f)]['intensity']))
	
	mz_interval = 50
	rt_interval = 5
	j=0
	outpath = "Data/Images"
	outfile = open(join(outpath,'metadata.json'),'a')
	for i in range(len(df['Sequence'])):
		if os.path.exists(datapath+'Images/'+filename+'-'+str(i+1)+'.png'):
			print(str(i+1)+' of '+str(len(df['Sequence']))+' was already created')
			continue
  
		if df['Retention time'][i]-rt_interval < min(df['Retention time'])+wash_out or df['Retention time'][i]+rt_interval > max(df['Retention time']) or df['m/z'][i]-mz_interval < min(df['m/z']) or df['m/z'][i]+mz_interval > max(df['m/z']):
			j+=1
			print(str(i+1)+' of '+str(len(df['Sequence']))+' was out of bounds')
			continue
  
		interval = {
				'mz' : {'min':df['m/z'][i]-mz_interval,'max':df['m/z'][i]+mz_interval},
				'rt' : {'min':df['Retention time'][i]-rt_interval,'max':df['Retention time'][i]+rt_interval}
			}
		# print(interval)
			  
		# Define the intervals for the given resolution
		x_d = (float(interval['mz']['max']) - float(interval['mz']['min']))/resolution['x']
		y_d = (float(interval['rt']['max']) - float(interval['rt']['min']))/resolution['y']
		# Create the initial array.
		# elements are given by (x,y)
		ms1_array = {}
		# print('Intervals',x_d,y_d)
		# Collect inverval statistics.
		stats = { 
			'x' : {},#x-axis, m/z
			'y' : [],#y-axis rt
			'clashed_cells':0,#Number of pixels in which there are more than one value
		}
		  
		# Get sorted list of scan ids.
		scan_ids = []
		for scan_id in mzml['ms1']:
			scan_ids.append(int(scan_id))
		intensitys = []
		for intensity in mzml['ms1']:
			intensitys.append(int(intensity))
		 
 
		for scan_id in sorted(scan_ids):
			scan_id = str(scan_id)
			# Get the intervals
			scan_time = float(mzml['ms1'][scan_id]['scan_time'])
			if scan_time < interval['rt']['min'] or scan_time > interval['rt']['max']:
				continue
			stats['y'].append(scan_time)
		# Calculate the y axis. 
			y_n = int((scan_time - interval['rt']['min'])/y_d)
			# print (scan_time,y_n)
			l = 0
			for mz_elem in mzml['ms1'][scan_id]['mz']:
				if mz_elem < interval['mz']['min'] or mz_elem > interval['mz']['max']:
					continue
				stats['x'][mz_elem] = 0
				x_n = int((mz_elem - interval['mz']['min'])/x_d)
				_key = (x_n,y_n)
				# Current strategy for collapsing the intensity values is taking their logs
				intensity_val = math.log(mzml['ms1'][scan_id]['intensity'][l] / maxint)
				try:
					ms1_array[_key].append(intensity_val)
				except KeyError:
					ms1_array[_key] = [intensity_val]
				l+=1
 
		# Create the final image.
		image = []
		for y_i in range(0,resolution['y']):
			row = []
			for x_i in range(0,resolution['x']):
				_key = (x_i,y_i)
				try:
					# Current strategy for normalizing intensity is mean.
					intensity = np.mean(ms1_array[_key])
					if len(ms1_array[_key])>1:
						stats['clashed_cells']+=1
				except KeyError:
					intensity = 0.0
				row.append(intensity)
			image.append(row)
		image = image[::-1]
		image = np.ma.masked_equal(image,0)
		colMap = cm.jet
		colMap.set_bad('black')
		# Save image
		fig = plt.figure()
		fig.set_size_inches((3,3))
		ax = plt.Axes(fig, [0., 0., 1., 1.])
		ax.set_axis_off()
		fig.add_axes(ax)
		plt.set_cmap('hot')
		ax.imshow(image, aspect='equal',cmap = colMap)#,vmin = 5, vmax = 16)
		plt.savefig(datapath+'Images/'+filename+'-'+str(i+1)+'.png')
		print(str(i+1)+' of '+str(len(df['Sequence']))+' subimages created')
 
		new_metadata = {}
		new_metadata.update({"image" : filename+'-'+str(i+1)})
		for ele in df.columns[1:]:
			if str(df[ele][i]) == 'nan' or str(df[ele][i]) == ' ' or ";" in str(df[ele][i]):
				continue
			else:
				new_metadata.update({str(ele) : str(df[ele][i])})
		# full = {filename+'-'+str(i+1) : new_metadata}
		outfile.write(json.dumps(new_metadata)+'\n')
	outfile.close()
	print(str(j)+' Images were out of bounds')

if __name__ == '__main__':

	import sys
	inputs = sys.argv[1]
	datapath = '/data/ProteomeToolsRaw/' #Server datapath
	if not os.path.exists(datapath+"Images"):
		os.mkdir(datapath+'Images')


	if inputs == 'all':
		URL = sys.argv[2]+"/"
		ftp = FTP('ftp.pride.ebi.ac.uk')
		ftp.login('anonymous')
		ftp.retrbinary('RETR ' + URL + 'README.txt' ,open(datapath+'intermediary.txt', 'wb').write)
		ftp.close()
		print('Downloading intermediary file(s)')
		df = pd.read_csv(datapath+'ntermediary.txt',sep='\t')
		os.remove(datapath+'intermediary.txt')
		urls = []
		for f in df.loc[df['TYPE'] == 'RAW',]['URI']:
			urls.append('https://www' + f[15:-4])
	elif inputs == 'all-PT':
		URL2019 = '/pride/data/archive/2019/05/PXD010595/'
		URL2017 = '/pride/data/archive/2017/02/PXD004732/'
		ftp = FTP('ftp.pride.ebi.ac.uk')
		ftp.login('anonymous')
		ftp.retrbinary('RETR ' + URL2019 + 'README.txt' ,open(datapath+'intermediary2019.txt', 'wb').write)
		ftp.retrbinary('RETR ' + URL2017 +  'README.txt' ,open(datapath+'intermediary2017.txt', 'wb').write)
		ftp.close()
		#extract the urls we need and remove the intermediary files
		print('Downloading intermediary files')
		df1 = pd.read_csv(datapath+'intermediary2019.txt',sep='\t')
		df2 = pd.read_csv(datapath+'intermediary2017.txt',sep='\t') 
		os.remove(datapath+'intermediary2019.txt')
		os.remove(datapath+'intermediary2017.txt')
		urls = []
		for f in df1.loc[df1['TYPE'] == 'RAW',]['URI']:
			urls.append('https://www' + f[15:-4])
		for f in df2.loc[df2['TYPE'] == 'RAW',]['URI']:
			urls.append('https://www' + f[15:-4])
	else: urls = [inputs]
	
	for f in urls:
		year = f[41:45]
		filename = f[59:]
		if not os.path.exists(datapath+filename):
			os.mkdir(datapath+filename)

		download(url = f)
		formatFile(filename = filename)
		internalmzML(filename = filename)

		wash_out = 8
		interval = {
				'mz' : {'min':360,'max':1250},
				'rt' : {'min':wash_out,'max':60} #~8 min washout
			}
		resolution = {'x':500,'y':300}
		
		full_image(interval,resolution,filename = filename,show=False)
		
		# shutil.rmtree(datapath+'Images')	#For testing metadata
		# os.mkdir(datapath+'Images')		#For testing metadata

		resolution = {'x':100,'y':100}
		sub_images(wash_out,resolution,filename = filename)

# python proteomeTools.py all-PT
# python proteomeTools.py https://www.ebi.ac.uk/pride/data/archive/2019/05/PXD010595/01974c_BC1-TUM_missing_first_3_01_01-ETD-1h-R4
# python proteomeTools.py https://www.ebi.ac.uk/pride/data/archive/2019/05/PXD010595/01974c_BA1-TUM_missing_first_1_01_01-ETD-1h-R4
# python proteomeTools.py https://www.ebi.ac.uk/pride/data/archive/2019/05/PXD010595/02208a_GE7-TUM_second_addon_55_01_01-ETD-1h-R1
# python proteomeTools.py all /pride/data/archive/2019/05/PXD010595

