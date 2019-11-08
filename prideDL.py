if __name__ == '__main__':
	import shutil
	import pandas as pd
	import csv
	import matplotlib.cm as cm
	import matplotlib.pyplot as plt
	import matplotlib as mpl
	from ftplib import FTP
	import re
	from datetime import datetime
	from pathlib import Path
	import glob
	import requests
	import json
	import sys		
	import math
	import numpy as np
	from zipfile import ZipFile
	from bs4 import BeautifulSoup
	import subprocess
	import os 
	from os.path import join

def download():
	#Start Process
	start = datetime.now()
	raws = raws.replace(' ','%20')
	#Check if Raw file exists
	if os.path.exists(datapath+filename+'/file.raw') or os.path.exists(datapath+filename+'/mzML.json') or os.path.exists(datapath+filename+'/file.mzML'):
		print('raw or parsed file exists')
	else:
		print('downloading raw file')
		os.system('wget -q --show-progress -O '+datapath+filename+'/file.raw'+' -c '+url[:-10]+raws+'.raw')
		end1 = datetime.now()
		diff1 = end1 - start
		print('Rawfile downloaded \ntime: '+str(diff1))
		
def formatFile():
	#Check whether the docker file is implemented or not
	output = subprocess.check_output('docker image ls',shell = True)
	if 'thermorawparser' not in str(output):
		os.chdir('..')
		os.chdir('ThermoRawFileParser/')
		subprocess.run('docker build --no-cache -t thermorawparser .', shell= True)
		os.chdir('..')
		os.chdir('MassSpecPipeline/')

	if os.path.exists(datapath+filename+'/file.mzML') or os.path.exists(datapath+filename+'/mzML.json'):
		print('mzML file exists')
	else:
		print('Formatting file to mzML')
		subprocess.run('docker run -v \"'+datapath[:-1]+':/data_input\" -i -t thermorawparser mono bin/x64/Debug/ThermoRawFileParser.exe -i=/data_input/'+filename+'/file.raw -o=/data_input/'+filename+'/ -f=1 -m=1', shell=True)
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
	try:
		#mass to charge (m/z)
		mz = spectrum['m/z array']
	except:
		import pprint
		pprint.pprint(spectrum)
		quit()
	#ion intensity
	intensity = spectrum['intensity array']
	return {'scan_time':scan_time,'intensity':intensity.tolist(),'mz':mz.tolist()}

def internalmzML():
	if os.path.exists(datapath+filename+'/mzML.json'):
		print('mzML data already extracted')
	else:
		print('Extracting data from mzML')
		data = mzml.MzML(datapath+filename+'/file.mzML')

		#Extracted data
		extracted = {'ms1':{}}
		i= 0 
		#Extract the necessary data from spectra
		for spectrum in data:
			
			#i+=1
			#if i%1000 == 0 :
			#	print(i)
		
			if spectrum['ms level'] != 1:
				continue
			#Scan id
			scan_id = int(spectrum['id'].split('scan=')[1])
			
			#Deal with ms level 1 spectra
			ms1_spectrum = process_ms1(spectrum)
			extracted['ms1'][scan_id] = {'mz':ms1_spectrum['mz'],'intensity':ms1_spectrum['intensity'],'scan_time':ms1_spectrum['scan_time']}

		f = open(datapath+filename+'/mzML.json','w')
		f.write(json.dumps(extracted))
		f.close()
		# os.remove(datapath+filename+'/file.mzml')

def full_image(interval,resolution,show=False):
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

def sub_images(resolution):
	print('Creating subimages')

	df = pd.read_csv(datapath+filename+'/'+pepfile)
	mzml = json.load(open(datapath+filename+'/mzML.json'))
	
	maxint = 0 #GLOBAL MAXIMA
	for f in mzml['ms1']:
		maxint = max(maxint,max(mzml['ms1'][str(f)]['intensity']))
	
	mz_interval = 50 #INTERVALS
	rt_interval = 5  #INTERVALS
	j=0
	outpath = datapath+'Images'
	if not os.path.exists(outpath):
		os.mkdir(outpath)

	outfile = open(join(outpath,'metadata.json'),'a')
	for i in range(len(df['Sequence'])):
		if os.path.exists(datapath+'Images/'+filename+'-'+str(i+1)+'.png'):
			#print(str(i+1)+' of '+str(len(df['Sequence']))+' was already created')
			continue
  
		if df['Retention time'][i]-rt_interval < min(df['Retention time'])+wash_out or df['Retention time'][i]+rt_interval > max(df['Retention time']) or df['m/z'][i]-mz_interval < min(df['m/z']) or df['m/z'][i]+mz_interval > max(df['m/z']):
			j+=1
			#print(str(i+1)+' of '+str(len(df['Sequence']))+' was out of bounds')
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
		fig.set_size_inches((2,2))
		ax = plt.Axes(fig, [0., 0., 1., 1.])
		ax.set_axis_off()
		fig.add_axes(ax)
		plt.set_cmap('hot')
		ax.imshow(image, aspect='equal',cmap = colMap)#,vmin = 5, vmax = 16)
		plt.savefig(datapath+'Images/'+filename+'-'+str(i+1)+'.png')
		plt.close(fig)
		#print('{0}\r'.format(str(i+1)+' of '+str(len(df['Sequence']))),) #PRINT CREATED IMAGE
		print("Progress {:2.1%}".format(i / len(df['Sequence'])), end="\r")

		new_metadata = {}
		new_metadata.update({"image" : filename+'-'+str(i+1)})
		for ele in df.columns[1:]:
			if str(df[ele][i]) == 'nan' or str(df[ele][i]) == ' ' or ";" in str(df[ele][i]):
				continue
			else:
				new_metadata.update({str(ele) : str(df[ele][i])})

		outfile.write(json.dumps(new_metadata)+'\n')

	outfile.close()
	print(str(j)+' Images were out of bounds')

def validated_input(prompt, valid_values):
    valid_input = False
    while not valid_input:
        value = input(prompt + ' ' + '/'.join(valid_values)+"\n")
        valid_input = value.lower() in valid_values
    return value

if __name__ == '__main__':

	accession = sys.argv[1]
	pepfile = input("What's the name of the MaxQuant output file?\n")

	datapath = '/data/ProteomeToolsRaw/' #Server datapath
	url  = 'https://www.ebi.ac.uk/pride/archive/projects/'+accession+'/files'
	html = requests.get(url).text
	soup = BeautifulSoup(html,'html.parser')
	for div in soup.find_all('div', {'class': 'grid_6 omega'}):
		url = div.find('a')['href']
		break

	os.system('wget -q --show-progress -O '+datapath+'readme.txt '+url+'/README.txt')
	df = pd.read_csv(datapath+'readme.txt',sep='\t')
	os.remove(datapath+'readme.txt')
	searchfiles = df.loc[df['TYPE'] == 'SEARCH',]['URI']

	for zips in searchfiles:
		zips = zips.replace(' ','%20')

		os.system('wget -q --show-progress -O '+datapath+'file.zip'+' -c '+zips)
		with ZipFile(datapath+'file.zip','r') as zipped:
			ziplist = zipped.namelist()
		for a in ziplist:
			if pepfile in str(a):
				subprocess.run('unzip -j '+datapath+'file.zip '+a+' -d '+datapath+'/',shell = True)
				break

		df2 = pd.read_csv(datapath+pepfile,sep='\t')
		df2 = df2.loc[df2['Sequence'] != ' ',]
		rawfiles = np.unique(df2['Raw file'])
		print(rawfiles)

		for raws in rawfiles:
			filename = raws
			
			if not os.path.exists(datapath+filename+'/file.zip'): #Move or rm zip.file
				os.system('mv '+datapath+'file.zip '+datapath+filename+'/file.zip')
			else:
				os.system('rm '+datapath+'file.zip')
			
			df3 = df2.loc[df2['Raw file'] == raws,]
			pd.DataFrame.to_csv(df3,datapath+pepfile)
			if not os.path.exists(datapath+filename+'/'+pepfile): #Move or rm txt.file
				os.system('mv '+datapath+pepfile+' '+datapath+filename+'/'+pepfile)
			else:
				os.system('rm '+datapath+pepfile)
			if filename == "01625b_GA1-TUM_first_pool_1_01_01-2xIT_2xHCD-1h-R1":
				continue
			print('\nfile: '+filename)

			download()
			formatFile()
			internalmzML()

			wash_out = 8
			interval = {
					'mz' : {'min':360,'max':1250},
					'rt' : {'min':wash_out,'max':60}
				}
			resolution = {'x':500,'y':300}
			# full_image(interval,resolution,show=False)

			resolution = {'x':100,'y':100}
			sub_images(resolution)
		
# python3 prideDL.py PXD004732
# python3 prideDL.py PXD010595