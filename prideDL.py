if __name__ == '__main__':
	import shutil
	import pandas as pd
	import csv
	import matplotlib.cm as cm
	import matplotlib.pyplot as plt
	import bisect
	import matplotlib as mpl
	from ftplib import FTP
	import re
	from datetime import datetime
	from pathlib import Path
	import glob
	import pickle
	import requests
	import json
	import sys		
	import math
	import numpy as np
	from zipfile import ZipFile
	from bs4 import BeautifulSoup
	from pyteomics import mzml,mzid,mgf
	import subprocess
	import os 
	from os.path import join

def download(file):
	#Start Process
	# file = file.replace(' ','%20') #URL handling
	
	#Check if Raw file exists
	if not (os.path.exists(datapath+filename+'/file.raw') or os.path.exists(datapath+filename+'/mzML.json') or os.path.exists(datapath+filename+'/file.mzML')):
		print('downloading raw file         ', end = '\r')
		os.system('wget -q --show-progress -O '+datapath+filename+'/file.raw'+' -c '+url+'/'+raws+'.raw')


def formatFile():
	#Check whether the docker file is implemented or not
	dockerls = subprocess.check_output('docker image ls',shell = True)
	if not 'thermorawparser' in str(dockerls):
		os.chdir('..')
		os.chdir('ThermoRawFileParser/')
		subprocess.run('docker build --no-cache -t thermorawparser .', shell= True)
		os.chdir('..')
		os.chdir('MassSpecPipeline/')

	if not (os.path.exists(datapath+filename+'/file.mzML') or os.path.exists(datapath+filename+'/mzML.json')):
		print('Formatting file to mzML         ', end = '\r')
		subprocess.run('docker run -v \"'+datapath[:-1]+':/data_input\" -i -t thermorawparser mono bin/x64/Debug/ThermoRawFileParser.exe -i=/data_input/'+filename+'/file.raw -o=/data_input/'+filename+'/ -f=1 -m=1', shell=True)
		os.remove(datapath+filename+'/file-metadata.txt')
		# os.remove(datapath+filename+'/file.raw')
		
	
def process_ms1(spectrum):
	#Scan information
	scan_info = spectrum['scanList']
	#Time
	scan_time = scan_info['scan'][0]['scan start time']
	mz = spectrum['m/z array'] 
	#ion intensity
	intensity = spectrum['intensity array']
	return {'scan_time':scan_time,'intensity':intensity.tolist(),'mz':mz.tolist()}


def internalmzML():
	if not os.path.exists(datapath+filename+'/mzML.json'):
		print('Extracting data from mzML         ', end = '\r')
		data = mzml.MzML(datapath+filename+'/file.mzML')

		#Extracted data
		extracted = {'ms1':{}}
		i= 0 
		#Extract the necessary data from spectra
		for spectrum in data:
		
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


def createImages(resolution,subimage_interval):

	print('Preparing data for image creation', end = '\r')
	mzml = json.load(open(datapath+filename+'/mzML.json'))
	
	mzlistlist = []
	rtlist = []
	for f in mzml['ms1']:
		mzlistlist.append(mzml['ms1'][f]['mz'])
		rtlist.append(mzml['ms1'][f]['scan_time'])
	mzlist = np.unique(sorted([item for sublist in mzlistlist for item in sublist]))
	
	wash_out = 8 #8 minutes of washout of the instrument (proetometools)
	interval = {
		'mz' : {'min':min(mzlist),'max':max(mzlist)},
		'rt' : {'min':wash_out,'max':max(rtlist)}
	}
	
	# Define the intervals for the given resolution
	mz_bin = (float(interval['mz']['max']) - float(interval['mz']['min']))/resolution['x']
	rt_bin = (float(interval['rt']['max']) - float(interval['rt']['min']))/resolution['y']
	
	if not os.path.exists(datapath+filename+'/'+str(resolution['x'])+'x'+str(resolution['y'])+'.txt'):
		#Create an empty array for layer use
		ms1_array = {}
		
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

			# Calculate the y axis. 
			y_n = int((scan_time - interval['rt']['min'])/rt_bin)
			# print (scan_time,y_n)
			l = 0
			
			for mz_elem in mzml['ms1'][scan_id]['mz']:
				if mz_elem < interval['mz']['min'] or mz_elem > interval['mz']['max']:
					continue
				x_n = int((mz_elem - interval['mz']['min'])/mz_bin)
				_key = (x_n,y_n)
				# Current strategy for collapsing the intensity values is taking their logs
				intensity_val = math.log(mzml['ms1'][scan_id]['intensity'][l])
				try:
					ms1_array[_key].append(intensity_val)
				except KeyError:
					ms1_array[_key] = [intensity_val]

				l+=1

		# Create the final image.
		run_i = 0				#For printing purposes
		nonzero_counter = 0		#How many pixels have non-zero values
		total_datapoints = 0	#How many datapoints does the file contain.
		image = []
		for y_i in range(0,resolution['y']):
			run_i+=1
			print("Creating full image: {:2.1%}                  ".format(run_i / resolution['y']), end = '\r') #Print how far we are	
			row = []
			for x_i in range(0,resolution['x']):
				_key = (x_i,y_i)
				try:
					intensity = np.mean(ms1_array[_key]) #Current strategy for normalizing intensity is mean.
					total_datapoints+=(len(ms1_array[_key]))
					nonzero_counter+=1
				except KeyError:
					intensity = 0.0
				row.append(intensity)
			image.append(row)
		print('Saving image files          ', end = '\r')

		imagedata = [image, nonzero_counter, total_datapoints]
		#Save as txt file
		with open(datapath+filename+'/'+str(resolution['x'])+'x'+str(resolution['y'])+'.txt', "wb") as pa:
			pickle.dump(imagedata, pa)

		#Creating the full image
		fullimage = image[::-1]
		fullimage = np.ma.masked_equal(fullimage,0)
		
		#Setup colormap
		colMap = cm.jet
		colMap.set_bad('darkblue')
		
		plt.imshow(fullimage,cmap=colMap,extent = [interval['mz']['min'], interval['mz']['max'], interval['rt']['min'], interval['rt']['max']],aspect = 'auto',vmax = 16,vmin = 6)
		plt.tight_layout()
		plt.xlabel('m/z', fontsize=12)
		plt.ylabel('Retention time - Minutes', fontsize=12)
		plt.axis([interval['mz']['min'], interval['mz']['max'], interval['rt']['min'], interval['rt']['max']])
		plt.tight_layout()

		print('Full image saved         ', end = '\r')
		if not os.path.exists(datapath+filename+'/'+str(resolution['x'])+'x'+str(resolution['y'])+'.png'):
			plt.savefig(datapath+filename+'/'+str(resolution['x'])+'x'+str(resolution['y'])+'.png')		
	
	else: #If the image data exists, just recall it instead of making it
		print('Loading image data         ', end = '\r')
		with open(datapath+filename+'/'+str(resolution['x'])+'x'+str(resolution['y'])+'.txt', "rb") as pa:
			imagedata = pickle.load(pa)
			image 			= imagedata[0]
			nonzero_counter = imagedata[1]
			total_datapoints= imagedata[2]

	#Create the sub-images
	#figuring out all of the mz and rt intervals 
	value = interval['mz']['min']
	mzrangelist = [value]
	for i in range(int(resolution['x'])):
		value+=mz_bin
		mzrangelist.append(value)

	value = interval['rt']['min']
	rtrangelist = [value]
	for i in range(int(resolution['y'])):
		value+=rt_bin
		rtrangelist.append(value)
	
	j=0 #Statistics about outlying images
	
	imgpath = datapath+'Images'
	if not os.path.exists(imgpath):
		os.mkdir(imgpath)
	outfile = open(imgpath+'/metadata.json','a') #The metadata file
	i = 0
	for index, rows in df2.iterrows():
		i+=1
		print("Creating subimages: {:2.1%}               ".format(i / len(df2['Sequence'])), end = '\r') #Print how far we are

		if rows['Retention time']-subimage_interval['rt'] < interval['rt']['min'] or rows['Retention time']+subimage_interval['rt'] > interval['rt']['max'] or rows['m/z']-subimage_interval['mz'] < interval['mz']['min'] or rows['m/z']+subimage_interval['mz']> interval['mz']['max']:
			j+=1 #Check if this image can be created in our range or not
			continue

		if os.path.exists(datapath+'Images/'+filename+'-'+str(i)+'.png'):
			continue

		mzlower = get_lower_bound(mzrangelist,rows['m/z'] - subimage_interval['mz']) 
		mzupper = get_lower_bound(mzrangelist,rows['m/z'] + subimage_interval['mz'])
		rtlower = get_lower_bound(rtrangelist,rows['Retention time'] - subimage_interval['rt'])
		rtupper = get_lower_bound(rtrangelist,rows['Retention time'] + subimage_interval['rt'])
		
		subimage = []
		k = 0
		for lines in image:
			if k < rtupper and k > rtlower:
				subimage.append(lines[mzlower:mzupper])
			k+=1

		subimage = subimage[::-1]
		subimage = np.ma.masked_equal(subimage,0)

		colMap = cm.jet
		colMap.set_bad('darkblue')

		# Save image
		fig = plt.figure()
		fig.set_size_inches((mzupper - mzlower)/200,(rtupper - rtlower)/200)
		ax = plt.Axes(fig, [0., 0., 1, 1])
		ax.set_axis_off()
		fig.add_axes(ax)
		plt.set_cmap('hot')
		ax.imshow(subimage, aspect='equal',cmap = colMap,vmin = 5, vmax = 16)
		plt.savefig(datapath+'Images/'+filename+'-'+str(i)+'.png')
		plt.close(fig)

		new_metadata = {}
		new_metadata.update({"image" : filename+'-'+str(i)})
		for ele in df2.columns[1:]:
			if str(rows[ele]) == 'nan' or str(rows[ele]) == ' ' or ";" in str(rows[ele]):
				continue
			else:
				new_metadata.update({str(ele) : str(rows[ele])})
		outfile.write(json.dumps(new_metadata)+'\n')
	outfile.close()
	
	print('Calculating end statistics:           ', end = '\r')
	
	mzlist_inrange = [i for i in mzlist if i > interval['mz']['min'] and i < interval['mz']['max']]
	rtlist_inrange = [i for i in rtlist if i > interval['rt']['min'] and i < interval['rt']['max']]

	outfile = open(datapath+'end_statistics.json','a')
	end_stats = {}
	end_stats['accession']			= accession
	end_stats['filename']			= filename	
	end_stats['unique mz'] 			= len(mzlist_inrange)
	end_stats['unique rt'] 			= len(rtlist_inrange)
	end_stats['datapoints'] 		= total_datapoints
	end_stats['data per pixel'] 	= total_datapoints / nonzero_counter 
	end_stats['Out of bounds']		= j
	outfile.write(json.dumps(end_stats)+'\n')
	outfile.close()
	print('Done!                                 ')


def get_lower_bound(haystack, needle):

    idx = bisect.bisect(haystack, needle)
    if idx > 0 and idx < len(haystack):
        return idx
    else:
        raise ValueError(f"{needle} is out of bounds of {haystack}")


if __name__ == '__main__':

	accession = sys.argv[1] #Get the accession number 
	pepfile = sys.argv[2]	#Get the name of the maxquant file

	# datapath = '/data/ProteomeToolsRaw/' #Server datapath
	datapath = 'Data/' #Server datapath
	
	url  = 'https://www.ebi.ac.uk/pride/archive/projects/'+accession+'/files'	
	html = requests.get(url).text			  #Webscraping the pride database
	soup = BeautifulSoup(html,'html.parser')									

	for div in soup.find_all('div', {'class': 'grid_6 omega'}):
		url = div.find('a')['href'] #Get the FTP link! 
		break

	os.system('wget -q --show-progress -O '+datapath+'readme.txt '+url+'/README.txt')
	df = pd.read_csv(datapath+'readme.txt',sep='\t')
	os.remove(datapath+'readme.txt')
	searchfiles = df.loc[df['TYPE'] == 'SEARCH',]['URI']

	for zips in searchfiles:			#For loop- for going through all the search files
		zips = zips.replace(' ','%20') 	#URL handling

		os.system('wget -q -O '+datapath+'file.zip'+' -c '+zips) #Get the Zip file
		with ZipFile(datapath+'file.zip','r') as zipped:
			ziplist = zipped.namelist() #Get a list of all contents of it

		for a in ziplist:
			if pepfile in str(a): 		#Fine the peptide file and extract it
				os.system('unzip -o -qq -j '+datapath+'file.zip '+a+' -d '+datapath)
				break
		
		df = pd.read_csv(datapath+pepfile,sep='\t', low_memory=False) #Read in file
		df = df.loc[df['Sequence'] != ' ',] 		#Remove empty sequences
		rawfiles = np.unique(df['Raw file'])		#A list containing all the different raw files this search file has data on

		for raws in rawfiles:
			filename = raws 
	
			if not os.path.exists(datapath+filename):	#Make the file directory if it doesnt exist
				os.system('mkdir '+datapath+filename)	

			if not os.path.exists(datapath+filename+'/file.zip'): #Move or rm zip.file
				os.system('mv '+datapath+'file.zip '+datapath+filename+'/file.zip')
			else:
				os.system('rm '+datapath+'file.zip')
			
			df2 = df.loc[df['Raw file'] == raws,] 		#Take only the part of the data that we need for this raw file
			pd.DataFrame.to_csv(df,datapath+pepfile)	#save this to file for moving

			if not os.path.exists(datapath+filename+'/'+pepfile): #Move or rm txt.file
				os.system('mv '+datapath+pepfile+' '+datapath+filename+'/'+pepfile)
			else:
				os.system('rm '+datapath+pepfile)

			if filename == "01625b_GA1-TUM_first_pool_1_01_01-2xIT_2xHCD-1h-R1": #This cannot be converted to mzml for some reason. So we skip it
				continue

			print('\nfile: '+filename) #Print what file we're working on

			download(raws)
			formatFile()
			internalmzML()

			resolution = {'x':1500,'y':1000}
			subimage_interval  = {'mz':50,'rt':5}
			createImages(resolution,subimage_interval)
		
# python3 prideDL.py PXD004732 allPeptides.txt
# python3 prideDL.py PXD010595 allPeptides.txt