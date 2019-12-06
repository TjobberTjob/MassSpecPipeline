if __name__ == '__main__':
	import shutil
	import pandas as pd
	import csv
	import matplotlib.cm as cm
	import matplotlib.pyplot as plt
	import bisect
	import matplotlib as mpl
	import pickle
	import requests
	import json
	import sys		
	import math
	import numpy as np
	import subprocess
	import os
	from zipfile import ZipFile
	from bs4 import BeautifulSoup
	from pyteomics import mzml,mzid,mgf
 

def zipfile_finder(accession, path, metapath):
	url = 'http://ftp.pride.ebi.ac.uk/pride/data/archive/'+accession

	#Download readme file
	os.system('wget -q -O '+path+accession[-9:]+'-readme.txt '+url+'/README.txt')

	#Handle and remove readme file
	df = pd.read_csv(path+accession[-9:]+'-readme.txt',sep='\t')
	os.remove(path+accession[-9:]+'-readme.txt')

	searchfiles = df.loc[df['TYPE'] == 'SEARCH',]['URI']
	return searchfiles, url


def rawfile_finder(zipfile, path, maxquant_file):
	#Handle spaces in urls
	zipfile = zipfile.replace(' ','%20')
	zipfilename = zipfile[63:]

	#Download zip file
	os.system('wget -q --show-progress -O '+path+zipfilename+' '+zipfile)

	#Get a list of files with directories from zip file
	with ZipFile(path+zipfilename,'r') as zipped:
		ziplist = zipped.namelist()

	#Extract the peptide file from the zipfile
	for a in ziplist:
		if maxquant_file in a:
			with ZipFile(path+zipfilename) as z:
				with z.open(a) as zf, open(path+zipfilename+'-'+maxquant_file, 'wb') as zfg:
					shutil.copyfileobj(zf, zfg)
					break
		else:
			continue

	#Go through the maxquant output file and get all the raw files
	df = pd.read_csv(path+zipfilename+'-'+maxquant_file,sep='\t', low_memory=False)
	df = df.loc[df['Sequence'] != ' ',] #Remove empty sequences 	
	rawfiles = np.unique(df['Raw file'])
	
	return rawfiles, df, zipfilename


def filehandling(filename, zipfilename, path, maxquant_file, df, url):
	filepath = path+filename+'/'
	#Make the file directory if it doesnt exist
	if not os.path.exists(filepath):	
		os.mkdir(filepath)

	#Move or rm zip.file
	if not os.path.exists(filepath+zipfilename): 
		shutil.copyfile(path+zipfilename, filepath+'file.zip')

	#Move or rm txt.file
	if not os.path.exists(filepath+zipfilename+'-'+maxquant_file): 
		df2 = df.loc[df['Raw file'] == filename,]
		pd.DataFrame.to_csv(df,filepath+maxquant_file)

	#Download the raw file
	if not (os.path.exists(filepath+'/file.raw') or os.path.exists(filepath+'/file.mzML') or os.path.exists(filepath+'/mzML.json')):
		os.system('wget -q --show-progress -O '+filepath+'/file.raw -c '+url+'/'+filename+'.raw')
	return df2, filepath


def formatFile(filename, path, filepath):
	print('Formatting file to mzML               ', end = '\r')
	#Check whether the docker file is implemented or not
	if not (os.path.exists(filepath+'file.mzML') or os.path.exists(filepath+'mzML.json')):
		dockerls = subprocess.check_output('docker image ls',shell = True)
		if not 'thermorawparser' in str(dockerls):
			os.chdir('..')
			try:
				os.system('git clone https://github.com/compomics/ThermoRawFileParser.git')
			except Exception:
				pass
			os.chdir('ThermoRawFileParser/')
			os.system('docker build --no-cache -t thermorawparser .')
			os.chdir('../MassSpecPipeline/')

		print('docker run -v \"'+os.getcwd()+path[:-1]+':/data_input\" -i -t thermorawparser mono bin/x64/Debug/ThermoRawFileParser.exe -i=/data_input/'+filename+'/file.raw -o=/data_input/'+filename+'/ -f=1 -m=1')#, shell=True)		
		quit()
		os.remove(filepath+'file-metadata.txt')
		os.remove(path+filename+'/file.raw')


def process_ms1(spectrum):
	#Scan information
	scan_info = spectrum['scanList']
	#Time
	scan_time = scan_info['scan'][0]['scan start time']
	mz = spectrum['m/z array'] 
	#ion intensity
	intensity = spectrum['intensity array']
	return {'scan_time':scan_time,'intensity':intensity.tolist(),'mz':mz.tolist()}


def internalmzML(path):
	#Extract the data from the mzml, if we havnt already
	if not os.path.exists(path+'mzML.json'):
		print('Extracting data from mzML                     ', end = '\r')
		data = mzml.MzML(path+'file.mzML')

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

		f = open(path+'mzML.json','w')
		f.write(json.dumps(extracted))
		f.close()
		# os.remove(path+'file.mzml')


def get_lower_bound(haystack, needle):
	idx = bisect.bisect(haystack, needle)
	if idx > 0 and idx < len(haystack):
		return idx
	else:
		raise ValueError(f"{needle} is out of bounds of {haystack}")


def createImages(filename, path, filepath, metapath, resolution, subimage_interval, df, savepng):

	print('Preparing data for image creation              ', end = '\r')
	mzml = json.load(open(filepath+'/mzML.json'))
	
	mzlistlist = []
	rtlist  = []
	intlist = []
	for f in mzml['ms1']:
		mzlistlist.append(mzml['ms1'][f]['mz'])
		rtlist.append(mzml['ms1'][f]['scan_time'])
		intlist.append(mzml['ms1'][f]['intensity'])
	mzlist = np.unique(sorted([item for sublist in mzlistlist for item in sublist]))
	intlist = [item for sublist in intlist for item in sublist]
	lowbound = math.log(np.percentile(intlist,0.5))
	highbound = math.log(np.percentile(intlist,99.5))

	wash_out = 8 #8 minutes of washout of the instrument (proetometools)
	interval = {
		'mz' : {'min':min(mzlist),'max':max(mzlist)},
		'rt' : {'min':wash_out,'max':max(rtlist)}
	}
	
	# Define the intervals for the given resolution
	mz_bin = (float(interval['mz']['max']) - float(interval['mz']['min']))/resolution['x']
	rt_bin = (float(interval['rt']['max']) - float(interval['rt']['min']))/resolution['y']
	
	if not os.path.exists(filepath+'/'+str(resolution['x'])+'x'+str(resolution['y'])+'.txt'):
		#Create an empty array for layer use
		ms1_array = {}
		
		# Get sorted list of scan ids.
		scan_ids = [int(scan_id) for scan_id in mzml['ms1']]

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
			print("Creating full image: {:2.1%}                                     ".format(run_i / resolution['y']), end = '\r') #Print how far we are	
			row = []
			for x_i in range(0,resolution['x']):
				_key = (x_i,y_i)
				try:
					meanintensity = np.mean(ms1_array[_key]) #Current strategy for normalizing intensity is mean.
					minintensity  = min(ms1_array[_key]) #Current strategy for normalizing intensity is mean.
					maxintensity  = max(ms1_array[_key]) #Current strategy for normalizing intensity is mean.
					inputpoints   = len(ms1_array[_key]) #Amount of inputs into this array
					pixelpoint 	  = [meanintensity,minintensity,maxintensity,inputpoints]
					total_datapoints+=(len(ms1_array[_key]))
					nonzero_counter+=1
				except KeyError:
					pixelpoint = [0,0,0,0]
				row.append(pixelpoint)
			image.append(row)
		print('Saving image files                            ', end = '\r')

		imagedata = [image, nonzero_counter, total_datapoints]
		#Save as txt file
		with open(filepath+str(resolution['x'])+'x'+str(resolution['y'])+'.txt', "wb") as pa:
			pickle.dump(imagedata, pa)
		
		#######Creating the full image (.png)#######
		if savepng:
			for ll in range(4):
				fullimage = [[y[ll] for y in x] for x in image]
				if ll == 0:
					titleofplot = 'Mean'
				elif ll == 1:
					titleofplot = 'Min'
				elif ll == 2:
					titleofplot = 'Max'
				elif ll == 3:
					titleofplot = 'Collapsed'

				if not os.path.exists(filepath+str(resolution['x'])+'x'+str(resolution['y'])+'-'+titleofplot+'.png'):
					#Set color of missing data
					fullimage = np.ma.masked_equal(fullimage,0)
					#Setup colormap
					colMap = cm.jet
					colMap.set_bad('darkblue')
					if titleofplot != 'Collapsed':
						plt.imshow(fullimage,cmap=colMap,extent = [interval['mz']['min'], interval['mz']['max'], interval['rt']['min'], interval['rt']['max']],aspect = 'auto', vmin = lowbound, vmax = highbound)
					else:
						plt.imshow(fullimage,cmap=colMap,extent = [interval['mz']['min'], interval['mz']['max'], interval['rt']['min'], interval['rt']['max']],aspect = 'auto')
					plt.tight_layout()
					plt.title(titleofplot)
					plt.xlabel('m/z', fontsize=12)
					plt.ylabel('Retention time - Minutes', fontsize=12)
					plt.axis([interval['mz']['min'], interval['mz']['max'], interval['rt']['min'], interval['rt']['max']])
					# plt.colorbar(extend = 'both')
					plt.tight_layout()
					plt.savefig(filepath+str(resolution['x'])+'x'+str(resolution['y'])+'-'+titleofplot+'.png')
					plt.close()
		#############################################

	else: #If the image data exists, just recall it instead of making it
		print('Loading image data                              ', end = '\r')
		with open(filepath+str(resolution['x'])+'x'+str(resolution['y'])+'.txt', "rb") as pa:
			imagedata 		= pickle.load(pa)
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

	imgpath = path+'images/'
	if not os.path.exists(imgpath):
		os.mkdir(imgpath)

	if not os.path.exists(metapath):
		os.mkdir(metapath)	
	outfile = open(metapath+'subimage.json','a') #The metadata file
	
	i = 0
	outbound = 0
	inbound  = 0
	for index, rows in df.iterrows():
		i+=1
		if i % int(df.shape[0]/40) == 0:
			print("Creating subimages: {:2.1%}                                  ".format(i / df.shape[0]), end = '\r') #Print how far we are

		if rows['Retention time']-subimage_interval['rt'] < interval['rt']['min'] or rows['Retention time']+subimage_interval['rt'] > interval['rt']['max'] or rows['m/z']-subimage_interval['mz'] < interval['mz']['min'] or rows['m/z']+subimage_interval['mz']> interval['mz']['max']:
			outbound+=1 #Check if this image can be created in our range or not
			continue
		inbound+=1
		if os.path.exists(imgpath+filename+'-'+str(i)+'.png'):
			continue
		
		onemzlength = len(mzrangelist)/(max(mzrangelist)-min(mzrangelist))
		onertlength = len(rtrangelist)/(max(rtrangelist)-min(rtrangelist))

		mzlower = int(get_lower_bound(mzrangelist,rows['m/z']) - onemzlength*subimage_interval['mz'])
		mzupper = int(get_lower_bound(mzrangelist,rows['m/z']) + onemzlength*subimage_interval['mz']) 
		rtlower = int(get_lower_bound(rtrangelist,rows['Retention time']) - onertlength*subimage_interval['rt']) 
		rtupper = int(get_lower_bound(rtrangelist,rows['Retention time']) + onertlength*subimage_interval['rt'])
		subimage = []
		k = 0
		for lines in image:
			if rtlower < k < rtupper:
				subimage.append(lines[mzlower:mzupper])
			k+=1
			if k > rtupper:
				break

		#Save image as txt file
		with open(imgpath+filename+'-'+str(i)+'.txt', 'wb') as pa:
			pickle.dump(subimage, pa)

		########Create image########
		if savepng:
			#get the mean values from the imagedata
			newimage = [[y[0] for y in x] for x in subimage]
			#set color of missing data
			newimage = np.ma.masked_equal(newimage,0)

			colMap = cm.jet
			colMap.set_bad('darkblue')

			fig = plt.figure()
			fig.set_size_inches(2,2)#(mzupper - mzlower)/100,(rtupper - rtlower)/100)
			ax = plt.Axes(fig, [0., 0., 1, 1])
			ax.set_axis_off()
			fig.add_axes(ax)
			plt.set_cmap('hot')
			ax.imshow(newimage, aspect='equal',cmap = colMap, vmin = lowbound, vmax = highbound)
			plt.savefig(imgpath+filename+'-'+str(i)+'.png')
			plt.close()
		###########################

		new_metadata = {}
		new_metadata.update({"image" : filename+'-'+str(i)})
		for ele in df.columns[1:]:
			if str(rows[ele]) == 'nan' or str(rows[ele]) == ' ' or ";" in str(rows[ele]):
				continue
			else:
				new_metadata.update({str(ele) : str(rows[ele])})
		outfile.write(json.dumps(new_metadata)+'\n')
	outfile.close()
	
	print('Calculating end statistics:                            ', end = '\r')
	
	mzlist_inrange = [i for i in mzlist if i > interval['mz']['min'] and i < interval['mz']['max']]
	rtlist_inrange = [i for i in rtlist if i > interval['rt']['min'] and i < interval['rt']['max']]

	end_stats = {}
	end_stats['accession']			= accession
	end_stats['filename']			= filename	
	end_stats['unique mz'] 			= len(mzlist_inrange)
	end_stats['unique rt'] 			= len(rtlist_inrange)
	end_stats['datapoints'] 		= total_datapoints
	end_stats['data per pixel'] 	= total_datapoints / nonzero_counter 
	end_stats['In bounds']			= inbound	
	end_stats['Out of bounds']		= outbound

	outfile = open(metapath+'sub_statistics.json','a')
	outfile.write(json.dumps(end_stats)+'\n')
	outfile.close()
	print('Done!                                               ')


def combined(accession, maxquant_file, path, metapath):
	#Find all zip files
	output      = zipfile_finder(accession = accession, path = datapath, metapath = metapath)
	searchfiles = output[0]
	url 	    = output[1]

	for zips in searchfiles:
		#Find all raw files in the zip file
		output   = rawfile_finder(zipfile = zips, path = datapath, maxquant_file = pepfile)
		rawfiles = output[0]
		df		 = output[1]
		zipfilename = output[2]

		for raws in rawfiles:
			filename = raws 

			#Skip this special case. Something wrong
			if filename == "01625b_GA1-TUM_first_pool_1_01_01-2xIT_2xHCD-1h-R1": 
				continue

			print('\nfile: '+accession+'/'+filename) 
			print('downloading raw file                  ', end = '\r')
			output   = filehandling(filename = filename, zipfilename = zipfilename, path = datapath, maxquant_file = pepfile, df = df, url = url)
			df2 	 = output[0]
			filepath = output[1]

			formatFile(filename = filename, path = datapath, filepath = filepath)
			internalmzML(path = filepath)

			#Set the resolution for the large image, and the intervals for the smaller ones
			reso   	 = {'x':1250,'y':1000}
			interval = {'mz':10,'rt':2}
			createImages(filename = filename, path = datapath, filepath = filepath, metapath = metapath,resolution = reso, subimage_interval = interval, df = df2, savepng = False)

		os.remove(datapath+zipfilename)
		os.remove(datapath+zipfilename+'-'+pepfile)


if __name__ == '__main__':
	#Path to data
	datapath = '/data/ProteomeToolsRaw/' #Server datapath
	# datapath = 'Data/' 
	metapath = datapath+'metadata/'

	#Assigning accession number and maxquant output file name
	pepfile = 'allPeptides.txt'
	input = sys.argv[1]
	if str(input) == 'accessions' or str(input) == 'accessions_filtered':
		for line in reversed(list(open(metapath+sys.argv[1]+'.json'))):
			data  = json.loads(line)
			accession = data['accession']
			try:
				combined(accession, maxquant_file = pepfile, path = datapath, metapath = metapath)
			except KeyboardInterrupt:
				print('Problem occured with: '+accession+'. unable to proceed at this time')
				pass
	else:
		accession = input
		with open(metapath+'accessions.txt', "rb") as pa:
			pride_accessions = pickle.load(pa)
		for a in pride_accessions:
			if accession in a:
				accession = a
				break
		try:
			combined(accession, maxquant_file = pepfile, path = datapath, metapath = metapath)
		except KeyboardInterrupt:
			print('Problem occured with: '+accession+'. unable to proceed at this time')
			pass
	
# python3 pdownload.py PXD004732
# python3 pdownload.py PXD010595
# python3 pdownload.py accessions_filtered