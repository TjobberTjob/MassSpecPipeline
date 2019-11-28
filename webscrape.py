if __name__ == '__main__':
	import requests
	import json
	from zipfile import ZipFile
	from os.path import join
	import time
	import re
	import random
	from bs4 import BeautifulSoup
	import os
	from urllib import parse
	import pickle
	import codecs
	import sys


def get_accessions(path):
	all_accessions = []
	accessions	 = 'accessions.txt'
	if os.path.exists(path+accessions):
		rm = validated_input('File already exists, overwrite?',('y','n'))
		if rm == 'y':
			os.remove(path+accessions)

	url  = 'http://ftp.pride.ebi.ac.uk/pride/data/archive/'
	page = requests.get(url).text	
	soup = BeautifulSoup(page, 'html.parser')
	#Level 1 - Getting the years
	for level1 in soup.find_all('a', href = True):
		if '20' not in level1['href']:
			continue #To avoid the non years
		page_2 = requests.get(url+level1['href']).text
		soup_2 = BeautifulSoup(page_2, 'html.parser')
		#Level 2 - year subfolders
		for level2 in soup_2.find_all('a', href = True):
			try:
				int(level2['href'][0:2]) #Only look at numerics...
			except Exception:
				continue

			print('Getting accessions from '+level1['href']+level2['href'], end = '\r')

			page_3 = requests.get(url+level1['href']+level2['href']).text
			soup_3 = BeautifulSoup(page_3, 'html.parser')

			#Level 3- actual accession numbers.
			for level3 in soup_3.find_all('a', href = True):
				accession = level3['href'].replace('/','')
				if len(accession) == len('PRD000000'):
					all_accessions.append(accession)

	with open(path+accessions, "wb") as pa:
		pickle.dump(all_accessions, pa)


def accessions_metadata(path):
	metadata = 'accessions.json'
	accessions	 = 'accessions.txt'
	if os.path.exists(path+metadata):
		overwrite = validated_input('Metadata already exists, wanna overwrite?',('y','n'))
		if overwrite == 'y':
			os.system('rm '+path+metadata)
		else:
			quit()

	with open(path+accessions, "rb") as pa:
		pride_accessions = pickle.load(pa) #Loading the data
	i = 0
	outfile = open(join(path,metadata),'a')

	for f in pride_accessions:
		i   += 1
		url  = 'https://www.ebi.ac.uk/pride/archive/projects/'+f
		html = requests.get(url).text
		soup = BeautifulSoup(html,'html.parser')
		my_dict = {}
		my_dict['accession'] = f
		for div in soup.find_all('div', {'class': 'grid_7 right-column'}):
			for div2 in div.find_all('div', {'class': 'grid_12'}): 
				try:
					name  = div2.find('h5').text 
				except Exception:
					continue
				try:
					value = div2.find('a').text 
				except Exception:
					value = "Not available"
				my_dict[name] = value 

		for div in soup.find_all('div', {'class': 'grid_16 left-column'}):
			plist = []
			for p in div.find_all('p'):
				plist.append(p.text.strip())
			my_dict['Submission Date'] = plist[-2]
			my_dict['Publication Date'] = plist[-1]

		url  = 'https://www.ebi.ac.uk/pride/archive/projects/'+f+'/files'
		html = requests.get(url).text
		soup = BeautifulSoup(html,'html.parser')

		filetypes = []
		for div in soup.find_all('div', {'class': 'grid_23 clearfix file-list'}):
			for h5 in div.find_all('h5'):
				name  = h5.text[re.search(' ',h5.text).span()[1]:]
				value = h5.text[:re.search(' ',h5.text).span()[0]]
				my_dict[name] = value

			for ta in  div.find_all('table'):
				try:
					filetypes.append(ta.find('td').text.strip()[re.search('\.',ta.find('td').text.strip()).span()[1]:])
				except Exception:
					print('f')
			my_dict['filetypes'] = filetypes

		
		url = 'https://www.ebi.ac.uk/pride/archive/projects/'+f+'/files'
		page = requests.get(url).text
		my_dict['maxquant'] = 'maxquant' in page.lower()

		#Check for allpeptides.txt		
		if my_dict['maxquant'] == True and '.zip' in my_dict['filetypes']:
			for div in soup.find_all('div', {'class': 'grid_6 omega'}):
				url = div.find('a')['href'] #Update URL with FTP link
				break
			if my_dict['maxquant'] == True and '.zip' in my_dict['filetypes']:
				os.system('wget -q -O '+metapath+'readme.txt '+url+'/README.txt')
			df = pd.read_csv(metapath+'readme.txt',sep='\t')
			os.remove(metapath+'readme.txt')
			searchfiles = df.loc[df['TYPE'] == 'SEARCH',]['URI']
			os.system('wget -q -O '+metapath+'file.zip '+searchfiles[0])

			with ZipFile(path+'file.zip','r') as zipped:
				ziplist = zipped.namelist()

			for xx in ziplist:
				if 'allPeptides.txt' in xx:
					my_dict['allpeptides'] == True
					break
		else:
			my_dict['allpeptides'] == False

		print('Progress {:2.1%}'.format(i / len(pride_accessions)), end='\r')

		outfile.write(json.dumps(my_dict)+'\n')
	outfile.close()


def validated_input(prompt, valid_values):
	valid_input = False
	while not valid_input:
		value	   = input(prompt + ' | ' + ' / '.join(valid_values)+"\n")
		valid_input = value in valid_values
	return value


if __name__ == '__main__':
	datapath = 'Data/'
	# datapath = '/data/ProteomeToolsRaw/'
	metapath = datapath+'metadata/'

	cmd = sys.argv[1]
	
	#Download pride accession numbers.
	if cmd == 'accessions':
		get_accessions(path = metapath)

	#Get Metadata for all accession numbers
	if cmd == 'metadata':
		accessions_metadata(path = metapath)

	#python3 webscrape.py accessions
	#python3 webscrape.py metadata
