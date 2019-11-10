import requests
import json
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

accessions     = 'pride_accessions.txt'
all_accessions = []
def get_accessions():

	if os.path.exists(datapath+accessions):
		rm = validated_input('File already exists, overwrite?',('y','n'))
		if rm == 'y':
			os.system('rm '+accessions)

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
			except AttributeError :
				continue

			print('Going into '+level1['href']+level2['href'])

			
			page_3 = requests.get(url+level1['href']+level2['href']).text
			soup_3 = BeautifulSoup(page_3, 'html.parser')
			i = 0
			#Level 3- actual accession numbers.
			for level3 in soup_3.find_all('a', href = True):
				accession = level3['href'].replace('/','')
				if len(accession) == len('PXD000000'):
					all_accessions.append(accession)

	with open(datapath+accessions, "wb") as pa:
		pickle.dump(all_accessions, pa)

def accessions_metadata():

	metadata = 'accession_metadata.json'
	if os.path.exists(datapath+metadata):
		overwrite = validated_input('Metadata already exists, wanna overwrite?',('y','n'))
		if overwrite == 'y':
			os.system('rm '+datapath+metadata)
		else:
			quit()

	with open(datapath+accessions, "rb") as pa:
		pride_accessions = pickle.load(pa) #Loading the data
	i = 0
	outfile = open(join(datapath,metadata),'a')

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
				except AttributeError :
					continue
				try:
					value = div2.find('a').text 
				except AttributeError :
					value = "Not available"
				my_dict[name] = value 

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
				except AttributeError:
					print('f')
			my_dict['filetypes'] = filetypes

		url  = 'https://www.ebi.ac.uk/pride/ws/archive/project/'
		page = requests.get(url).text
		my_dict['maxquant'] = 'maxquant' in page.lower()
		print('Progress {:2.1%}'.format(i / len(pride_accessions)), end='\r')

		outfile.write(json.dumps(my_dict)+'\n')
	outfile.close()

def filtering():

	try:
		os.system('rm '+path+'accession_filtered.json')
	except Exception:
		print('no filtered version exist')

	outfile = open(path+'accession_filtered.json','w')
	lines = set()
	for line in open(path+'accession_metadata.json','r'):
		data = json.loads(line)

		if data['maxquant'] == True: #Add filter here
			outfilewrite(line)
			lines.add(line)
	outfile.close

def validated_input(prompt, valid_values):
	valid_input = False
	while not valid_input:
		value       = input(prompt + ' | ' + ' / '.join(valid_values)+"\n")
		valid_input = value in valid_values
	return value
 
if __name__ == '__main__':
	datapath = '/data/ProteomeToolsRaw/'

	cmd = sys.argv[1]

	#Download pride accession numbers.
	if cmd == 'accession':
		get_accessions()

	if cmd == 'metadata':
		accessions_metadata()

	if cmd == 'filter':
		filtering()
