import requests
import json
from os.path import join
import time
import random
from bs4 import BeautifulSoup
import os
import pickle
import codecs
import sys

accessions = 'pride_accessions.txt'
all_accessions = []
def get_accessions():

	if os.path.exists(datapath+accessions):
		rm = validated_input('File already exists, overwrite?',('y','n'))
		if rm == 'y':
			os.system('rm '+accessions)

	url = 'http://ftp.pride.ebi.ac.uk/pride/data/archive/'
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

def validated_input(prompt, valid_values):
	valid_input = False
	while not valid_input:
		value = input(prompt + ' | ' + ' / '.join(valid_values)+"\n")
		valid_input = value in valid_values
	return value

 
if __name__ == '__main__':
	datapath = 'data/'

	with open(datapath+accessions, "rb") as pa:
		b = pickle.load(pa)
	print(len(b))
	
	cmd = sys.argv[1]

	#Download pride accession numbers.
	if cmd == 'accession':
		get_accessions()