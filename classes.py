if 'import' == 'import':
	import json
	import os
	from glob import glob
	import sys
	import shutil
	import random
	import re
	import numpy as np
	from collections import defaultdict

def validated_input(prompt, valid_values):
	valid_input = False
	while not valid_input:
		value = input(prompt + ' | ' + ' / '.join(valid_values)+"\n")
		valid_input = value in valid_values
	return value

def numbered_input(prompt, min, max):
	valid_input = False
	while valid_input == False:
		value = input(prompt+'\n')
		valid_input = int(value) >= min and int(value) <= max
	return value

#Pathfinding
datapath = "/data/ProteomeToolsRaw/Images/"
trainpath = datapath+'training/'
valpath = datapath+'validation/'

#Reset images

dirs = glob(datapath+"/*/")
udirs = np.unique(dirs)

if dirs != [] and (dirs[0] == trainpath[:-1] or dirs[0] == valpath[:-1]):
	dirs = [os.path.dirname(p) for p in glob.glob(datapath+"/*/*/*")]
	udirs = np.unique(dirs) 

if len(udirs) != 0:
	reset = validated_input('Do you want to reset image folders?', ('y','n'))
	if reset == 'y':
		os.system('find '+str(datapath)+' -mindepth 2 -name \"*.png\" -exec mv -t '+ str(datapath)+ ' {} +')
		os.system("rm -rf "+trainpath)
		os.system("rm -rf "+valpath)


#Move images
def classifyImages(classes):
	split = validated_input('Do you wanna split the data into training and validation?', ('y','n'))
	if split == "yes" or split == "y":
		splitratio = numbered_input("What should the validation % be?",0,100)

		if not os.path.exists(trainpath):
			os.mkdir(trainpath)
		if not os.path.exists(valpath):
			os.mkdir(valpath)

		imgdata = {}
		
		#Preparing metadata
		print("Preparing Metadata")
		for line in open(datapath+'metadata_filtered.json'):
			data = json.loads(line)

			names = data['image']+".png"
			imgdata[names] = data[imClass]
	
			if not os.path.exists(trainpath+data[imClass]):
				os.mkdir(trainpath+data[imClass])
			if not os.path.exists(valpath+data[imClass]):
				os.mkdir(valpath+data[imClass])

		#CREATING TRAINING DAT
		print("Sorting into training data")
		imgadata = defaultdict(list)
		for k, v in imgdata.items():
			imgadata[v].append(k)
		for f in imgadata:
			for g in imgadata[f]:
				try:
					shutil.move(datapath+g, trainpath+f+"/"+g)
				except Exception:
					pass
		
		#CREATING VALIDATION DATA
		print("Sorting into Validation data")
		for f in imgadata:
			splits = round(len(imgadata[f])*(int(splitratio)/100))
			mlist = random.sample(imgadata[f],k=splits)
			for g in mlist:
				try:
					shutil.move(trainpath+f+"/"+g, valpath+f+"/"+g)
				except Exception:
					pass

	elif split == "no" or split == "n":
		for line in open(datapath+'metadata_filtered.json'):
			data = json.loads(line)
			if not os.path.exists(datapath+data[imClass]):
				os.mkdir(datapath+data[imClass])
			try:
				shutil.move(datapath+data['image']+".png", datapath+data[imClass]+"/")
			except Exception:
				pass
	else: quit()


if __name__ == '__main__':


	i = 1
	while i == 1:
		for line in open(datapath+'metadata_filtered.json'):
			data = json.loads(line)
			i = i+1
	distlist = list(data.keys())
	distlist.append('reset')
	imClass = validated_input('What do you want to classify based on?',distlist)

	if imClass == 'reset':
		quit()
	else:
		classifyImages(classes = imClass)


