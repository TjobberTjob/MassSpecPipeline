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


def resetImage(path, trainpath, valpath):
	os.system('find '+str(path)+' -mindepth 2 -name \"*.png\" -exec mv -t '+ str(path)+ ' {} +')
	try:
		shutil.rmtree(trainpath)
		shutil.rmtree(valpath)
	except Exception:
		pass


def classifyImages(path, trainpath, valpath, imgClass):
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
		for line in open(path+'metadata_filtered.json'):
			data = json.loads(line)

			names = data['image']+".png"
			imgdata[names] = data[imgClass]
	
			if not os.path.exists(trainpath+data[imgClass]):
				os.mkdir(trainpath+data[imgClass])
			if not os.path.exists(valpath+data[imgClass]):
				os.mkdir(valpath+data[imgClass])

		#CREATING TRAINING DAT
		print("Sorting into training data")
		imgadata = defaultdict(list)
		for k, v in imgdata.items():
			imgadata[v].append(k)
		for f in imgadata:
			for g in imgadata[f]:
				try:
					shutil.move(path+g, trainpath+f+"/"+g)
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
		for line in open(path+'metadata_filtered.json'):
			data = json.loads(line)
			if not os.path.exists(path+data[imgClass]):
				os.mkdir(path+data[imgClass])
			try:
				shutil.move(path+data['image']+".png", path+data[imgClass]+"/")
			except Exception:
				pass
	else: quit()


if __name__ == '__main__':

	datapath = 'Data/Images/'
	# datapath = "/data/ProteomeToolsRaw/Images/"
	trainpath = datapath+'training/'
	valpath = datapath+'validation/'

	for line in open(datapath+'metadata_filtered.json'):
		data = json.loads(line)
		break
	distlist = list(data.keys())
	distlist.append('reset')
	Class = validated_input('What do you want to classify based on?',distlist)

	if Class == 'reset':
		resetImage(path = datapath, trainpath = trainpath, valpath = valpath)
	else:
		classifyImages(path = datapath, trainpath = trainpath, valpath = valpath, imgClass = Class)


