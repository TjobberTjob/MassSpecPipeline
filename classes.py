if 'import' == 'import':
	import json
	import os
	from glob import glob
	import sys
	import shutil
	import random
	import re
	import pickle
	import numpy as np
	from collections import defaultdict
	import pandas as pd

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


def classifyImages(path, trainpath, valpath, metapath, imgClass):
	splitratio = numbered_input("What should the training % be?",0,100)
	
	if not os.path.exists(trainpath):
		os.mkdir(trainpath)
	if not os.path.exists(valpath):
		os.mkdir(valpath)

	
	if classorreg == 'classification':
		imgdata = {}
		#Preparing metadata
		print("Preparing Metadata")
		for line in open(metapath+'subimage_filtered.json'):
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

	elif classorreg == 'regression':
		imgdata = []
		#Preparing metadata
		print("Preparing Metadata")
		for line in open(metapath+'subimage_filtered.json'):
			data = json.loads(line)

			name = data['image']+".png"
			Class = data[imgClass]
			imgdata.append([name,Class])
		splits = round(len(imgdata)*(int(splitratio)/100))
		random.shuffle(imgdata)
		trainlist   = imgdata[0:splits]
		vallist 	= imgdata[splits:]
		

		print('Moving training files')
		a = []
		for f in trainlist:
			a.append(f[0])

		for f in a:
			try:
				# os.system('mv '+path+f+' '+trainpath+f)
				shutil.move(path+f, trainpath+f)
			except Exception:
				pass
			df = pd.DataFrame(trainlist, columns = ['image','class'])
			with open(trainpath+'data.txt', "wb") as pa:
				pickle.dump(df, pa)
		
		print('Moving validation files')
		a = []
		for f in vallist:
			a.append(f[0])
		for f in a:
			try:
				# os.system('mv '+path+f' '+valpath+f)
				shutil.move(path+f, valpath+f)
			except Exception:
				pass
			df = pd.DataFrame(vallist, columns = ['image','class'])
			with open(valpath+'data.txt', "wb") as pa:
				pickle.dump(df, pa)
	else: quit()


if __name__ == '__main__':
	classorreg = sys.argv[1]

	datapath  = 'Data/'
	imagepath = 'Data/images/'
	# datapath  = "/data/ProteomeToolsRaw/"
	# imagepath  = "/data/ProteomeToolsRaw/images/"
	trainpath = imagepath+'training/'
	valpath   = imagepath+'validation/'
	metapath  = datapath+'metadata/'

	for line in open(metapath+'subimage_filtered.json'):
		data = json.loads(line)
		break
	distlist = list(data.keys())
	distlist.append('reset')
	Class = validated_input('What do you want to classify based on?',distlist)

	if Class == 'reset':
		resetImage(path = imagepath, trainpath = trainpath, valpath = valpath)
	else:
		classifyImages(path = imagepath, trainpath = trainpath, valpath = valpath,  metapath = metapath, imgClass = Class)


