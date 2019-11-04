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
	print("Do you want to reset image folders? y/n")
	reset = input()
	if reset == "yes" or reset == "y":
		os.system('find '+str(datapath)+' -mindepth 2 -name \"*.png\" -exec mv -t '+ str(datapath)+ ' {} +')
		os.system("rm -rf "+trainpath)
		os.system("rm -rf "+valpath)


#Move images
def classifyImages(classes):
	print("Do you wanna split the data into training and validation? y/n")
	split = input()
	if split == "yes" or split == "y":

		print("What should the validation % be?")
		splitratio = input()

		if not os.path.exists(trainpath):
			os.mkdir(trainpath)
		if not os.path.exists(valpath):
			os.mkdir(valpath)

		imgdata = {}
		
		#Preparing metadata
		print("Preparing Metadata")
		for line in open(datapath+'metadata_cleaned.json'):
			try:
				data = json.loads(line)
			except Exception:
				print(" ")
			print(data['m/z'])
			print(float(data['m/z']) > 424 and float(data['m/z']) < 425)
			quit()
			if float(data['m/z']) > 424 and float(data['m/z']) < 425:
				continue
			else:
				names = data['image']+".png"
				imgdata[names] = data[imClass]
				
				if not os.path.exists(trainpath+data[imClass]):
					os.mkdir(trainpath+data[imClass])
				
				if not os.path.exists(valpath+data[imClass]):
					os.mkdir(valpath+data[imClass])

		#CREATING TRAINING DATA
		print("Sorting into training data")
		imgadata = defaultdict(list)
		for k, v in imgdata.items():
			imgadata[v].append(k)
		for f in imgadata:
			for g in imgadata[f]:
				shutil.move(datapath+g, trainpath+f+"/"+g)
		
		#CREATING VALIDATION DATA
		print("Sorting into Validation data")
		for f in imgadata:
			splits = round(len(imgadata[f])*(int(splitratio)/100))
			mlist = random.sample(imgadata[f],k=splits)
			for g in mlist:
				shutil.move(trainpath+f+"/"+g, valpath+f+"/"+g)

	elif split == "no" or split == "n":
		for line in open(datapath+'metadata_cleaned.json'):
			data = json.loads(line)
			if not os.path.exists(datapath+data[imClass]):
				os.mkdir(datapath+data[imClass])
			shutil.move(datapath+data['image']+".png", datpath+data[imClass]+"/")
	else: quit()


if __name__ == '__main__':

	imClass = sys.argv[1]

	if imClass == "reset":
		quit()
	else:
		classifyImages(classes = imClass)


