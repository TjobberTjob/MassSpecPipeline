import json
import os
import glob
import sys
import shutil
import random
import re
datapath = "/data/ProteomeToolsRaw/Images/"
trainpath = datapath+'training/'
valpath = datapath+'validation/'

if not os.path.exists(datapath):
	os.mkdir(datapath)

#Reset images
dirs = [os.path.dirname(p) for p in glob.glob(datapath+"/*/*")]
udirs = [] 
for x in dirs:  
	if x not in udirs: 
		udirs.append(x)


if dirs != [] and (dirs[0] == trainpath[:-1] or dirs[0] == valpath[:-1]):
	dirs = [os.path.dirname(p) for p in glob.glob(datapath+"/*/*/*")]
	udirs = [] 
	for x in dirs:  
		if x not in udirs:
			udirs.append(x)


if not udirs == [ ]:
	print("Do you want to reset image folders? y/n")
	reset = input()
	if reset == "yes" or reset == "y":
		for files in udirs:
			images = [f for f in glob.glob(files + "**/*.png", recursive=True)]
			for imgs in images:
				os.system("mv "+imgs+" "+datapath)
			shutil.rmtree(files) 	
	else: quit()

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

		for line in open(datapath+'metadata.json'):
			data = json.loads(line)
			if not os.path.exists(trainpath+data[imClass]):
				os.mkdir(trainpath+data[imClass])
			os.system("mv "+datapath+data['image']+".png "+trainpath+data[imClass]+"/")	

		dirs = [os.path.dirname(p) for p in glob.glob(trainpath+"/*/*")]
		udirs = [] 
		for x in dirs:  
			if x not in udirs: 
				udirs.append(x)

		for f in udirs:
			folderclass = f[[m.start() for m in re.finditer('/', f)][-1]+1:]
			if not os.path.exists(valpath+str(folderclass)):
				os.mkdir(valpath+str(folderclass))
			mlist = os.listdir(f+"/")
			print(mlist)
			splits = round(len(os.listdir(f))*(int(splitratio)/100))
			print(splits)
			mlist = random.choices(mlist,k=splits)
			if mlist == [[]]:
				continue
			else:
				for image in mlist:
					os.system("mv "+trainpath+str(folderclass)+"/"+str(image)+" "+valpath+str(folderclass))

	else:
		for line in open(datapath+'metadata.json'):
			data = json.loads(line)
			if not os.path.exists(datapath+data[imClass]):
				os.mkdir(datapath+data[imClass])
			os.system("mv "+datapath+data['image']+".png "+datapath+data[imClass]+"/")	


if __name__ == '__main__':

	imClass = sys.argv[1]

	if imClass == "reset":
		quit()
	else:
		classifyImages(classes = imClass)


