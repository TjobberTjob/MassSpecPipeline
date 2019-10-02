# Used to move .raw and .zip files from directory into a folder of same name
import os 
import glob
datapath = "/data/ProteomeToolsRaw/"
for files in [f for f in glob.glob(datapath + "**/*.raw", recursive=True)]:
	if not os.path.exists(files[:-4]+"/"):
		os.mkdir(files[:-4]+"/")
	os.system("mv "+files+" "+files[:-4])

for files in [f for f in glob.glob(datapath + "**/*.zip", recursive=True)]:
	if not os.path.exists(files[:-4]+"/"):
		os.mkdir(files[:-4]+"/")
	os.system("mv "+files+" "+files[:-4])

