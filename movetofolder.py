import os 
import glob
datapath = "data/ProteomeToolsRaw/"
for files in [f for f in glob.glob(datapath + "**/*", recursive=True)]:
	if not os.path.exists(files[:-4]+"/"):
		os.mkdir(files[:-4]+"/")
	os.system("mv "+files+" "+files[:-4])

