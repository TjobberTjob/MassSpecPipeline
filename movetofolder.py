import os 
import glob
for files in [f for f in glob.glob(path + "**/*.raw", recursive=True)]:
	if not os.path.exists(datapath+files[:-4]):
		os.mkdir(datapath+files[:-4])
	os.system("mv "+datapath+"/"+files+" "files[:-4])

