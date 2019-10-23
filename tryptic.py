datapath = "/data/ProteomeToolsRaw/tryptic/"
import os 
import glob
files = [os.path.dirname(p) for p in glob.glob(datapath+"/*/*.zip")]
print(files)