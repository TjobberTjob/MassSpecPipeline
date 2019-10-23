datapath = "/data/ProteomeToolsRaw/tryptic/"

import glob
files = [os.path.dirname(p) for p in glob.glob(datapath+"/*/*")]
print(files)