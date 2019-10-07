import json
import os
import glob
import sys
import shutil
import random
import re

datapath = "/data/ProteomeToolsRaw/Images/"
i = 0
for line in open(datapath+'metadata.json'):
	i+=1
	data = json.loads(line)
	if data['image'] == "01974C_BA1-TUM_missing_first_1_01_01-DDA-1h-R4-930.png":
		print(i)
