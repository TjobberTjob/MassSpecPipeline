import json
import os
import glob
import sys
import shutil
import random
import re

i = 0
j=0
datapath = "/data/ProteomeToolsRaw/Images/"
for line in open(datapath+'metadata.json'):
	j = j+1
	if j%10000 == 0:
		print(j)
	data = json.loads(line)
	imname = data['image']
	if not os.path.exists(datapath+imname):
		i = i+1
		print("doesnt exist")
print(i)