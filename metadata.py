import json
import os
import glob
import sys
import shutil
import random
import re

datapath = "/data/ProteomeToolsRaw/Images/"
i = 0
while i < 2:
	i+=1
	for line in open(datapath+'metadata.json'):
		data = json.loads(line)
		print(data)
