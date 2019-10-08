i = 0
datapath = "/data/ProteomeToolsRaw/Images/"
for line in open(datapath+'metadata.json'):
	data = json.loads(line)
	imname = data['image']
	if not os.path.exists(datapath+imname):
		i = i+1
		print("doesnt exist")
print(i)