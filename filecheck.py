i = 0
for line in open(datapath+'metadata.json'):
	data = json.loads(line)
	imname = data['image']
	if not os.path.exists("/data/ProteomeToolsRaw/Images/"+imname):
		i = i+1
		print("doesnt exist")
print(i)