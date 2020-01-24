import gzip
import json

# for line in open(f'/data/ProteomeToolsRaw/PXD010595/02101a_GA8-TUM_proteo_TMT_8_01_01-ETD-1h-R1/mzML.json', 'r'):
#     data = json.loads(line)
#     with gzip.GzipFile('/data/ProteomeToolsRaw/PXD010595/02101a_GA8-TUM_proteo_TMT_8_01_01-ETD-1h-R1/mzML2.json', 'w') as fout:
#         fout.write(json.dumps(data).encode('utf-8'))
i = 0
while i < 10:
    with gzip.GzipFile('/data/ProteomeToolsRaw/PXD010595/02101a_GA8-TUM_proteo_TMT_8_01_01-ETD-1h-R1/mzML2.json', 'r') as fin:
        data = json.loads(fin.read().decode('utf-8'))
    print(data)
    i+=1
