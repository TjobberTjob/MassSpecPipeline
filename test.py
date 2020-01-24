import gzip
import json

for line in open(f'{path}{str(file)}.json', 'r'):
    data = json.loads(line)
    with gzip.GzipFile('/data/ProteomeToolsRaw/PXD010595/02101a_GA8-TUM_proteo_TMT_8_01_01-ETD-1h-R1/mzML2.json', 'w') as fout:
        fout.write(json.dumps(data).encode('utf-8'))