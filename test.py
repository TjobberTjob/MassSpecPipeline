import gzip
import json
import os

for f in os.listdir('/data/ProteomeToolsRaw/'):
    if os.path.isdir(f'/data/ProteomeToolsRaw/{f}') and f[0:3] == 'PXD' or f[0:3] == 'PRD':
        print(f'/data/ProteomeToolsRaw/{f}')
        for g in os.listdir(f'/data/ProteomeToolsRaw/{f}'):
            path = f'/data/ProteomeToolsRaw/{f}/{g}/'
            print(path)
            data = json.load(open(f'{path}mzML.json'))

            with gzip.GzipFile(f'{path}mzML2.json', 'w') as fout:
                fout.write(json.dumps(data).encode('utf-8'))
            quit()