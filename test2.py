import gzip
import json
import os
i = 0

for f in os.listdir('/data/ProteomeToolsRaw/'):
    i += 1
    if os.path.isdir(f'/data/ProteomeToolsRaw/{f}') and f[0:3] == 'PXD' or f[0:3] == 'PRD':
        for g in os.listdir(f'/data/ProteomeToolsRaw/{f}'):
            if g == 'Animal_20150123_ASK_RAN_TAM_C':
                print(f'{i} / {str(len(os.listdir("/data/ProteomeToolsRaw/")))}')

