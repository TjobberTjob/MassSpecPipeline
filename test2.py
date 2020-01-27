import gzip
import json
import os
i = 0
a = [f'/data/ProteomeToolsRaw/{f}/{g}/allPeptides.txt' for f in os.listdir('/data/ProteomeToolsRaw/') for g in os.listdir(f'/data/ProteomeToolsRaw/{f}') if os.path.exists(f'/data/ProteomeToolsRaw/{f}/{g}/allPeptides.txt')]
for f in os.listdir('/data/ProteomeToolsRaw/'):
    if os.path.isdir(f'/data/ProteomeToolsRaw/{f}') and f[0:3] == 'PXD' or f[0:3] == 'PRD':
        for g in os.listdir(f'/data/ProteomeToolsRaw/{f}'):
            i += 1
            if g == 'SH_mitoD_7_190615_run1':
                print(f'{i} / {str(len(a))}')

