import gzip
import json
import os
i = 0
a = [f'/data/ProteomeToolsRaw/{f}/{g}/allPeptides.txt' for f in os.listdir('/data/ProteomeToolsRaw/') for g in os.listdir(f'/data/ProteomeToolsRaw/{f}') if os.path.exists(f'/data/ProteomeToolsRaw/{f}/{g}/allPeptides.txt')]
print(a)
quit()
for f in os.listdir('/data/ProteomeToolsRaw/'):
    i += 1
    if os.path.isdir(f'/data/ProteomeToolsRaw/{f}') and f[0:3] == 'PXD' or f[0:3] == 'PRD':
        for g in os.listdir(f'/data/ProteomeToolsRaw/{f}'):
            if g == '20160704_UCL_Mouse_IgG_30min_2a':
                print(f'{i} / {str(len(os.listdir("/data/ProteomeToolsRaw/")))}')

