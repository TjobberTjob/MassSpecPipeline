import gzip
import json
import os

for f in os.listdir('/data/ProteomeToolsRaw/'):
    if os.path.isdir(f) and f[0:3] == 'PXD' or f[0:3] == 'PRD':
        print(f)
        quit()