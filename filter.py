import glob
import json
import os
import numpy as np
import random
import glob
import sys
from collections import defaultdict, Counter


def filter(path, file):
    if file == 'subimage':
        if sys.argv[2] == 'combine':
            outfile = open(f'{path}subimage.json', 'a')
            for imagejson in glob.glob(f'{datapath}subimage-*.json'):
                for line in open(imagejson, 'r'):
                    outfile.write(json.dumps(line) + '\n')
                os.remove(imagejson)
            outfile.close()
            quit()

        if os.path.exists(f'{path}subimage_filtered.json'):
            print('Removing old filtered version')
            os.remove(f'{path}{str(file)}_filtered.json')

        lines_seen = set()
        outfile = open(f'{path}{str(file)}_filtered.json', 'w')

        if sys.argv[2] == 'Sequence':
            seen = [json.loads(line)['Sequence'] for line in open(f'{path}subimage.json') if 'Sequence' in json.loads(line)]
            Seen = np.unique(seen)
            a = {}
            for f in Seen:
                a[str(f)] = seen.count(f)
            seen = [f[0] for f in Counter(a).most_common(10)]

            amounts = [len(seen[f]) for f in Seen]
            Seen = defaultdict(list)
            for f in seen:
                random.shuffle(seen[f])
                Seen[f] = seen[f][0:min(amounts)]

            for line in open(f'{path}{str(file)}.json', 'r'):
                data = json.loads(line)

                if 'size' in data and data['size'] == [166, 66, 4] and 'Sequence' in data and data['Sequence'] in Seen \
                        and line not in lines_seen:
                    data['Seq_class'] = Seen.index(data['Sequence'])
                    lines_seen.add(line)
                    outfile.write(json.dumps(data) + '\n')
            outfile.close()

        elif sys.argv[2] == 'PTM':
            seen = defaultdict(list)
            for line in open(f'{path}subimage.json'):
                if 'Modifications' in json.loads(line):
                    if not json.loads(line)['Modifications'] == 'Oxidation (M)':
                        seen[0].append(json.loads(line)['image'])
                    else:
                        seen[1].append(json.loads(line)['image'])

            amounts = [len(seen[f]) for f in seen]
            Seen = defaultdict(list)
            for f in seen:
                random.shuffle(seen[f])
                Seen[f] = seen[f][0:min(amounts)]

            for line in open(f'{path}{str(file)}.json', 'r'):
                data = json.loads(line)
                if 'size' in data and data['size'] == [166, 66, 4] and 'Modifications' in data and line not in lines_seen:
                    if data['image'] in Seen[0]:
                        data['Modi_class'] = 0
                        lines_seen.add(line)
                        outfile.write(json.dumps(data) + '\n')
                    elif data['image'] in Seen[1]:
                        data['Modi_class'] = 1
                        lines_seen.add(line)
                        outfile.write(json.dumps(data) + '\n')
            outfile.close()

        elif sys.argv[2] == 'Charge':
            seen = defaultdict(list)
            for line in open(f'{path}subimage.json'):
                if 'Charge' in json.loads(line):
                    charge = json.loads(line)['Charge']
                    seen[charge].append(json.loads(line)['image'])

            amounts = [len(seen[f]) for f in seen]
            Seen = defaultdict(list)
            for f in seen:
                random.shuffle(seen[f])
                Seen[f] = seen[f][0:min(amounts)]

            for line in open(f'{path}{str(file)}.json', 'r'):
                data = json.loads(line)

                if 'size' in data and data['size'] == [166, 66, 4] and 'Charge' in data and data['Charge'] in Seen \
                        and line not in lines_seen:
                    lines_seen.add(line)
                    outfile.write(json.dumps(data) + '\n')
            outfile.close()


    elif file == 'accessions':
        if os.path.exists(f'{path}{str(file)}_filtered.json'):
            print('Removing old filtered version')
            os.remove(f'{path}{str(file)}_filtered.json')

        lines_seen = set()
        outfile = open(f'{path}{str(file)}_filtered.json', 'w')

        for line in open(f'{path}{str(file)}.json', 'r'):
            data = json.loads(line)
            if 'allpeptides' in data and data['allpeptides'] and 'filetypes' in data and 'raw' in data['filetypes'] and line not in lines_seen:  ### FILTER HERE ###
                outfile.write(line)
                lines_seen.add(line)

                outfile.write(json.dumps(data) + '\n')
        outfile.close()
    else:
        print('no bueno')
        quit()


if __name__ == '__main__':
    # Read datapath from config file
    with open('config.json') as json_file:
        data = json.load(json_file)

    datapath = f'{data["path"]}metadata/'

    filetofilter = sys.argv[1]
    filter(datapath, filetofilter)


# python3 filehandler.py filter accessions
# python3 filehandler.py filter subimage
# python3 filehandler.py move /data/ProteomeToolsRaw/images/ /home/tochr15/images/ tochr15@yeast.imada.sdu.dk
# python3 filehandler.py move /data/ProteomeToolsRaw/metadata/ /home/tochr15/metadata/ tochr15@yeast.imada.sdu.dk
