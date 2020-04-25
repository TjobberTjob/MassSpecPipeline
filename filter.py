import glob
import json
import os
import numpy as np
import random
import glob
import sys
from collections import defaultdict, Counter


def mostcommon(x, amount, path):
    seen = [json.loads(line)[x] for line in open(f'{path}subimage.json') if x in json.loads(line)]
    Seen = np.unique(seen)
    a = {}
    for f in Seen:
        a[str(f)] = seen.count(f)
    seen = [f[0] for f in Counter(a).most_common(amount)]

    amounts = [len(seen[f]) for f in Seen]
    Seen = defaultdict(list)
    for f in seen:
        random.shuffle(seen[f])
        Seen[f] = seen[f][0:min(amounts)]

    return Seen


def filter(path, file):
    if file == 'subimage':
        if sys.argv[2] == 'combine':

            if os.path.exists(f'{path}subimage.json'):
                allimgs = [json.loads(line)['image'] for line in open(f'{path}subimage.json') if
                           'image' in json.loads(line)]
            else:
                allimgs = []

            with open(f'{path}subimage.json', 'a') as outfile:
                for imagejson in glob.glob(f'{path}subimage-*.json'):
                    for line in open(imagejson, 'r'):
                        data = json.loads(line)
                        if data['image'] not in allimgs:
                            outfile.write(json.dumps(data) + '\n')
                    # os.remove(imagejson)
            quit()

        lines_seen = set()
        outfile = open(f'{path}{str(file)}_filtered.json', 'w')

        getsizes = [json.loads(lines)['size'] for lines in open(f'{path}subimage.json') if 'size' in json.loads(lines)]
        mostcommonsize = np.unique(getsizes)
        a = {}
        for f in mostcommonsize:
            a[str(f)] = getsizes.count(f)
        for f in Counter(a).most_common(2):
            f2 = str(f[0]).replace(" ", "")[1:-1].split(',')
            if len(f2) == 3:
                ms1size = f[0]
                break

        getscores = [float(json.loads(lines)['Score']) for lines in open(f'{path}subimage.json') if
                     'Score' in json.loads(lines)
                     and 'size' in json.loads(lines) and str(json.loads(lines)['size']) == ms1size]
        getabovehere = np.percentile(getscores, 60)

        if sys.argv[2] == 'Sequence':
            Seen = mostcommon('Sequence', 10)

            for line in open(f'{path}{str(file)}.json', 'r'):
                data = json.loads(line)

                if 'size' in data and data['size'] == str(ms1size) and 'Sequence' in data and data[
                    'Sequence'] in Seen \
                        and line not in lines_seen:
                    data['Seq_class'] = Seen.index(data['Sequence'])
                    lines_seen.add(line)
                    outfile.write(json.dumps(data) + '\n')
            outfile.close()


        elif sys.argv[2] == 'Length':

            seen = [json.loads(line)['Length'] for line in open(f'{path}subimage.json') if
                    'Length' in json.loads(line)]
            Seen = np.unique(seen)
            Seen = sorted(Seen)
            for line in open(f'{path}{str(file)}.json', 'r'):
                data = json.loads(line)

                if 'size' in data and data['size'] == str(ms1size) and 'Length' in data and \
                        data['Length'] in Seen and 'Score' in data and float(data['Score']) > getabovehere:
                    data['Length_class'] = Seen.index(data['Length'])
                    lines_seen.add(line)
                    outfile.write(json.dumps(data) + '\n')
            outfile.close()


        elif sys.argv[2] == 'PTM':
            seen = defaultdict(list)
            for line in open(f'{path}subimage.json'):
                if 'Modifications' in json.loads(line):
                    if not json.loads(line)['Modifications'] == 'Unmodified':
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
                if 'size' in data and data['size'] == str(
                        ms1size) and 'Modifications' in data and line not in lines_seen:
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
                data = json.loads(line)
                if 'Charge' in data and 'size' in data and str(data['size']) == ms1size and 'Score' in data \
                        and float(data['Score']) > getabovehere:
                    charge = json.loads(line)['Charge']
                    seen[charge].append(json.loads(line)['image'])
            amounts = defaultdict(list)
            for f in seen:
                amounts[f] = len(seen[f])
            minamount = min(f for f in amounts.values() if f >= (0.1 * sum(amounts.values())))

            Seen = defaultdict(list)
            for f in seen:
                if len(seen[f]) >= minamount:
                    random.shuffle(seen[f])
                    Seen[f] = seen[f][0:minamount]


            i = 0
            for line in open(f'{path}{str(file)}.json', 'r'):
                data = json.loads(line)

                if 'charge' in data:
                    charge = data['charge']
                    if data['image'] in Seen[charge]:
                        lines_seen.add(line)
                        outfile.write(json.dumps(data) + '\n')
                        i += 1
            outfile.close()
            print(f'Length of filtered file: {i}')


    elif file == 'accessions':
        if os.path.exists(f'{path}{str(file)}_filtered.json'):
            print('Removing old filtered version')
            os.remove(f'{path}{str(file)}_filtered.json')

        lines_seen = set()
        outfile = open(f'{path}{str(file)}_filtered.json', 'w')

        for line in open(f'{path}{str(file)}.json', 'r'):
            data = json.loads(line)
            if 'allpeptides' in data and data['allpeptides'] and 'filetypes' in data and 'raw' in data[
                'filetypes'] and line not in lines_seen:
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
