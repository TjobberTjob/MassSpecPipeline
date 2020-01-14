import requests
import json
from zipfile import ZipFile
from os.path import join
import re
from bs4 import BeautifulSoup
import os
import pickle
import sys


def get_accessions(path):
    all_accessions = []
    accessions = 'accessions.txt'

    url = 'http://ftp.pride.ebi.ac.uk/pride/data/archive/'
    page = requests.get(url).text
    soup = BeautifulSoup(page, 'html.parser')
    # Level 1 - Getting the years
    for level1 in soup.find_all('a', href=True):
        if '20' not in level1['href']:
            continue  # To avoid the non years
        page_2 = requests.get(f'{url}{level1["href"]}').text
        soup_2 = BeautifulSoup(page_2, 'html.parser')
        # Level 2 - year subfolders
        for level2 in soup_2.find_all('a', href=True):
            try:
                int(level2['href'][0:2])  # Only look at numerics...
            except:
                continue
            print(f'Getting accessions from {level1.text}{level2.text}', end='\r')
            page_3 = requests.get(f'{url}{level1["href"]}{level2["href"]}').text
            soup_3 = BeautifulSoup(page_3, 'html.parser')

            # Level 3- actual accession numbers.
            for level3 in soup_3.find_all('a', href=True):
                accession = level3['href'].replace('/', '')
                if len(accession) == len('PRD000000'):
                    all_accessions.append(accession)

    with open(f'{path}{accessions}', "wb") as pa:
        pickle.dump(all_accessions, pa)


def accessions_metadata(file_list, path):
    metafile = 'accessions.json'
    outfile = open(join(metapath, metafile), 'a')
    for i, f in enumerate(file_list):
        try:
            print('Progress {:2.1%}                                       '.format(i / len(file_list)), end='\r')
            api = f'https://www.ebi.ac.uk/pride/ws/archive/project/{f}'
            apijson = requests.get(api).json()

            metafile = {}
            maxquant = False
            for g in apijson:
                if ':' not in str(apijson[g]):
                    metafile[g] = apijson[g]
                if 'maxquant' in str(apijson[g]).lower():
                    maxquant = True
            metafile['maxquant'] = maxquant

            files = f'https://www.ebi.ac.uk/pride/ws/archive/file/list/project/{f}'
            filesjson = requests.get(files).json()
            filetypes = []
            for f in filesjson['list']:
                filetype = f['fileName'][re.search('\.', f['fileName']).span()[1]:]
                if filetype not in filetypes:
                    filetypes.append(filetype)
            metafile['filetypes'] = filetypes

            if metafile['maxquant'] and 'zip' in metafile['filetypes']:
                try:
                    for f in filesjson['list']:
                        filetype = f['fileName'][re.search('\.', f['fileName']).span()[1]:]
                        if f['fileType'] == 'SEARCH' and filetype == 'zip':
                            zipfile = f['downloadLink']
                            break
                    os.system(f'wget -q -O {path}file.zip {zipfile}')

                    with ZipFile(f'{path}file.zip', 'r') as zipped:
                        ziplist = zipped.namelist()
                    os.remove(f'{path}file.zip')

                    for xx in ziplist:
                        if 'allPeptides.txt' in xx:
                            metafile['allpeptides'] = True
                            break
                except:
                    metafile['allpeptides'] = False
            else:
                metafile['allpeptides'] = False

            outfile.write(json.dumps(metafile) + '\n')
        except:
            pass


def update_metadata(mpath):
    metadata = 'accessions.json'
    accessions = 'accessions.txt'

    if os.path.exists(f'{mpath}{accessions}'):
        os.remove(f'{mpath}{accessions}')
    get_accessions(path=mpath)

    with open(f'{mpath}{accessions}', "rb") as pa:
        pride_accessions = pickle.load(pa)  # Loading the data

    all_accessions = [f for f in pride_accessions]
    had_accessions = [json.loads(line)['accession'] for line in open(f'{mpath}{metadata}') if
                      'accession' in json.loads(line)]
    missing_accessions = [f for f in all_accessions if f not in had_accessions]

    accessions_metadata(missing_accessions, mpath)


def validated_input(prompt, valid_values):
    valid_input = False
    while not valid_input:
        value = input(prompt + ' | ' + ' / '.join(valid_values) + "\n")
        valid_input = value in valid_values
    return value


if __name__ == '__main__':
    # Read datapath from config file
    with open('config.json') as json_file:
        data = json.load(json_file)

    datapath = data['path']
    metapath = f'{datapath}metadata/'
    if not os.path.exists(metapath):
        os.mkdir(metapath)

    cmd = sys.argv[1]

    # Download pride accession numbers.
    if cmd == 'accessions':
        accessions = 'accessions.txt'
        if os.path.exists(f'{metapath}{accessions}'):
            rm = validated_input('File already exists, overwrite?', ('y', 'n'))
            if rm == 'y':
                os.remove(f'{metapath}{accessions}')
            if rm == 'n':
                quit()
        get_accessions(path=metapath)

    # Get Metadata for all accession numbers
    if cmd == 'metadata':
        metadata = 'accessions.json'
        accessions = 'accessions.txt'

        if os.path.exists(f'{metapath}{metadata}'):
            overwrite = validated_input('Metadata already exists, wanna overwrite?', ('y', 'n'))
            if overwrite == 'y':
                os.remove(f'{metapath}{metadata}')
            else:
                quit()

        with open(f'{metapath}{accessions}', "rb") as pa:
            pride_accessions = pickle.load(pa)  # Loading the data

        accessions_metadata(pride_accessions, metapath)

    if cmd == 'update':
        update_metadata(metapath)

# python3 pwebscrape.py accessions
# python3 pwebscrape.py metadata
