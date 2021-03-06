import bisect
import glob
import gzip
import json
import math
import os
import pickle
import shutil
import subprocess
import sys
import time
from multiprocessing.dummy import Pool as ThreadPool
from pathlib import Path
from zipfile import ZipFile
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from pyteomics import mzml


def get_lower_bound(haystack, needle):
    idx = bisect.bisect(haystack, needle)
    if idx > 0 and idx < len(haystack):
        return idx
    else:
        raise ValueError(f"{needle} is out of bounds of {haystack}")


def pride_file_finder(accnr, path):
    url = f'https://www.ebi.ac.uk/pride/ws/archive/file/list/project/{accnr}/'
    serverIssue = False
    privacyIssue = False
    try:
        urljson = requests.get(url).json()
        zipfiles = []
        rawfiles = []

        # If zipfiles have the same name as rawfiles and we have the allpeptides, dont download
        for jsonelem in urljson['list']:
            filetype = jsonelem['fileName'].split('.')[-1]
            if (jsonelem['fileType'] == 'SEARCH' or jsonelem['fileType'] == 'OTHER') and filetype == 'zip' and \
                    jsonelem['downloadLink'][-9:] != "Fasta.zip":
                zipfiles.append(jsonelem['downloadLink'])
            if jsonelem['fileType'] == 'RAW' and filetype == 'raw':
                rawfiles.append(jsonelem['downloadLink'])
    except:
        if not multiprocessing:
            print("API connection issue")
        if 'status' in urljson.keys() and urljson['status'] == 'UNAUTHORIZED':
            privacyIssue = True
        else:
            serverIssue = True
        return privacyIssue, serverIssue, []

    allCheck = ['allPeptides.txt' in os.listdir(f'{path}{accnr}/{files}/') for files in os.listdir(f'{path}{accnr}/')
                if len(os.listdir(f'{path}{accnr}/')) == len(rawfiles)]
    if False in allCheck or allCheck == []:
        haveallMQF = False
    else:
        haveallMQF = True

    return zipfiles, rawfiles, haveallMQF


def zipfile_downloader(zipfile, path, maxquant_file):
    # Handle spaces in urls
    zipfileurl = zipfile.translate(str.maketrans({"(": r"\(", ")": r"\)", " ": r"%20"}))
    zipfilename = zipfile.replace(' ', '-')[63:].replace('(', '-').replace(')', '-')

    # Download zip file
    if os.path.exists(f'{path}{zipfilename}'):
        os.remove(f'{path}{zipfilename}')
    if not multiprocessing:
        print('Downloading zip file                                                    ', end='\r')
    while not os.path.exists(f'{path}{zipfilename}'):
        os.system(f'wget -q -O {path}{zipfilename} {zipfileurl}')
        if os.stat(f'{path}{zipfilename}').st_size == 0:
            os.remove(f'{path}{zipfilename}')
        time.sleep(5)

    # os.system(f'curl {zipfileurl} --output {path}{zipfilename}')

    # Get a list of files with directories from zip file
    with ZipFile(f'{path}{zipfilename}', 'r') as zipped:
        ziplist = zipped.namelist()

    # Extract the peptide file from the zipfile
    for a in ziplist:
        if maxquant_file in a:
            with ZipFile(f'{path}{zipfilename}') as z:
                with z.open(a) as zf, open(f'{path}{zipfilename[:-4]}-{maxquant_file}', 'wb') as zfg:
                    shutil.copyfileobj(zf, zfg)
                break

    # Go through the maxquant output file and get all the raw files
    df = pd.read_csv(f'{path}{zipfilename[:-4]}-{maxquant_file}', sep='\t', low_memory=False)
    df = df.loc[df['Sequence'] != ' ',]  # Remove empty sequences
    rawfiles = np.unique(df['Raw file'])

    os.remove(f'{path}{zipfilename}')  # Remove useless zipfile
    os.remove(f'{path}{zipfilename[:-4]}-{maxquant_file}')
    return rawfiles, df


def file_handler(accnr, filename, path, maxquant_file, df, rawfiles):
    filepath = f'{path}{accnr}/{filename}/'
    # Make the file directory if it doesnt exist
    if not os.path.exists(filepath):
        os.mkdir(filepath)

    # Check if filespecific allPeptides.txt exists
    df2 = df.loc[df['Raw file'] == filename,]
    if not os.path.exists(f'{filepath}{maxquant_file}'):
        pd.DataFrame.to_csv(df2, f'{filepath}{maxquant_file}')

    # Download the raw file
    if not multiprocessing:
        print('Downloading raw file                                                    ', end='\r')

    if not (os.path.exists(f'{filepath}file.mzML') or os.path.exists(f'{filepath}mzML.json')):
        for fileraw in rawfiles:
            if filename in fileraw:
                os.system(f'wget -q -c -O {filepath}file.raw {fileraw}')
                # os.system(f'curl {fileraw} --output {filepath}file.raw')
                break

    return df2, filepath


def file_formatter(accnr, filename, path, filepath):
    if not multiprocessing:
        print('Formatting file to mzML										', end='\r')

    with open('config.json') as json_file:
        config = json.load(json_file)
    formatusing = config['format_software']

    if not (os.path.exists(f'{filepath}file.mzML') or os.path.exists(f'{filepath}mzML.json')):
        if formatusing == 'conda':
            if path[0] == '/':
                relpath = path
            else:
                relpath = f'{os.getcwd()}{path}'

            condalist = subprocess.check_output('conda list', shell=True)
            if not 'thermorawfileparser' in str(condalist):
                os.system('conda install -c bioconda thermorawfileparser')

            condalist = subprocess.check_output('conda list', shell=True)
            if not 'thermorawfileparser' in str(condalist):
                print('Conda issue. Cannot install thermorawfileparser, try re-installing conda')
                quit()

            os.system(
                f'mono /opt/conda/bin/ThermoRawFileParser.exe -i={relpath}{accnr}/{filename}/file.raw -o={relpath}{accnr}/{filename}/ -f=1 -m=1 >/dev/null 2>&1')

        elif formatusing == 'docker':
            # Check whether the docker file is implemented or not
            dockerls = subprocess.check_output('docker image ls', shell=True)
            try:
                if not 'thermorawparser' in str(dockerls):
                    if not os.path.exists(f'{Path(os.getcwd()).parent}/ThermoRawFileParser'):
                        os.mkdir(f'{Path(os.getcwd()).parent}/ThermoRawFileParser')
                        os.system(
                            f'git clone https://github.com/compomics/ThermoRawFileParser.git {Path(os.getcwd()).parent}/ThermoRawFileParser')
                        os.system(
                            'cd .. && cd ThermoRawFileParser/ && docker build --no-cache -t thermorawparser . && cd '
                            '../MassSpecPipeline/')
                    else:
                        os.system(
                            'cd .. && cd ThermoRawFileParser/ && docker build --no-cache -t thermorawparser . && cd '
                            '../MassSpecPipeline/')
            except:
                print('Docker issue. Cannot install thermorawfileparser, try re-installing docker')
                quit()
                return

            if not 'thermorawparser' in str(dockerls):
                print('Docker issue. Cannot install thermorawfileparser, try re-installing docker')
                quit()

            if path[0] == '/':
                relpath = path[:-1]
            else:
                relpath = f'{os.getcwd()}{path[:-1]}'  # Either gives path as root path or have data as a sub folder to the one the code is in

            os.system(f'chmod -R a+rwx {path}*')
            if os.path.exists(f'{filepath}file.raw'):
                os.system(f'docker run -v "{relpath}:/data_input" -i -t thermorawparser mono '
                          f'bin/x64/Debug/ThermoRawFileParser.exe -i=/data_input/{accnr}/{filename}/file.raw -o=/data_input/{accnr}/{filename}/ -f=1 -m=1')

        else:
            print('format software not specified in config.json. Choices are "conda" or "docker"')
            quit()

        if os.path.exists(f'{filepath}file-metadata.txt'):
            os.remove(f'{filepath}file-metadata.txt')
        if os.path.exists(f'{filepath}file.raw'):
            os.remove(f'{filepath}file.raw')


def process_ms1(spectrum):
    # Scan information
    scan_info = spectrum['scanList']
    # Time
    scan_time = scan_info['scan'][0]['scan start time']
    mz = spectrum['m/z array']
    # ion intensity
    intensity = spectrum['intensity array']
    return {'scan_time': scan_time, 'intensity': intensity.tolist(), 'mz': mz.tolist()}


def process_ms2(spectrum):
    # Fish out the precursors.
    try:
        precursors = spectrum['precursorList']
        if precursors['count'] != 1:
            if not multiprocessing:
                print("Number of precursors different than 1, not designed for that")
            quit()
        ion = precursors['precursor'][0]['selectedIonList']
        if ion['count'] != 1:
            if not multiprocessing:
                print("More then one selected ions, not designed for that")
            quit()

        ion = ion['selectedIon'][0]['selected ion m/z']
        ms1_scan = int(precursors['precursor'][0]['spectrumRef'].split('scan=')[1])

        # Fish out the scan index
        scan_index = spectrum['index']

        scan_info = spectrum['scanList']
        # m/z and intensity arrays
        mz = spectrum['m/z array']
        intensity = spectrum['intensity array']
    except:
        scan_index = []
        ion = []
        ms1_scan = []
        mz = []
        intensity = []

    return {'scan_index': scan_index, 'precursor_scan': ms1_scan, 'precursor_ion': ion, 'm/z': mz, 'rt': intensity}


def extract_from_mzml(path):
    # Extract the data from the mzml, if we havnt already
    if not os.path.exists(f'{path}mzML.json'):
        if not multiprocessing:
            print('Extracting data from mzML                                                    ', end='\r')
        data = mzml.MzML(f'{path}file.mzML')

        # Extracted data
        extracted = {'ms1': {}, 'ms2': {}}
        # Extract the necessary data from spectra
        for spectrum in data:
            if spectrum['ms level'] == 1:
                # Scan id
                scan_id = int(spectrum['id'].split('scan=')[1])

                # Deal with ms level 1 spectra
                ms1_spectrum = process_ms1(spectrum)
                extracted['ms1'][scan_id] = {'mz': ms1_spectrum['mz'],
                                             'intensity': ms1_spectrum['intensity'],
                                             'scan_time': ms1_spectrum['scan_time']}

            elif spectrum['ms level'] == 2:
                # Scan id
                scan_id = int(spectrum['id'].split('scan=')[1])

                # Deal with ms level 1 spectra
                ms2_spectrum = process_ms2(spectrum)
                extracted['ms2'][scan_id] = {'scan_index': ms2_spectrum['scan_index'],
                                             'precursor_scan': ms2_spectrum['precursor_scan'],
                                             'precursor_ion': ms2_spectrum['precursor_ion'],
                                             'm/z_array': [mz for mz in ms2_spectrum['m/z']],
                                             'rt_array': [rt for rt in ms2_spectrum['rt']]}

            else:
                pass

        with gzip.GzipFile(f'{path}mzML.json', 'w') as fout:
            fout.write(json.dumps(extracted).encode('utf-8'))
        fout.close()
        # os.remove(f'{path}file.mzML')


def image_preparameters(filepath):
    if not multiprocessing:
        print('Preparing parameter for image creation                                                    ', end='\r')
    with gzip.GzipFile(f'{filepath}mzML.json', 'r') as fin:
        mzml = json.loads(fin.read().decode('utf-8'))

    mzlist = np.unique(sorted([item for f in mzml['ms1'] for item in mzml['ms1'][f]['mz']]))
    rtlist = [mzml['ms1'][f]['scan_time'] for f in mzml['ms1']]
    intlist = [item for f in mzml['ms1'] for item in mzml['ms1'][f]['intensity']]

    lowbound = math.log(np.percentile(intlist, 0.5))
    highbound = math.log(np.percentile(intlist, 99.5))

    interval = {
        'mz': {'min': min(mzlist), 'max': max(mzlist)},
        'rt': {'min': min(rtlist), 'max': max(rtlist)}
    }

    # Define the intervals for the given resolution
    with open('config.json') as json_file:
        config = json.load(json_file)
    mz_bin = float(config['mz_bin'])
    rt_bin = float(config['rt_bin'])

    resolution = {'x': int((max(mzlist) - min(mzlist)) / mz_bin), 'y': int((max(rtlist) - min(rtlist)) / rt_bin)}

    return mzml, [lowbound, highbound], interval, [mz_bin, rt_bin], resolution


def full_png_image(image, filepath, resolution, interval, lowbound, highbound):
    listnames = ['Mean', 'Min', 'Max', 'Collapsed']
    for i in range(4):
        fullimage = [[y[i] for y in x] for x in image]
        fullimage.reverse()
        titleofplot = listnames[i]

        if not os.path.exists(f'{filepath}{str(resolution["x"])}x{str(resolution["y"])}-{titleofplot}.png'):
            # Set color of missing data
            fullimage = np.ma.masked_equal(fullimage, 0)
            # Setup colormap
            colMap = cm.jet
            colMap.set_bad('darkblue')
            if titleofplot != 'Collapsed':
                plt.imshow(fullimage, cmap=colMap,
                           extent=[interval['mz']['min'], interval['mz']['max'], interval['rt']['min'],
                                   interval['rt']['max']], aspect='auto', vmin=lowbound, vmax=highbound)
            else:
                plt.imshow(fullimage, cmap=colMap,
                           extent=[interval['mz']['min'], interval['mz']['max'], interval['rt']['min'],
                                   interval['rt']['max']], aspect='auto')
            plt.tight_layout()
            plt.title(titleofplot)
            plt.xlabel('m/z', fontsize=12)
            plt.ylabel('Retention time - Minutes', fontsize=12)
            plt.axis([interval['mz']['min'], interval['mz']['max'], interval['rt']['min'], interval['rt']['max']])
            if titleofplot == 'Collapsed':
                plt.colorbar(extend='both')
            plt.tight_layout()
            plt.savefig(f'{filepath}{str(resolution["x"])}x{str(resolution["y"])}-{titleofplot}.png')
            plt.close()


def full_txt_image(mzmlfile, interval, bins, resolution, filepath, bounds, savepng):
    mz_bin = bins[0]
    rt_bin = bins[1]
    # Create an empty array for layer use
    ms1_array = {}

    # Get sorted list of scan ids.
    scan_ids = [int(scan_id) for scan_id in mzmlfile['ms1']]

    for scan_id in sorted(scan_ids):
        scan_id = str(scan_id)

        # Get the intervals
        scan_time = float(mzmlfile['ms1'][scan_id]['scan_time'])
        if scan_time < interval['rt']['min'] or scan_time > interval['rt']['max']:
            continue

        # Calculate the y axis.
        y_n = int((scan_time - interval['rt']['min']) / rt_bin)
        for index, mz_elem in enumerate(mzmlfile['ms1'][scan_id]['mz']):
            if mz_elem < interval['mz']['min'] or mz_elem > interval['mz']['max']:
                continue
            x_n = int((mz_elem - interval['mz']['min']) / mz_bin)
            _key = (x_n, y_n)
            # Current strategy for collapsing the intensity values is taking their logs
            intensity_val = math.log(mzmlfile['ms1'][scan_id]['intensity'][index])
            try:
                ms1_array[_key].append(intensity_val)
            except KeyError:
                ms1_array[_key] = [intensity_val]

    # Create the final image.
    image = []
    for y_i in range(0, resolution['y']):
        if y_i % 25 == 0:
            if not multiprocessing:
                print('Creating full image: {:2.1%}                                                    '.format(
                    y_i / resolution['y']), end='\r')  # Print how far we are
        row = []
        for x_i in range(0, resolution['x']):
            _key = (x_i, y_i)
            try:
                meanintensity = np.mean(ms1_array[_key])  # Current strategy for normalizing intensity is mean.
                minintensity = min(ms1_array[_key])
                maxintensity = max(ms1_array[_key])
                inputpoints = len(ms1_array[_key])  # Amount of inputs into this point
                pixelpoint = [meanintensity, minintensity, maxintensity, inputpoints]
            except KeyError:
                pixelpoint = [0, 0, 0, 0]
            row.append(pixelpoint)
        image.append(row)
    if not multiprocessing:
        print('Saving image files                                                    ', end='\r')

    imagedata = [image, interval, bins, resolution, bounds]

    # Save as txt file
    with open(f'{filepath}{str(resolution["x"])}x{str(resolution["y"])}.txt', "wb") as pa:
        pickle.dump(imagedata, pa)
    lowbound = bounds[0]
    highbound = bounds[1]

    if savepng:  # save full image to png
        full_png_image(image, filepath, resolution, interval, lowbound, highbound)
    return image


def sub_png_image(subimage, imgpath, filename, index, lowbound, highbound):
    newimage = [[y[0] for y in x] for x in subimage]
    newimage.reverse()
    newimage = np.ma.masked_equal(newimage, 0)

    colMap = cm.jet
    colMap.set_bad('darkblue')

    fig = plt.figure()
    fig.set_size_inches(2, 2)  # (mzupper - mzlower)/100,(rtupper - rtlower)/100)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    plt.set_cmap('hot')
    ax.imshow(newimage, aspect='equal', cmap=colMap, vmin=lowbound, vmax=highbound)
    plt.savefig(f'{imgpath}{filename}-{str(index + 1)}.png')
    plt.close()


def sub_txt_image(accnr, interval, bins, image, bounds, resolution, mzmlfile, path, mpath, df, subimage_interval,
                  filename, multiprocessing,
                  savepng):
    lowbound = bounds[0]
    highbound = bounds[1]
    mz_bin = bins[0]
    rt_bin = bins[1]
    mzrangelist = [interval['mz']['min'] + i * mz_bin for i in range(int(resolution['x']))]
    rtrangelist = [interval['rt']['min'] + i * rt_bin for i in range(int(resolution['y']))]

    imgpath = f'{path}images/'
    if not os.path.exists(imgpath):
        os.mkdir(imgpath)

    if not os.path.exists(mpath):
        os.mkdir(mpath)

    outfile = open(f'{mpath}subimage-{accnr}.json', 'a')

    df.reset_index(drop=True, inplace=True)
    for index, rows in df.iterrows():
        if int((index + 1) / int(df.shape[0]) * 100) % 5 == 0:
            if not multiprocessing:
                print(f'Creating subimages: {int(((index + 1) / df.shape[0]) * 100)}%         ', end='\r')  # Print how far we are

        # if not 450 < rows['m/z'] < 455:  # Add a filter here to increase speed and reduce storage
        #     continue

        if rows['Retention time'] - subimage_interval['rt'] < interval['rt']['min'] or rows['Retention time'] + \
                subimage_interval['rt'] > interval['rt']['max'] or rows['m/z'] - subimage_interval['mz'] < \
                interval['mz']['min'] or rows['m/z'] + subimage_interval['mz'] > interval['mz']['max']:
            continue

        # Find the closest and highest intensity scan in msms scan list
        scannumbers = rows['MSMS Scan Numbers'].split(';')
        ms2list = [[f, mzmlfile['ms2'][f]['precursor_ion'], max(mzmlfile['ms2'][f]['rt_array'])] for f in scannumbers]
        sortedms2list = sorted(ms2list, key=lambda x: x[2], reverse=True)
        closestmz = [f[1] for f in sortedms2list][(np.abs(np.array([f[1] for f in sortedms2list]) - rows['m/z'])).argmin()]
        ms2scan = [f for f in sortedms2list if f[1] == closestmz][0][0]

        if os.path.exists(f'{imgpath}{accnr}-{filename}-{ms2scan}.json'):
            continue

        mzlen = int(subimage_interval['mz'] / mz_bin)
        rtlen = int(subimage_interval['rt'] / rt_bin)

        mzlower = int(get_lower_bound(mzrangelist, rows['m/z']) - mzlen)
        mzupper = int(get_lower_bound(mzrangelist, rows['m/z']) + mzlen)
        rtlower = int(get_lower_bound(rtrangelist, rows['Retention time']) - rtlen)
        rtupper = int(get_lower_bound(rtrangelist, rows['Retention time']) + rtlen)

        ms1info = [lines[mzlower:mzupper] for lines in image[rtlower:rtupper]]
        ms1size = str([f for f in np.array(ms1info).shape])

        ms2info = [mzmlfile['ms2'][ms2scan]['m/z_array'], np.log([mzmlfile['ms2'][ms2scan]['rt_array']]).tolist()[0]]
        ms2info = [[ms2info[0][index], ms2info[1][index]] for index in range(len(ms2info[0]))]
        ms2size = str([f for f in np.array(ms2info).shape])

        fullsubimage = {'ms1': ms1info, 'ms2': ms2info}

        # Save image as json file
        imageoutfile = open(f'{imgpath}{accnr}-{filename}-{ms2scan}.json', 'w')
        imageoutfile.write(json.dumps(fullsubimage))

        if savepng:  # save subimages to png
            sub_png_image(ms1info, imgpath, filename, index, lowbound, highbound)

        new_metadata = {}
        new_metadata['image'] = f'{accnr}-{filename}-{ms2scan}.json'
        new_metadata['accession'] = accnr
        new_metadata['ms1size'] = ms1size
        new_metadata['ms2size'] = ms2size
        for ele in df.columns:
            if str(rows[ele]) == 'nan' or str(rows[ele]) == ' ' or ";" in str(rows[ele]):
                continue
            else:
                new_metadata[str(ele)] = str(rows[ele])
        outfile.write(json.dumps(new_metadata) + '\n')
    outfile.close()


def local_file_main(path, filename, mpath):
    global zipfile, rawfile
    maxquant_file = 'allPeptides.txt'
    filepath = f'{sysinput}{filename}/'
    accnr = sysinput.split('/')[-2:-1][0]
    print(f'\nfile: {accnr}/{filename}')

    for file in os.listdir(f'{filepath}'):
        if file.endswith('.zip'):
            zipfile = file
    for file in os.listdir(f'{filepath}'):
        if file.endswith('.raw'):
            rawfile = file
    if os.path.exists(f'{filepath}{maxquant_file}'):
        allPep = True
    else:
        allPep = False

    if 'rawfile' in locals() and (allPep or 'zipfile' in locals()):
        if not os.path.exists(f'{filepath}{maxquant_file}'):
            # Get a list of files with directories from zip file
            with ZipFile(f'{filepath}{zipfile}', 'r') as zipped:
                ziplist = zipped.namelist()

            # Extract the peptide file from the zipfile
            for a in ziplist:
                if maxquant_file in a:
                    with ZipFile(f'{filepath}{zipfile}') as z:
                        with z.open(a) as zf, open(f'{filepath}allPeptides.txt', 'wb') as zfg:
                            shutil.copyfileobj(zf, zfg)
                        break

            df = pd.read_csv(f'{filepath}{maxquant_file}', sep='\t', low_memory=False)
            df = df.loc[df['Sequence'] != ' ',]  # Remove empty sequences
            df = df.loc[df['Raw file'] == rawfile,]
            pd.DataFrame.to_csv(df, f'{filepath}{maxquant_file}')
        else:
            df = pd.read_csv(f'{filepath}{maxquant_file}', sep='\t', low_memory=False)

        main_p2(accnr, filename, path, mpath, filepath, df)
    else:
        try:
            os.system('rm /data/ProteomeToolsRaw/*.*')
        except:
            pass
        if not multiprocessing:
            print(f'Necessary files dont exist in {f}')
        quit()


def main_p2(accnr, filename, path, mpath, filepath, df2, multiprocessing):
    file_formatter(accnr, filename, path, filepath)
    extract_from_mzml(filepath)

    imagefiles = glob.glob(f'{path}{accnr}/{filename}/*x*.txt')
    binsmatch = False
    if not imagefiles == []:
        for imagefile in imagefiles:
            if not multiprocessing:
                print(f'Fetching image files                                               ', end='\r')
            with open(imagefile, "rb") as pa:
                imagedata = pickle.load(pa)
            image = imagedata[0]
            interval = imagedata[1]
            bins = imagedata[2]
            resolution = imagedata[3]
            bounds = imagedata[4]

            with open('config.json') as configjson:
                config = json.load(configjson)
            mz_bin = float(config['mz_bin'])
            rt_bin = float(config['rt_bin'])

            if bins[0] == mz_bin and bins[1] == rt_bin:
                with gzip.GzipFile(f'{filepath}mzML.json', 'r') as fin:
                    mzml = json.loads(fin.read().decode('utf-8'))
                binsmatch = True
                break

    with open('config.json') as configjson:
        config = json.load(configjson)
    savepng = config['save_images_as_png'] == 'True'

    if not binsmatch:
        output = image_preparameters(filepath)
        mzml = output[0]
        bounds = output[1]
        interval = output[2]
        bins = output[3]
        resolution = output[4]

        output = full_txt_image(mzml, interval, bins, resolution, filepath, bounds, savepng)
        image = output[0]

    subimage_interval = {'mz': config['mz_interval'], 'rt': config['rt_interval']}
    sub_txt_image(accnr, interval, bins, image, bounds, resolution, mzml, path, mpath, df2, subimage_interval, filename,
                  multiprocessing,
                  savepng)


def main_p1(accnr, maxquant_file, path, mpath, multiprocessing):
    if not multiprocessing:
        print(f'\nAccessions: {accnr}')

    if not os.path.exists(f'{path}{accnr}/'):
        os.mkdir(f'{path}{accnr}/')

    # Find all zip files
    output = pride_file_finder(accnr, path)
    allZip = output[0]
    allRaw = output[1]
    haveallMQF = output[2]
    if isinstance(allZip, bool) or isinstance(allRaw, bool):
        serverissue = allZip
        restrictedissue = allRaw

    if "serverissue" in globals() and serverissue:
        print(f'{accnr}: ✔ - 0/{len(allRaw)} PRIDE Servers cannot be reached')
        return
    if "restrictedissue" in globals() and restrictedissue:
        print(f'{accnr}: ✔ - 0/{len(allRaw)} Restricted PRIDE project')
        return

    # skip files if skip-incomplete or acquire_only_new is true
    if haveallMQF:
        if acquire_only_new:
            print(f'{accnr}: ✔ - {len(allRaw)}/{len(allRaw)} Rawfiles extracted')
            return
    else:
        if skip_incomplete:
            print(f'Accession: {accnr}: ✔ - Skipping')
            return

    for zips in reversed(allZip):
        try:  # TRY ALL ZIPS
            if not haveallMQF:  # If we dont have all needed files, we need to get them from the API
                output = zipfile_downloader(zips, path, maxquant_file)
                rawfiles = output[0]
                df = output[1]

                knownrawfiles = []
                for filenames in rawfiles:
                    for onlinenames in allRaw:
                        if filenames in onlinenames:
                            knownrawfiles.append(filenames)

                for filename in knownrawfiles:
                    try:  # TRY ALL RAWS IN ZIP
                        if not multiprocessing:
                            print(f'file: {accnr}/{filename}                                               ')

                        output = file_handler(accnr, filename, path, pepfile, df, allRaw)
                        df2 = output[0]
                        filepath = output[1]

                        main_p2(accnr, filename, path, mpath, filepath, df2, multiprocessing)
                        if not multiprocessing:
                            print(f'{filename}: ✔                         ')

                    except Exception as error:
                        exc_type, exc_obj, exc_tb = sys.exc_info()
                        if not multiprocessing:
                            if errormessages:
                                print(f'Rawfile error. {filename}: ✖ | Error Class: {exc_type} |'
                                      f' Error: {error} | Line: {exc_tb.tb_lineno}')
                                del (exc_type, exc_obj, exc_tb)
                            else:
                                print(
                                    f'Rawfile error. {filename}: ✖')
                        pass

            else:  # If we have all needed files, we dont need to get them from the API
                rawfiles = os.listdir(f'{path}{accnr}')
                knownrawfiles = []
                for filenames in rawfiles:
                    for onlinenames in allRaw:
                        if filenames in onlinenames:
                            knownrawfiles.append(filenames)

                for filename in knownrawfiles:
                    try:  # TRY ALL RAWS IN ZIP
                        if not multiprocessing:
                            print(f'\nfile: {accnr}/{filename}                                               ')

                        # set path and import allpeptides
                        filepath = f'{path}{accnr}/{filename}/'
                        df2 = pd.read_csv(f'{filepath}{maxquant_file}', sep=',', low_memory=False)

                        # Download raw file
                        if not (os.path.exists(f'{filepath}file.mzML') or os.path.exists(f'{filepath}mzML.json')):
                            for rawfiles in allRaw:
                                if filename in rawfiles:
                                    os.system(f'wget -q -c -O {filepath}file.raw {rawfiles}')
                                    # os.system(f'curl {fileraw} --output {filepath}file.raw')
                                    break

                        main_p2(accnr, filename, path, mpath, filepath, df2, multiprocessing)
                        if not multiprocessing:
                            print(f'{filename}: ✔                         ')

                    except Exception as error:
                        exc_type, exc_obj, exc_tb = sys.exc_info()
                        if not multiprocessing:
                            if errormessages:
                                print(f'Rawfile error. {filename}: ✖ | Error Class: {exc_type} |'
                                      f' Error: {error} | Line: {exc_tb.tb_lineno}')
                                del (exc_type, exc_obj, exc_tb)
                            else:
                                print(
                                    f'Rawfile error. {filename}: ✖')
                        pass

        except Exception as error:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            if not multiprocessing:
                if errormessages:
                    print(f'Zipfile error. {zips.split("/")[-1]}: ✖ | Error Class: {exc_type} |'
                          f' Error: {error} | Line: {exc_tb.tb_lineno}')  # 'issue occoured, going to next zipfile')
                    del (exc_type, exc_obj, exc_tb)
                else:
                    print(f'Zipfile error. {zips.split("/")[-1]}: ✖')

            if os.path.exists(f'{path}{zips.replace(" ", "-")[63:].replace("(", "-").replace(")", "-")}'):
                os.remove(f'{path}{zips.replace(" ", "-")[63:].replace("(", "-").replace(")", "-")}')
            pass

    allCheck = [files for files in os.listdir(f'{path}{accnr}/') if 'mzML.json' in os.listdir(f'{path}{accnr}/{files}')]
    if len(allCheck) == len(allRaw):
        print(f'{accnr}: ✔ - {len(allCheck)}/{len(allRaw)} All rawfiles extracted')
    else:
        if len(allCheck) != len(allRaw):
            print(f'{accnr}: ✖ - {len(allCheck)}/{len(allRaw)} Some or all MaxQuant file names do not '
                  f'match PRIDE file names')
        else:
            print(f'{accnr}: ✖ - {len(allCheck)}/{len(allRaw)} Rawfiles extracted')


if __name__ == '__main__':
    # Read datapath from config file
    with open('config.json') as json_file:
        data = json.load(json_file)

    datapath = data['path']
    if not os.path.exists(datapath):
        print('Path given in config.json doesn\'t exist. Please specify working path.')
        quit()

    metapath = f'{datapath}metadata/'
    acquire_only_new = data['acquire_only_new'] == 'True'
    multi = data['multi_processing'] == 'True'
    nr_processes = data['nr_processes']
    errormessages = data['receive_error_messages'] == 'True'
    skip_incomplete = False

    # Assigning accession number and maxquant output file name
    pepfile = 'allPeptides.txt'
    sysinput = sys.argv[1]

    # Options#
    if str(sysinput) == 'reset':  # Reset files and folders if you want to remake all images in another setting
        if os.path.exists(f'{datapath}images/'):
            shutil.rmtree(f'{datapath}images/')
        for imagejson in glob.glob(f'{metapath}subimage*.json'):
            os.remove(imagejson)

    elif str(sysinput)[0] == '/':  # For local fine purposes.
        dirsinpath = os.listdir(sysinput)
        for f in dirsinpath:
            local_file_main(datapath, f, metapath)

    elif str(sysinput) == 'complete':  # For re-creating images from already downloaded and parsed files
        listofowned = [f for f in os.listdir(datapath) if
                       os.path.isdir(f'{datapath}{f}') and f[0:3] == 'PRD' or f[0:3] == 'PXD']
        skip_incomplete = True
        for accession in listofowned:
            if multi:

                multiprocessing = True
                accessions = [(f, pepfile, datapath, metapath, multiprocessing) for f in listofowned]
                pool = ThreadPool(nr_processes)
                pool.starmap(main_p1, accessions)

            else:
                multiprocessing = False
                main_p1(str(accession), pepfile, datapath, metapath, multiprocessing)

    elif str(sysinput) == 'pride' or str(sysinput) == 'pridefiltered':  # Going through the metadata
        if str(sysinput) == 'pride':
            pridefile = 'accessions'
        else:
            pridefile = 'accessions_filtered'

        if multi:
            multiprocessing = True
            accessions = [(json.loads(linez)['accession'], pepfile, datapath, metapath, multiprocessing)
                          for linez in reversed(list(open(f'{metapath}{pridefile}.json'))) if
                          'accession' in json.loads(linez) and json.loads(linez)["maxquant"]]
            pool = ThreadPool(nr_processes)
            pool.starmap(main_p1, accessions)

        else:
            multiprocessing = False
            for line in reversed(list(open(f'{metapath}{pridefile}.json'))):
                data = json.loads(line)
                accession = data['accession']
                main_p1(str(accession), pepfile, datapath, metapath, multiprocessing)

    elif str(sysinput)[0:3] == 'PRD' or str(sysinput)[0:3] == 'PXD':  # For single accessions usage
        accession = sysinput
        multiprocessing = False
        main_p1(str(accession), pepfile, datapath, metapath, multiprocessing)

    else:
        print('Input not recognized. Check readme file for all possible inputs')
        quit()

# python3 extractor.py PXD004732
# python3 extractor.py PXD010595
# python3 extractor.py pridefiltered
# python3 extractor.py owned
# python3 extractor.py /mnt/c/Users/TobiaGC/Dropbox/Universitet/CompBiomed/Speciale/MassSpecPipeline/Data/PXD010595/


# Seq_class (4)  val_loss: 0.0092 - val_accuracy: 0.9775
# Seq_class (10) val_loss: 0.7285 - val_accuracy: 0.8244
# m/z val_mse: 4000
# Length val_accuracy: 0.5160
