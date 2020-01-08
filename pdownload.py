import bisect
import json
import math
import os
import pickle
import re
import shutil
import subprocess
import sys
import time
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


def filefinder(accnr):
    url = 'https://www.ebi.ac.uk/pride/ws/archive/file/list/project/' + accnr
    urljson = requests.get(url).json()
    zipfiles = []
    rawfiles = []

    for f in urljson['list']:
        filetype = f['fileName'][re.search('\.', f['fileName']).span()[1]:]
        if f['fileType'] == 'SEARCH' and filetype == 'zip':
            zipfiles.append(f['downloadLink'])
        if f['fileType'] == 'RAW' and filetype == 'raw':
            rawfiles.append(f['downloadLink'])

    return zipfiles, rawfiles


def zipfile_downloader(zipfile, path, maxquant_file):
    # Handle spaces in urls
    zipfile = zipfile.replace(' ', '%20')
    zipfilename = zipfile[63:]

    # Download zip file
    try:
        os.remove(path + zipfilename)
    except:
        pass
    os.system('wget -q --show-progress -O ' + path + zipfilename + ' ' + zipfile)

    # Get a list of files with directories from zip file
    with ZipFile(path + zipfilename, 'r') as zipped:
        ziplist = zipped.namelist()

    # Extract the peptide file from the zipfile
    for a in ziplist:
        if maxquant_file in a:
            with ZipFile(path + zipfilename) as z:
                with z.open(a) as zf, open(path + zipfilename[:-4] + '-' + maxquant_file, 'wb') as zfg:
                    shutil.copyfileobj(zf, zfg)
                    break
        else:
            continue

    # Go through the maxquant output file and get all the raw files
    df = pd.read_csv(path + zipfilename[:-4] + '-' + maxquant_file, sep='\t', low_memory=False)
    df = df.loc[df['Sequence'] != ' ',]  # Remove empty sequences
    rawfiles = np.unique(df['Raw file'])

    return rawfiles, df, zipfilename


def filehandling(accnr, filename, zipfilename, path, maxquant_file, df, rawfiles):
    accessionpath = path + accnr + '/'
    filepath = accessionpath + filename + '/'
    # Make the file directory if it doesnt exist
    if not os.path.exists(accessionpath):
        os.mkdir(accessionpath)
    if not os.path.exists(filepath):
        os.mkdir(filepath)

    # Check if old version exists
    if os.path.exists(path + filename + '/'):
        for files in os.listdir(path + filename + '/'):
            if not os.path.exists(filepath + files):
                shutil.move(path + filename + '/' + files, filepath + files)
        shutil.rmtree(path + filename + '/')

    # Move or rm zip.file
    if os.path.exists(filepath + 'file.zip'):
        if os.path.getsize(path + zipfilename) > os.path.getsize(filepath + 'file.zip'):
            shutil.copyfile(path + zipfilename, filepath + 'file.zip')
    else:
        shutil.copyfile(path + zipfilename, filepath + 'file.zip')

    # Check if filespecific allPeptides.txt exists
    if not os.path.exists(filepath + maxquant_file):
        df2 = df.loc[df['Raw file'] == filename,]
        pd.DataFrame.to_csv(df2, filepath + maxquant_file)
    else:
        df2 = pd.read_csv(filepath + maxquant_file)

    # Download the raw file
    print('Downloading raw file                                                    ', end='\r')
    if not (os.path.exists(filepath + 'file.mzML') or os.path.exists(filepath + 'mzML.json')):
        if os.path.exists(filepath + 'file.raw'):
            if os.path.getsize(filepath + 'file.raw') == 0:  # If this is an empty file with nothing in it, remove it
                # (causes problems with download)
                os.remove(filepath + 'file.raw')
        for f in rawfiles:
            if filename in f or len(rawfiles) == 1:
                os.system('wget -q --show-progress -O ' + filepath + '/file.raw -c ' + f)
                break

    return df2, filepath


def formatFile(accnr, filename, path, filepath):
    print('Formatting file to mzML										', end='\r')
    # Check whether the docker file is implemented or not
    if not (os.path.exists(filepath + 'file.mzML') or os.path.exists(filepath + 'mzML.json')):
        dockerls = subprocess.check_output('docker image ls', shell=True)
        if not 'thermorawparser' in str(dockerls):
            try:
                os.system(
                    'cd .. && git clone https://github.com/compomics/ThermoRawFileParser.git && cd MassSpecPipeline/')
            except:
                pass
            os.system('cd .. && cd ThermoRawFileParser/ && docker build --no-cache -t thermorawparser . && cd '
                      '../MassSpecPipeline/')

        if path[0] == '/':
            relpath = path[:-1]
        else:
            relpath = os.getcwd() + path[:-1]
        os.system('chmod -R a+rwx ' + path + '*')
        os.system('docker run -v \"' + relpath + ':/data_input\" -i -t thermorawparser mono '
                                                 'bin/x64/Debug/ThermoRawFileParser.exe -i=/data_input/' + accnr +
                  '/' + filename + '/file.raw -o=/data_input/' + accnr + '/' + filename + '/ -f=1 -m=1')
        os.remove(filepath + 'file-metadata.txt')
        os.remove(filepath + 'file.raw')


def process_ms1(spectrum):
    # Scan information
    scan_info = spectrum['scanList']
    # Time
    scan_time = scan_info['scan'][0]['scan start time']
    mz = spectrum['m/z array']
    # ion intensity
    intensity = spectrum['intensity array']
    return {'scan_time': scan_time, 'intensity': intensity.tolist(), 'mz': mz.tolist()}


def internalmzML(path):
    # Extract the data from the mzml, if we havnt already
    if not os.path.exists(path + 'mzML.json'):
        print('Extracting data from mzML                                                    ', end='\r')
        data = mzml.MzML(path + 'file.mzML')

        # Extracted data
        extracted = {'ms1': {}}
        i = 0
        # Extract the necessary data from spectra
        for spectrum in data:

            if spectrum['ms level'] != 1:
                continue
            # Scan id
            scan_id = int(spectrum['id'].split('scan=')[1])

            # Deal with ms level 1 spectra
            ms1_spectrum = process_ms1(spectrum)
            extracted['ms1'][scan_id] = {'mz': ms1_spectrum['mz'], 'intensity': ms1_spectrum['intensity'],
                                         'scan_time': ms1_spectrum['scan_time']}

        f = open(path + 'mzML.json', 'w')
        f.write(json.dumps(extracted))
        f.close()
        os.remove(path + 'file.mzML')


def preparameters(filepath):
    print('Preparing parameter for image creation                                                    ', end='\r')
    mzml = json.load(open(filepath + '/mzML.json'))

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
    # mz_bin = (float(interval['mz']['max']) - float(interval['mz']['min']))/resolution['x']
    mz_bin = 0.75
    # rt_bin = (float(interval['rt']['max']) - float(interval['rt']['min']))/resolution['y']
    rt_bin = 0.06
    resolution = {'x': int((max(mzlist) - min(mzlist)) / mz_bin), 'y': int((max(rtlist) - min(rtlist)) / rt_bin)}

    return mzml, [mzlist, rtlist, intlist], [lowbound, highbound], interval, [mz_bin, rt_bin], resolution


def fullpng(image, filepath, resolution, interval, lowbound, highbound):
    listnames = ['Mean', 'Min', 'Max', 'Collapsed']
    for i in range(4):
        fullimage = [[y[i] for y in x] for x in image]
        fullimage.reverse()
        titleofplot = listnames[i]

        if not os.path.exists(
                filepath + str(resolution['x']) + 'x' + str(resolution['y']) + '-' + titleofplot + '.png'):
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
            plt.savefig(filepath + str(resolution['x']) + 'x' + str(resolution['y']) + '-' + titleofplot + '.png')
            plt.close()


def fullimg(mzmlfile, interval, bins, resolution, filepath, bounds, savepng):
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
    nonzero_counter = 0  # How many pixels have non-zero values
    total_datapoints = 0  # How many datapoints does the file contain.
    image = []
    for y_i in range(0, resolution['y']):
        if y_i % 25 == 0:
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
                total_datapoints += (len(ms1_array[_key]))
                nonzero_counter += 1
            except KeyError:
                pixelpoint = [0, 0, 0, 0]
            row.append(pixelpoint)
        image.append(row)
    print('Saving image files                                                    ', end='\r')
    # image.reverse()
    imagedata = [image, nonzero_counter, total_datapoints]
    # Save as txt file
    with open(filepath + str(resolution['x']) + 'x' + str(resolution['y']) + '.txt', "wb") as pa:
        pickle.dump(imagedata, pa)
    lowbound = bounds[0]
    highbound = bounds[1]

    if savepng:  # save full image to png
        fullpng(image, filepath, resolution, interval, lowbound, highbound)

    return imagedata


def subpng(subimage, imgpath, filename, index, lowbound, highbound):
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
    plt.savefig(imgpath + filename + '-' + str(index + 1) + '.png')
    plt.close()


def subimgs(interval, bins, resolution, path, df, subimage_interval, filename, image, bounds, savepng):
    mz_bin = bins[0]
    rt_bin = bins[1]
    mzrangelist = [interval['mz']['min'] + i * mz_bin for i in range(int(resolution['x']))]
    rtrangelist = [interval['rt']['min'] + i * rt_bin for i in range(int(resolution['y']))]

    imgpath = path + 'images/'
    if not os.path.exists(imgpath):
        os.mkdir(imgpath)

    metapath = path + 'metadata/'
    if not os.path.exists(metapath):
        os.mkdir(metapath)

    lowbound = bounds[0]
    highbound = bounds[1]

    outfile = open(metapath + 'subimage.json', 'a')  # The metadata file

    outbound = 0
    inbound = 0
    inmzbound = 0
    df.reset_index(drop=True, inplace=True)
    for index, rows in df.iterrows():
        if (index + 1) % int(df.shape[0] / 40) == 0:
            print('Creating subimages: {:2.1%}                                                    '.format(
                (index + 1) / df.shape[0]), end='\r')  # Print how far we are

        if rows['Retention time'] - subimage_interval['rt'] < interval['rt']['min'] or rows['Retention time'] + \
                subimage_interval['rt'] > interval['rt']['max'] or rows['m/z'] - subimage_interval['mz'] < \
                interval['mz']['min'] or rows['m/z'] + subimage_interval['mz'] > interval['mz']['max']:
            outbound += 1  # Check if this image can be created in our range or not
            continue
        inbound += 1
        if os.path.exists(imgpath + filename + '-' + str(index + 1) + '.png'):
            continue

        if not 450 < rows['m/z'] < 455:  # Filter
            continue
        inmzbound += 1

        mzlen = int(subimage_interval['mz'] / mz_bin)
        rtlen = int(subimage_interval['rt'] / rt_bin)

        mzlower = int(get_lower_bound(mzrangelist, rows['m/z']) - mzlen)
        mzupper = int(get_lower_bound(mzrangelist, rows['m/z']) + mzlen)
        rtlower = int(get_lower_bound(rtrangelist, rows['Retention time']) - rtlen)
        rtupper = int(get_lower_bound(rtrangelist, rows['Retention time']) + rtlen)

        subimage = [lines[mzlower:mzupper] for lines in image[rtlower:rtupper]]
        subimage2 = np.array(subimage)

        # Save image as txt file
        with open(imgpath + filename + '-' + str(index + 1) + '.txt', 'wb') as imagefile:
            pickle.dump(subimage, imagefile)

        if savepng:  # save subimages to png
            subpng(subimage, imgpath, filename, index, lowbound, highbound)

        new_metadata = {}
        new_metadata['image'] = filename + '-' + str(index + 1)
        new_metadata['accession'] = accession
        new_metadata['size'] = subimage2.shape
        for ele in df.columns:
            if str(rows[ele]) == 'nan' or str(rows[ele]) == ' ' or ";" in str(rows[ele]):
                continue
            else:
                new_metadata[str(ele)] = str(rows[ele])
        outfile.write(json.dumps(new_metadata) + '\n')
    outfile.close()

    return [inbound, outbound, inmzbound], metapath


def endstats(inputlists, interval, accnr, filename, total_datapoints, nonzero_counter, inorout, mpath):
    print('Calculating end statistics:                                                    ', end='\r')
    mzlist = inputlists[0]
    rtlist = inputlists[1]

    mzlist_inrange = [i for i in mzlist if interval['mz']['min'] < i < interval['mz']['max']]
    rtlist_inrange = [i for i in rtlist if interval['rt']['min'] < i < interval['rt']['max']]

    inbound = inorout[0]
    outbound = inorout[1]
    inmzbound = inorout[2]

    end_stats = {}
    end_stats['accession'] = accnr
    end_stats['filename'] = filename
    end_stats['unique mz'] = len(mzlist_inrange)
    end_stats['unique rt'] = len(rtlist_inrange)
    end_stats['datapoints'] = total_datapoints
    end_stats['data per pixel'] = total_datapoints / nonzero_counter
    end_stats['In bounds'] = inbound
    end_stats['Out of bounds'] = outbound
    end_stats['in mz range'] = inmzbound

    outfile = open(mpath + 'sub_statistics.json', 'a')
    outfile.write(json.dumps(end_stats) + '\n')
    outfile.close()
    print('Done!                                                    ')


def combined(accnr, maxquant_file, path):
    # Find all zip files
    output = filefinder(accnr)
    allZip = output[0]
    allRaw = output[1]

    for zips in reversed(allZip):
        # finds raw files for this zip file
        output = zipfile_downloader(zips, path, maxquant_file)
        rawfiles = output[0]
        df = output[1]
        zipfilename = output[2]

        for raws in rawfiles:
            filename = str(raws)

            # Skip this special case. Something wrong... dont know, dont care
            not_working = ['01625b_GA1-TUM_first_pool_1_01_01-2xIT_2xHCD-1h-R1', '01790a_BE1-TUM_second_pool_71_01_01-3xHCD-1h-R1', '01709a_GD2-TUM_first_pool_110_01_01-2xIT_2xHCD-1h-R1'']
            if filename in not_working:
                continue

            print('\nfile: ' + accnr + '/' + filename)
            output = filehandling(accnr, filename, zipfilename, path, pepfile, df, allRaw)
            df2 = output[0]
            filepath = output[1]

            formatFile(accnr, filename, path, filepath)
            internalmzML(filepath)

            output = preparameters(filepath)
            mzml = output[0]
            inputlists = output[1]
            bounds = output[2]
            interval = output[3]
            bins = output[4]
            resolution = output[5]

            # Make the image
            if not os.path.exists(filepath + str(resolution['x']) + 'x' + str(resolution['y']) + '.txt'):
                output = fullimg(mzml, interval, bins, resolution, filepath, bounds, savepng=False)
                image = output[0]
                nonzero_counter = output[1]
                total_datapoints = output[2]
            # Retrieve if exist already
            else:
                print('Loading image data                                                    ', end='\r')
                with open(filepath + str(resolution['x']) + 'x' + str(resolution['y']) + '.txt', "rb") as pa:
                    output = pickle.load(pa)
                image = output[0]
                nonzero_counter = output[1]
                total_datapoints = output[2]

            subimage_interval = {'mz': 25, 'rt': 5}
            output = subimgs(interval, bins, resolution, path, df2, subimage_interval, filename, image, bounds,
                             savepng=False)
            inorout = output[0]
            metapath = output[1]

            endstats(inputlists, interval, accnr, filename, total_datapoints, nonzero_counter, inorout, metapath)

        os.remove(datapath + zipfilename)
        os.remove(datapath + zipfilename[:-4] + '-' + pepfile)


if __name__ == '__main__':
    # Read datapath from config file
    with open('config.json') as json_file:
        data = json.load(json_file)

    datapath = data['path']
    metapath = datapath + 'metadata/'

    try:
        os.system('rm ' + datapath + '*.*')
    except:
        pass

    # Assigning accession number and maxquant output file name
    pepfile = 'allPeptides.txt'
    sysinput = sys.argv[1]
    if str(sysinput) == 'reset':
        try:
            shutil.rmtree(datapath + 'images')
            os.remove(metapath + 'subimage.json')
            os.remove(metapath + 'subimage_filtered.json')
            os.remove(metapath + 'sub_statistics.json')
        except:
            pass
    elif str(sysinput) == 'owned':
        listofowned = [f for f in os.listdir(datapath) if os.path.isdir(datapath+f) and f[0:3] == 'PRD' or f[0:3] == 'PXD']
        for accession in listofowned:
            try:
                combined(str(accession), pepfile, datapath)
            except KeyboardInterrupt:
                print('Problem occured with: ' + accession + '. unable to proceed at this time')
                pass
    elif str(sysinput) == 'accessions' or str(sysinput) == 'accessions_filtered':
        for line in reversed(list(open(metapath + sys.argv[1] + '.json'))):
            data = json.loads(line)
            accession = data['accession']
            try:
                combined(str(accession), pepfile, datapath)
            except KeyboardInterrupt:
                print('Problem occured with: ' + accession + '. unable to proceed at this time')
                pass
    else:
        accession = sysinput
        combined(str(accession), pepfile, datapath)


# python3 pdownload.py PXD004732
# python3 pdownload.py PXD010595
# python3 pdownload.py accessions_filtered
