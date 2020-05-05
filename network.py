import json
import os
import random
import sys
import matplotlib.pyplot as plt
import numpy as np
from keras.engine.saving import load_model
from simplejson import loads
from classes.datagenerator import DataGenerator
from classes.model import Network_Model


def create_network_file(length_ms2, filetouse):
    print('Creating MS2 data structure')
    outfile = open(f'{metapath}subimage_filtered_network.json', 'w')

    if length_ms2 == 'max':
        length_ms2 = min([int(brokenlines.split('[')[1].split(',')[0]) for line in open(f'{metapath}{filetouse}') for
                          brokenlines in line.split(', "') if 'ms2size' in brokenlines and 'both' in line])

    for line in open(f'{metapath}{filetouse}', 'r'):
        linesplit = line.split(', "')
        for parts in linesplit:
            if 'ms2size' in parts and 'both' in line:
                size = int(parts.split('[')[1].split(',')[0])
        if 'both' in line and size > length_ms2:
            data = loads(line)
            outfile.write(json.dumps(data) + '\n')
    outfile.close()

    return length_ms2


def data_fetcher(path, filetouse, imageclass, test_accessions):
    print('Data fetching')

    for lines in open(f'{path}{filetouse}', 'r'):
        data = loads(lines)
        ms1size = eval(data['ms1size'])
        ms2size = eval(data['ms2size'])
        break

    accs = [acc.split(', "')[1][-10:-1] for acc in open(f'{path}{filetouse}', 'r')]
    accs = np.unique(accs)
    random.shuffle(accs)
    test_accs = accs[0:test_accessions]

    testfiles = []
    trainvalfiles = []
    for acc in open(f'{path}{filetouse}'):
        accession = acc.split(', "')[1][-10:-1]
        name = acc.split(', "')[0][11:-1]
        if accession in test_accs:
            testfiles.append(name)
        else:
            trainvalfiles.append(name)

    with open('config.json') as json_file:
        config = json.load(json_file)['networkattributes']
    splitratio = config['training_percentage'] / 100
    random.shuffle(trainvalfiles)
    splits = round(len(trainvalfiles) * float(splitratio))
    trainfiles = trainvalfiles[:splits]
    validationfiles = trainvalfiles[splits:]

    partition = {'train': trainfiles, 'validation': validationfiles, 'test': testfiles}

    labels = {}
    for line in open(f'{path}{filetouse}'):
        data = loads(line)
        name = data['image']
        label = data[imageclass]

        labels[name] = label

    for f in partition:
        print(f'Datapoint in {f}: {len(partition[f])}')

    return partition, labels, ms1size, ms2size


def history_plot(metric, metapath, imageclass):
    plt.plot(history.history[f'{metric}'])
    plt.plot(history.history[f'val_{metric}'])
    plt.title(f'Plot of {metric} over time')
    plt.ylabel(f'{metric}')
    plt.xlabel('Epoch')
    plt.legend(['Training', 'Validation'], loc='upper left')
    plt.savefig(f'{metapath}{imageclass.replace("/", "")}.png')


def nnmodel(ms1size, ms2size, n_channels, length_ms2, classification, n_classes, imageclass, metapath, patience,
            which_ms_touse):
    model_network = Network_Model(which_ms_touse, classification, n_classes, ms1size, ms2size, n_channels, length_ms2,
                                  metapath,
                                  imageclass, patience)
    model = model_network.get_network()
    callbacks_list = model_network.get_callbacks()
    # print(model.summary())
    # plot_model(model, to_file="model.png")

    return model, callbacks_list


if __name__ == '__main__':
    # Read datapath from config file
    with open('config.json') as json_file:
        data = json.load(json_file)

    datapath = data['path']
    metapath = f'{datapath}metadata/'
    imagepath = f'{datapath}images/'

    with open('config.json') as json_file:
        config = json.load(json_file)['networkattributes']
    test_accessions = config['test_accessions']
    n_channels = config['n_channels']
    batch_size = config['batch_size']
    epochs = config['epochs']
    patience = config['early_stopping']
    setseed = config['set_seed'] == 'True'
    which_ms_touse = config['which_ms_touse']
    length_ms2 = config['length_ms2']

    if os.path.exists(f'{metapath}subimage_filtered.json'):
        filetouse = 'subimage_filtered.json'
    else:
        print('No filtered datafile exists')
        quit()

    # Creating network files
    if which_ms_touse == 'both' or which_ms_touse == 'ms2':
        length_ms2 = create_network_file(length_ms2, filetouse)

    if setseed:
        random.seed(1)
        os.environ['PYTHONHASHSEED'] = str(1)
        np.random.seed(1)
        # set_seed(1)

    classification = sys.argv[1].lower()
    if not (classification == 'c' or classification == 'r'):
        print('classification or regression problem not input correctly.')
        quit()
    classification = classification == 'c'

    if classification:
        imageclass = f'{sys.argv[2]}_class'
    else:
        imageclass = sys.argv[2]

    if (which_ms_touse == 'both' or which_ms_touse == 'ms2') and os.path.exists(
            f'{metapath}subimage_filtered_network.json'):
        filetouse = 'subimage_filtered_network.json'
    elif os.path.exists(f'{metapath}subimage_filtered.json'):
        filetouse = 'subimage_filtered.json'

    nameofclass = imageclass.replace('/', '')
    output = data_fetcher(metapath, filetouse, imageclass, test_accessions)
    partition = output[0]
    labels = output[1]
    ms1size = output[2]
    ms2size = output[3]

    if classification:
        classes = [f.split('"')[-2] for line in open(f'{metapath}{filetouse}', 'r') for f in line.split(', "') if
                   imageclass in f]
        n_classes = len(np.unique(classes))
    else:
        n_classes = 1

    params = {'ms1size': (ms1size[0], ms1size[1]),
              'ms2size': (ms2size[0], ms2size[1]),
              'batch_size': batch_size,
              'n_classes': n_classes,
              'n_channels': n_channels,
              'shuffle': True,
              'which_ms_touse': which_ms_touse,
              'length_ms2': length_ms2,
              'classification': classification,
              'imagepath': imagepath,
              }

    training_generator = DataGenerator(imagepath, partition['train'], labels, **params)
    validation_generator = DataGenerator(imagepath, partition['validation'], labels, **params)

    output = nnmodel(ms1size, ms2size, n_channels, length_ms2, classification, n_classes, nameofclass, metapath,
                     patience,
                     which_ms_touse)
    model = output[0]
    callbacks_list = output[1]
    history = model.fit_generator(generator=training_generator, validation_data=validation_generator, epochs=epochs,
                                  callbacks=callbacks_list)

    # PLOT HISTORY
    if classification:
        history_plot('accuracy', metapath, imageclass)
    else:
        history_plot('r2', metapath, imageclass)

    # TEST MODEL
    if len(partition['test']) > 0:
        model = load_model(f'{metapath}Best-{imageclass}.h5')
        test_generator = DataGenerator(imagepath, partition['test'], labels, **params)
        testaccuracy = model.evaluate_generator(test_generator)
        if classification:
            print(f'Accuracy on test data. Loss: {testaccuracy[0]}. Accuracy: {testaccuracy[1]}')
        else:
            print(f'Accuracy on test data. Loss: {testaccuracy[0]}. R^2: {testaccuracy[1]}')

# python3 network.py r m/z
# python3 network.py c Length_class
# python3 network.py c Seq_class
# python3 network.py c Modi_class
