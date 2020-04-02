# Disclaimer! Currently only works for MS files analyzed by MaxQuant.

TOOL is a connector between MS and AI, downloading, handling, formatting and extracting data from raw MS files into singular image files with metadata.
Before we go into examples of usage, first step is to change config.json file

The only MUST CHANGE is the path. Change this to the desired path. Everything else is explained here but not necessary to change.
path: Where all data gets put to - Need for large storage (end path with "/")
mz_bin: Size of m/z bins used for large image - Defaults to 0.75
rt_bin: Size of retention time bins used for large image - Defaults to 0.06
mz_interval: m/z interval used for subimages around the peptides - Defaults to +-25
rt_interval: rt interval used for subimages around the peptides - Defaults to +-5
acquire_only_new: A boolean value for skipping PRIDE projects if you have already created all images - Defaults to "True"
skip_incomplete: A boolean value for skipping PRIDE projects where all necessary files arent in local directories - Defaults to "False"
multithread: A Boolean value for running the extraction multithreaded (reduces error messages)- Defaults to "False"
nr_threads: Amount of threads used for multi threading - Defaults to 4
filterbroken: A method to skip PRIDE projects that has problems with them - Defaults to "False"
formatsoftware: What software is used for formatting from .raw to .mzml - Possibilies are "conda" or "docker" so the used method must be installed - Defaults to "conda"
### networkattributes ###
batch_size: Size of batches for stochastic gradient descent - Defaults to 124
epochs: Max epochs used for neural network training - Defaults to 100
n_channels: Number of information channels used in neural network training - Defaults to 4
test_accessions: Number of accession projects set aside to testing (not validation) - Defaults to 1
early_stopping: Number of epochs with no improvement before stopping - Defaults to 10}}
#NOTE: Values given as "True" and "False" needs to be in quotes.

Usage:
Part - PRIDE Webscraper:
Starting of we present a method to getting information on all PRIDE projects.
To begin with we need the accessions numbers for all PRIDE projects. This is gotten by using the scraper.py using the input "accessions".
This will create, if it doesnt already exists, the metadata folder and create a txt file with all accession numbers from pride
Example:
python3 scraper.py accessions

Next we need to get metadata for all accessions in the accession list.
This is also done using scraper.py using the input "metadata".
This is a lengthy process as it needs to downloads a lot of files from PRIDE database.
For this reason we have released the newest version as of April 2020 on the website.
Example:
python3 scraper metadata

After running or downloading the metadata, you can utilize the update input to scraper.py aswell.
This will update the metadata file with all new pride accession numbers.
Example:
python3 scraper update


Part - Metadata handling:
This is an extra part we added to help with metadata handling.
It can create a filtered version of either PRIDE metadata or subimage metadata. Two possible use cases.
(1): PRIDE metadata filter.
This will add a filtered version of PRIDE metadata file, which removes all non maxquant entries, and entries without raw data.
It also removes any duplicates, if any occoured by error.
Example:
python3 filehandler.py accessions

(2): Subimage metadata filter.
This will add a filtered version of the subimage metadata file created by the extractor.py script.
This is used as a precurser to neural networks or as extra data or class creation.
By default it will add a class called Modi_class that's 1 if the peptide is Oxidized, and 0 otherwise.
If other filters or classes needs to be put in, it has to be done manually in the code, as there are too many possibilities to add them all.
Example:
python3 filehandler.py subimage


Part - Extracting:
This will extract everything from either local or PRIDE database raw files.
This has five use possible use cases:
(1): A single pride accession. This will download all files needed one at a time and run through untill everything is done.
If a zipfile fails it will proceed to next zipfile.
Example:
python3 extractor.py PXD000000

(2): Usage for all accession in a filtered or unfiltered version of the metadata. Goes through all accessions,
skips if there's a problem with some of them (and there will be)
Example:
python3 extractor.py accessions OR python3 extractor.py accessions_filtered

(3): Usage for all accession currently in your local directory that has all informations.
Used if you change a parameter for subimages and you need to rerun everything faster.
This will only run if every allPeptides.txt file exists in needed folders.
Example:
python3 extractor.py owned

(4): Used for local rawfiles instead of pride accessions. If you have a local rawfile you want peptides extracted from.
Has to have raw file and maxquant or zipfile in the folder.
Example:
python3 extractor.py /path/to/folder/


(5): Used to reset everything regarding the subimages
This is used if you need to remove all subimages, metadata on subimages, substatistics on files.
Often only needed if you change subimage parameters and want everything to stay consistent within files.
If you dont reset you will have multiple of the same name in the subimage metadata file.
Example:
python3 extractor.py reset


Part - Neural networks:
This will automatically use subimage_filtered.json if it exists otherwise it will use subimage.json as it's data source.
It needs to be given three arguments:
(1) Whether its a classification problem or a regression problem. C for classification and R for regression.
(2) The class or variable the network needs to predict. can be continuous like m/z or a class like length of sequence.
(3) The percentage of training data, given as a fraction (fx 0.5 for 50%).
Example:
python3 network.py C Length 0.8 OR python3 network.py R m/z 0.75



