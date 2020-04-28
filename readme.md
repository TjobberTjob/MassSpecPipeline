# Mass Spectrometry Pipeline.

Code aimed at facilitating downloading mass spectrometry datasets for repurposing.

#### Disclaimer! Currently only works for MS files analyzed by MaxQuant. Furthermore it only works on Ubuntu / Linux / Bash and python version 3.5 or newer
##### git clone https://github.com/TjobberTjob/MassSpecPipeline

## Requirements.

This repository was coded and tested using python 3.6.8.

Required python packages:

`pip install pyteomics`

`pip install matplotlib`

`pip install keras`

`pip install pydot`

`pip install pandas`

`pip install numpy`

`pip install simplejson`

`pip install lxml`

`tensorflow (CPU or GPU)`

Optional conda installations:

`conda install -c bioconda thermorawfileparser`                              

`conda install -y tensorflow-gpu=2.1.0`                                 

`conda install -y pytorch=1.4.0 torchvision cudatoolkit=10.1 -c pytorch`

`conda install -c anaconda graphviz`

Optional docker installations:
`git clone https://github.com/compomics/ThermoRawFileParser.git`       

`docker build --no-cache -t thermorawparser /ThermoRawFileParser`         


### config file
#### extraction attributes:
TOOL is a connector between MS and AI, downloading, handling, formatting and extracting data from raw MS files into singular image files with metadata.
Before we go into examples of usage, first step is to change config.json file

The only MUST CHANGE is the path. Change this to the desired path. Everything else is explained here but not necessary to change.

path: Where all data gets put to - Need for large storage (end path with "/")

mz_bin: Size of m/z bins used for large image - Defaults to 0.75

rt_bin: Size of retention time bins used for large image - Defaults to 0.06

mz_interval: m/z interval used for subimages around the peptides - Defaults to +-25

rt_interval: rt interval used for subimages around the peptides - Defaults to +-5

acquire_only_new: A boolean value for skipping PRIDE projects if you have already created all images - Defaults to "True"

multithread: A Boolean value for running the extraction multithreaded (reduces error messages)- Defaults to "True"

nr_threads: Amount of threads used for multi threading - Defaults to 10

formatsoftware: What software is used for formatting from .raw to .mzml - Possibilies are "conda" or "docker" so the used method must be installed - Defaults to "conda"

Errormessage: If files fail, whether to just skip them or skip them and get the error message. Used for debugging - Defaults to "False"

savepng: Whether or not to save images as png files or only txt - Dafults to "False"

filterscore: When filtering subimage metadata, do you want to have a score amount that the peptide has to achieve, and what is the percentile of all score you want to have above. - Defaults to ["False", 55],

filteramount: When filtering subimage metadata, do you want to have an amount of datapoints pr class. The filter filters off percentages of total datapoints in either class (so 10(%) means that each class has to have 10% of all datapoints or it wont be incorporated) - Defaults to ["False", 50], 

topxamount: When you want only x classes. This filteres away the classes that arent in the most abundant classes. (fx only 10 most abundant sequences) - Defaults to "max" that is no top limit, else write a numerical value

#### networkattributes:
batch_size: Size of batches for stochastic gradient descent - Defaults to 64

epochs: Max epochs used for neural network training - Defaults to 100

n_channels: Number of information channels used in neural network training - Defaults to 4

test_accessions: Number of accession projects set aside to testing (not validation) - Defaults to 3

training_percentage: The amount of data going to training vs validating - Defaults to 80 (80% train - 20% validation)

early_stopping: Number of epochs with no improvement before stopping - Defaults to 10

setseed: Makes sure that it always initializes files in the same training and test groupings. Defaults to "True"

MS: Whether to use just ms1, ms2 or ms1 and ms2 data for the neural networks - Defaults to "ms1" ("both" and "ms2" are the other option)

lenms2: The length of the MS2 data used. This is needed because of the size of the neural network for MS2 data. 
 Has two options, a numerical value where all data points with less data are filtered away, or "max" where it takes the
 n highest values from all ms2 data, where n is the shortest ms2 data in a single datapoint. Only used if MS is "both" - Defualts to "max"

#### NOTE: Values given as "True" and "False" needs to be in quotes.

### Usage:
#### Part - PRIDE Webscraper:
Starting of we present a method to getting information on all PRIDE projects.
To begin with we need the accessions numbers for all PRIDE projects. This is gotten by using the scraper.py using the input "accessions".
This will create, if it doesnt already exists, the metadata folder and create a txt file with all accession numbers from pride
Example:
python scraper.py accessions

Next we need to get metadata for all accessions in the accession list.
This is also done using scraper.py using the input "metadata".
This is a lengthy process as it needs to downloads a lot of files from PRIDE database.
For this reason we have released the newest version as of April 2020 on the website.
Example:
python scraper metadata

After running or downloading the metadata, you can utilize the update input to scraper.py aswell.
This will update the metadata file with all new pride accession numbers.
Example:
python scraper update


####Part - Metadata handling:
This is an extra part we added to help with metadata handling.
It can create a filtered version of either PRIDE metadata or subimage metadata. Two possible use cases.

(1): PRIDE metadata filter.
This will add a filtered version of PRIDE metadata file, which removes all non maxquant entries, and entries without raw data.
It also removes any duplicates, if any occoured by error.

Example:

python filter.py accessions

(2): Subimage metadata filter.

This will add a filtered version of the subimage metadata file created by the extractor.py script.
This is used as a precurser to neural networks or as extra data or class creation.
There are a few defaults added at the time of publishing: (all of these sort away the data with wrong filesize, which in a necessary step)

PTM - Gives a class called Modi_class with 0 for no modification and 1 for any modification

Length - Gives a class called Length_class

Sequence - Sort only the x most abundant sequences from the data and creates Sequence_class

Charge - Filters away only the 80 lowest percentiles of scores. meaning the data become cleaner (used when you have multi millions of subimages) 

Example:

python filter.py subimage PTM/Length/Sequence/Charge


####Part - Extracting:
This will extract everything from either local or PRIDE database raw files.
This has five use possible use cases:

(1): A single pride accession. This will download all files needed one at a time and run through untill everything is done.
If a zipfile fails it will proceed to next zipfile.

Example:

python extractor.py PXD000000

(2): Usage for all accession in a filtered or unfiltered version of the metadata. Goes through all accessions,
skips if there's a problem with some of them (and there will be)

Example:

python extractor.py pride OR python3 extractor.py pridefiltered

(3): Usage for all accession currently in your local directory that has all informations.
Used if you change a parameter for subimages and you need to rerun everything faster.
This will only run if every allPeptides.txt file exists in needed folders.

Example:

python extractor.py owned

(4): Used for local rawfiles instead of pride accessions. If you have a local rawfile you want peptides extracted from.
Has to have raw file and maxquant or zipfile in the folder.

Example:

python extractor.py /path/to/folder/


(5): Used to reset everything regarding the subimages
This is used if you need to remove all subimages, metadata on subimages, substatistics on files.
Often only needed if you change subimage parameters and want everything to stay consistent within files.
If you dont reset you will have multiple of the same name in the subimage metadata file.

Example:

python extractor.py reset


Part - Neural networks:
This will automatically use subimage_filtered.json if it exists otherwise it will use subimage.json as it's data source.
It needs to be given three arguments:

(1) Whether its a classification problem or a regression problem. C for classification and R for regression.

(2) The class or variable the network needs to predict. can be continuous like m/z or a class like length of sequence.

Example:

python network.py c Length & python network.py r m/z 



