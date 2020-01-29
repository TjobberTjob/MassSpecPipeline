## Mass Spectrometry Data Pipeline.

This pipeline is designed to facilitate downloading datasets available in PRIDE for re-analyses of MS-based proteomics datasets.

Using this pipeline you can:

* Obtain the list of currently held PRIDE accessions.



##1. Identifying datasets - PRIDE metadata mining

In order to identify suitable datasets for analysis it is necessary to analyse the 

#1.1 Obtaining a list of currently held PRIDE accessions.

Type the following command from the project root directory:

`python3 scraper.py accessions`

The full list of accessions is downloaded into {project_root}/Data/metadata/accessions.txt and stored as a pickle file that you can open by:

`pickle.load(open('Data/metadata/accessions.txt','rb'))`
