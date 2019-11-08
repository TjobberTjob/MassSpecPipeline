Usage of pipeline for MaxQuant datasets ONLY:
python proteomeTools.py all-PT -> Downloads all the data from proteomeTools datasets

python proteomeTools.py https://www.ebi.ac.uk/pride/data/archive/2019/05/PXD010595/01974c_BC1-TUM_missing_first_3_01_01-ETD-1h-R4 -> Downloads the data from this one link. Do not specify filetype (e.g .raw at the end)

python proteomeTools.py all /pride/data/archive/2019/05/PXD010595 -> Downloads everything from this pride archive

After this is done move to classes.py
python classes.py Sequence -> Moves all images into folders classified by their sequence
After this youll be asked if you want to classify them into train and validation dataset, answer as needed. 
If you answer "yes" youll be asked to give the validation %. train will be 1-(val%) 

If there is already folders existing with classifications, youll be prompted to reset folders. Not doing this will stop the code