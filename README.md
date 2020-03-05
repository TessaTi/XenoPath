# Xenopath tool v. 0.1

Xenopath is a tool annotating  the xenobiotic degradation potential of microrganism from metagenomic data.  


# Installation 







# Download and Create the reference database

Xenopath tool rely on Swissprot database and annotation data. The following script will directly downlod the newest version of swissprot and create the database index. 

xenopath_makedb.py -h
usage: xenopath_makedb.py [-h] --database DATABASE --output OUTPUT
                          [--threads THREADS]


example usage: 
python xenopath_makedb.py --database swissprot --output ./swissprot_db --threads 5

