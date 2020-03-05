# Xenopath tool v. 0.1

Xenopath is a tool annotating  the xenobiotic degradation potential of microrganism from metagenomic data.  


# Installation 






# Download and Create the reference databases

Xenopath tool rely on Swissprot database and annotation data for reads-Open Reading Frames annotation, and integrate Kaiju software for taxonomic annotation. 

### 1. DB for ORFs annotation

The following script will directly downlod the newest version of swissprot and create the database index.

python xenopath_makedb.py -h
usage: xenopath_makedb.py [-h] --database DATABASE --output OUTPUT
                          [--threads THREADS]


example usage: 
python xenopath_makedb.py --database swissprot --output ./swissprot_db --threads 5

### 2. DB for taxonomic annnotation (Kaiju: https://github.com/bioinformatics-centre/kaiju)


Can donwload the taxonomic database that best suits your data. 

kaiju-makedb -s <DB> # see the link above for information
    
example: 
    kaiju-makedb -s refseq 
    


### 