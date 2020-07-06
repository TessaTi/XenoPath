# Xenopath tool v. 0.1

Xenopath is a tool annotating the xenobiotic degradation potential of microrganism from metagenomic data.  


# Installation 
Install the software requirements listed in: requirements.txt
Ideally you would have anaconda/miniconda installed such that you can create a new environment with all the requirement sowtware 
by: 

conda env create --file xenopath_requirements.txt



# Getting started: Download and Create the reference databases

Xenopath tool relies on Swissprot database and annotation data for reads-Open Reading Frames annotation, and integrate Kaiju software for taxonomic annotation. 

### 1. DB for ORFs annotation
Download XenoPathDB: 


example usage: 
python xenopath_makedb.py --database swissprot --output ./swissprot_db --threads 5

### 2. DB for taxonomic annnotation (Kaiju: https://github.com/bioinformatics-centre/kaiju)


Can donwload the taxonomic database that best suits your data. 

kaiju-makedb -s <DB> # see the link above for information
    
example: 
    kaiju-makedb -s refseq 
    


### Configuration file and sample file
XenoPath requires a configuration file (example:) and a metadata file with the sample information (example: ): 
The configuration file accept the following parameters: 
[kaiju_setting]
node = /home/user/node.dmp
names = /home/user/names.dmp
fmi =  /home/user/kaiju_db.fmi
a = None 
s = None
E = None
p = None 

[diamond_setting]
db_diamond = /home/user/swissprot.dmd 
mode = sensitive

[general_setting]
threads = 1

The metadata file accepts the following column names: 'samples','file','typereads','groups'.
Sample: the sample name
File: the path to the file
Typereads: paired, unpaired, trimmomatic (paired and unpaired reads)
Groups: control/case/non-specified
Example: 
    'samples','file','typereads','groups'
           S1  file/location/ single 
