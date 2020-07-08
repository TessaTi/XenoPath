# XenoPath<sub>xs</sub> v. 0.1

Xenopath<sub>xs</sub> is a tool annotating the xenobiotic degradation potential of microorganism from metagenomic and metatranscriptomic data, accepting short-reads paired-end, single-end as well as the paired and unpaired output from [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).
The tool accepts cleaned reads, we provide a tutorial for the reads pre-processing steps (here).

XenoPath<sub>xs</sub> relies on other tools:
 - [Kaiju](https://github.com/bioinformatics-centre/kaiju)
 - [Diamond](https://github.com/bbuchfink/diamond)
 - [Minpath](https://omics.informatics.indiana.edu/MinPath/)




# 1. Getting started: 

### Installation 
software requirements listed in: xenopath_env.yml\
Ideally you can create a new environment in conda with all the requirements:\

conda env create -f xenopath_env.yml\
conda activate xenopath_env\

### Download and Create the reference databases

#### DB for taxonomic annnotation [Kaiju](https://github.com/bioinformatics-centre/kaiju)

Donwloading the taxonomic database that best suits your data.\
example:\
    kaiju-makedb -s refseq 
    
Download data from this repository:\
- swissprot_06-29_xeno.dmnd\
- swissprot_xenobiotic_degra_2020-06-29.csv\
- map2ec # for Minpath

# 2. Configuration file and sample file
XenoPath<sub>xs</sub> requires a configuration file (example:config.cfg) and a metadata file with the sample information (example: metadata.cfg ): 
The configuration file accept the following parameters: 

\[kaiju_setting]\
kaijuscript = /home/teresa/prova\
node = /home/user/Kaiju_data/kaijudb/nodes.dmp\
names = /home/user/Kaiju_data/kaijudb/names.dmp\
fmi =  /home/user/Kaiju_data/kaijudb/kaiju_db.fmi\
a = None\
s = None\
E = None\
p = None

\[diamond_setting]\
db_diamond = /home/user/pipeline/swissprot_06-29_xeno.dmnd\
mode = sensitive

\[general_setting]\
threads = 5\
folderout = /home/user/pipeline/xenopath_20200703-170138\
XenoPathfolder = /home/user/pipeline\
annotab = /home/user/data/swissprot_xenobiotic_degra_2020-06-29.csv\
typereads = 'paired'

\[minpath_setting]\
minpathec = /home/user/Minpath_tool/MinPath/data/map2ec\
minpathpy = /home/user/Progetti/Minpath_tool/MinPath

The metadata file accepts the following column names: 'samples','file','typereads','groups'.
    
example :\
samples,file,typereads,pair,groups\
S010,/home/user/pipeline/metacent_mock/S010_R1.fastq,single,r1,O\
Y2,/home/user/pipeline/metacent_mock/Y2_R1.fastq,single,r1,Y

Sample: the sample name
File: the path to the file
Typereads: paired, unpaired, trimmomatic (paired and unpaired reads)
Groups: control/case/non-specified
Example: 
'samples','file','typereads','groups'
S1  file/location/ single 

# 3. Running XenoPath 
    
### command line:
    python xeno_path.py --config_file config.cfg --metadata metadata.txt
    
### filter tables based on number of reads or relative abundance: 
    python filter_table.py --dataframe file.txt --samples 10 --x 10 --p 5 --out ./
    
    --samples: dataset number of samples\
    --x: minimum number of reads\
    --p: in at least --p sample\
    
### output files
The output folder contains intermediate data and the following tables summarizing the information for the dataset submitted: 
 - <strong>samples metaphenotype</strong>
1. all_ectaxa.csv #ec taxa read count 
1. all_path_abundance.csv #  pathway abundance
1. all_pathtaxa_coverage.csv # pathway coverage 
    
- <strong>full taxonomy at phylum/family/genus/species level</strong>
1. phylum_fulltaxa.tsv  
1. family_fulltaxa.tsv
1. genus_fulltaxa.tsv
1. species_fulltaxa.tsv
    
- <strong>taxonomy of xenobiotic subset</strong>
1. phylum_xenotaxa.tsv
1. family_xenotaxa.tsv       
1. genus_fulltaxa.tsv
1. species_xenotaxa.tsv
   

    
### visualization 
    
    --to be completed
    
    
