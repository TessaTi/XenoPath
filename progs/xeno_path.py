#!/usr/bin/env python
#################
# XENOPATH
#################
import time
from datetime import datetime,date
import argparse
import subprocess
import pandas as pd
import sys, os,csv

try:
    import ConfigParser as configparser
except ImportError:
    import configparser


from taxa_run import pairedfun,singlefun,trimmomaticfun

def readconfigfile(configfile):
    config = configparser.ConfigParser()
    config.read(configfile)
   
    #possible arguments from config file
    #node,names,fmi,a,s,E,x,X,p='','','','','','','','',''
    #db_diamond,mode = '',''
    settings_kj = ['node','names','fmi','a','s','E','p']
    settings_dm = ['mode','db_diamond']
    #general_setting = ['threads','folderout','typereads','annotab','XenoPathfolder']
    #general_setting = ['threads','folderout','typereads','database']
    general_setting = ['threads','folderout','typereads','annotab','XenoPathfolder']
    minpath_setting = ['minpathec','minpathpy']
    #minpath_setting = ['minpathec']
    kaiju_setting,diamond_setting = '',''
    l, ld, g = {},{},{}
    for k,v in config.items():
        if k == 'kaiju_setting':
            for setting in range(0,len(settings_kj)):
                var = config.get(k,settings_kj[setting])
                if var == None:
                    pass
                else:
                    l[settings_kj[setting]] = var
                #print(k,v,i, config.get(k,i))
            kaiju_setting = (str(' '.join([ ('-'+k+' '+v) if v != 'None' else '' for k,v in l.items() ])))
        elif k == 'diamond_setting':
            for setting in range(0,len(settings_dm)):
                var = config.get(k,settings_dm[setting])
                if var == None:
                    pass
                else:
                    ld[settings_dm[setting]] = var
            diamond_setting = (str(' '.join([ ('-'+k+' '+v) if v != 'None' else '' for k,v in ld.items() ])))
        elif k == 'general_setting':
            for setting in range(0,len(general_setting)):
                var = config.get(k,general_setting[setting])
                if var == None:
                    pass
                else:
                    g[general_setting[setting]] = var
        elif k == 'minpath_setting':
            for setting in range(0,len(minpath_setting)):
                var = config.get(k,minpath_setting[setting])
                if var == None:
                    pass
                else:
                    g[minpath_setting[setting]] = var

    ###set error if readtype is not specified
    try:
        if g['typereads'] != None and g['typereads'] in ['paired', 'unpaired','trimomatic']:
            print('Reading configuration...')
    except ValueError:
        print("wrong configuration! Check typereads information")
    
    return(l,ld,g) #return dict

def xenotool(configfile,samplefile):

    ##########################################
    # KAIJU, DIAMOND AND GENERAL CONFIGURATION
    ##########################################
    kaiju_config, diamond_config, general_config = readconfigfile(configfile)
    print('kaiju congif', kaiju_config)
    print('diamond congif', diamond_config)

    #threads,minpath_map,minpath_folder = general_config['threads'],general_config['minpathec'],general_config['minpathpy']
    minpathec = general_config['minpathec'] #new
    database = diamond_config['db_diamond']
    tabfile = general_config['annotab']
    threads  = general_config['threads']
    folderout = general_config['folderout']
    
    typereads = general_config['typereads']
    print('general config', general_config)
    print('threads', threads)
    print('typereads', typereads)
    ########## ########## ########## ##########
    ########## OPEN METADATA file as pandas
    # open file with forward and reverse reads
    F = pd.read_csv(samplefile, header = 0)
    columns = ['samples','file','typereads','pair','groups']
    try:
        if list(F.columns) in columns:
            None
    except ValueError:
        print(list(F.columns),columns)
        print("wrong metadata fields! Try again...")
        sys.exit(1) 

    samples = list(set(F['samples'].values.tolist())) # list of samples
    groups = F['groups'].values.tolist()   # group for each sample
    print('table',F.head)
    print(samples)
    print(groups)

    ##create output directory in the current folder if not specified
    pathout = ''
    print('folderout',folderout)
    if (folderout == 'None') or (folderout =="''"):
        timestr = time.strftime("%Y%m%d-%H%M%S")
        cwd = os.getcwd()
        pathout = os.path.join(cwd,'xenopath'+'_'+timestr)
        print('create folder', pathout)

        try:
            if not os.path.exists(pathout):
                print('creating output directory: ', pathout)
                os.mkdir(pathout)
        except ValueError:
            print("Folder already exist! Check folder...")
            sys.exit()

    else:
        print('pathout already exist', pathout, folderout)
        pathout = folderout

    # GET READ TYPE INFORMATION

    #tabfile,scriptfolder = general_config['annotab'], general_config['XenoPathfolder']
    scriptfolder = general_config['XenoPathfolder']
    try:
        if typereads != None:
            print('typereads check', typereads)
            if typereads == "'paired'":
                #pairedfun(samples,F,kaiju_config,diamond_config,threads,tabfile,scriptfolder,general_config,pathout)
                pairedfun(samples,F,kaiju_config,diamond_config,threads,tabfile,scriptfolder,general_config,pathout,database) #new tabfile 
            elif typereads == "'single'":
                singlefun(samples,F,kaiju_config,diamond_config,threads,tabfile,scriptfolder,general_config,pathout,database)
                #singlefun(samples,F,kaiju_config,diamond_config,threads,tabfile,scriptfolder,general_config,pathout)
                
            elif typereads == "'trimmomatic'":
                trimmomaticfun(samples,F,kaiju_config,diamond_config,threads,tabfile,scriptfolder,general_config,pathout,database)
                #trimmomaticfun(samples,F,kaiju_config,diamond_config,threads,tabfile,scriptfolder,general_config,pathout)
                
    except ValueError:
        print("typereads information missing! Try again...")
        sys.exit(1) 

    print('end',time.time())

def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x

def Main():
    parser = argparse.ArgumentParser(description='python XenoPath.py --config_file config.cfg \
    --metadata metadata.csv ')
    parser.add_argument('--config_file',required=True, type=extant_file,help='1. configuration setting file')
    parser.add_argument('--metadata',required=True, type=extant_file,help='2. samples fastq location, read type and group info')

    args = parser.parse_args()

    xenotool(args.config_file,args.metadata)

if __name__ == '__main__':
    begin_time = datetime.now()
    print(begin_time)
    Main()
    print(datetime.now())
    print(datetime.now() - begin_time)
