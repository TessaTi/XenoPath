#!/usr/bin/env python
import sys
import argparse
import pandas as pd
import numpy as np
from datetime import datetime,date

def reformatting_level(report,out,type):
    
    # open the report table
    table= pd.read_csv(report,header=0, delimiter='\t',keep_default_na=False)
    print(table.loc[table['taxon_name'].str.contains('NA'), 'taxon_name'])
    print(table[table.isna().any(axis=1)])
    table = table.fillna('unclassified;unclassified;unclassified;unclassified;')
    
    table = table.replace(to_replace='NA', value ='unclassified')    
    table['file'] = table['file'].str.split('/').str[-1].str.split('.out').str[0]
    table.loc[table['taxon_name'].str.contains('Viruses'),'taxon_name'] = 'Viruses;Viruses;Viruses;Viruses;' 
    table.loc[table['taxon_name'].str.contains('viral'),'taxon_name'] = 'Viruses;Viruses;Viruses;Viruses;'
    table.loc[table['taxon_name'] ==  'unclassified', 'taxon_name' ] = 'unclassified;unclassified;unclassified;unclassified;'    
    table.loc[table['taxon_name'] ==  'NA', 'taxon_name' ] = 'unclassified;unclassified;unclassified;unclassified;'
    #change column name as samples name 
    level =['phylum', 'family','genus','species']

    for l in range(0,len(level)):
        phylum = table.copy(deep=True)
        phylum['taxon_id'] = phylum['taxon_id'].astype(str).str.split('.', expand = True)[0]
        phylum['taxon_name'] = phylum['taxon_name'].astype(str)
        
        phylum['taxon_name'] = phylum['taxon_name'].str.split(';').str[l] #level table
        phylum = phylum.fillna('unclassified')
        phylum.loc[phylum['taxon_name'] ==  'NA', 'taxon_name' ] = 'unclassified'
        print(phylum[phylum.isna().any(axis=1)])
        dout = phylum.groupby(['file','taxon_name'])['reads'].agg('sum').reset_index(name ='read_count')#, as_index=False).size().reset_index(name='campione')
        # row taxa column samples
        d= dout.pivot(index='taxon_name', columns='file', values='read_count')
        d = d.fillna(0)
        
        d.to_csv(out+'/'+level[l]+'_'+type+'.tsv')
        del (phylum)
        
def Main(): 
    
    parser = argparse.ArgumentParser(description='python ...')
    
    parser.add_argument('--reportable',help='1. table')
    parser.add_argument('--out',help='2. output')
    parser.add_argument('--type',help='3. full taxa or subsetted seno taxonomy')
    
    args = parser.parse_args()
    
    if len(sys.argv) != 7:
        print(len(sys.argv))
        print("You haven't specified any arguments. Use -h to get more details on how to use this command.")
        parser.print_help()
        sys.exit()
    reformatting_level(args.reportable,args.out,args.type)

if __name__ == '__main__':
    begin_time = datetime.now()
    print(begin_time)
    Main()
    print(datetime.now())
    print(datetime.now() - begin_time)


