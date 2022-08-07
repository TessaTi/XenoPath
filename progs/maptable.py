#!/usr/bin/env python
'''
Output: 
1. file of EC taxa abundances
2. file of pathway coverage by taxa
3. file with saple pathway abundances
4. file with ec!taxa|pathway information 
'''
import sys
import argparse
import csv
import pandas as pd
import time
from datetime import datetime

#########################
# open detail file as list
#########################
def opendetail(detail):
    D = {}
    for line in detail: 
        line = line.strip()
        if line.startswith('path'):
            path = line.split('# ')[1]
            x = D.get(path,[])
            D[path] = x
        else: 
            ec = line.split(' ')[0]
            x = D.get(path,[])
            x.append(ec)
            D[path] = x
            
    lval = []
    for k,v in D.items(): 
        for i in v: 
            lval.append(i)
    return(D,lval)


def xeno_mock(s,detailfile,mock_xenobiotic,maptaxainfo,pathname,out): 

    ##add pathway name 
    Dpath = pd.read_csv(pathname,sep ='\t', header=0) #open path list with pathway description
    Dpath.columns = ['MAP','Kegg_pathway']
    Dpath['MAP'] = Dpath['MAP'].apply(str)

    
    ##open deatil and find the path
    detail = open(detailfile,'r').readlines()
    D,Dvalues = opendetail(detail)
  
    ##open the csv file read and ec
    Mock_xenobiotic = pd.read_csv(mock_xenobiotic, sep='\t')
    Mock_xenobiotic.columns = ['qseqid', 'EC_number'] #qseqid == read
    #subset to ec in minpath
    Mock_xenobiotic = Mock_xenobiotic[Mock_xenobiotic['EC_number'].isin(Dvalues)]

    ##open functional annotation and taxa
    Maptaxainfo = pd.read_csv(maptaxainfo, header=0, sep = '\t')
    
    #intersection with map and pathway name for each read
    M = pd.merge(Mock_xenobiotic,Maptaxainfo, how = 'inner', on=['qseqid','EC_number'])

    # ec taxa output file
    M['taxa_name'] = M['taxa_name'].str.split(";").str[-2].str.lstrip(' ')#.str.join(';') #take only the gene name
    #first output save
    M  = M[['EC_number', 'taxa_name']]
    M = M.replace(to_replace='NA', value ='unclassified')
    M = M.fillna('unclassified')
    md = M.copy()
    M['EC_number'] = M['EC_number'].apply(str)
    M['taxa_name'] = M['taxa_name'].apply(str)
    M['EC_number_taxa'] = M[['EC_number', 'taxa_name']].apply(lambda x: '|'.join(x), axis=1)
    M = M.drop(['EC_number', 'taxa_name'], axis=1)
    M = M.groupby(['EC_number_taxa']).size().reset_index(name='counts')
    print('total ec taxa name', M['counts'].sum(), M.shape)
    M['CPM'] = ((M['counts']/(M.shape[0]) * 1000000))
    M = M.round({'CPM': 2})
    M.to_csv(out+'/'+s+'_ectaxa.csv',index = False,sep="\t",header=True,columns=["EC_number_taxa","counts","CPM"])

    #pathway - taxa coverage 
    result = [x for x in md[['EC_number','taxa_name']].values.tolist()] #based on the ec number duplicate the line
    l =[]
    for i in result:
        ec =i[0]
        for k,v in D.items():
            for val in v:
                if ec == val:
                    l.append((ec,i[1],k.split('00')[1]))

    #create daframe from list of pathways obtained from minpath
    df = pd.DataFrame(l, columns =['EC_number', 'taxa_name', 'MAP'])

    #merge path with pathway name
    df1 = pd.merge(df,Dpath,on='MAP')

    #pathway coverage taxa
    path_coverage = df1.groupby(['MAP','Kegg_pathway','taxa_name']).size().reset_index(name='counts')

    #path_coverage = path_coverage[path_coverage.counts > 5]
    path_coverage['coverage']= path_coverage['counts'] / path_coverage.groupby(['MAP','Kegg_pathway'])["counts"].transform("sum")

    
    path_coverage['path_taxa'] = path_coverage[['MAP','Kegg_pathway','taxa_name']].apply(lambda x: '|'.join(x), axis=1)
    path_coverage.to_csv(out+'/'+s+'_pathtaxa_coverage.csv',index = False,sep="\t",header=True, 
                         columns=['path_taxa','counts','coverage'])

    ## pathway abundance
    p_coverage = df1.groupby(['MAP','Kegg_pathway']).size().reset_index(name='counts')
    p_coverage['MAP_Kegg_pathway'] = p_coverage[['MAP','Kegg_pathway']].apply(lambda x: '|'.join(x), axis=1)
    p_coverage['CPM'] = ((p_coverage['counts']/(p_coverage.shape[0]) * 1000000))
    p_coverage.to_csv(out+'/'+s+'_path_abundance.csv',index = False,sep="\t",header=True,columns=['MAP_Kegg_pathway','counts','CPM'])

    
    df1 = df1.astype(str)
  
    df1['EC_number_taxa_pathway'] = df1[['EC_number', 'taxa_name','MAP','Kegg_pathway']].apply(lambda x: '|'.join(x), axis=1)
    df1 = df1.drop(['EC_number', 'taxa_name','MAP','Kegg_pathway'], axis=1)

    df1.to_csv(out+'/'+s+'_ectaxapath.csv',index = False,sep="\t",header=True,columns=["EC_number_taxa_pathway"])
    
def Main(): 
    
    parser = argparse.ArgumentParser(description='python ..')    
    parser.add_argument('--samplesname',help='1. sample name')
    parser.add_argument('--detail',help='2. detail from minpath')
    parser.add_argument("--read_ec",help="3. read ec")
    parser.add_argument("--xeno",help="4. fun and taxa annotation")
    parser.add_argument("--pathname",help="5. pathway name")
    parser.add_argument("--output_folder",help="6. output folder")
 
    args = parser.parse_args()

    if len(sys.argv) != 13:
        print("You haven't specified any arguments. Use -h to get more details on how to use this command.")
        parser.print_help()
        sys.exit()
            
    xeno_mock(args.samplesname,args.detail,args.read_ec,args.xeno,args.pathname,args.output_folder)

if __name__ == '__main__':
    begin_time = datetime.now()
    print(begin_time)
    Main()
    print(datetime.now())
    print(datetime.now() - begin_time)

