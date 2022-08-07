#!/usr/bin/env python
import pandas as pd
import argparse
import os 
import sys
from datetime import datetime,date


def taxonomymap(highest,taxa_folder, D, sample, out): 
    print('get fun 2')
  
    print('save xeno classification')
    ### for each read get the taxa name 
    df = pd.read_csv(os.path.join(taxa_folder,sample+'_names.tsv'), header=None, sep='\n')
    tt = df[0].str.split('\t', expand=True)
    #tt = pd.read_csv(os.path.join(taxa_folder,sample+'_names.tsv'), sep="\t", header=None, engine='python')

    tt.columns=["status", "qseqid", "ncbi_taxon_id", "score","sequences_matched_taxon","accession_number_sequence_matched","fragment_seq_matched","taxa_name"]

    mm = pd.merge(tt, highest, on = "qseqid", how='inner')
    
    fout = os.path.join(out,sample+'_xeno.csv')
    print('I am saving the output',fout)
    mm.to_csv(fout, header=True, sep='\t',index = False,index_label = False)
    n_annotated = mm.shape[0]
    print('sample',sample,'numbereads',n_annotated)
    #D = D.append({'sample':sample,'numbereads':n_annotated}, ignore_index=True)
    
    #return(D)

def mapswissprot(sample,input_folder,taxa_folder,swissprot, D, out):
    #sname,input_folder,taxonomy,Dictiomary_number_reads
    print('load data, get fun')
    
    #open diamond file
    f1 = os.path.join(input_folder,sample+'.out') #read 1 
   
    #open file
    r1 = pd.read_csv(f1, delimiter="\t", header=None, engine='python') #open file 
    r1.columns = ['qseqid','qlen', 'sseqid', 'sallseqid', 'slen', 'qstart', 'qend', 'sstart','send','qseq', 'sseq', 'evalue','bitscore','score', 'length', 'pident', 'nident','mismatch', 'positive', 'gapopen', 'gaps', 'ppos', 'qframe', 'btop', 'stitle', 'salltitles', 'qcovhsp', 'qcovhsp']
    '''
    ## r2 read 
    if os.path.exists(input_folder,sample+'_j.out'):
        f2 = os.path.join(input_folder,sample+'_j.out') #read 2
        r2 = pd.read_csv(f2, delimiter="\t", header=None, engine='python') #open file 
        r2.columns = ['qseqid','qlen', 'sseqid', 'sallseqid', 'slen', 'qstart', 'qend', 'sstart','send','qseq', 'sseq', 'evalue','bitscore','score', 'length', 'pident', 'nident','mismatch', 'positive', 'gapopen', 'gaps', 'ppos', 'qframe', 'btop', 'stitle', 'salltitles', 'qcovhsp', 'qcovhsp']
    
        r0 = pd.concat([r1, r2], axis=0)
    else:
    ''' 
    r0 = r1.copy(deep = True)
    r0['coverage_query'] = abs(round(100*((r0.qend - r0.qstart)/r0.qlen),2))     
    #REPLACE EMPTY STRING
    r0 = r0.replace(r'^\s*$', 'NA', regex=True)
    r0 = r0.replace(r'^Viruses*$', 'NA; NA; NA; NA; NA;', regex=True)
 
    # coverage and percentage identity
    r1subset = r0[(r0.pident >= 40.0) & (r0.coverage_query >= 80.0 )].copy()
    r1subset['sseqid'] = r1subset['sseqid'].str.split("|").str[1]
    res = pd.merge(swissprot, r1subset, on='sseqid')

    # drop empty lines with not ec number nor ko information
    resnonull = res.dropna(subset=['EC_number', 'Cross_reference_KO'])

    #order from highest to lowest
    df_out = resnonull.sort_values(['pident', 'coverage_query'], ascending=[False,False]).drop_duplicates(['qseqid']).reset_index(drop=True)
    
    df_out['qseqid'] = df_out['qseqid'].str.split("/").str[0]
    
    # take max value again after sorting 
    highest = df_out.sort_values(['pident', 'coverage_query'], ascending=[False,False]).drop_duplicates(['qseqid']).reset_index(drop=True)
    

    taxonomymap(highest,taxa_folder,D,sample, out)
    
    #return(D)
    
       

def funannot(sample,input_folder,taxonomy,swissprotab,out): 
    # create empty dictionary for the output 
    Dictiomary_number_reads = pd.DataFrame(columns=['sample', 'numbereads'])
    # open swissprot
    # compare the query with the swissprot data 
    print('get fun')
    swissprot = pd.read_csv(swissprotab, delimiter=",", header=0, engine='python')
    swissprot.columns = ['MAP','EC_number','sseqid', 'Entry name', 'Status', 'Protein_names', 'Gene_names',
                         'Organism','Organism ID', 'Length', 'Cross_reference_KEGG', 'Cross_reference_BioCyc', 
                         'Cross_reference_Pfam', 'Cross_reference_eggNOG', 'Cross_reference_STRING',
                         'Cross_reference_IntAct', 'Gene_ontology_cellular_component', 
                         'Gene_ontology_biological_process', 'Gene_ontology_molecular_function',
                         'Cross_reference_KO', 'Cross_reference_RefSeq', 'Cross_reference_EMBL',
                         'Cross_reference_CCDS','Pathway','Cross_reference_UniPathway',
                         'Cross-reference_Reactome','Kegg_pathway']
    # open samplesname
    ##fname = open(os.path.abspath(samplesname),'r')
    ##for line in fname:
        ##sample=line.strip()
        # parse r1 and r2 results from diamond 

    Dictiomary_number_reads = mapswissprot(sample,input_folder,taxonomy,swissprot,Dictiomary_number_reads,out)
    return(Dictiomary_number_reads)
    

def Main(): 
    
    parser = argparse.ArgumentParser(description='python ...')
    parser.add_argument('--samplesname',help='1. file with sample name')
    parser.add_argument("--diamond_folder",help="2. input folder diamond")
    parser.add_argument("--taxa_folder",help="3. taxa information")
    parser.add_argument("--swissprot_tab",help="4. tab swissprot")
    parser.add_argument("--output_folder",help="5. output folder")
 
    args = parser.parse_args()
   
    if len(sys.argv) != 11:
        print("You haven't specified any arguments. Use -h to get more details on how to use this command.")
        parser.print_help()
        sys.exit()

    funannot(args.samplesname,args.diamond_folder,args.taxa_folder,args.swissprot_tab,args.output_folder)
    
if __name__ == '__main__':
    begin_time = datetime.now()
    print(begin_time)
    Main()
    print(datetime.now())
    print(datetime.now() - begin_time)
