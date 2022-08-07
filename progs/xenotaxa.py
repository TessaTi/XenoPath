#!/usr/bin/env python
import sys
import argparse
import pandas as pd
import os
from datetime import datetime,date
def xeno_taxonomy(s,xenoread_ec,taxa_kaiju,output):
    
    '''
    C       NB501204:67:HFTWHBGX5:1:11101:10000:1398        1678    255     28026,547043,566552,1150460,    WP_003834224.1,WP_004220774.1,WP_033500442.1,WP_126016489.1,     
    NAAIAGARASEQTNPTDEGYGVDHAVVQDAAGKLYDVFASNTPEGRKRLA,
    '''
    X = pd.read_csv(xenoread_ec, sep = '\t',header = None) #xeno red 
    X.columns = ['read', 'EC_number']
    T = pd.read_csv(taxa_kaiju,sep='\t', header=None,names= ['status','read','score1','score2','score3','score4','score5']) #full taxonoy for that sample
    #T.columns = ['status','read','score1','score2','score3','score4','score5']

    xt = T[T["read"].isin(X['read'])]
    print('read xeno shape', X.shape)
    print('full read shape', T.shape)
    print('samples', s, xt.shape)

    xt.to_csv(output+'/'+s+'_xeno.txt',sep ='\t', header = False, index=False )
    
def Main(): 
    
    parser = argparse.ArgumentParser(description='python ...')
    
    parser.add_argument('--samplesname',help='1. sample name')
    parser.add_argument('--xenoread_ec',help='2. xeno reads')
    parser.add_argument('--taxa_kaiju',help='3. taxa')
    parser.add_argument('--out',help='4. output')
    
    args = parser.parse_args()
    
    if len(sys.argv) != 9:
        print("You haven't specified any arguments. Use -h to get more details on how to use this command.")
        parser.print_help()
        sys.exit()
    xeno_taxonomy(args.samplesname,args.xenoread_ec,args.taxa_kaiju,args.out)

if __name__ == '__main__':
    begin_time = datetime.now()
    print(begin_time)
    Main()
    print(datetime.now())
    print(datetime.now() - begin_time)    
