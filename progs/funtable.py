#!/usr/bin/env python
import sys
import pandas as pd
import argparse
import os
from datetime import datetime,date
import time

def outables(output):
    fmt = ['_path_abundance.csv','_ectaxa.csv','_pathtaxa_coverage.csv']
    for x in fmt:
        frames = []
        ff = [os.path.join(output,'tables',i) for i in os.listdir(os.path.join(output,'tables')) if i.endswith(x)]
        
        
        for filename in ff:
            print(filename)
            sample = filename.split('/')[-1].split('_')[0]
            df = pd.read_csv(filename, index_col=0, header=0,sep= '\t')
            
            if 'CPM' in list(df.columns):
                df.drop('CPM', axis=1, inplace=True)
                df=df.rename(columns = {'counts':sample})
            elif 'coverage' in list(df.columns) and 'counts' in list(df.columns):
                df.drop('counts', axis=1, inplace=True)
                df=df.rename(columns = {'coverage':sample})
            
            frames.append(df)
        frame = pd.concat(frames, axis=1,sort=False)
            
        frame.fillna(0, inplace=True)
        frame.to_csv(output+'/'+'all'+x, sep = '\t',index=True)
    
    
def Main():

    parser = argparse.ArgumentParser(description='python ..')
    parser.add_argument("--out",help="1. folder")

    args = parser.parse_args()

    if len(sys.argv) != 3:
        print("You haven't specified any arguments. Use -h to get more details on how to use this command.")
        parser.print_help()
        sys.exit()

    outables(args.out)

if __name__ == '__main__':
    begin_time = datetime.now()
    print(begin_time)
    Main()
    print(datetime.now())
    print(datetime.now() - begin_time)
