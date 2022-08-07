#!/usr/bin/env python
'''
subset ec-taxa based on min number fo count and dataset composition
first row, index row
'''
import pandas as pd
import sys
import argparse
from datetime import datetime,date

def filterdf(dataframe,samples,x,p,out):
    
    df = pd.read_csv(dataframe, sep = '\t', header= 0,index_col=0)
    df = df.astype(int)
    print(df)
    print(df.shape)
    if p == 0:
        print(df[df > int(x)].count(axis = 1))
        df1 = df[df >= int(x)] 
    else:
        dff = df[df >= int(x)].count(axis = 1) 
        #print(dff)
        dff1 = dff[dff >= float(p)*int(samples)]
        # drop row not in dff1 
        df1 = df[df.index.isin(dff1.index)].dropna()
        
    print(df1)
    print(df1.shape)
    df1.to_csv(out+'filtered_table'+'_x_'+x+'_p_'+p+'.csv')


def Main():

    parser = argparse.ArgumentParser(description='python filter_table.py --dataframe file.txt --samples 2 --x 10 --p 1 --out ./')
    parser.add_argument('--dataframe',help='1. table')
    parser.add_argument("--samples",help="2. number of samples")
    parser.add_argument("--x",help="3. number of count/relative ab")
    parser.add_argument("--p",help="4. percentage of sample")
    parser.add_argument("--out",help="5. output folder")

    args = parser.parse_args()

    if len(sys.argv) != 11:
        print("You haven't specified any arguments. Use -h to get more details on how to use this command.")
        parser.print_help()
        sys.exit()

    filterdf(args.dataframe,args.samples,args.x,args.p,args.out)

if __name__ == '__main__':
    begin_time = datetime.now()
    print(begin_time)
    Main()
    print(datetime.now())
    print(datetime.now() - begin_time)


