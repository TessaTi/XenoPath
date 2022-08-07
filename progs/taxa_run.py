#!/usr/bin/env python
import pandas as pd 
import os 
import sys 
import csv 
import subprocess 
import time 

#########################################
#### READS PAIRED ####
#########################################
def kaijufun(kaijuarg):
    cl = subprocess.Popen([("kaiju %s") %(str(kaijuarg))],shell = True)
    cl.wait() 

### taxa name 
def kaijuadd(kaijusett):
    clname = subprocess.Popen([("kaiju-addTaxonNames %s") %(str(kaijusett))],shell=True)
    clname.wait()

def diamondfun(diamondarg):
    dl = subprocess.Popen([("diamond blastx %s") %(str(diamondarg))],shell = True)
    dl.wait()

#merge annotation
def funannot(sett):
    print(sett)
    cl = subprocess.Popen(sett)
    cl.wait()

def read_ec_minpath(sample,funtax):
    print('get read ec..')
    real = pd.read_csv(os.path.join(funtax,sample+'_xeno.csv'), delimiter ='\t', header = 0)
    real = real.drop_duplicates('qseqid') # all the xenobiotic real reads
    real = real[['qseqid','EC_number']]
    real.to_csv(funtax+'/'+sample+'_read_ec.csv', header=False, sep='\t', index=False)

def minpathrun(sett):
    cl = subprocess.Popen(sett)
    cl.wait()

def maketable(sett):
    cl = subprocess.Popen(sett)
    cl.wait()    

def kaijutables(sett):
    cl = subprocess.Popen([("kaiju2table %s") %(str(sett))],shell = True)
    cl.wait()

def reformatable(setlevel):
    cl = subprocess.Popen(setlevel)
    cl.wait()

def getxenotaxa(settxeno):
    cl = subprocess.Popen(settxeno)
    cl.wait()

def outaxa(output,kaiju_config, taxaout,scriptfolder):
    D = {'fulltaxa':['.out','kaiju.tsv'],'xenotaxa':['.txt','kaiju_xeno.tsv']}
    # taxa output full report taxonomy for all samples
    for k,v in D.items():
        if not os.path.exists(os.path.join(output,v[1])):
            settable = ' -t '+kaiju_config['node']+' -n '+kaiju_config['names']+' -r '+' species '+' -l '+ 'phylum,family,genus,species'+' -o '+os.path.join(output,v[1])+' '+os.path.join(taxaout,'*'+v[0])+' -v'
            #print('settable',settable)
            kaijutables(settable)
            #reformat table at level
            #setlevel = ['python',os.path.join(scriptfolder,'level.py'), '--reportable',os.path.join(output,v[1]),'--out',output,'--type',k]
            setlevel = ['python', os.path.join(scriptfolder,'level.py'), '--reportable',os.path.join(output,v[1]),'--out',output,'--type',k]
            #print(setlevel)
            reformatable(setlevel)
 
def outables(sett):
    cl = subprocess.Popen(sett)
    cl.wait()

def pairedfun(samples,F,kaiju_config,diamond_config,threads,tabfile,scriptfolder,general_config,output,database):
#def pairedfun(samples,F,kaiju_config,diamond_config,threads,tabfile,scriptfolder,general_config,output): 
   
    taxaout = os.path.join(output,'taxonomy')
    funout = os.path.join(output,'diamond')
    funtax = os.path.join(output,'classificationxp')
    minpathfolder = os.path.join(output,'minpathfolder')
    tables = os.path.join(output,'tables')
   
    try: 
        if not os.path.exists(taxaout):
            os.mkdir(taxaout)
        if not os.path.exists(funout):
            os.mkdir(funout)
        if not os.path.exists(funtax):
            os.mkdir(funtax)
        if not os.path.exists(minpathfolder):
            os.mkdir(minpathfolder)
        if not os.path.exists(tables):
            os.mkdir(tables)
   
    except ValueError: 
       sys.exit('taxa/fun folder already exist')
    
    taxalevel= ['phylum','family','genus','species']	

    D ={} # dictionary with number of reads annotated as xenobiotic
    for s in samples: 
        #print('sample',s)
        #print('table', F.columns)
        i = F[(F['samples'] == s) & (F['pair'] == 'r1')]['file'].values[0]
        #print('i', i)
        j = F[(F['samples'] == s) & (F['pair'] == 'r2')]['file'].values[0]
        # run kaiju 
        print('taxonomic classification starts..',time.time())
        
        # arguments 
        kaijuarg =' -z '+ str(threads) + ' -t '+kaiju_config['node']+' -f '+kaiju_config['fmi']+' -i '+str(i)+' -j '+str(j)+' -v '+'-o '+os.path.join(taxaout,s+'.out').strip("'")
        
        #print(kaijuarg)                                
        if not os.path.exists(os.path.join(taxaout,s+'.out')):
            kaijufun(kaijuarg) 
        print('taxonomic classification ends..',time.time())

        
        kaijuadd_sett =  ' -t '+kaiju_config['node']+' -n '+kaiju_config['names']+' -r superkingdom,phylum,family,genus,species '+' -v '+ ' -i '+os.path.join(taxaout,s+'.out')+' -o '+os.path.join(taxaout,s+'_'+'names'+'.tsv')
        if not os.path.exists(os.path.join(taxaout,s+'_names.tsv')):
            kaijuadd(kaijuadd_sett)
        print('end..',time.time())
    #########################################
    #### READS FUNCTIONAL CLASSIFICATION ####
    #########################################
        
        mode = ''
        if 'mode' in diamond_config.keys() or diamond_config['mode'] != None:
            mode = ' --'+diamond_config['mode']
        diamond_sett_i = ' -p '+ str(threads) +' --db '+ database+ ' -q '+ i +' -o ' + os.path.join(funout,s+'_i.out') +' -f 6 qseqid qlen sseqid sallseqid slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qframe btop stitle salltitles qcovhsp qcovhsp' + mode  
        diamond_sett_j = ' -p '+ str(threads) +' --db '+database+ ' -q '+ j +' -o ' + os.path.join(funout,s+'_j.out') +' -f 6 qseqid qlen sseqid sallseqid slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qframe btop stitle salltitles qcovhsp qcovhsp'  
    
        #print(diamond_sett_i)
        #print(diamond_sett_j)
        print('xenobiotic classification starts..',time.time())
        
        if not os.path.exists(os.path.join(funout,s+'_i.out')):
            diamondfun(diamond_sett_i)
            diamondfun(diamond_sett_j)
        print('end..',time.time())
        
        #concat frame
        mlist = [os.path.join(funout,s+'_i.out'), os.path.join(funout,s+'_j.out')]
        rf = []
        for f in mlist:
            if os.path.exists(f):
                rf.append('exist')
            else: 
                rf.append('none')
        try:
            if len(list(set(rf))) == 1 and list(set(rf))[0] == 'exist':                                                           
                #concatentate file 
                if not os.path.exists(os.path.join(funout,s+'.out')):
                    print(mlist)
                    cmd = ['cat',mlist[0],mlist[1]]
                    #print(mlist[0],mlist[1])
                    #print(cmd)
                    with open((funout+'/'+s+'.out'), "w+") as fin:
                        sb = subprocess.Popen(cmd, stdout = fin)
                        sb.wait()
                    fin.close()
        except Exception as e:
            print >> sys.stderr, "does not exist"
            print >> sys.stderr, "Exception: %s" % str(e)
            sys.exit(1)
         
           
        #for fx in mlist:
        #    os.remove(fx)        
    #########################################
    #### READS FUNCTIONAL CLASSIFICATION ####
    #########################################
        #p = 'getfun.py'
        p = str(os.path.join(scriptfolder,'getfun.py'))
        #annotsett = ['python',p,'--samplesname',s,'--diamond_folder',funout,'--taxa_folder',taxaout,'--swissprot_tab',database+'/swissprot_xenobiotic_degra_2020-06-29.csv','--output_folder',funtax]
        annotsett = ['python',p,'--samplesname',s,'--diamond_folder',funout,'--taxa_folder',taxaout,'--swissprot_tab',tabfile,'--output_folder',funtax]
        #print(annotsett)
        if not os.path.exists(os.path.join(funtax,s+'_xeno.csv')):
            funannot(annotsett)
       
    #########################################
    #### READS EC  ####
    #########################################
        read_ec_minpath(s,funtax)
        print('pathway..',time.time())
        #settminpath = ['python',os.path.join(general_config['minpathpy'],'MinPath1.4.py'),'-any',os.path.join(funtax,s+'_read_ec.csv'),'-map',gene$
        
        settminpath = ['python',os.path.join(general_config['minpathpy'],'MinPath1.4.py'),'-any',os.path.join(funtax,s+'_read_ec.csv'),'-map',general_config['minpathec'],'-report',os.path.join(minpathfolder,s+'.report'),'-details',os.path.join(minpathfolder,s+'.details')]
        #print(settminpath)
        if not os.path.exists(os.path.join(minpathfolder,s+'.details')):
            minpathrun(settminpath)
        print('end...',time.time())
    #########################################
    #### OUTPUT  ####
    #########################################
        settout = ['python',os.path.join(scriptfolder,'maptable.py'),'--samplesname',s,'--detail',os.path.join(minpathfolder,s+'.details'),'--read_ec',os.path.join(funtax,s+'_read_ec.csv'),'--xeno',os.path.join(funtax,s+'_xeno.csv'),'--pathname',general_config['XenoPathfolder']+'/database/xenobiotic_degradation_pathways.txt','--output_folder',tables]
        
        #settout = ['python',os.path.join(scriptfolder,'maptable.py'),'--samplesname',s,'--detail',os.path.join(minpathfolder,s+'.details'),'--read_ec',os.path.join(funtax,s+'_read_ec.csv'),'--xeno',os.path.join(funtax,s+'_xeno.csv'),'--pathname',os.path.join(scriptfolder,'data','xenobiotic_degradation_pathways.txt'),'--output_folder',tables]                
        print('saving XenoPath output..')
        #print(settout)
        if not os.path.exists(os.path.join(tables,s+'_ectaxapath.csv')):
            maketable(settout)
        
    #########################################
    #### extract xeno out 
    #########################################
        #settxeno = ['python','xenotaxa.py','--samplesname',s,'--xenoread_ec',os.path.join(funtax,s+'_read_ec.csv'),'--taxa_kaiju',os.path.join(taxaout,s+'.out'),'--out',taxaout]
        #'python',os.path.join(scriptfolder,'xenotaxa.py')
        settxeno = ['xenotaxa.py','--samplesname',s,'--xenoread_ec',os.path.join(funtax,s+'_read_ec.csv'),'--taxa_kaiju',os.path.join(taxaout,s+'.out'),'--out',taxaout]
        #print(settxeno)
        if not os.path.exists(os.path.join(taxaout,s+'_xeno.txt')):
            getxenotaxa(settxeno)
  
    #get all samples table
    print('saving dataset tables..')
    outaxa(output,kaiju_config,taxaout,scriptfolder)
    #outaxa(output,kaiju_config,taxaout)
    settablex = ['funtable.py','--out',output]
    #print('settablex')
    outables(settablex)
   
    return

#def singlefun(samples,F,kaiju_config,diamond_config,threads,general_config,output,database):
def singlefun(samples,F,kaiju_config,diamond_config,threads,tabfile,scriptfolder,general_config,output,database):
    
    taxaout = os.path.join(output,'taxonomy')
    funout = os.path.join(output,'diamond')
    funtax = os.path.join(output,'classificationxp')
    minpathfolder = os.path.join(output,'minpathfolder')
    tables = os.path.join(output,'tables')

    
    try:
        if not os.path.exists(taxaout):
            os.mkdir(taxaout)
        if not os.path.exists(funout):
            os.mkdir(funout)
        if not os.path.exists(funtax):
            os.mkdir(funtax)
        if not os.path.exists(minpathfolder):
            os.mkdir(minpathfolder)
        if not os.path.exists(tables):
            os.mkdir(tables)

    except ValueError:
       sys.exit('taxa/fun folder already exist')

    taxalevel= ['phylum','family','genus','species']

    D ={} # dictionary with number of reads annotated as xenobiotic
    for s in samples:
        #print('sample',s)
        #print('table', F, F.columns)
        i = F[(F['samples'] == s) & (F['typereads'] == 'single')]['file'].values[0]
        #print('i', i)
        # run kaiju
        print('taxonomic classification starts',time.time())
        
        #script_kaiju_paired = os.path.join(foldescript, script_kaiju_single.sh)
        # arguments
        kaijuarg = ' -z '+ str(threads)+' -t '+kaiju_config['node']+' -f '+kaiju_config['fmi']+' -i '+str(i)+' -v '+'-o '+os.path.join(taxaout,s+'.out').strip("'")
        #print(kaijuarg)
        if not os.path.exists(os.path.join(taxaout,s+'.out')):
            kaijufun(kaijuarg) 
        print('end...',time.time())
        print('start add taxa name',time.time())
        kaijuadd_sett =  ' -t '+kaiju_config['node']+' -n '+kaiju_config['names']+' -r superkingdom,phylum,family,genus,species '+' -v '+ ' -i '+os.path.join(taxaout,s+'.out')+' -o '+os.path.join(taxaout,s+'_'+'names'+'.tsv')
        #print(kaijuadd_sett)
        if not os.path.exists(os.path.join(taxaout,s+'_names.tsv')):
            kaijuadd(kaijuadd_sett)
        print('end...', time.time())
    #########################################
    #### READS FUNCTIONAL CLASSIFICATION ####
    #########################################
        
        mode = ''
        if 'mode' in diamond_config.keys() or diamond_config['mode'] != None:
            mode = ' --'+diamond_config['mode']
        diamond_sett_i = ' -p '+ str(threads) +' --db '+ database+ ' -q '+ i +' -o ' + os.path.join(funout,s+'.out') +' -f 6 qseqid qlen sseqid sallseqid slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qframe btop stitle salltitles qcovhsp qcovhsp' + mode
        
        #print(diamond_sett_i)

        print('xenobiotic classification starts',time.time())
        
        if not os.path.exists(os.path.join(funout,s+'.out')):
            diamondfun(diamond_sett_i)
        print('end...',time.time())

    #########################################
    #### READS FUNCTIONAL CLASSIFICATION ####
    #########################################
        p = 'getfun.py'
        #p = str(os.path.join(scriptfolder,'getfun.py'))
        #annotsett = ['python',p,'--samplesname',s,'--diamond_folder',funout,'--taxa_folder',taxaout,'--swissprot_tab',database+'/swissprot_xenobiotic_degra_2020-06-29.csv','--output_folder',funtax]
        annotsett = [p,'--samplesname',s,'--diamond_folder',funout,'--taxa_folder',taxaout,'--swissprot_tab',tabfile,'--output_folder',funtax]
        #print('annotsett)
        if not os.path.exists(os.path.join(funtax,s+'_xeno.csv')):
            funannot(annotsett)

    #########################################
    #### READS EC  ####
    #########################################
        read_ec_minpath(s,funtax)
        #settminpath = ['python','MinPath1.4.py','-any',os.path.join(funtax,s+'_read_ec.csv'),'-map',database+'/map2ec','-report',os.path.join(minpathfolder,s+'.report'),'-details',os.path.join(minpathfolder,s+'.details')]
        
        settminpath = ['python',os.path.join(general_config['minpathpy'],'MinPath1.4.py'),'-any',os.path.join(funtax,s+'_read_ec.csv'),'-map',general_config['minpathec'],'-report',os.path.join(minpathfolder,s+'.report'),'-details',os.path.join(minpathfolder,s+'.details')]
        print(settminpath)
        if not os.path.exists(os.path.join(minpathfolder,s+'.details')):
            minpathrun(settminpath)

    #########################################
    #### OUTPUT  ####
    #########################################
        #settout = ['python',os.path.join(scriptfolder,'maptable.py'),'--samplesname',s,'--detail',os.path.join(minpathfolder,s+'.details'),'--read_ec',os.path.join(funtax,s+'_read_ec.csv'),'--xeno',os.path.join(funtax,s+'_xeno.csv'),'--pathname',os.path.join(scriptfolder,'data','xenobiotic_degradation_pathways.txt'),'--output_folder',tables]
        settout = ['python',os.path.join(scriptfolder,'maptable.py'),'--samplesname',s,'--detail',os.path.join(minpathfolder,s+'.details'),'--read_ec',os.path.join(funtax,s+'_read_ec.csv'),'--xeno',os.path.join(funtax,s+'_xeno.csv'),'--pathname',general_config['XenoPathfolder']+'/database/xenobiotic_degradation_pathways.txt','--output_folder',tables]
        print('saving XenoPath output..')
        #print(settout)
        if not os.path.exists(os.path.join(tables,s+'_ectaxapath.csv')):
            maketable(settout)
        

    #########################################
    #### extract xeno out
    #########################################
        settxeno = ['xenotaxa.py','--samplesname',s,'--xenoread_ec',os.path.join(funtax,s+'_read_ec.csv'),'--taxa_kaiju',os.path.join(taxaout,s+'.out'),'--out',taxaout]
        #settxeno = ['python',os.path.join(scriptfolder,'xenotaxa.py'),'--samplesname',s,'--xenoread_ec',os.path.join(funtax,s+'_read_ec.csv'),'--taxa_kaiju',os.path.join(taxaout,s+'.out'),'--out',taxaout]
        #print(settxeno)
        if not os.path.exists(os.path.join(taxaout,s+'_xeno.txt')):
            getxenotaxa(settxeno)
    
    
    print('saving dataset tables..')
    #outaxa(output,kaiju_config,taxaout)
    outaxa(output,kaiju_config,taxaout,scriptfolder)
    #settablex = ['python', os.path.join(scriptfolder,'funtable.py'),'--out',output]
    settablex = ['funtable.py','--out',output]
    outables(settablex)    
    return

def trimmomaticfun(samples,F,kaiju_config,diamond_config,threads,tabfile,scriptfolder,general_config,output,database):
#def trimmomaticfun(samples,F,kaiju_config,diamond_config,threads,tabfile,general_config,output,database):    
    taxaout = os.path.join(output,'taxonomy')
    funout = os.path.join(output,'diamond')
    funtax = os.path.join(output,'classificationxp')
    minpathfolder = os.path.join(output,'minpathfolder')
    tables = os.path.join(output,'tables')

    
    try:
        if not os.path.exists(taxaout):
            os.mkdir(taxaout)
        if not os.path.exists(funout):
            os.mkdir(funout)
        if not os.path.exists(funtax):
            os.mkdir(funtax)
        if not os.path.exists(minpathfolder):
            os.mkdir(minpathfolder)
        if not os.path.exists(tables):
            os.mkdir(tables)

    except ValueError:
       sys.exit('taxa/fun folder already exist')

    taxalevel= ['phylum','family','genus','species']

    D ={} # dictionary with number of reads annotated as xenobiotic
    for s in samples:
        #print('sample',s)
        #print('table', F.head, F.columns)
        i = F[(F['samples'] == s) & (F['typereads'] == 'paired') & (F['pair'] == 'r1')]['file'].values[0]
        j = F[(F['samples'] == s) & (F['typereads'] == 'paired') & (F['pair'] == 'r2')]['file'].values[0]
        ii = F[(F['samples'] == s) & (F['typereads'] == 'unpaired') & (F['pair'] == 'r1')]['file'].values[0]
        jj = F[(F['samples'] == s) & (F['typereads'] == 'unpaired') & (F['pair'] == 'r2')]['file'].values[0]
        #print('i', i,j,ii,jj)
        # run kaiju
        print('taxonomic classification starts',time.time())
        
        #run kaiju for paired end reads
        kaijuarg= '-z '+ str(threads) +' -t '+kaiju_config['node']+' -f '+kaiju_config['fmi']+' -i '+str(i)+' -j '+str(j)+' -v '+'-o '+os.path.join(taxaout,s+'_p.out').strip("'")
        #print(kaijuarg)
        if not os.path.exists(os.path.join(taxaout,s+'_p.out')):
            kaijufun(kaijuarg) 
        #run for single-and
           
        kaijuarg = '-z '+ str(threads)+' -t '+kaiju_config['node']+' -f '+kaiju_config['fmi']+' -i '+str(ii)+' -v '+'-o '+os.path.join(taxaout,s+'_ur1.out')
        #print(kaijuarg)
        if not os.path.exists(os.path.join(taxaout,s+'_ur1.out')):
            kaijufun(kaijuarg) 
            kaijuarg = '-z '+ str(threads)+' -t '+kaiju_config['node']+' -f '+kaiju_config['fmi']+' -i '+str(jj)+' -v '+'-o '+os.path.join(taxaout,s+'_ur2.out')
            if not os.path.exists(os.path.join(taxaout,s+'_ur2.out')):
                kaijufun(kaijuarg) 
        print('end...',time.time())
        
        #concat frame
        mlist = [os.path.join(taxaout,s+'_ur1.out'), os.path.join(taxaout,s+'_ur2.out'),os.path.join(taxaout,s+'_p.out')]
    
        rf = []
        for f in mlist:
            if os.path.exists(f):
                rf.append('exist')
            else:
                rf.append('none')
        #print('list rf', rf, list(set(rf)))
        try:
            if len(list(set(rf))) == 1 and list(set(rf))[0] == 'exist':
                #concatentate file
                if not os.path.exists(os.path.join(taxaout,s+'.out')):
                    #print(mlist)
                    cmd = ['cat',mlist[0],mlist[1],mlist[2]]
                    #print(mlist[0],mlist[1],mlist[2])
                    #print(cmd)
                    with open((taxaout+'/'+s+'.out'), "w+") as fin:
                        sb = subprocess.Popen(cmd, stdout = fin)
                        sb.wait()
                    fin.close()
        except Exception as e:
            print >> sys.stderr, "does not exist"
            print >> sys.stderr, "Exception: %s" % str(e)
            sys.exit(1)
         
        #for fx in mlist:
        #    os.remove(fx)
        
        kaijuadd_sett =  ' -t '+kaiju_config['node']+' -n '+kaiju_config['names']+' -r superkingdom,phylum,family,genus,species '+' -v '+ ' -i '+os.path.join(taxaout,s+'.out')+' -o '+os.path.join(taxaout,s+'_'+'names'+'.tsv')
        if not os.path.exists(os.path.join(taxaout,s+'_names.tsv')):
            kaijuadd(kaijuadd_sett)
    #########################################
    #### READS FUNCTIONAL CLASSIFICATION ####
    #########################################
        
        mode = ''
        if 'mode' in diamond_config.keys() or diamond_config['mode'] != None:
            mode = ' --'+diamond_config['mode']
        diamond_sett_i = ' -p '+ str(threads) +' --db '+ database + ' -q '+ i +' -o ' + os.path.join(funout,s+'_i.out') +' -f 6 qseqid qlen sseqid sallseqid slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qframe btop stitle salltitles qcovhsp qcovhsp' + mode
        diamond_sett_j = ' -p '+ str(threads) +' --db '+ database + ' -q '+ j +' -o ' + os.path.join(funout,s+'_j.out') +' -f 6 qseqid qlen sseqid sallseqid slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qframe btop stitle salltitles qcovhsp qcovhsp'
        diamond_sett_ii = ' -p '+ str(threads) +' --db '+database + ' -q '+ ii +' -o '+ os.path.join(funout,s+'_ii.out') +' -f 6 qseqid qlen sseqid sallseqid slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qframe btop stitle salltitles qcovhsp qcovhsp'
        diamond_sett_jj = ' -p '+ str(threads) +' --db '+ database + ' -q '+ jj +' -o ' + os.path.join(funout,s+'_jj.out') +' -f 6 qseqid qlen sseqid sallseqid slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qframe btop stitle salltitles qcovhsp qcovhsp'

        #print(diamond_sett_i)
        #print(diamond_sett_j)
        #print(diamond_sett_ii)
        #print(diamond_sett_jj)
        print('xenobiotic classification starts',time.time())
      
        if not os.path.exists(os.path.join(funout,s+'.out')):
            diamondfun(diamond_sett_i)
            diamondfun(diamond_sett_j)
            diamondfun(diamond_sett_ii)
            diamondfun(diamond_sett_jj)
        print('end...',time.time())
        
        mlist = [os.path.join(funout,s+'_i.out'), os.path.join(funout,s+'_j.out'),os.path.join(funout,s+'_ii.out'),os.path.join(funout,s+'_jj.out')]
        
        rf = []
        for f in mlist:
            if os.path.exists(f):
                rf.append('exist')
            else:
                rf.append('none')
        if len(list(set(rf))) == 1 and list(set(rf))[0] == 'exist':
            #concatentate file
            if not os.path.exists(os.path.join(funout,s+'.out')):
                #print(mlist)
                cmd = ['cat',mlist[0],mlist[1],mlist[2],mlist[3]]
                #print(mlist[0],mlist[1],mlist[2],mlist[3])
                #print(cmd)
                with open((funout+'/'+s+'.out'), "w+") as fin:
                    sb = subprocess.Popen(cmd, stdout = fin)
                    sb.wait()
                fin.close()
        else:
            sys.exit()
        #for fx in mlist:
        #    os.remove(fx)
        
    #########################################
    #### READS FUNCTIONAL CLASSIFICATION ####
    #########################################
        p = str(os.path.join(scriptfolder,'getfun.py'))
        #p = 'getfun.py'
        #annotsett = ['python',p,'--samplesname',s,'--diamond_folder',funout,'--taxa_folder',taxaout,'--swissprot_tab',database+'/swissprot_xenobiotic_degra_2020-06-29.csv','--output_folder',funtax]
        annotsett = ['python',p,'--samplesname',s,'--diamond_folder',funout,'--taxa_folder',taxaout,'--swissprot_tab',tabfile,'--output_folder',funtax]
        #print(annotsett)
        if not os.path.exists(os.path.join(funtax,s+'_xeno.csv')):
            funannot(annotsett)

    #########################################
    #### READS EC  ####
    #########################################
        read_ec_minpath(s,funtax)
        #settminpath = ['python','MinPath1.4.py','-any',os.path.join(funtax,s+'_read_ec.csv'),'-map',database+'/map2ec','-report',os.path.join(minpathfolder,s+'.report'),'-details',os.path.join(minpathfolder,s+'.details')]
        settminpath = ['python',os.path.join(general_config['minpathpy'],'MinPath1.4.py'),'-any',os.path.join(funtax,s+'_read_ec.csv'),'-map',general_config['minpathec'],'-report',os.path.join(minpathfolder,s+'.report'),'-details',os.path.join(minpathfolder,s+'.details')]
        #print(settminpath)
        if not os.path.exists(os.path.join(minpathfolder,s+'.details')):
            minpathrun(settminpath)

    #########################################
    #### OUTPUT  ####
    #########################################
        #settout = ['python',os.path.join(scriptfolder,'maptable.py'),'--samplesname',s,'--detail',os.path.join(minpathfolder,s+'.details'),'--read_ec',os.path.join(funtax,s+'_read_ec.csv'),'--xeno',os.path.join(funtax,s+'_xeno.csv'),'--pathname',os.path.join(scriptfolder,'database','xenobiotic_degradation_pathways.txt'),'--output_folder',tables]
        #settout = ['python','maptable.py','--samplesname',s,'--detail',os.path.join(minpathfolder,s+'.details'),'--read_ec',os.path.join(funtax,s+'_read_ec.csv'),'--xeno',os.path.join(funtax,s+'_xeno.csv'),'--pathname',database+'/xenobiotic_degradation_pathways.txt','--output_folder',tables]
        settout = ['python',os.path.join(scriptfolder,'maptable.py'),'--samplesname',s,'--detail',os.path.join(minpathfolder,s+'.details'),'--read_ec',os.path.join(funtax,s+'_read_ec.csv'),'--xeno',os.path.join(funtax,s+'_xeno.csv'),'--pathname',general_config['XenoPathfolder']+'/database/xenobiotic_degradation_pathways.txt','--output_folder',tables]
        print('saving XenoPath output..')
        #print(settout)
        if not os.path.exists(os.path.join(tables,s+'_ectaxapath.csv')):
            maketable(settout)
        

    #########################################
    #### extract xeno out
    #########################################
        #settxeno = ['python','xenotaxa.py','--samplesname',s,'--xenoread_ec',os.path.join(funtax,s+'_read_ec.csv'),'--taxa_kaiju',os.path.join(taxaout,s+'.out'),'--out',taxaout]
        settxeno = ['python',os.path.join(scriptfolder,'xenotaxa.py'),'--samplesname',s,'--xenoread_ec',os.path.join(funtax,s+'_read_ec.csv'),'--taxa_kaiju',os.path.join(taxaout,s+'.out'),'--out',taxaout]
        #settout = ['python',os.path.join(scriptfolder,'maptable.py'),'--samplesname',s,'--detail',os.path.join(minpathfolder,s+'.details'),'--read_ec',os.path.join(funtax,s+'_read_ec.csv'),'--xeno',os.path.join(funtax,s+'_xeno.csv'),'--pathname',general_config['XenoPathfolder']+'/database/xenobiotic_degradation_pathways.txt','--output_folder',tables]
        #print(settxeno)
        if not os.path.exists(os.path.join(taxaout,s+'_xeno.txt')):
            getxenotaxa(settxeno)

    #get taxa
    #outaxa(output,kaiju_config,taxaout)
    #settablex = ['python', 'funtable.py','--out',output]
    outaxa(output,kaiju_config,taxaout,scriptfolder)    
    settablex = ['funtable.py','--out',output]
    outables(settablex)


    return
