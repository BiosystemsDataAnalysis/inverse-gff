
#%% 

# script that creates the ranges of the non-coding parts, i.e. those parts that do not overlap with gene definitions
# fvdk, 2024

#%%

from itertools import compress
import subprocess
import sys
import multiprocessing

import resource


_packagesInstalled = True

try:
    import pandas as pd
except ImportError as e:
    print("pandas is not installed, please install before running")
    _packagesInstalled = False

try:
    import numpy as np
except ImportError as e:
    print("numpy is not installed, please install before running")
    _packagesInstalled = False


try:
    import psutil
except ImportError as e:
    print("psutil is not installed, please install before running")
    _packagesInstalled = False

# if not all packages were install, quit
if not(_packagesInstalled):
    sys.exit(1)


import argparse
from os.path import exists

# check if script is running in debug/interactive mode or not
interactive_mode = hasattr(sys, 'ps1')


#%%

parser = argparse.ArgumentParser()
parser.add_argument("gff",help="The input gff file")
parser.add_argument("out",help="The output gff file")
parser.add_argument("-t","--type",help="Look in the type column (default info column)",default=False,action="store_true")
parser.add_argument("-k","--key",help="The key to look for",type=str,default="transcript")
parser.add_argument("-c","--case",help="Do a case sensitive key match (default case insensitive)",default=False,action="store_true")
parser.add_argument("-m","--maxt",help="Set the maximum number of threads (default 1)",type=int,default=1)
parser.add_argument("-l","--lim",help="Set the memory limit (default 0.8)",type=float,default=0.8)

# only run this part when run in interactive/debug mode
if interactive_mode:
    sys.argv = ['','--gff','Mus_musculus.GRCm38.102.gff3','--out','test.gff3'] 
    
# parse the arguments
args = parser.parse_args()


def check_arguments():
    result = True    
    if len(args.gff)>0 and (not exists(args.gff)):
        print('The gff file cannot be found')
        result = False    

    return result

if not check_arguments():
    exit(1)

# find run lengths of zeros and ones, based on https://gist.github.com/akTwelve/dc0bbbf26fb14493898fc74cd2aa7f74
def rle(vector):
    # determine differences of vectors shifted one position from each other to determine locations
    changes = np.nonzero(np.append(vector,0) != np.append(0,vector))[0]
    # subtract the subsequent locations from each other to determine length
    # shift the changes vector by one and subtract
    rle_arr = (np.append(changes,changes[-1])-np.append(0,changes))[0:-1]
    # determine if there is a last stretch of ones
    remaining = len(vector) - rle_arr.sum()
    # append it to running length vector
    rle_arr = np.append(rle_arr,remaining)

    return rle_arr, changes


# %%


def create_chromosome_entries(chromosome, mpdict, sema, datablock:pd.DataFrame, arrayvec:np.array, columns,chr_region):

        try:
        
            print("adding non coding ranges for chromosome {0} ({1}-{2})".format(chromosome,int(chr_region.iloc[0].start),int(chr_region.iloc[0].stop)))

            for _,feature in datablock.iterrows():
                # remove one for zero based
                _start = int((feature.stop > feature.start) * feature.start + (feature.stop < feature.start) * feature.stop)-1
                # also include the last item 
                _stop =  int((feature.stop > feature.start) * feature.stop + (feature.stop < feature.start) * feature.start)
                # make range
                _range = [i for i in range(_start,_stop)]
                # set range to zero
                arrayvec[_range] = 0
                
            # find run lengths of array
            _,_pos = rle(arrayvec)
            # reshape array to start,stop matrix                
            _pos_start_stop = np.reshape(_pos,(len(_pos)//2,2))

            # begin with empty list
            row_list = []
            # add chromosome info to export
            # row_list.append(chr_region.to_dict('records')[0])
            _cnt = 1
            # loop over all sequences and create entry
            for _row in range(_pos_start_stop.shape[0]):
                s = _pos_start_stop[_row,0] + 1 # correct for 1 index again
                e = _pos_start_stop[_row,1]
                dict = {'id':chromosome,
                        'or':'MS',
                        'tp':'NC',
                        'start':s, 
                        'stop':e,
                        'nn1':'.',
                        'strand':'.',
                        'nn3':'.',
                        'info':'ID=NC:NC_{0}_{1:06d}'.format(chromosome,_cnt)} 
                _cnt = _cnt+1
                row_list.append(dict)

            # make it a dataframe  
            df = pd.DataFrame(row_list, columns=columns)
            # store result in dictionary
            mpdict[chromosome] = df
            # release semaphore
        except MemoryError as me:
            print("Memory error calculating chromosome/block {0}".format(chromosome))
            mpdict[chromosome] = -1
        
        sema.release()


## function to load gene definition file from gtf
def create_inverse_definition_file(gene_definition_filename):

    f_ = open(gene_definition_filename,'r')
    # first check the headers
    headers = []
    for l in f_:
        if l.startswith("#"):
            headers.append(l)
        else:
            break
    f_.close()

    print("importing gff file")
    # % read gene definitions
    gene_def = pd.read_csv(gene_definition_filename, sep='\t', skiprows=len(headers), header=None, low_memory=False)
    columns = ['id', 'or', 'tp', 'start', 'stop', 'nn1', 'strand', 'nn3', 'info']
    gene_def.columns = columns
    gene_def.drop(gene_def.index[-1], inplace=True)

    print("remove empty entries")
    # only select non na entries from table
    gene_def = gene_def.loc[~gene_def.tp.isna()].copy()

    # get unique chromosomes 
    chromosomes = gene_def.loc[gene_def.tp.str.upper()=='CHROMOSOME'].id.unique()

    # set number of threads
    concurrency = min(len(chromosomes),args.maxt)
    # create number of semaphores
    sema = multiprocessing.Semaphore(concurrency)
    # store the results
    block_results = multiprocessing.Manager().dict()
    # create vector for running processes
    mp_loop = []    


    for chromosome in chromosomes:
        # find unique chromosome entry
        chr_region = gene_def.loc[ (gene_def.id==chromosome) & (gene_def.tp.str.upper()=='CHROMOSOME')]
        if chr_region.shape[0]>1:
            print("skip chromosome {0}, too many ranges defined".format(chromosome))
            continue
        # skip empty chromosome definitions
        if chr_region.shape[0]==0:
            print("skip chromosome {0}, no range defined".format(chromosome))
            continue
        # arrayvec is now zero based
        arrayvec = np.ones((int(chr_region.iloc[0].stop-chr_region.iloc[0].start),1))
        # loop over all other features that were defined        
        if args.type:
            # if type selection
            if args.case:
                other_features = gene_def.loc[ (gene_def.id==chromosome) & (gene_def.tp == (args.key))]        
            else:
                other_features = gene_def.loc[ (gene_def.id==chromosome) & (gene_def.tp.str.upper() == (args.key.upper()))]        
        else:
            # case insensitive search for transcript in info part
            other_features = gene_def.loc[ (gene_def.id==chromosome) & (gene_def['info'].astype('str').str.contains(args.key,case=args.case))]
                    
        sema.acquire()
        _block_function = multiprocessing.Process(target=create_chromosome_entries,args=(chromosome,block_results,sema,other_features,arrayvec,columns,chr_region))
        mp_loop.append(_block_function)
        _block_function.start()

    # wait until all blocks have been processed
    for _p in mp_loop:
        _p.join()


    # check for memory errors and redo missing blocks, single thread
    for chromosome in block_results.keys():
        if type(block_results[chromosome]) is int:
            if block_results[chromosome]==-1:                
                print("redoing block {0}".format(chromosome) )
                chr_region = gene_def.loc[ (gene_def.id==chromosome) & (gene_def.tp.str.upper()=='CHROMOSOME')]
                if chr_region.shape[0]>1:
                    print("skip chromosome {0}, too many ranges defined".format(chromosome))
                    continue
                # skip empty chromosome definitions
                if chr_region.shape[0]==0:
                    print("skip chromosome {0}, no range defined".format(chromosome))
                    continue
                # arrayvec is now zero based
                arrayvec = np.ones((int(chr_region.iloc[0].stop-chr_region.iloc[0].start),1))
                # loop over all other features that were defined        
                if args.type:
                # if type selection
                    if args.case:
                        other_features = gene_def.loc[ (gene_def.id==chromosome) & (gene_def.tp == (args.key))]        
                    else:
                        other_features = gene_def.loc[ (gene_def.id==chromosome) & (gene_def.tp.str.upper() == (args.key.upper()))]        
                else:
                    # case insensitive search for transcript in info part
                    other_features = gene_def.loc[ (gene_def.id==chromosome) & (gene_def['info'].astype('str').str.contains(args.key,case=args.case))]
                
                sema = multiprocessing.Semaphore(1)
                sema.acquire()
                _block_function = multiprocessing.Process(target=create_chromosome_entries,args=(chromosome,block_results,sema,other_features,arrayvec,columns,chr_region))
                _block_function.start()
                # wait for function to finish
                _block_function.join()
                if type(block_results[chromosome]) is int:
                    if block_results[chromosome] == -1:
                        print("Not enough memory to process chromosome/block {0}, exiting ... ".format(chromosome))
                        sys.exit(1)


                

    # combine the blocks to a single dataframe
    df_inverse = block_results[chromosomes[0]]
    for _c in range(1,len(chromosomes)):
        df_inverse = pd.concat([df_inverse,block_results[chromosomes[_c]]],axis=0)        
                                
    # nake integers integers
    df_inverse['start']=df_inverse['start'].astype('int')
    df_inverse['stop']=df_inverse['stop'].astype('int')

    return df_inverse, headers


def limit_memory(maxperc=0.8): 
    m = psutil.virtual_memory()
    maxsize = int(m.total*maxperc)
    soft, hard = resource.getrlimit(resource.RLIMIT_AS) 
    resource.setrlimit(resource.RLIMIT_AS, (maxsize, hard)) 


# %% calling the main routine

limit_memory(args.lim)

igff2,_hdrs = create_inverse_definition_file(args.gff)

# write the regions to a file
print("writing to output file {0}".format(args.out))

filename = args.out
with open(filename,'w') as outfile:
    outfile.writelines(_hdrs)
    outfile.close()

igff2.to_csv(filename,header=None,index=None,mode="a",sep="\t")

# %%
