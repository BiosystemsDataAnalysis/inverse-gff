
#%% 

# script that creates the ranges of the non-coding parts, i.e. those parts that do not overlap with gene definitions
# fvdk, 2024

#%%

import importlib.metadata
from itertools import compress
import subprocess
import sys

# check if packages are installed, if not install them
def check_python_installation():
    installed = [x.name for x in importlib.metadata.distributions()]
    required = ['pandas','typed-argument-parser', 'numpy']

    ireq =  [not(i in installed) for i in required]
    not_installed = list(compress(required,ireq))
    
    return set(not_installed)

missing = check_python_installation()

if missing:
     print('trying to install missing packages')
     python = sys.executable
     subprocess.check_call([python, '-m', 'pip', 'install', *missing], stdout=subprocess.DEVNULL)


# import the packages

import pandas as pd
import numpy as np
from tap import Tap
from os.path import exists

# check if script is running in debug/interactive mode or not
interactive_mode = hasattr(sys, 'ps1')


#%%

class myArguments(Tap):
    gff: str # the input gff file
    out: str # the output gff file
    type:bool=False # set the default column to look for to type (instead of info)
    key:str='transcript' # set the default string to look for
    case:bool=False # set the search to case sensitive, default case insensitve


# only run this part when run in interactive/debug mode
if interactive_mode:
    sys.argv = ['','--gff','Mus_musculus.GRCm38.102.gff3','--out','test.gff3'] 
    
# parse the arguments
args = myArguments().parse_args()


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

    # initialize output array
    df_noncoding = pd.DataFrame(columns=columns)    

    # get chromosome regions
    chromosomes = gene_def.id.unique()
    for chromosome in chromosomes:
        # find unique chromosome entry
        chr_region = gene_def.loc[ (gene_def.id==chromosome) & (gene_def.tp=='chromosome')]
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
        
        for _,feature in other_features.iterrows():
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

        row_list = []
        # add chromosome info to export
        # row_list.append(chr_region.to_dict('records')[0])
        _cnt = 1
        print("adding non coding ranges for chromosome {0} ({1}-{2})".format(chromosome,int(chr_region.iloc[0].start),int(chr_region.iloc[0].stop)))
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
        # append if necessary
        if df_noncoding.shape[0]>0:
            df_noncoding = pd.concat([df_noncoding,df],ignore_index=True)
        else:
            df_noncoding = df.copy()
            
        
    # nake integers integers
    df_noncoding['start']=df_noncoding['start'].astype('int')
    df_noncoding['stop']=df_noncoding['stop'].astype('int')

    return df_noncoding, headers

# %% calling the main routine

igff2,_hdrs = create_inverse_definition_file(args.gff)

# write the regions to a file
print("writing to output file {0}".format(args.out))

filename = args.out
with open(filename,'w') as outfile:
    outfile.writelines(_hdrs)
    outfile.close()

igff2.to_csv(filename,header=None,index=None,mode="a",sep="\t")

# %%
