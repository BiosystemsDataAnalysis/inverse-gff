# inverse-gff

python script to create ranges that are not related to any gene


# required packages (automatically installed)

    pandas
    numpy
    typed-argument-parser
    
# usage

    usage: inverse_gtf.py [-h] [-t] [-k KEY] [-c] [-m MAXT] gff out

    positional arguments:
        gff                   The input gff file
        out                   The output gff file

    optional arguments:
        -h, --help            show this help message and exit
        -t, --type            Look in the type column (default info column)
        -k KEY, --key KEY     The key to look for
        -c, --case            Do a case sensitive key match (default case insensitive)
        -m MAXT, --maxt MAXT  Set the maximum number of threads (default 1)
    
# examples:

## case sensitive search in the info column for the word gene
        
python inverse-gff.py Mus_musculus.GRCm39.111.gff3 out.gff3 --key gene --case

## example 2, case insensitive search in the type (3rd) column for the word gene
        
python inverse_gtf.py Mus_musculus.GRCm39.111.gff3 test.111.gene.gff3 --key gene --type

## example 3, case insentive search for transcript in the info column using 10 cores
python inverse_gtf.py Mus_musculus.GRCm39.111.gff3 test.111.gene.gff3 -m 10