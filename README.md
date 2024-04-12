# inverse-gff

python script to create ranges that are not related to any gene


# required packages (automatically installed)

    pandas
    numpy
    typed-argument-parser
    
# usage

    python inverse-gff.py --gff file.gff3 --out out.gff3 [--key gene] [--case] [--type]

    --key : a string to look for in the column 
    --type: use the type column to look for the key (default a the info column is selected)
    --case: do a case sensitive search (default the search is case insensitive)
    