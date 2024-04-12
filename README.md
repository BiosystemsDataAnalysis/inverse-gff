# inverse-gff

python script to create ranges that are not related to any gene


# required packages (automatically installed)

    pandas
    numpy
    typed-argument-parser
    
# usage

    python inverse-gff.py --gff file.gff3 --out out.gff3 [--key gene] [--case] [--type]

    --key : a string to look for in the column 
    --type: (boolean, default False) use the type column to look for the key (default a the info column is selected)
    --case: (boolean, default False) do a case sensitive search (default the search is case insensitive)
    

# examples:

## case sensitive search in the info column for the word gene
        
python inverse-gff.py --gff Mus_musculus.GRCm39.111.gff3  --out out.gff3 --key gene --case

## example 2, case insensitive search in the type (3rd column) for the word gene
        
python inverse_gtf.py --gff Mus_musculus.GRCm39.111.gff3 --out test.111.gene.gff3 --key gene --type