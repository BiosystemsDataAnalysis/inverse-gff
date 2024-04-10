# inverse-gff

python script to create ranges that are not related to any gene


# required packages

    pandas
    numpy
    typed-argument-parser
    
# usage

    python inverse-gff.py --gff *file.gff3* --out *out.gff3* [--gene]

    __--gene__ option to switch between gene definitions (i.e. type == 'gene') or transcripts. Transcripts are default
