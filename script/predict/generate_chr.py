import sys
f1 = sys.argv[1]
f2 = sys.argv[2]
chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",]
chromosomes_xy = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22","chrX","chrY"]
with open(f2,'w') as out:
    if f1 == 'chr1-22':
        for chr in chromosomes:
            print(chr, file=out)
    if f1 == 'chr1-22XY':
        for chr in chromosomes_xy:
            print(chr, file=out)
    if f1 in chromosomes_xy:
        print(f1, file=out)