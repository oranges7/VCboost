import shlex
import os
import sys
from argparse import ArgumentParser, SUPPRESS
from collections import defaultdict
from sys import exit, stderr
from subprocess import check_output, PIPE, Popen
from subprocess import PIPE


def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=stderr, bufsize=8388608):
    return Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize, universal_newlines=True)


def FiterHeteSnpPhasing(ctgName, args):
    vcf_fn = args.vcf_fn

    contig_name = ctgName
    split_folder = args.split_folder
    variant_dict = defaultdict(str)
    qual_set = defaultdict(int)

    header = []


    unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (vcf_fn)))
    for row in unzip_process.stdout:
        row = row.rstrip()
        if row[0] == '#':
            header.append(row + '\n')
            continue
        columns = row.strip().split()
        ctg_name = columns[0]
        if contig_name and contig_name != ctg_name:
            continue
        pos = int(columns[1])
        ref_base = columns[3]
        alt_base = columns[4]
        genotype = columns[9].split(':')[0].replace('|', '/')

        if len(ref_base) == 1 and len(alt_base) == 1:
            if genotype == '0/1' or genotype == '1/0':
                variant_dict[pos] = row
                qual = float(columns[5])
                qual_set[pos] = qual


    print('[INFO] Total heterozygous SNP positions selected: {}: {}'.format(contig_name, len(variant_dict)))

    f = open(os.path.join(split_folder, '{}.vcf'.format(contig_name)), 'w')
    f.write(''.join(header))
    for key, row in sorted(variant_dict.items(), key=lambda x: x[0]):
        f.write(row + '\n')
    f.close()


def main():
    parser = ArgumentParser(description="Select heterozygous snp candidates for WhatsHap phasing")
    parser.add_argument('--split_folder', type=str, default=None,
                        help="Path to directory that stores small bed region for raw alignment. (default: %(default)s)")

    parser.add_argument('--vcf_fn', type=str, default=None,
                        help="Path of the input vcf file. (default: %(default)s)")
    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)
    chr_list = list(range(1, 23)) + ['X', 'Y']
    for num in chr_list:
        contig_name = 'chr' + str(num)
        FiterHeteSnpPhasing(contig_name, args)


if __name__ == "__main__":
    main()
