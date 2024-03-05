import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser, SUPPRESS


def calculate_metrics(confusion_matrix):
    num_classes = confusion_matrix.shape[0]
    precision = np.zeros(num_classes)
    recall = np.zeros(num_classes)
    f1 = np.zeros(num_classes)

    for i in range(num_classes):
        tp = confusion_matrix[i, i]
        fp = np.sum(confusion_matrix[:, i]) - tp
        fn = np.sum(confusion_matrix[i, :]) - tp

        precision[i] = tp / (tp + fp)
        recall[i] = tp / (tp + fn)
        f1[i] = 2 * (precision[i] * recall[i]) / (precision[i] + recall[i])

    return precision, recall, f1

def run(args):
    d={}
    combined_array = np.empty((0, 4))
    all_file = args.in_path
    output = args.out_path
    a = args.threshold
    or_file = args.original_file
    with open(all_file) as f:
        for i in f:
            single_file = i.strip()
            single_array = pd.read_csv(single_file+'.txt', sep=' ', header=None)
            combined_array = np.concatenate([combined_array, single_array], axis=0)

    chr = combined_array[:, 0].astype(str)
    pos = combined_array[:, 1].astype(int)
    pred_labels_one = combined_array[:, 2:].astype(float)
    class_labels = [0, 1]

    print(f"threshold: {a}")
    pred_labels = np.where(pred_labels_one[:, 1] > a, 1, 0)
    fin_file = output + '/vc_boost.vcf'
    for index, x in enumerate(pred_labels):
        if x == 0:
            d[str(chr[index]) + ' ' + str(pos[index])]=1
    add_pos = args.add_file
    with open(add_pos, 'r') as add:
        for line in add:
            atom = line.strip().split(' ')
            d[str(atom[0]) + ' ' + str(atom[1])] = 1
    with open(or_file, 'r') as vcf, open(fin_file, 'w') as fl:
        for line in vcf:
            if line[0] == '#':
                print(line.strip(), file=fl)
                continue
            atom = line.strip().split('\t')
            if args.conting != 'chr1-22':
                if atom[0] == 'chrX' or atom[0] == 'chrY':
                    continue
            if atom[0] + ' ' + atom[1] not in d:
                print(line.strip(), file=fl)

def main():
    parser = ArgumentParser(
        description="此文件旨在将结果文件合并，并输出到VCF文件中")
    parser.add_argument('--in_path', '-i', type=str, default=None,
                        help="输入bam文件的路径. (default: %(default)s)")
    parser.add_argument('--threshold', '-t', type=float, default=0.02,
                        help="输入阈值. (default: %(default)s)")
    parser.add_argument('--out_path', '-o', type=str, default=None,
                        help="输出文件的路径. (default: %(default)s)")
    parser.add_argument('--original_file', '-f', type=str, default=None,
                        help="原始文件. (default: %(default)s)")
    parser.add_argument('--add_file', '-a', type=str, default=None,
                        help="添加文件. (default: %(default)s)")
    parser.add_argument('--conting', '-c', type=str, default=None,
                        help="conting. (default: %(default)s)")
    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    run(args)


if __name__ == '__main__':
    main()

