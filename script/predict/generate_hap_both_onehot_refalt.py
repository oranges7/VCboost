import sys, pysam, time, os, copy, argparse, subprocess, random, re, datetime
import numpy as np
from subprocess import Popen, PIPE, STDOUT

from argparse import ArgumentParser, SUPPRESS
from sys import exit, stderr



def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=stderr, bufsize=8388608):
    process = subprocess.Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize,
                               universal_newlines=True)

    # 等待命令执行完成
    stdout, stderr = process.communicate()

    # 检查命令是否成功执行
    if process.returncode == 0:
        # 返回命令的标准输出
        return stdout
    else:
        # 返回错误信息
        return stderr


def msa(seq_list, ref, mincov, maxcov):
    ref1 = ref.replace('N', '_').replace('R', '-').replace('Y', '-').replace('W', '-').replace('K', '-').replace('M',
                                                                                                                 '-').replace(
        'S', '-').replace('B', '-')

    mapping = {'A': 1, 'G': 4, 'T': 3, 'C': 2, '-': 0}
    # rev_base_map = {0: 'A', 1: 'G', 2: 'T', 3: 'C', 4: '-'}
    np.random.seed(812)
    sample = list(seq_list.keys())
    if len(sample) < mincov:
        return np.zeros((72, 5)), np.zeros((72, 5))
    if len(sample) > maxcov:
        sample = random.sample(sample, min(len(sample), maxcov))

    sample = sorted(sample)

    fa_tmp_file = ''.join(['>%s_SEQ\n%s\n' % (read_name, seq_list[read_name]) for read_name in sample])

    fa_tmp_file += '>ref_SEQ\n%s' % ref1

    gap_penalty = 1.0
    msa_process = Popen(['muscle', '-quiet', '-gapopen', '%.1f' % gap_penalty, '-maxiters', '1', '-diags1'],
                        stdout=PIPE, stdin=PIPE, stderr=PIPE)
    hap_file = msa_process.communicate(input=fa_tmp_file.encode('utf-8'))
    # print(hap_file)
    if len(hap_file) == 0:
        print('hapfile length 0')

    tmp = hap_file[0].decode('utf-8')[1:].replace('\n', '').split('>')
    # print(tmp)
    zz_0 = []
    for seq in tmp:
        p1, p2 = seq.split('_SEQ')
        if p1 != 'ref':
            zz_0.append(p2)
        else:
            ref_real_0 = p2

    try:
        ref_real_0_mat = np.eye(5)[[mapping[x] for x in ref_real_0]]
    except UnboundLocalError:
        return np.zeros((72, 5)), np.zeros((72, 5))

    mat = np.array([[mapping[c] for c in x] for x in zz_0], dtype=int)
    h0_mat = np.sum(np.eye(5)[mat], axis=0).astype(np.float32)
    alt_mat = (h0_mat / (np.sum(h0_mat, axis=1)[:, np.newaxis]))

    alt_mat *= 100
    ref_real_0_mat = ref_real_0_mat[:128, :]
    ref_real_0_mat *= 100
    if ref_real_0_mat.shape[0] < 128:
        rows_to_add = 128 - ref_real_0_mat.shape[0]
        zeros_matrix = np.zeros((rows_to_add, ref_real_0_mat.shape[1]))
        ref_real_0_mat = np.vstack((ref_real_0_mat, zeros_matrix))
    if alt_mat.shape[0] < 128:
        rows_to_add = 128 - alt_mat.shape[0]
        zeros_matrix = np.zeros((rows_to_add, alt_mat.shape[1]))

        alt_mat = np.vstack((alt_mat, zeros_matrix))

    return ref_real_0_mat[:72, :], alt_mat[:72, :]


def gen_hap_data(d):
    mapping = {'A': 25, 'G': 50, 'T': 75, 'C': 100, '-': 0}
    hap_0_data = np.zeros((80, 33))
    for index, i in enumerate(d[0].values()):
        if index >= 40:
            continue
        if len(i) < 33:
            continue
        hap_0_line = np.array([mapping[x] for x in i])
        hap_0_data[index, :] = hap_0_line
    hap_1_data = np.zeros((40, 33))
    for index, i in enumerate(d[1].values()):
        if index >= 40:
            continue
        if len(i) < 33:
            continue
        hap_1_line = np.array([mapping[x] for x in i])
        hap_1_data[index, :] = hap_1_line
    return hap_0_data, hap_1_data


def _normalize_mq(x):
    return int(100 * min(x, 60) / 60)


def gen_data_one_hap(d):
    mapping = {'A': 25, 'G': 50, 'T': 75, 'C': 100, '-': 0}
    hap_0_data = np.zeros((80, 33))
    hap_1_data = np.zeros((80, 33))
    mq_data = np.zeros((80, 33))
    count = 0
    for index, i in enumerate(d.values()):
        if index >= 80:
            continue
        # print(len(i[0]))
        if len(i[0]) < 33:
            continue
        hap_0_line = np.array([mapping[x] for x in i[0]])
        hap_0_data[index, :] = hap_0_line
        if i[1] == 1:
            hap_1_line = np.tile(np.array([100]), (1, 33))
        else:
            hap_1_line = np.tile(np.array([0]), (1, 33))
        hap_1_data[index, :] = hap_1_line
        mq_line = np.tile(np.array([_normalize_mq(i[2])]), (1, 33))
        mq_data[index, :] = mq_line
        count += 1
    return hap_0_data, hap_1_data, count, mq_data


def gen_data(d, ref):
    ref1 = ref.replace('N', '-').replace('R', '-').replace('Y', '-').replace('W', '-').replace('K', '-').replace('M',
                                                                                                                 '-').replace(
        'S', '-').replace('B', '-')
    mapping = {'A': 25, 'G': 50, 'T': 75, 'C': 100, '-': 0}
    ref_data = np.array([[mapping[x] for x in ref1]])

    h0, h1, count, mq = gen_data_one_hap(d)
    if count < 80:
        rows_to_add = 80 - count
        zero_rows = np.zeros((rows_to_add, ref_data.shape[1]))
        ref_data = np.vstack([np.tile(ref_data, (count, 1)), zero_rows])
    else:
        ref_data = np.tile(ref_data, (80, 1))


    hap_data = np.dstack([ref_data, h0, h1, mq])
    return hap_data


def run(cp0, fastafile, samfile, pileup_file,miss_output, hap_reads_0, hap_reads_1):
    mincov = 1
    maxcov = 100

    flag = 0x4 | 0x100 | 0x200 | 0x400 | 0x800
    chrom = cp0.strip().split(' ')[0]
    v_pos = int(cp0.strip().split(' ')[1])
    # info = [cp0.strip().split(' ')[2].replace(',', '-').replace('N', '-').replace('R', '-').replace('Y', '-').replace('W', '-').replace('K', '-').replace('M', '-').replace('S', '-').replace('B', '-'),cp0.strip().split(' ')[3].replace(',', '-').replace('N', '').replace('R', '-').replace('Y', '-').replace('W', '-').replace('K', '-').replace('M', '-').replace('S', '-').replace('B', '-')]
    info1 = cp0.strip().split(' ')[2].replace('N', '-').replace('R', '-').replace('Y', '-').replace('W', '-').replace(
        'K', '-').replace('M', '-').replace('S', '-').replace('B', '-')
    info2 = cp0.strip().split(' ')[3].replace('N', '-').replace('R', '-').replace('Y', '-').replace('W', '-').replace(
        'K', '-').replace('M', '-').replace('S', '-').replace('B', '-')
    ref = fastafile.fetch(chrom, v_pos - 1, v_pos + 200)
    if ',' in info1:
        info1 = info1.split(',')[0]
    if ',' in info2:
        info2 = info2.split(',')[0]

    ref = info2 + ref[len(info1):]

    ref = ref[:160]

    d = {0: {}, 1: {}}
    in_run = 0
    for pcol in samfile.pileup(chrom, v_pos - 1, v_pos, min_base_quality=12, min_mapping_quality=28, flag_filter=flag,
                               truncate=True):
        for pread in pcol.pileups:
            dt = pread.alignment.query_sequence[max(0, pread.query_position_or_next):pread.query_position_or_next + 160]
            mq = pread.alignment.mapping_quality

            if pread.alignment.qname in hap_reads_0:
                d[0][pread.alignment.qname] = dt
            elif pread.alignment.qname in hap_reads_1:
                d[1][pread.alignment.qname] = dt
        if len(d[0]) == 0:
            ref_data_0, data_0 = np.zeros((72, 5)), np.zeros((72, 5))
        else:
            ref_data_0, data_0 = msa(d[0], ref, 2, 100)
        if len(d[1]) == 0:
            ref_data_1, data_1 = np.zeros((72, 5)), np.zeros((72, 5))
        else:
            ref_data_1, data_1 = msa(d[1], ref, 2, 100)
        data = np.dstack([ref_data_0, data_0, ref_data_1, data_1])
        # if data_0 is None or data_1 is None:
        #     continue
        in_run = 1
        s = chrom + ' ' + str(v_pos) + ' ' + ' '.join(str(x) for x in data.reshape(-1).astype(np.int16)) + '\n'
        pileup_file.write(s)
    if in_run == 0:
        s = chrom + ' ' + str(v_pos) + '\n'
        miss_output.write(s)


def gen(args):
    fb = args.pos_path
    num_line_result = subprocess_popen(["wc", "-l", fb])
    num_line = int(num_line_result.split()[0])

    sam_path = args.bam_path
    ref_path = args.ref_path

    pileup_file = open(args.out_path + '/' + fb.split('/')[-1] + '.pileup', "w")
    miss_output = open(args.out_path + '/' + 'miss' + '.pos', "a")
    samfile = pysam.Samfile(sam_path, "rb")
    fastafile = pysam.FastaFile(ref_path)
    kf = open(fb, 'r')
    chrom = args.chromosome
    length = fastafile.get_reference_length(chrom)
    hap_dict = {1: [], 2: []}
    for pread in samfile.fetch(chrom, 0, length + 1):
        if pread.has_tag('HP'):
            hap_dict[pread.get_tag('HP')].append(pread.qname)

    hap_reads_0 = set(hap_dict[1])
    hap_reads_1 = set(hap_dict[2])
    for chr_pos in kf:
        run(chr_pos, fastafile, samfile, pileup_file,miss_output, hap_reads_0, hap_reads_1)



def main():
    parser = ArgumentParser(
        description="The purpose of this file is to try to use the snp site in the vcf file to phase the phasing bam file, and extract the hap information as the feature of the pf model.")
    parser.add_argument('--bam_path', '-b', type=str, default=None,
                        help="Enter the path to the bam file. (default: %(default)s)")
    parser.add_argument('--chromosome', '-c', type=str, default=None,
                        help="chromosome. (default: %(default)s)")
    parser.add_argument('--ref_path', '-r', type=str,default=None,
                        help="Input the path of the REF file. (default: %(default)s)")
    parser.add_argument('--pos_path', '-p', type=str, default=None,
                        help="The path of the input sites file. (default: %(default)s)")
    parser.add_argument('--out_path', '-o', type=str, default=None,
                        help="The path of the output file. (default: %(default)s)")
    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    gen(args)


if __name__ == '__main__':
    main()
