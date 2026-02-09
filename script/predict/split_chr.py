import sys
f1 = sys.argv[1]
output_folder = sys.argv[2]
xy=sys.argv[3]
chromosomes = ['chr' + str(i) for i in range(1, 23)]
chromosome_dict = {chromosome: -1 for chromosome in chromosomes}
file_list=[]
if xy == 'chr1-22':
    with open(f1,'r') as all:
        for chr_n in all:
            if chr_n[0]=='#':
                continue
            chr_pos=chr_n.strip().split('\t')
            if chr_pos[0] not in chromosomes:
                continue
            chr = chr_pos[0]
            chromosome_dict[chr] += 1

            num = chromosome_dict[chr]//5000
            out_file = output_folder+'/'+chr_pos[0]+'_'+str(num)
            if out_file not in file_list:
                file_list.append(out_file)
                with open(out_file,'w') as chr_file:
                    print(chr_pos[0], chr_pos[1], chr_pos[3], chr_pos[4],file=chr_file)
            else:
                with open(out_file,'a') as chr_file:
                    print(chr_pos[0], chr_pos[1], chr_pos[3], chr_pos[4],file=chr_file)
else:
    chromosomes.append('chrX')
    chromosomes.append('chrY')
    chromosome_dict = {chromosome: -1 for chromosome in chromosomes}
    with open(f1, 'r') as all:
        for chr_n in all:
            if chr_n[0] == '#':
                continue
            chr_pos = chr_n.strip().split('\t')
            chr = chr_pos[0]
            if chr not in chromosomes:
                continue
            chromosome_dict[chr] += 1
            num = chromosome_dict[chr] // 5000

            out_file = output_folder + '/' + chr_pos[0] + '_' + str(num)
            if out_file not in file_list:
                file_list.append(out_file)
                with open(out_file, 'w') as chr_file:
                    print(chr_pos[0], chr_pos[1], chr_pos[3], chr_pos[4], file=chr_file)
            else:
                with open(out_file, 'a') as chr_file:
                    print(chr_pos[0], chr_pos[1], chr_pos[3], chr_pos[4], file=chr_file)