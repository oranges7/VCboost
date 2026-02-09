import sys
import os

merge = sys.argv[1]
benchmark = sys.argv[2]
mode = sys.argv[3]
file_path, full_file_name = os.path.split(merge)

d={}
g={}
with open(merge, 'r') as m,open(benchmark,'r') as ben,open(file_path+'/f_'+mode+'.txt','w') as ff,open(file_path+'/p_'+mode+'.txt','w') as p:
    for a in m:
        if a[0] == '#':
            continue
        c = []
        b = a.strip().split('\t')

        c.append(b[3])
        c.append(b[4])
        c.append(b[5])
        c.append(b[6])
        d[b[0] + ' ' + str(b[1])] = c
    for e in ben:
        if e[0] == '#':
            continue
        f = []
        h = e.strip().split('\t')

        f.append(h[3])
        f.append(h[4])

        g[h[0] + ' ' + str(h[1])] = f

    for key in d:
        if len(d[key][0]) != len(d[key][1]):

            if key not in g:
                if key.strip().split(' ')[0] == 'chrX':
                    continue
                if key.strip().split(' ')[0] == 'chrY':
                    continue
                print(key, d[key][0], d[key][1], d[key][2], d[key][3], file=ff)
            if key in g:
                if d[key][0] == g[key][0] and d[key][1] == g[key][1]:
                    print(key, g[key][0], g[key][1], d[key][2], d[key][3], file=p)
                else:
                    print(key, d[key][0], d[key][1], d[key][2], d[key][3], g[key][0], g[key][1], file=ff)
