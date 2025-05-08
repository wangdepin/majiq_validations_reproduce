import sys

if len(sys.argv) != 3:
	print ('usage: %s input_file output_file' % sys.argv)
in_name = sys.argv[1]
out_name = sys.argv[2]


in1 = open(in_name)
out1 = open(out_name, 'w+')
out1.write('#Created using: %s\n' % ' '.join(sys.argv))

for ll in in1.readlines():
    try:
        if ll.startswith('#'): continue
        grp1 = []
        grp2 = []
        tab=ll.strip().split()
        n = len(tab)
        if n == 0 : continue
        i = 0
        while(i<n):
            if tab[i] == '-o':
                i+=1
                outname = tab[i].split('/')[-1]
            elif tab[i] == '-grp1':
                i+=1
                while(not tab[i].startswith('-') and i<n):
                    exp = '.'.join(tab[i].split('/')[-1].split('.')[:-1])
                    grp1.append(exp)
                    i+=1
            elif tab[i] == '-grp2':
                i+=1    
                print (i, n, tab)
                while(i<n and not tab[i].startswith('-')):
                    exp = '.'.join(tab[i].split('/')[-1].split('.')[:-1])
                    grp2.append(exp)
                    i+=1 
            elif tab[i] == '-n' or tab[i] == '--names':
                i+=1
                nn1 = tab[i]
                nn2 = tab[i+1]
            else:
                i+=1
        s = "%s\t%s\t%s\t%s\t%s" % (nn1, nn2, outname, ','.join(grp1), ','.join(grp2))
        #print(s)
        out1.write('%s\n' %s)
    except Exception as e:
        print(e)
        print(tab)
        raise
in1.close()
out1.close()

