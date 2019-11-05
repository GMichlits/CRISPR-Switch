__author__ = 'georg.michlits'

mapped_file = open('Ulimap_v3.txt','r')

counts = {}
counts['guides'] = {}
counts['exp_rev'] = {}
for line in mapped_file:
    column = line.rstrip('\n').split('\t')
    #print(column)
    guide = column[0]
    exp = column[7]
    if not guide in counts['guides']:
        counts['guides'][guide] = 0
    counts['guides'][guide] += 1
    if not exp in counts['exp_rev']:
        counts['exp_rev'][exp] = 0
    counts['exp_rev'][exp] += 1

outfile = open('basic_checks_full.txt','w')
for guide in counts['guides']:
    outfile.write(guide + '\t' + str(counts['guides'][guide]) + '\n')
outfile.write('\n')
for exp in counts['exp_rev']:
    outfile.write(exp + '\t' + str(counts['exp_rev'][exp]) + '\n')
outfile.write('\n')
