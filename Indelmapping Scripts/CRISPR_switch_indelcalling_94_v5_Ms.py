__author__ = 'georg.michlits'

mapfile = open(input('input file name e.g. testmapp.txt: '),'r')
outname = input('enter outfile name (no suffix) e.g.testmapped1: ')
indeldict_file = open(outname + '_indeldict.txt','w')
#indices = open('ind_scar_94.txt','r')
scarinfofile = open('scar_info_CRISPR-switch_vali_v3.txt','r')
#Properties = open('ZUB_ext3.txt','r')

guideProperty = {}
#i = 0
#for line in Properties:
#    if i > 0:
#        column = line.rstrip('\n').split('\t')
#        guide = column[0]
#        cut_orientation = column[25]
#        guideProperty[guide] = cut_orientation
#    i += 1
guideProperty['Line1_1'] = 'as'
guideProperty['Line1_2'] = 's'
guideProperty['ROSA26_1'] = 'as'
guideProperty['ROSA26_2'] = 's'
guideProperty['ROSA26_3'] = 'as'
guideProperty['Olfr10_1'] = 'as'
guideProperty['Olfr52_1'] = 'as'
guideProperty['Olfr44_1'] = 'as'
guideProperty['Nup93_1'] = 's'
guideProperty['Traip_7'] = 's'
guideProperty['VEGFA_1'] = 'as'
guideProperty['IGDCC3_x'] = 'as'
guideProperty['LOC116437_x'] = 's'
guideProperty['VEGFA_3'] = 's'
guideProperty['MAX_x'] = 'as'
guideProperty['COMDA_x'] = 's'
guideProperty['EMX1_x'] = 's'
guideProperty['HCN1_x'] = 'as'
guideProperty['MFAP1_x'] = 's'
guideProperty['FANCF_x'] = 's'
guideProperty['LINC00971_x'] = 'as'
guideProperty['SNX1_x'] = 'as'

#this file loads the output from scarmapping_v2.py.
#It determines the deletion size 0bp del is wt.
#insertions, substitiones and exon junct. will follow in a later version.
#finally it gives for every guide del from -30bp to +30bp in every condition as output.

#Easy structure:

#def function for indel mapping

#make dic[guide] of
#    1) Hap Dipl
#        2) 2,4,8,18 days
#            3) -20del to +20del
#            -count
#

#
#read input file
#    for every guide:
#        determine indel and store in a dict of indices#
#
#go though dictonary and print





#def functions for indel mapping
def revcomp(DNA):
    upDNA = DNA.upper()
    replacement1 = upDNA.replace('A', 't')
    replacement2 = replacement1.replace('T', 'a')
    replacement3 = replacement2.replace('C', 'g')
    replacement4 = replacement3.replace('G', 'c')
    complimentary_seq3_5 = replacement4
    complimentary_seq5_3 = complimentary_seq3_5[::-1]
    return(complimentary_seq5_3.upper())
assert revcomp('ACTNNNg') == 'CNNNAGT'

def indelmapping(REFAmplicon,Amplicon,cutsite,direction_read):
    REFAmplicon = REFAmplicon.upper()
    Amplicon = Amplicon.upper()
    l1 = len(REFAmplicon)
    l2 = len(Amplicon)
    l = min(l1,l2)
    pos_start = 0
    if pos_start < 0:
        pos_start = 0
    pos_end = l
    if pos_end > l:
        pos_end = l
    MM = 0
    for pos in range(pos_start,pos_end):
        bp_REF = REFAmplicon[pos]
        bp_map = Amplicon[pos]
        #print(REFAmplicon[0:pos+1])
        #print(Amplicon[0:pos+1])
        if bp_map == 'N' or bp_map == 'X':
            bp_map = bp_REF
        if not bp_REF == bp_map:
            #print('missmatch')
            #print('position ' + str(pos))
            #print('cutsite ' + str(cutsite))
            #print('read orientation ' + direction_read)
            #if direction_read == 'RV':
            #    print('mutation position ' + str(int(cutsite)-int(pos)))
            #else:
            #    print('mutation position ' + str(int(pos)-int(cutsite)))
            MM += 1
            Subscore = 0
            sub_check = 0
            for move_l in range(1,11):                          #check if 1bp sub Sub
                if pos+move_l < l:
                    sub_check +=1
                    if REFAmplicon[pos+move_l] == Amplicon[pos+move_l]:
                        Subscore += 1
            if sub_check >= 6 and Subscore/sub_check > 8/11:
                if direction_read == 'RV':
                    return ('Sub:' + str(1) + ':' + str(int(cutsite)-int(pos)))
                return ('Sub:' + str(1) + ':' + str(int(pos)-int(cutsite)+1))
            else:                             #now we check for del
                for shift in range(1,22):               #check fo del
                    n=0
                    #print(REFAmplicon[pos+shift:])
                    #print(Amplicon[pos:])
                    for move_l in range(1,7):
                        if pos + move_l < l2 and pos+move_l+shift < l1:
                            n += 1
                            if not REFAmplicon[pos+move_l+shift] == Amplicon[pos+move_l]:
                                break
                        if n == 6:
                            #print('shift is ' + str(shift))
                            if direction_read == 'RV':
                                #calculate bp correction overlap (for equal bp at ends of shift that over result in wrong pos)
                                step_continue = 'yes'
                                correction_bp_overlap = 0
                                for step_corr in range(1,shift+1):
                                    if step_continue == 'yes':
                                        if REFAmplicon[pos-step_corr] == REFAmplicon[pos-step_corr+shift]:
                                            correction_bp_overlap +=1
                                        else:
                                            step_continue = 'no'
                                return('del:' + str(shift) + ':' + str((shift-int(cutsite)+int(pos)-correction_bp_overlap)*-1))#-shift))
                            return('del:' + str(shift) + ':' + str(int(pos)-int(cutsite)))
                for shift in range(1,22):               #check fo ins
                    n=0
                    for move_l in range(1,7):
                        if pos+move_l+shift < l2 and pos+move_l < l1:
                            n += 1
                            if not REFAmplicon[pos+move_l] == Amplicon[pos+move_l+shift]:
                                break
                        if n == 6:
                            if direction_read == 'RV':
                                #calculate bp correction overlap (for equal bp at ends of shift that over result in wrong pos)
                                step_continue = 'yes'
                                correction_bp_overlap = 0
                                for step_corr in range(1,shift+1):
                                    if step_continue == 'yes':
                                        if Amplicon[pos-step_corr] == Amplicon[pos-step_corr+shift]:
                                            correction_bp_overlap +=1
                                        else:
                                            step_continue = 'no'
                                return('ins:' + str(shift) + ':' + str(int(cutsite)-int(pos)+correction_bp_overlap))
                            return('ins:' + str(shift) + ':' + str(int(pos)-int(cutsite)))
            return('NA:0:0')
    if MM == 0:
        return('wt:0:0')


assert(indelmapping('tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGTTTTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGTTTTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','54','FW')) == 'wt:0:0'
assert(indelmapping('tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGTTTTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGTTTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','54','FW')) == 'del:1:2'
assert(indelmapping('tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGTTTTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGTTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','54','FW'))== 'del:2:1'
assert(indelmapping('tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGTTTTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','54','FW'))== 'del:3:0'
assert(indelmapping('tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGTTTTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','54','FW'))== 'del:4:-1'
assert(indelmapping('tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGTTTTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGTTTTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','54','FW')) == 'wt:0:0'
assert(indelmapping('tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGTTTTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGT.TTTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','54','FW')) == 'ins:1:0'
assert(indelmapping('tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGTTTTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGT..TTTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','54','FW')) == 'ins:2:0'
assert(indelmapping('tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGTTTTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGT...TTTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','54','FW')) == 'ins:3:0'
assert(indelmapping('tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGTTTTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGTT....TTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','54','FW')) == 'ins:4:1'
assert(indelmapping('tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGTTTTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGT.TTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','54','FW')) == 'Sub:1:1'
assert(indelmapping('tcttctagGATGCTCACAACTGTATTCCTGAGCTGGACAATGAGACAGCCATGTTTTCTGTCTACGATGGACATGGAGgtaactttaggagatcatattggtagtagtgtaagacccct','CTACCTCCCAGACGAGCCTCACCCTCCATTCTATGAGGTGTATCGGAACAGTGAGTCGGTGACCCCCAATCCACGGTCCCCACTTGAGGACTATTCCCTCCACATCATTGACCTTCACACTGGAG','54','FW')) == 'NA:0:0'

#scarinfo[guide] = (cut,frame)
scarinfo_dict = {}

for line in scarinfofile:                   #collect all information regarding the cut site
    if not line.startswith('\t'):
        column = line.rstrip('\n').split('\t')
        guide = column[0].replace('CTRL_','').split('_')[2]+'_' + column[0].replace('CTRL_','').split('_')[3]
        Amplicon = column[1]
        FWprimer = column[2]
        RVprimer = column[3]
        length = column[4]
        cut = column[5]
        frame = column[6]
        scarinfo_dict[guide] = (cut, frame)
scarinfofile.seek(0)

#make dic[guide] of
#    1) Hap Dipl
#        2) 2,4,8,18 days
#            3) -20del to +20del
#            -count
indelDICT = {}
for line in scarinfofile:                   #collect guides in DICT
    if not line.startswith('\t'):
        column = line.rstrip('\n').split('\t')
        guide = column[0].replace('CTRL_','').split('_')[2]+'_' + column[0].replace('CTRL_','').split('_')[3]
        indelDICT[guide] = {}
exp_type = ['d2+4OH','d10+4OH','d10-4OH']

exp_rev_dict = {
'R1':'Cas9Dox',
'R2':'Cas9ERT',
'R3':'U6Dox',
'R5':'Switch',
'R6':'KO-later',
'R7':'wt',
'R8':'noDNA'}
exp_FW_dict = {
'F1':'Ms_d2_+4OH',
'F2':'Ms_d10_+4OH',
'F3':'Ms_d10_-4OH'}

time = ['d2+4OH','d10+4OH','d10-4OH']
for guide in indelDICT:
    for barcode in exp_rev_dict:
        exp_rev = exp_rev_dict[barcode]
        indelDICT[guide][exp_rev] = {}
        for d in exp_FW_dict:
            exp_fw = exp_FW_dict[d]
            indelDICT[guide][exp_rev][exp_fw] = {}

#read input file
#    for every guide:
#        determine indel and store in a dict of indices#

count_mapping_events = 0

for line in mapfile:
    count_mapping_events += 1
    if count_mapping_events%1000000 == 0:
        print(str(count_mapping_events/1000000) + 'million reads processed')
    column = line.rstrip('\n').split('\t')
    guide = column[0]
    REF_Amplicon = column[1]
    Amplicon = column[2]
    length_REFamplicon = int(column[3])
    cutsite = int(column[4])
    frame = column[5]
    direction_read = column[6]
    FW_RV_indices = column[7]
    if not 'na' in FW_RV_indices and not 'F4' in FW_RV_indices and not 'F5' in FW_RV_indices and not 'R4' in FW_RV_indices and not 'R9' in FW_RV_indices:
        exp_fw = exp_FW_dict[FW_RV_indices.split(':')[0]]
        exp_rev = exp_rev_dict[FW_RV_indices.split(':')[1]]
        mutation = indelmapping(REF_Amplicon,Amplicon,cutsite,direction_read)
        if guideProperty[guide] == 'as': #invert the pos label of the mutation if guide in antisense direction
            mutation = ':'.join([mutation.split(':')[0],mutation.split(':')[1],str(int(mutation.split(':')[2])*-1)])
        mutation_type = mutation.split(':')[0]
        mutation_size = int(mutation.split(':')[1])
        mutation_pos = int(mutation.split(':')[2])
        if mutation not in indelDICT[guide][exp_rev][exp_fw]:
            indelDICT[guide][exp_rev][exp_fw][mutation] = 0
        indelDICT[guide][exp_rev][exp_fw][mutation] += 1

#indelDICT[guide][exp_rv][exp_fw][mutation_type][mutation_size]=count

indeldict_file.write('guide' + '\t' + 'exp_rev' + '\t' + 'exp_fw' + '\t' + 'mutation' + '\t' + 'count')
for guide in sorted(indelDICT):
    for exp_rev in sorted(indelDICT[guide]):
        for exp_fw in sorted(indelDICT[guide][exp_rev]):
            for mutation in sorted(indelDICT[guide][exp_rev][exp_fw]):
                count = str(indelDICT[guide][exp_rev][exp_fw][mutation])
                indeldict_file.write('\n' + guide + '\t' + exp_rev + '\t' + exp_fw + '\t' + mutation
                                     + '\t' + str(count))

print('count mapping events: ' + str(count_mapping_events))
