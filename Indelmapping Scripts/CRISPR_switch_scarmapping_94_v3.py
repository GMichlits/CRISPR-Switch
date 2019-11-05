__author__ = 'georg.michlits'

#Script identifies the guide based on the primer sequences used from PE file.
#It looks for the primer in the vicinity of the start of the read (+-5bp) and then matches the guide.
#next it extends or truncates the read so that it matches the REF Amplicon where the position of cut etc are known.
#finally it writes an output file, with Amplicon and matching inforamtion + index information

########## version 2

######### version 3
#uses scarinfofile_v2 which contains manually corrected scar infos. Sequences found in ctrl that differe from REF BL6 genome
#probably that are variants from an312 genome 4 SNPs were corrected in Rhox2f, Ints2, Traip3 and Ints10, plus a 15insertion
#that appears to be the wt variant in an312.

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


scarinfofile = open('scar_info_CRISPR-switch_vali_v3.txt','r')
PEsamfile = open('Uli_scar.sam','r') # here this is not a PE file name stayed for historic reasons
outfile = open('Ulimap_v3.txt','w')
indices = open('BC_Switch_ind_v2.txt','r')

dictIND = {}
for line in indices:
    col = line.rstrip('\n').split('\t')
    bc_FW = col[1].rstrip(' ')
    exp_FW = col[0]
    bc_RV = revcomp(col[3])
    exp_RV = col[2]
    dictIND[bc_FW] = exp_FW
    dictIND[bc_RV] = exp_RV
print(str(dictIND))

primDICT = {}
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
        primDICT[FWprimer] = {'guide':guide,                  #write a dictionary that holds that information
                                'Amplicon':Amplicon,
                                'length':length,
                                'cut':cut,
                                'frame':frame,
                                'direction_read':'FW'}
        primDICT[RVprimer] = {'guide':guide,                  #write a dictionary that holds that information
                                'Amplicon':Amplicon,
                                'length':length,
                                'cut':cut,
                                'frame':frame,
                                'direction_read':'RV'}

def MMmax(dna1,dna2,maxMM):
    c=0
    for i in range(0,len(dna1)):
        if not dna1[i] == dna2[i]:
            c = c + 1
            if c > maxMM:
                return False
                break
    return True

def searchprimer(DNA,primDICT):
    for primer in primDICT:
        first5prime = primer[0:5]
        for a in range(0,5):
            seed5bp_DNA = DNA[a:a+5]
            if seed5bp_DNA == first5prime:
                if MMmax(primer,DNA[a:a+len(primer)],1) == True:
                    return((primer,a,0))
    for primer in primDICT:
        for b in range(0,5):
            first5prime = primer[0+b:5+b]
            seed5bp_DNA = DNA[0:5]
            if seed5bp_DNA == first5prime:
                if MMmax(primer[b:],DNA[0:len(primer)-b],1) == True:
                    return((primer,0,b))
    return('not found',0,0)
assert searchprimer('xxxxTTTCGGCTGAGCTGTxCTCGxxxxx',['TTTCGGCTGAGCTGTACTCG'])[0] == 'TTTCGGCTGAGCTGTACTCG'
assert searchprimer('GGCTGAGCTGTxCTCGxxxxx',['TTTCGGCTGAGCTGTACTCG'])[0] == 'TTTCGGCTGAGCTGTACTCG'

nf = 0
l = 0
fw_primer_found = 0
rv_primer_found = 0

for line in PEsamfile:
    if not line == '\n':
        column = line.rstrip('\n').split('\t')
        BC_fw = column[1].split(':')[2][0:10]
        read = column[0]
        BC_rev = column[2].split(':')[2][0:10]
        if BC_fw in dictIND:
            exp_fw = dictIND[BC_fw]
        else:
            exp_fw = 'na'
        if BC_rev in dictIND:
            exp_rev = dictIND[BC_rev]
        else:
            exp_rev = 'na'
        l += 1
        #print(str(BC_rev) + '\t' + exp_rev)
        if l % 100000 == 0:
            print(str(l/1000000) + 'million')
        primer_tul = searchprimer(read,primDICT)
        if primer_tul[0] == 'not found':
            f = 'not found'
            #print('no primer found')
            nf += 1
        else:
            direction_read = primDICT[primer_tul[0]]['direction_read']
            if direction_read == 'FW':
                fw_primer_found +=1
                Amplicon = read[primer_tul[1]:]
                correct_truncated_primer = 'xxxxx'
                Amplicon = correct_truncated_primer[0:primer_tul[2]]+Amplicon
            if direction_read == 'RV':
                rv_primer_found +=1
                Amplicon = read[primer_tul[1]:]
                correct_truncated_primer = 'xxxxx'
                Amplicon = correct_truncated_primer[0:primer_tul[2]]+Amplicon
            guide = primDICT[primer_tul[0]]['guide']
            REFAmplicon = primDICT[primer_tul[0]]['Amplicon']
            REF_len = primDICT[primer_tul[0]]['length']
            cut = primDICT[primer_tul[0]]['cut']
            frame = primDICT[primer_tul[0]]['frame']
            if direction_read == 'RV':
                REFAmplicon = revcomp(REFAmplicon)
                cut = str(int(REF_len)-int(cut))
            #if primer_tul[0][5:15] in Amplicon:    ###this condition must always be met since above is else: (from 'not found')
            outfile.write(guide+'\t'+REFAmplicon+'\t'+Amplicon+'\t'+REF_len+'\t'+cut+'\t'+frame+'\t'+direction_read +'\t'+exp_fw+':'+exp_rev+'\n')
