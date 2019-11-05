__author__ = 'georg.michlits'

dict_file = open('v3_Ms_indeldict.txt','r')
output_file = open('Ms_summary_analysis_filtered_indels_V3.txt','w')
#output_file = open('Ms_indel_analysis_txt','w')

threshold = 0.005

mut_dict = {}
total_count = {}
for i,line in enumerate(dict_file):
    if i > 0:
        column = line.rstrip('\n').split('\t')
        guide = column[0]
        exp_rev = column[1]
        exp_fw = column[2]
        mutation = column[3]
        count = int(column[4])
        if not exp_rev in mut_dict:

            mut_dict[exp_rev] = {}
            total_count[exp_rev] = {}
        if not guide in mut_dict[exp_rev]:

            mut_dict[exp_rev][guide] = {}
            total_count[exp_rev][guide] = {}
        if not exp_fw in mut_dict[exp_rev][guide]:

            mut_dict[exp_rev][guide][exp_fw] = {}
            total_count[exp_rev][guide][exp_fw] = 0
        total_count[exp_rev][guide][exp_fw] += count
        mut_dict[exp_rev][guide][exp_fw][mutation] = count

output_file.write('experiment\tguide\twt_d10+4OH\tSub_d10+4OH\tins_d10+4OH\tdel_d10+4OH\tNA_d10+4OH'
                  + '\twt_d10-4OH\tSub_d10-4OH\tins_d10-4OH\tdel_d10-4OH\tNA_d10-4OH'
                  + '\twt_d2+4OH\tSub_d2+4OH\tins_d2+4OH\tdel_d2+4OH\tNA_d2+4OH')

mut_classes = ['wt','Sub','ins','del','NA']
for exp_rev in sorted(mut_dict):
    for guide in sorted(mut_dict[exp_rev]):
        output_file.write('\n' + exp_rev + '\t' + guide)
        for exp_fw in sorted(mut_dict[exp_rev][guide]):
            count_dict = {'wt':0, 'Sub':0, 'ins':0, 'del':0, 'NA':0}
            total_reads_sample = total_count[exp_rev][guide][exp_fw]
            for mutation in mut_dict[exp_rev][guide][exp_fw]:
                if float(mutation.split(':')[2]) < 21 and float(mutation.split(':')[2]) > -21:
                    for element in count_dict:
                        if element in mutation:
                            count = mut_dict[exp_rev][guide][exp_fw][mutation]
                            if count/total_reads_sample >= threshold:
                                if element == 'wt':
                                    count_dict[element] += count
                                elif element == 'Sub':
                                    mutation_size = int(mutation.split(':')[1])
                                    mutation_pos = int(mutation.split(':')[2])
                                    if 5 > mutation_pos > -5:
                                        count_dict[element] += count
                                elif element == 'ins':
                                    mutation_size = int(mutation.split(':')[1])
                                    mutation_pos = int(mutation.split(':')[2])
                                    if -5 < mutation_pos < 5:
                                        count_dict[element] += count
                                elif element == 'del':
                                    add_to_count = 'no'
                                    mutation_size = int(mutation.split(':')[1])
                                    mutation_pos = int(mutation.split(':')[2])
                                    if mutation_pos < 0:
                                        if mutation_pos + mutation_size >= 0:
                                            add_to_count = 'yes'
                                        rel_pos = mutation_pos+mutation_size
                                    elif mutation_pos >= 0:
                                        rel_pos = mutation_pos
                                    if add_to_count == 'yes':
                                        count_dict[element] += count
                                    else:
                                        if -5 < rel_pos < 5:
                                            count_dict[element] += count
                                elif element == 'NA':
                                    count_dict[element] += count
            for mut_cat in mut_classes:
                sum_count = count_dict[mut_cat]
                output_file.write('\t' + str(sum_count))