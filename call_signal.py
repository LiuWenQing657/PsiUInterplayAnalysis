
###############################
# Statistic method claiming
###############################

# cutoff: set p value as cutoff
# deletion ratio report: report treated deletion ratio
# multi T report: report one psiU during multi T, deletion ratio is the biggest one
# site report: report multi T region


###############################
# packages import
###############################

import pandas as pd
import sys
import time
import scipy.stats as stats


###############################
# define functions
###############################

# define a function to calculate total counts and deletion of untreated and treated files and merge as one dataframe
def calc_merg(treated_file_calc, untreated_file_calc):

    # create a dictionary to store result
    result_dic_calc = {}

    # read treated file first
    treated_calc = open(treated_file_calc)
    first_line_calc = True

    for line_calc in treated_calc:
        # do not read the first line
        if first_line_calc:
            first_line_calc = False
            continue

        line_ls_calc = line_calc.split("\t")
        id_calc = line_ls_calc[0] + '_' + line_ls_calc[1]   # id of term
        base_calc = line_ls_calc [2]                        # ref base of term

        # only calculate T site
        if base_calc != 'T':
            continue

        del_count_calc = int(line_ls_calc[8])
        non_del_count_calc = int(line_ls_calc[3]) + int(line_ls_calc[4]) + int(line_ls_calc[5]) + int(line_ls_calc[6])
            # total count of term

        if del_count_calc == 1 and non_del_count_calc == 0:       # situation that site without coverage
            del_count_calc = 0

        result_dic_calc[id_calc] = [del_count_calc, non_del_count_calc, 0, 0]
            # last two 0 is left for untreated del_count and non-del-count

    treated_calc.close()

    # read untreated file then
    untreated_calc = open(untreated_file_calc)
    first_line_calc = True

    for line_calc in untreated_calc:
        # do not read the first line
        if first_line_calc:
            first_line_calc = False
            continue

        line_ls_calc = line_calc.split("\t")
        id_calc = line_ls_calc[0] + '_' + line_ls_calc[1]   # id of term
        base_calc = line_ls_calc [2]                        # ref base of term

        del_count_calc = int(line_ls_calc[8])
        non_del_count_calc = int(line_ls_calc[3]) + int(line_ls_calc[4]) + int(line_ls_calc[5]) + int(line_ls_calc[6])
            # total count of term
        if del_count_calc == 1 and non_del_count_calc == 0:  # situation that site without coverage
            del_count_calc = 0

        # add this information into dictionary if id is in dictionary already
        if result_dic_calc.__contains__(id_calc):
            result_dic_calc[id_calc][2] = del_count_calc
            result_dic_calc[id_calc][3] = non_del_count_calc

    untreated_calc.close()

    # LOG
    sys.stderr.write("\t%s\tFile merge: Done! \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    return result_dic_calc


# define a function to call signal and give a result table
def call_sig(treated_file_call, untreated_file_call, pvalue_cutoff_call, del_ratio_cutoff_call, gene_name_file_call):

    # get calculating result
    dic_tot_call = calc_merg(treated_file_call, untreated_file_call)
        # key. chr_name_index
        # 0. del_count_treated
        # 1. non_del_count_treated
        # 2. del_count_untreated
        # 3. non_del_count_untreated


    # call signal, as continous T
    start_site_call = ''                    # store starting site of T
    end_site_call = ''                      # store ending site of T
    T_bool_call = False                     # store whether there is start T or not
    T_cont_call = False                     # store whether there is continues T or not
    curr_del_treated_call = []              # store current deletion count of treated
    curr_non_del_treated_call = []          # store current non deletion count of treated
    curr_del_untreated_call = []            # store current deletion count of untreated
    curr_non_del_untreated_call = []        # store current non deletion count of treated

    output_chr_name_call = []               # store output chr name
    output_chr_index_call = []              # store output chr index
    output_del_ratio_call = []              # store output del ratio of treated result
    output_tot_count_call = []              # store output total count of treated result
    output_start_call = []                  # store start site, used in .bed file
    output_end_call = []                    # store end site, used in .bed file
    output_treated_del_count_call = []      # store the del count of treated result
    output_untreated_del_ratio_call = []    # store the del ratio of untreated result
    output_untreated_tot_count_call = []    # store the total count of untreated result
    output_untreated_del_count_call = []    # store the del count of untreated result
    output_pvalue_call = []                 # store the pvalue of sites

    last_key_call = 'NR_100000000'      # store last key


    for key_call in dic_tot_call:

        # give a judge whether T is continues or not
        if int(chr_extract(key_call)[1]) - int(chr_extract(last_key_call)[1]) == 1:
            T_cont_call = True
        else:
            T_cont_call = False

        # starting T
        if T_bool_call == False:
            T_bool_call = True
            start_site_call = key_call
            curr_del_treated_call.append(dic_tot_call[key_call][0])
            curr_non_del_treated_call.append(dic_tot_call[key_call][1])
            curr_del_untreated_call.append(dic_tot_call[key_call][2])
            curr_non_del_untreated_call.append(dic_tot_call[key_call][3])
        # extending T
        elif T_bool_call == True and T_cont_call == True:
            # if last RNA is end by T and next RNA is start by T, there will be problems, but not take into consideration now
            curr_del_treated_call.append(dic_tot_call[key_call][0])
            curr_non_del_treated_call.append(dic_tot_call[key_call][1])
            curr_del_untreated_call.append(dic_tot_call[key_call][2])
            curr_non_del_untreated_call.append(dic_tot_call[key_call][3])
        # not extending T
        elif T_bool_call == True and T_cont_call == False:

            # manipulating last T left
            end_site_call = last_key_call
            max_result_call = get_max(curr_del_treated_call, curr_non_del_treated_call, curr_del_untreated_call, curr_non_del_untreated_call)
                # max_result_call[0] total count of treated
                # max_result_call[1] max del ratio of treated
                # max_result_call[2] p value
                # max_result_call[3] del ratio difference
                # max_result_call[4] total count of untreated
                # max_result_call[5] del ratio of untreated
                # max_result_call[6] del count of treated
                # max_result_call[7] del count of untreated
            # p < 0.1 means there is signal, del ratio difference > cutoff, than call signal
            if max_result_call[2] < pvalue_cutoff_call and max_result_call[3] > del_ratio_cutoff_call:
                output_del_ratio_call.append(max_result_call[1])
                output_tot_count_call.append(max_result_call[0])
                output_chr_name_call.append(chr_extract(start_site_call)[0])
                output_start_call.append(int(chr_extract(start_site_call)[1]) - 1)
                output_end_call.append(int(chr_extract(end_site_call)[1]))
                # index, one T or continous T
                if start_site_call == end_site_call:
                    output_chr_index_call.append(chr_extract(start_site_call)[1])
                else:
                    output_chr_index_call.append(chr_extract(start_site_call)[1] + '-' + chr_extract(end_site_call)[1])
                output_pvalue_call.append(max_result_call[2])
                output_treated_del_count_call.append(max_result_call[6])
                output_untreated_del_count_call.append(max_result_call[7])
                output_untreated_tot_count_call.append(max_result_call[4])
                output_untreated_del_ratio_call.append(max_result_call[5])

            # initialize
            curr_del_treated_call = []
            curr_non_del_treated_call = []
            curr_del_untreated_call = []
            curr_non_del_untreated_call = []

            # manipulating current T
            start_site_call = key_call
            curr_del_treated_call.append(dic_tot_call[key_call][0])
            curr_non_del_treated_call.append(dic_tot_call[key_call][1])
            curr_del_untreated_call.append(dic_tot_call[key_call][2])
            curr_non_del_untreated_call.append(dic_tot_call[key_call][3])

        # store last_key_call
        last_key_call = key_call

    # create output dataframe
    if gene_name_file_call == '':
        output_df_call = pd.DataFrame({'chr_name': output_chr_name_call,
                                       "site": output_chr_index_call,
                                       "treated_total_counts": output_tot_count_call,
                                       "treated_deletion_counts": output_treated_del_count_call,
                                       "treated_deletion_ratio": output_del_ratio_call,
                                       "ctrl_total_counts": output_untreated_tot_count_call,
                                       "ctrl_deletion_counts": output_untreated_del_count_call,
                                       "ctrl_deletion_ratio": output_untreated_del_ratio_call,
                                       "p_value": output_pvalue_call})
    else:
        output_gene_name_call = get_gene_name(output_chr_name_call, gene_name_file_call)
        output_df_call = pd.DataFrame({'chr_name': output_chr_name_call,
                                       "site": output_chr_index_call,
                                       "treated_total_counts": output_tot_count_call,
                                       "treated_deletion_counts": output_treated_del_count_call,
                                       "treated_deletion_ratio": output_del_ratio_call,
                                       "ctrl_total_counts": output_untreated_tot_count_call,
                                       "ctrl_deletion_counts": output_untreated_del_count_call,
                                       "ctrl_deletion_ratio": output_untreated_del_ratio_call,
                                       "p_value": output_pvalue_call,
                                       "gene_name": output_gene_name_call})

    # create .bed file
    output_df_bed_call = pd.DataFrame({'chr_name': output_chr_name_call, "start": output_start_call, "end": output_end_call})

    return output_df_call, output_df_bed_call


# define a function to get gene name of every output id
def get_gene_name(gene_id_ls_ggn, gene_name_file_ggn):

    gene_name_dic_ggn = {}                  # store the gene name information in dictionary

    # get gene name information
    gene_name_ggn = open(gene_name_file_ggn, 'r')

    for line_ggn in gene_name_ggn:
        line_ls_ggn = line_ggn.split('\t')
        gene_name_dic_ggn[line_ls_ggn[0]] = line_ls_ggn[1][:-1]

    # get gene_name output list
    gene_name_ls_ggn = []

    for ii in range(len(gene_id_ls_ggn)):

        gene_id_ggn = gene_id_ls_ggn[ii].split('.')[0]
        try:
            gene_name_ls_ggn.append(gene_name_dic_ggn[gene_id_ggn])
        except:
            gene_name_ls_ggn.append('NONE')

    return gene_name_ls_ggn


# define a function to extract chr name and chr index from chr_name_index
def chr_extract(name_index_chr):

    list_chr = name_index_chr.split('_')

    # get index
    index_chr = list_chr[-1]

    # get name
    name_chr = ''
    for j in range(len(list_chr) - 1):
        name_chr = name_chr + list_chr[j] + '_'

    name_chr = name_chr[:-1]

    return name_chr, index_chr


# define a function to get max difference with its del ratio and total count
def get_max(del_treated_get, non_del_treated_get, del_untreated_get, non_del_untreated_get):

    del_ratio_treated_get = []

    for ii in range(len(del_treated_get)):
        if (del_treated_get[ii] + non_del_treated_get[ii]) == 0:
            del_ratio_treated_get.append(0)
        else:
            del_ratio_treated_get.append(del_treated_get[ii] / (del_treated_get[ii] + non_del_treated_get[ii]))

    max_treated_del_ratio_get = max(del_ratio_treated_get)

    for ii in range(len(del_treated_get)):
        if max_treated_del_ratio_get == del_ratio_treated_get[ii]:
            max_treated_count_get = del_treated_get[ii] + non_del_treated_get[ii]
            max_treated_del_count_get = del_treated_get[ii]
            max_untreated_count_get = del_untreated_get[ii] + non_del_untreated_get[ii]
            max_untreated_del_count_get = del_untreated_get[ii]

            if del_treated_get[ii] == 0 and non_del_treated_get[ii] == 0:       # treated with no coverage, cannot do test, return p = 1
                pvalue_get = 1
            elif del_untreated_get[ii] == 0  and non_del_untreated_get[ii] ==0:  # untreated with no coverage, cannot do test, return p = 1
                pvalue_get = 1
            else:
                obs_matrix_get = [[del_treated_get[ii], del_untreated_get[ii]],[non_del_treated_get[ii], non_del_untreated_get[ii]]]
                tot_sample_get = del_treated_get[ii] + del_untreated_get[ii] + non_del_treated_get[ii] + non_del_untreated_get[ii]
                min_num_get = min(del_treated_get[ii], del_untreated_get[ii], non_del_treated_get[ii], non_del_untreated_get[ii])
                # if total sample < 40 or min sample number <1, do fisher test
                if tot_sample_get < 40 or min_num_get < 1:
                    oddsratio_get, pvalue_get = stats.fisher_exact(obs_matrix_get, alternative = 'greater')
                # if total samples >= 40 and min sample number >=5, do chi2 test
                elif tot_sample_get >= 40 and min_num_get >=5:
                    chi2_get, pvalue_get, dof_get, expected_get = stats.chi2_contingency(obs_matrix_get, correction = False)
                # if total samples >= 40 and min sample <5 and >=1, do adjusted chi2 test
                else:
                    chi2_get, pvalue_get, dof_get, expected_get = stats.chi2_contingency(obs_matrix_get, correction = True)

            if (del_untreated_get[ii] + non_del_untreated_get[ii]) == 0:
                max_del_ratio_dif_get = max_treated_del_ratio_get
                max_untreated_del_ratio_get = 0
            else:
                max_untreated_del_ratio_get = del_untreated_get[ii] / (del_untreated_get[ii] + non_del_untreated_get[ii])
                max_del_ratio_dif_get = max_treated_del_ratio_get - max_untreated_del_ratio_get

            break

    return max_treated_count_get, max_treated_del_ratio_get, pvalue_get, max_del_ratio_dif_get, max_untreated_count_get, max_untreated_del_ratio_get, max_treated_del_count_get, max_untreated_del_count_get


# Convert seconds into days-hours-minutes-seconds format
def convert_seconds(tot_seconds):
    out_seconds = tot_seconds % 60
    tot_minutes = int(tot_seconds // 60)
    out_minutes = int(tot_minutes % 60)
    tot_hours = int(tot_minutes // 60)
    out_hours = int(tot_hours % 24)
    tot_days = int(tot_hours // 24)

    return ("%d days  %d hours  %d minutes  %.1f seconds." %(tot_days, out_hours, out_minutes, out_seconds))


###############################
# main program
###############################

if __name__ == "__main__":

    start_time = time.time()    # for the use of stat total time

    ###########################
    # input
    ###########################

    treated_file = ""                                       # input treated file
    untreated_file = ""                                     # input untreated file
    pvalue_cutoff = 0.05                                    # cutoff of pvalue, default = 0.1
    del_ratio_cutoff = 0.05                                 # cutoff of reads count, default = 10
    output_stat_file = ""                                   # output statistics result
    output_bed_file = ""                                    # output bed file
    gene_name_file = ""                                     # store the ID and gene name information

    for i in range(len(sys.argv)):
        input_para = sys.argv[i]
        if input_para == '-i':
            treated_file = sys.argv[i + 1]
        elif input_para == '-u':
            untreated_file = sys.argv[i + 1]
        elif input_para == '-o':
            output_stat_file = sys.argv[i + 1]
        elif input_para == '-b':
            output_bed_file = sys.argv[i + 1]
        elif input_para == '-pc':
            pvalue_cutoff = float(sys.argv[i + 1])
        elif input_para == '-rc':
            del_ratio_cutoff = int(sys.argv[i + 1])
        elif input_para == '-x':
            gene_name_file = sys.argv[i + 1]
        elif input_para == '-h' or input_para == '--help':
            print(
                '--------------------------------------------------Welcome to use PBS-seq signal calling--------------------------------------------------')
            print('Usage:')
            print('    python call-signal_2G_ver8.py [Options]* -u <untreated_in> -i <treated_in> -o <result_out>')
            print('')
            print('    <untreated_in>  Input untreated parse result (.bmat file)')
            print('    <treated_in>    Input treated parse result (.bmat file)')
            print('    <result_out>    Output signal calling result (.csv file)')
            print('')
            print('[Options]')
            print('')
            print('  Gene name:')
            print('    -x <gene_name>  import the ID and gene name information if want to get gene name information in the output file')
            print('')
            print('  Cutoff:')
            print('    -pc [FLOAT]     Pvalue cutoff, < [FLOAT] will be called, default =', pvalue_cutoff)
            print('    -rc [FLOAT]     Deletion ratio cutoff, difference between treated and untreated > [FLOAT] will be called, default =', del_ratio_cutoff)
            print('')
            print('  Output:')
            print('    -b              Output .bed file of called sites, default: do not output bed file')
            print('')
            print('  Other:')
            print('    -h / --help     help')
            exit(0)

    ###########################
    # processing & output
    ###########################

    # LOG
    sys.stderr.write("%s\tSignal Calling: Starts \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    call_sig_result = call_sig(treated_file, untreated_file, pvalue_cutoff, del_ratio_cutoff, gene_name_file)
        # call_sig_result[0]: statistic result
        # call_sig_result[1]: bed file

    # output
    call_sig_result[0].to_csv(output_stat_file, sep=',', index = False)
    if output_bed_file != "":
        call_sig_result[1].to_csv(output_bed_file, sep="\t", header = None, index = False)

    #LOG
    end_time = time.time()
    total_time = end_time - start_time

    sys.stderr.write("\n%s\tSignal Calling: Done! \n"
                     "    Total time used: %s \n"
                     % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), convert_seconds(total_time)))
