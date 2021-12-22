#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 12 09:55:45 2021
The Accuracy function requires pandas and the Get_Header function.
It takes the coriell ID as a string such as 'HG00445', the coriell's
tabix vcf and the thermo run vcf for that coriell.
Because the thermo vcfs use '/' instead of '|' to delineate GTs, we
require two dicts for the concordance part of the function.
A script using this function can then iterate through a dict of coriells with
their vcfs as tupules as shown below.

@author: howardwetsman
"""
import os
import time
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'


def Get_Header(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        count = 0
        for i in range(len(lines)):
            if '##' in lines[i]:
                count = count+1
    return count


def Accuracy(coriell, tabix_df, run_df):
    print(f'Working {coriell}')
    tabix_df = tabix_df.rename(columns={coriell: 'EXPECTED', 'ALT': 'T_ALT'})
    tabix_df[['#CHROM', 'POS', 'ID', 'REF', 'T_ALT', 'EXPECTED']]
    tabix_df['EXPECTED'] = tabix_df['EXPECTED'].map(expected_replace_dict)
    tabix_df.POS = tabix_df.POS.astype(str)
    tabix_pos = tabix_df[(tabix_df.EXPECTED != 'WT')]
    tabix_neg = tabix_df[(tabix_df.EXPECTED == 'WT')]

    run_df['#CHROM'] = run_df['#CHROM'].str[3:]
    run_df['FOUND'] = run_df[coriell].str[0:3]
    run_df = run_df.rename(columns={'ALT': 'R_ALT'})
    run_df = run_df[['#CHROM', 'POS', 'REF', 'R_ALT', 'FOUND']]
    run_df.FOUND = run_df.FOUND.map(found_replace_dict)
    run_df.POS = run_df.POS.astype(str)
    compare_pos = pd.merge(tabix_pos, run_df, how='outer', on=['#CHROM', 'POS'])  # ,'REF'])
    compare_pos = compare_pos[~compare_pos['EXPECTED'].isnull()]
    compare_neg = pd.merge(tabix_neg, run_df, how='outer', on=['#CHROM', 'POS'])  # ,'REF'])
    compare_neg = compare_neg[compare_neg.EXPECTED == 'WT']

    con_TP = compare_pos[compare_pos.EXPECTED == compare_pos.FOUND]
    con_tp = con_TP.shape[0]
    print('Concordant_True_Positives:', con_tp)
    call_TP = compare_pos[~compare_pos.FOUND.isnull()]
    call_tp = call_TP.shape[0]
    print('Call_True_Positives:', call_tp)

    combined_TP = call_TP.append(con_TP)
    combined_TP.sort_values(['#CHROM', 'POS']).drop_duplicates(subset=['#CHROM', 'POS'], keep=False)
    examine_TP = combined_TP[combined_TP['EXPECTED'] != combined_TP['FOUND']]

    FN = compare_pos[compare_pos.FOUND.isnull()]
    fn = FN.shape[0]
    print('False_Negatives:', fn)

    compare_neg = pd.merge(tabix_neg, run_df, how='outer', on=['#CHROM', 'POS'])  # ,'REF'])
    compare_neg = compare_neg[compare_neg.EXPECTED == 'WT']
    compare_neg.to_csv('compare_neg.csv')
    TN = compare_neg[compare_neg.FOUND.isnull()]
    tn = TN.shape[0]
    print('True_Negatives:', tn)

    FP = compare_neg[~compare_neg.FOUND.isnull()]
    fp = FP.shape[0]
    print('False_Positives:', fp)

    call_PPA = call_tp/(call_tp + fn)
    con_PPA = con_tp/(con_tp+fn)
    call_PPV = call_tp/(call_tp+fp)
    con_PPV = con_tp/(con_tp+fp)

    print('Call_PPA: ', call_PPA)
    print('Con_PPA: ', con_PPA)
    print('Call_PPV: ', call_PPV)
    print('Con_PPV: ', con_PPV)

    print(f'Printing output to Excel {path}{coriell}_for_Validation.xlsx')
    toc_content = {1: 'The sheets of this workbook contain all data found. Empty spaces should be interpretted that nothing was found.',
                   2: 'True Positives were those variants expected that were found in the run.',
                   3: 'False Positives were those variants seen in the run that were not expected. An empty sheet is indicative of no false positives.',
                   4: 'True Negatives were those variants for which the control was WT and were not seen in the run. It is normal to have empty space in the FOUND column.',
                   5: 'False Negatives were those variants for which the variant was not WT and were not seen in the run. The FOUND column is empty by definition.',
                   6: 'Concordance shows HET rather than specific genotype to avoid confounding with false phasing.',
                   7: '',
                   8: f'True Negatives: {tn}',
                   9: f'False Negatives: {fn}',
                   10: f'Call True Positives: {call_tp}',
                   11: f'Con True Positives: {con_tp}',
                   12: f'False Positives: {fp}',
                   13: f'Call PPA: {round(100*call_PPA,2)}%',
                   14: f'Call PPV: {round(100*call_PPV,2)}%',
                   15: f'Concordant PPA: {round(100*con_PPA,2)}%',
                   16: f'Concordant PPV: {round(100*con_PPV,2)}%'}
    TOC_df = pd.DataFrame.from_dict(toc_content, orient='index')
    with pd.ExcelWriter(f'{coriell}_for_Validation.xlsx') as writer:
        TOC_df.to_excel(writer, sheet_name='Summary - How_to_use_this_file', index=False)
        call_TP.to_excel(writer, sheet_name='Call_True_Positives', index=False)
        con_TP.to_excel(writer, sheet_name='Concordant_True_Positives', index=False)
        examine_TP.to_excel(writer, sheet_name='TPs_to_be_examined', index=False)
        FP.to_excel(writer, sheet_name='False_Positives', index=False)
        TN.to_excel(writer, sheet_name='True_Negatives', index=False)
        FN.to_excel(writer, sheet_name='False_Negatives', index=False)

    return call_tp, con_tp, fp, tn, fn, FN


# main
time0 = time.time()
client = "Lab"
panel = 'Cardio'
path = f'./{client}/{panel}/'


def Get_Samples(client, panel):
    files = os.listdir(f'./{client}/{panel}')
    pos_cont_files = [x for x in files if 'Pos_Control' in x]
    samples = [x.strip('_Pos_Control.vcf') for x in pos_cont_files]
    samples = [x.strip(f'{client}_{panel}_') for x in samples]
    return samples


def Get_Sample_Files(sample, client, panel):
    files = os.listdir(f'./{client}/{panel}')
    sample_files = [x for x in files if sample in x]
    return sample_files


def Make_File_Dict(client, panel):
    dict1 = {}
    files = os.listdir(f'./{client}/{panel}')
    pos_cont_files = [x for x in files if 'Pos_Control' in x]
    samples = Get_Samples(client, panel)
    print(pos_cont_files)
    run_files = [x for x in files if 'Pos_Control' not in x]
    print(run_files)
    for coriell_file in pos_cont_files:
        coriell = coriell_file.strip(f'{client}_{panel}_')
        coriell = coriell.strip('_Pos_Control.vcf')
        coriell_run_files = [x for x in run_files if coriell in x]
        print(coriell_run_files[0])
        dict1[coriell] = coriell_run_files[0]
    print('\n', dict1)
    return dict1


expected_replace_dict = {'0|1': 'HET', '1|0': 'HET', '1': 'HOM', '1|1': 'HOM',
                         '2|2': 'HOM2', '3|3': 'HOM3', '4|4': 'HOM4', '5|5': 'HOM5',
                         '0|2': 'HET2', '2|0': 'HET2', '0|3': 'HET3', '3|0': 'HET3',
                         '0|4': 'HET4', '4|0': 'HET4', '0|5': 'HET5', '5|0': 'HET5',
                         '1|2': 'CH12', '2|1': 'CH12', '1|3': 'CH13', '3|1': 'CH13',
                         '1|4': 'CH14', '4|1': 'CH14', '1|5': 'CH15', '5|1': 'CH15',
                         '2|3': 'CH23', '3|2': 'CH23', '2|4': 'CH24', '4|2': 'CH24',
                         '2|5': 'CH25', '5|2': 'CH25', '3|4': 'CH34', '4|3': 'CH34',
                         '3|5': 'CH35', '5|3': 'CH35', '4|5': 'CH45', '5|4': 'CH45',
                         '0': 'WT', '0|0': 'WT'}
found_replace_dict = {'0/1': 'HET', '1/0': 'HET', '1': 'HOM', '1/1': 'HOM',
                      '2/2': 'HOM2', '3/3': 'HOM3', '4/4': 'HOM4', '5/5': 'HOM5',
                      '0/2': 'HET2', '2/0': 'HET2', '0/3': 'HET3', '3/0': 'HET3',
                      '0/4': 'HET4', '4/0': 'HET4', '0/5': 'HET5', '5/0': 'HET5',
                      '1/2': 'CH12', '2/1': 'CH12', '1/3': 'CH13', '3/1': 'CH13',
                      '1/4': 'CH14', '4/1': 'CH14', '1/5': 'CH15', '5/1': 'CH15',
                      '2/3': 'CH23', '3/2': 'CH23', '2/4': 'CH24', '4/2': 'CH24',
                      '2/5': 'CH25', '5/2': 'CH25', '3/4': 'CH34', '4/3': 'CH34',
                      '3/5': 'CH35', '5/3': 'CH35', '4/5': 'CH45', '5/4': 'CH45',
                      '0': 'WT', '0/0': 'WT'}
# file_dict = {'HG03706': 'HG03706_Accuracy_1_2021-06-15_17.vcf',
#              'HG02476': 'HG02476_Accuracy_1_2021-06-15_17.vcf',
#              'NA19652': 'NA19652_Accuracy_1_2021-06-15_17.vcf',
#              'HG01941': 'HG01941_Accuracy_1_2021-06-15_17.vcf'}

nuc_chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
              '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
              '21', '22', 'X', 'Y']
c_call_tp = c_con_tp = c_fp = c_tn = c_fn = 0
total_FN_df = pd.DataFrame()

# eventual data dict
# {'sample':{'control':control_vcf},{'accuracy1':accuracy1_vcf},{'accuracy2':accuracy2_vcf}    }


def Get_Date_Run(file):
    datetime = file.split('.')[0].split(f'{client}_{panel}_')[1].split('_', 1)[1]
    date, run = datetime.split('_')
    return date, run


def Get_Accuracy1(sample_df):
    if sample_df.shape[0] == 1:
        return sample_df.index[0]
    else:
        df = sample_df[sample_df.date == sample_df.date.min()]
        if df.shape[0] == 1:
            file = df.index[0]
            # sample_dict['accuracy1'] = file
            print('\n', file)
        else:
            df1 = df[df.run == sd.run.min()]
            file = df1.index[0]
            # sample_dict['accuracy1'] = file
    return file


def Get_Precision(sample_df, accuracy1):
    sample_df.drop(accuracy1, inplace=True, axis=0)
    print(sample_df)
    if sample_df.shape[0] == 1:
        return sample_df.index[0]
    else:
        df = sample_df[sample_df.date == sample_df.date.min()]
        if df.shape[0] == 1:
            return df.index[0]

        else:
            df1 = df[df.run == sd.run.min()]
            return df1.index[0]


def Make_Sample_Dict(sample_files):
    sample_dict = {}
    sample_df = pd.DataFrame()
    pos_cont = [x for x in sample_files if 'Pos_Control' in x]
    sample_dict['control'] = pos_cont[0]
    run_files = [x for x in sample_files if 'Pos_Control' not in x]
    for file in run_files:
        print(file)
        files_dict = {}
        date, run = Get_Date_Run(file)
        date = pd.to_datetime(date, format='%Y-%m-%d')
        print(date)
        run = int(run)
        print(run)
        sample_df.loc[file, 'date'] = date
        sample_df.loc[file, 'run'] = run

        print('\n', sample_df)
    # get accuracy1
    accuracy1 = Get_Accuracy1(sample_df)
    sample_dict['accuracy1'] = accuracy1
    precision = Get_Precision(sample_df, accuracy1)
    sample_dict['precision'] = precision
    print(sample_dict)

    # iterate run files and order them by datetime
    print(run_files)


###############
#Program Start#
###############
# define list of samples
samples = Get_Samples(client, panel)
# iterate list of smaples and create sample dicts
for sample in samples:
    sample_files = Get_Sample_Files(sample, client, panel)
    print('\n', sample_files)
    # once have sample files go through them to create sample dict
    sample_dict = Make_Sample_Dict(sample_files)


file_dict = Make_File_Dict(client, panel)
print('\n', file_dict)


for key, value in file_dict.items():
    print(key, value)
    coriell = key
    tabix_coriell = f'{client}_{panel}_{coriell}_Pos_Control.vcf'
    run_vcf = f'{path}{value}'
    tabix_header = Get_Header(f'{path}{tabix_coriell}')
    tabix_df = pd.read_csv(f'{path}{tabix_coriell}', header=tabix_header, sep='\t')
    run_header = Get_Header(run_vcf)
    run_df = pd.read_csv(run_vcf, header=run_header, sep='\t')
    temp_df = run_df
    temp_df['SAMPLE'] = coriell
    temp_df = temp_df.rename(columns={coriell: 'EXPECTED'})
    # Steve_df = Steve_df.append(temp_df)

    call_tp, con_tp, fp, tn, fn, temp_df = Accuracy(coriell, tabix_df, run_df)

    total_FN_df = total_FN_df.append(temp_df)
    c_call_tp = c_call_tp + call_tp
    c_con_tp = c_con_tp + con_tp
    c_fp = c_fp + fp
    c_tn = c_tn + tn
    c_fn = c_fn + fn
    call_PPA = c_call_tp/(c_call_tp+c_fn)
    call_PPV = c_call_tp/(c_call_tp+c_fp)
    con_PPA = (c_con_tp/(c_con_tp+fn))
    con_PPV = (c_con_tp/(c_con_tp+fp))

    print('True Negative:', tn)
    print('False Negative:', fn)
    print('Call True Positive:', call_tp)
    print('Concordant True Positive:', con_tp)
    print('False Positive:', fp)
    print('Call PPA:', f'{round(100*call_PPA,2)}%')
    print('Call PPV:', f'{round(100*call_PPV,2)}%')
    print('Concordant PPA:', f'{round(100*con_PPA,2)}%')
    print('Concordant PPV:', f'{round(100*con_PPV,2)}%')
    print(total_FN_df.shape)
    print()
duplicate_df = pd.DataFrame()
duplicates = total_FN_df.pivot_table(index=['#CHROM', 'POS'], aggfunc='size')
indexlist = duplicates.index.tolist()
l = len(indexlist)
for i in range(l):
    t = indexlist[i]
    c, p = t[0], int(t[1])
    duplicate_df.loc[i, '#CHROM'] = c
    duplicate_df.loc[i, 'POS'] = p
    duplicate_df.loc[i, 'COUNT'] = int(duplicates[i])

duplicate_df['#CHROM'] = pd.Categorical(duplicate_df['#CHROM'],
                                        categories=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                                                    '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
                                                    '21', '22', 'X', 'Y'], ordered=True)
duplicate_df.sort_values(['#CHROM', 'POS'], inplace=True)
duplicate_df['COUNT'] = duplicate_df['COUNT'].astype(int)
duplicate_df['POS'] = duplicate_df['POS'].astype(int)
duplicate_df['POS'] = duplicate_df['POS'].astype(str)


files_df = pd.DataFrame.from_dict(file_dict, orient='index', columns=['File'])
summary_df = pd.DataFrame()
print('For run: ')
summary_df.loc['True Negative:', 'Accuracy_Measure'] = c_tn
summary_df.loc['False Negative:', 'Accuracy_Measure'] = c_fn
summary_df.loc['Call True Positive:', 'Accuracy_Measure'] = c_call_tp
summary_df.loc['Con True Positive:', 'Accuracy_Measure'] = c_con_tp
summary_df.loc['False Positive:', 'Accuracy_Measure'] = c_fp
summary_df.loc['Call PPA:', 'Accuracy_Measure'] = f'{round(100*call_PPA,2)}%'
summary_df.loc['Call PPV:', 'Accuracy_Measure'] = f'{round(100*call_PPV,2)}%'
summary_df.loc['Concordant PPA:', 'Accuracy_Measure'] = f'{round(100*con_PPA,2)}%'
summary_df.loc['Concordant PPV:', 'Accuracy_Measure'] = f'{round(100*con_PPV,2)}%'
print()
print('Writing last file for client')
with pd.ExcelWriter(f'{client}_{panel}_Summary_For_Accuracy.xlsx') as writer:
    duplicate_df.to_excel(writer, sheet_name='False_Negative_Run_Duplications', index=False)
    files_df.to_excel(writer, sheet_name='Files_in_Run')
    summary_df.to_excel(writer, sheet_name='Summary_Statistics')
print('Done')
print()


time1 = time.time()
print(f'Total = {time1-time0}')


#
