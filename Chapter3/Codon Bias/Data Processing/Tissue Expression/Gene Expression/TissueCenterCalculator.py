import pandas as pd
import numpy as np

cb = pd.read_csv("HumanTEcodonbias.csv", index_col = 1)
te = pd.read_csv("TE.csv")

ele = pd.read_csv("ElevatedGenes.csv")
enh = pd.read_csv("EnhancedGenes.csv")
enr = pd.read_csv("EnrichedGenes.csv")

tissues = list(te.columns)[1:]

sum_table = pd.DataFrame(columns = ['Tissue',
                                    'Elevated Transcripts',
                                    'Enhanced Transcripts',
                                    'Enriched Transcripts'])

ele_centers = pd.DataFrame(columns = cb.columns[2:])
enh_centers = pd.DataFrame(columns = cb.columns[2:])
enr_centers = pd.DataFrame(columns = cb.columns[2:])

for tis in tissues:
    ele_sub = cb[cb['Gene ID'].isin(ele['Gene ID'].where(ele['Tissue'] == tis))]
    enh_sub = cb[cb['Gene ID'].isin(enh['Gene ID'].where(enh['Tissue'] == tis))]
    enr_sub = cb[cb['Gene ID'].isin(enr['Gene ID'].where(enr['Tissue'] == tis))]

    sum_table = sum_table.append({'Tissue': tis,
                                  'Elevated Transcripts':len(ele_sub),
                                  'Enhanced Transcripts':len(enh_sub),
                                  'Enriched Transcripts':len(enr_sub)}, ignore_index = True)

    if len(ele_sub) != 0:
        ele_centers = ele_centers.append(ele_sub.mean(), ignore_index = True)
    if len(enh_sub) != 0:
        enh_centers = enh_centers.append(enh_sub.mean(), ignore_index = True)
    if len(enr_sub) != 0:
        enr_centers = enr_centers.append(enr_sub.mean(), ignore_index = True)

ele_centers.insert(0, column = 'Tissue', value = sum_table['Tissue'].where(sum_table['Elevated Transcripts'] != 0).dropna().reset_index(drop=True))
enh_centers.insert(0, column = 'Tissue', value = sum_table['Tissue'].where(sum_table['Enhanced Transcripts'] != 0).dropna().reset_index(drop=True))
enr_centers.insert(0, column = 'Tissue', value = sum_table['Tissue'].where(sum_table['Enriched Transcripts'] != 0).dropna().reset_index(drop=True))
    
ele_centers.to_csv("ElevatedCodonBiasCenters.csv", index = False)
enh_centers.to_csv("EnhancedCodonBiasCenters.csv", index = False)
enr_centers.to_csv("EnrichedCodonBiasCenters.csv", index = False)

sum_table.to_csv("TESummaryTranscriptTable.csv", index = False)
