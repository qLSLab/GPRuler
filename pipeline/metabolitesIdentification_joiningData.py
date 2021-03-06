import os
import ast
import sys
import pandas as pd
import genericLib as gL

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]

# setting input data
testModel = sys.argv[1]
if testModel == 'recon':
    ## Recon 3
    outputFileName = 'recon3'
    dfMetsInfer_classicMethod = 'recon3D_metabolites_wInferredIds'
    dfMetsInfer_fw_metacyc1 = 'recon3_mappingMetaCyc_empty'
    dfMetsInfer_fw_metacyc2 = 'recon3_mappingMetaCyc_100'
    dfMetsInfer_fw_metacyc3 = 'recon3_mappingMetaCyc_91_99'

    dfMetsInfer_fw_KC1 = 'recon3_mappingKeggC_empty'
    dfMetsInfer_fw_KC2 = 'recon3_mappingKeggC_100'
    dfMetsInfer_fw_KC3 = 'recon3_mappingKeggC_91_99'

    dfMetsInfer_fw_KG1 = 'recon3_mappingKeggG_empty'
    dfMetsInfer_fw_KG2 = 'recon3_mappingKeggG_100'
    dfMetsInfer_fw_KG3 = 'recon3_mappingKeggG_91_99'

    dfMetsInfer_fw_C1 = 'recon3_mappingChebi_empty'
    dfMetsInfer_fw_C2 = 'recon3_mappingChebi_100'
    dfMetsInfer_fw_C3 = 'recon3_mappingChebi_91_99'

elif testModel == 'y7':
    ## Yeast 7
    outputFileName = 'yeast7'
    dfMetsInfer_classicMethod = 'y7_metabolites_wInferredIds'

    dfMetsInfer_fw_metacyc1 = 'yeast7_mappingMetaCyc_empty'
    dfMetsInfer_fw_metacyc2 = 'yeast7_mappingMetaCyc_100'
    dfMetsInfer_fw_metacyc3 = 'yeast7_mappingMetaCyc_91_99'

    dfMetsInfer_fw_KC1 = 'yeast7_mappingKeggC_empty'
    dfMetsInfer_fw_KC2 = 'yeast7_mappingKeggC_100'
    dfMetsInfer_fw_KC3 = 'yeast7_mappingKeggC_91_99'

    dfMetsInfer_fw_KG1 = 'yeast7_mappingKeggG_empty'
    dfMetsInfer_fw_KG2 = 'yeast7_mappingKeggG_100'
    dfMetsInfer_fw_KG3 = 'yeast7_mappingKeggG_91_99'

    dfMetsInfer_fw_C1 = 'yeast7_mappingChebi_empty'
    dfMetsInfer_fw_C2 = 'yeast7_mappingChebi_100'
    dfMetsInfer_fw_C3 = 'yeast7_mappingChebi_91_99'
elif testModel == 'y8':
    ## Yeast 8
    outputFileName = 'yeast8'
    dfMetsInfer_classicMethod = 'y8_metabolites_wInferredIds'
    dfMetsInfer_fw_metacyc1 = 'yeast8_mappingMetaCyc_empty'
    dfMetsInfer_fw_metacyc2 = 'yeast8_mappingMetaCyc_100'
    dfMetsInfer_fw_metacyc3 = 'yeast8_mappingMetaCyc_91_99'

    dfMetsInfer_fw_KC1 = 'yeast8_mappingKeggC_empty'
    dfMetsInfer_fw_KC2 = 'yeast8_mappingKeggC_100'
    dfMetsInfer_fw_KC3 = 'yeast8_mappingKeggC_91_99'

    dfMetsInfer_fw_KG1 = 'yeast8_mappingKeggG_empty'
    dfMetsInfer_fw_KG2 = 'yeast8_mappingKeggG_100'
    dfMetsInfer_fw_KG3 = 'yeast8_mappingKeggG_91_99'

    dfMetsInfer_fw_C1 = 'yeast8_mappingChebi_empty'
    dfMetsInfer_fw_C2 = 'yeast8_mappingChebi_100'
    dfMetsInfer_fw_C3 = 'yeast8_mappingChebi_91_99'

elif testModel == 'hmr':
    ## HMRcore
    outputFileName = 'hmrCore'
    dfMetsInfer_classicMethod = 'hmrCore_metabolites_wInferredIds'

    dfMetsInfer_fw_metacyc1 = 'hmrCore_mappingMetaCyc_empty'
    dfMetsInfer_fw_metacyc2 = 'hmrCore_mappingMetaCyc_100'
    dfMetsInfer_fw_metacyc3 = 'hmrCore_mappingMetaCyc_91_99'

    dfMetsInfer_fw_KC1 = 'hmrCore_mappingKeggC_empty'
    dfMetsInfer_fw_KC2 = 'hmrCore_mappingKeggC_100'
    dfMetsInfer_fw_KC3 = 'hmrCore_mappingKeggC_91_99'

    dfMetsInfer_fw_KG1 = 'hmrCore_mappingKeggG_empty'
    dfMetsInfer_fw_KG2 = 'hmrCore_mappingKeggG_100'
    dfMetsInfer_fw_KG3 = 'hmrCore_mappingKeggG_91_99'

    dfMetsInfer_fw_C1 = 'hmrCore_mappingChebi_empty'
    dfMetsInfer_fw_C2 = 'hmrCore_mappingChebi_100'
    dfMetsInfer_fw_C3 = 'hmrCore_mappingChebi_91_99'

elif testModel == 'ownData':
    ## specify your input data
    outputFileName = ''
    dfMetsInfer_classicMethod = ''

    dfMetsInfer_fw_metacyc1 = ''
    dfMetsInfer_fw_metacyc2 = ''
    dfMetsInfer_fw_metacyc3 = ''

    dfMetsInfer_fw_KC1 = ''
    dfMetsInfer_fw_KC2 = ''
    dfMetsInfer_fw_KC3 = ''

    dfMetsInfer_fw_KG1 = ''
    dfMetsInfer_fw_KG2 = ''
    dfMetsInfer_fw_KG3 = ''

    dfMetsInfer_fw_C1 = ''
    dfMetsInfer_fw_C2 = ''
    dfMetsInfer_fw_C3 = ''


dMetMapping = {}

# joining outputs from MetaCyc
dfMetaCyc = pd.read_csv(os.path.join(RAWDIR,'metacyc_compounds_20201216152513.csv'), sep='\t', dtype = str)

df1 = pd.read_csv(os.path.join(OUTDIR, dfMetsInfer_fw_metacyc1 + '.tsv'), sep='\t', index_col=0)
df1['Matches'] = df1['Matches'].apply(ast.literal_eval)
df2 = pd.read_csv(os.path.join(OUTDIR, dfMetsInfer_fw_metacyc2 + '.tsv'), sep='\t', index_col=0)
df2['Matches'] = df2['Matches'].apply(ast.literal_eval)
df3 = pd.read_csv(os.path.join(OUTDIR, dfMetsInfer_fw_metacyc3 + '.tsv'), sep='\t', index_col=0)
df3['Matches'] = df3['Matches'].apply(ast.literal_eval)

dfAll = pd.concat([df1,df2,df3])

for row in dfAll.itertuples():
    dfIsolateMatches = dfMetaCyc[dfMetaCyc['Name'].isin(row.Matches)]
    if dfIsolateMatches.empty is True:
        lCorrespondences = []
    else:
        lCorrespondences = list(dfIsolateMatches['Id'].dropna()) + list(dfIsolateMatches['ChebiId'].dropna()) + list(dfIsolateMatches['KeggId'].dropna())
    if row.Name not in dMetMapping:
        dMetMapping[row.Name] = lCorrespondences
    else:
        dMetMapping[row.Name] += lCorrespondences

# joining outputs from KEGG
dfKeggC = pd.read_csv(os.path.join(RAWDIR,'kegg_compounds_20201216150802.csv'), sep='\t', dtype = str)
dfKeggG = pd.read_csv(os.path.join(RAWDIR,'kegg_glycans_20201216150802.csv'), sep='\t', dtype = str)
dfKeggC['Name'] = dfKeggC['Name'].apply(ast.literal_eval)
dfKeggC['ChebiId'] = dfKeggC['ChebiId'].apply(ast.literal_eval)
dfKeggG['Name'] = dfKeggG['Name'].apply(ast.literal_eval)
dfKeggG['ChebiId'] = dfKeggG['ChebiId'].apply(ast.literal_eval)
dfKeggC = dfKeggC.explode('Name')
dfKeggG = dfKeggG.explode('Name')

df1 = pd.read_csv(os.path.join(OUTDIR, dfMetsInfer_fw_KC1 + '.tsv'), sep='\t', index_col=0)
df1['Matches'] = df1['Matches'].apply(ast.literal_eval)
df2 = pd.read_csv(os.path.join(OUTDIR, dfMetsInfer_fw_KC2 + '.tsv'), sep='\t', index_col=0)
df2['Matches'] = df2['Matches'].apply(ast.literal_eval)
df3 = pd.read_csv(os.path.join(OUTDIR, dfMetsInfer_fw_KC3 + '.tsv'), sep='\t', index_col=0)
df3['Matches'] = df3['Matches'].apply(ast.literal_eval)

dfAll = pd.concat([df1,df2,df3])

for row in dfAll.itertuples():
    dfIsolateMatches1 = dfKeggC[dfKeggC['Name'].isin(row.Matches)]
    dfIsolateMatches2 = dfKeggC[dfKeggC['Id'].isin(row.Matches)]
    dfIsolateMatches = pd.concat([dfIsolateMatches1, dfIsolateMatches2])

    if dfIsolateMatches.empty is True:
        lCorrespondences = []
    else:
        lCorrespondences = list(dfIsolateMatches['Id'].dropna()) + [el for l in list(dfIsolateMatches['ChebiId'].dropna()) for el in l]
    if row.Name.strip() not in dMetMapping:
        dMetMapping[row.Name.strip()] = lCorrespondences
    else:
        dMetMapping[row.Name.strip()] += lCorrespondences


df1 = pd.read_csv(os.path.join(OUTDIR, dfMetsInfer_fw_KG1 + '.tsv'), sep='\t', index_col=0)
df1['Matches'] = df1['Matches'].apply(ast.literal_eval)
df2 = pd.read_csv(os.path.join(OUTDIR, dfMetsInfer_fw_KG2 + '.tsv'), sep='\t', index_col=0)
df2['Matches'] = df2['Matches'].apply(ast.literal_eval)
df3 = pd.read_csv(os.path.join(OUTDIR, dfMetsInfer_fw_KG3 + '.tsv'), sep='\t', index_col=0)
df3['Matches'] = df3['Matches'].apply(ast.literal_eval)

dfAll = pd.concat([df1,df2,df3])

for row in dfAll.itertuples():
    dfIsolateMatches1 = dfKeggG[dfKeggG['Name'].isin(row.Matches)]
    dfIsolateMatches2 = dfKeggG[dfKeggG['Id'].isin(row.Matches)]
    dfIsolateMatches = pd.concat([dfIsolateMatches1, dfIsolateMatches2])

    if dfIsolateMatches.empty is True:
        lCorrespondences = []
    else:
        lCorrespondences = list(dfIsolateMatches['Id'].dropna()) + [el for l in list(dfIsolateMatches['ChebiId'].dropna()) for el in l]
    if row.Name.strip() not in dMetMapping:
        dMetMapping[row.Name.strip()] = lCorrespondences
    else:
        dMetMapping[row.Name.strip()] += lCorrespondences


# joining outputs from ChEBI
dfChebiNames = pd.read_csv(os.path.join(RAWDIR, 'chebi_names_20201216.tsv.bz2'), sep = '\t', compression='bz2', dtype=str)
dfChebiUniprot = pd.read_csv(os.path.join(RAWDIR, 'chebi_uniprot_20201216.tsv.bz2'), sep = '\t', header=None, dtype=str, compression='bz2', names=['ID','NAME_wSymbols', 'NAME'])
dfChebiCompounds =  pd.read_csv(os.path.join(RAWDIR, 'chebi_compounds_20201216.tsv.bz2'), sep = '\t', compression='bz2', dtype=str)

df1 = pd.read_csv(os.path.join(OUTDIR, dfMetsInfer_fw_C1 + '.tsv'), sep='\t', index_col=0)
df1['Matches'] = df1['Matches'].apply(ast.literal_eval)
df2 = pd.read_csv(os.path.join(OUTDIR, dfMetsInfer_fw_C2 + '.tsv'), sep='\t', index_col=0)
df2['Matches'] = df2['Matches'].apply(ast.literal_eval)
df3 = pd.read_csv(os.path.join(OUTDIR, dfMetsInfer_fw_C3 + '.tsv'), sep='\t', index_col=0)
df3['Matches'] = df3['Matches'].apply(ast.literal_eval)

dfAll = pd.concat([df1,df2,df3])

for row in dfAll.itertuples():
    dfIsolateMatches = dfChebiNames[dfChebiNames['NAME'].isin(row.Matches)]
    if dfIsolateMatches.empty is True:
        lCorrespondences = []
    else:
        lCorrespondences = list(dfIsolateMatches['COMPOUND_ID'].dropna())
    if row.Name.strip() not in dMetMapping:
        dMetMapping[row.Name.strip()] = lCorrespondences
    else:
        dMetMapping[row.Name.strip()] += lCorrespondences

for row in dfAll.itertuples():
    dfIsolateMatches = dfChebiUniprot[dfChebiUniprot['NAME'].isin(row.Matches)]
    if dfIsolateMatches.empty is True:
        lCorrespondences = []
    else:
        lCorrespondences = list(dfIsolateMatches['ID'].dropna())
    if row.Name.strip() not in dMetMapping:
        dMetMapping[row.Name.strip()] = lCorrespondences
    else:
        dMetMapping[row.Name.strip()] += lCorrespondences

for row in dfAll.itertuples():
    dfIsolateMatches = dfChebiCompounds[dfChebiCompounds['NAME'].isin(row.Matches)]
    if dfIsolateMatches.empty is True:
        lCorrespondences = []
    else:
        lCorrespondences = list(dfIsolateMatches['ID'].dropna())
    if row.Name.strip() not in dMetMapping:
        dMetMapping[row.Name.strip()] = lCorrespondences
    else:
        dMetMapping[row.Name.strip()] += lCorrespondences

dMetMapping_woDuplicates = {}
for met in dMetMapping:
    dMetMapping_woDuplicates[met.strip()] = gL.unique(dMetMapping[met])

dfMatches = pd.DataFrame(dMetMapping_woDuplicates.items(), columns=['Name', 'Identifiers'])
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingFuzzy.tsv'), sep='\t', index = False)


## Joining FuzzyWuzzy and output from Step1
dfFuzzyMatches = pd.read_csv(os.path.join(OUTDIR, outputFileName + '_mappingFuzzy.tsv'), sep='\t')
dfFuzzyMatches['Identifiers_fuzzy'] = dfFuzzyMatches['Identifiers'].apply(ast.literal_eval)

dfClassicMatches = pd.read_csv(os.path.join(OUTDIR, dfMetsInfer_classicMethod + '.csv'), sep='\t')
dfClassicMatches['Identifiers_classic'] = dfClassicMatches['Identifiers'].apply(ast.literal_eval)
dfClassicMatches['Name'] = dfClassicMatches.Name.str.strip()

dfMerge = pd.merge(dfFuzzyMatches, dfClassicMatches, on = 'Name')
dfMerge['Identifiers'] = dfMerge.apply(lambda row: gL.unique(row['Identifiers_fuzzy'] + row['Identifiers_classic']), axis=1)

dfMerge_filter = dfMerge[['Name', 'Identifiers_fuzzy', 'Identifiers_classic', 'Identifiers']]
dfMerge_filter.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingFuzzyAndClassic.tsv'), sep='\t', index = False)
