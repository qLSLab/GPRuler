import pandas as pd
import sys
import genericLib as gL
import os
import time
import ast

start = time.time()

timeStamp = gL.getTimeStamp()

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]
FIGUREDIR = workingDirs[4]
DFDIR = workingDirs[5]

testModel = sys.argv[1]
if testModel == '1':
    ## Recon 3
    outputFileName = 'recon3'
    dfMetsInfer_classicMethod = 'recon3D_metabolites_wInferredIds_20201218172131'

    dfMetsInfer_fw_metacyc1 = 'recon3_mappingMetaCyc_empty_20210115103057'
    dfMetsInfer_fw_metacyc2 = 'recon3_mappingMetaCyc_100_20210115103057'
    dfMetsInfer_fw_metacyc3 = 'recon3_mappingMetaCyc_91_99_20210115103057_checked'

    dfMetsInfer_fw_KC1 = 'recon3_mappingKeggC_empty_20210115103057'
    dfMetsInfer_fw_KC2 = 'recon3_mappingKeggC_100_20210115103057'
    dfMetsInfer_fw_KC3 = 'recon3_mappingKeggC_91_99_20210115103057_checked'

    dfMetsInfer_fw_KG1 = 'recon3_mappingKeggG_empty_20210115103057'
    dfMetsInfer_fw_KG2 = 'recon3_mappingKeggG_100_20210115103057'
    dfMetsInfer_fw_KG3 = 'recon3_mappingKeggG_91_99_20210115103057_checked'

    dfMetsInfer_fw_C1 = 'recon3_mappingChebi_empty_20210115103057'
    dfMetsInfer_fw_C2 = 'recon3_mappingChebi_100_20210115103057'
    dfMetsInfer_fw_C3 = 'recon3_mappingChebi_91_99_20210115103057_checked'

elif testModel == '2':
    ## Yeast 7
    outputFileName = 'yeast7'
    dfMetsInfer_classicMethod = 'y7_metabolites_wInferredIds_20201219101924'

    dfMetsInfer_fw_metacyc1 = 'yeast7_mappingMetaCyc_empty_20210113143144'
    dfMetsInfer_fw_metacyc2 = 'yeast7_mappingMetaCyc_100_20210113143144'
    dfMetsInfer_fw_metacyc3 = 'yeast7_mappingMetaCyc_91_99_20210113143144_checked'

    dfMetsInfer_fw_KC1 = 'yeast7_mappingKeggC_empty_20210113143144'
    dfMetsInfer_fw_KC2 = 'yeast7_mappingKeggC_100_20210113143144'
    dfMetsInfer_fw_KC3 = 'yeast7_mappingKeggC_91_99_20210113143144_checked'

    dfMetsInfer_fw_KG1 = 'yeast7_mappingKeggG_empty_20210113143144'
    dfMetsInfer_fw_KG2 = 'yeast7_mappingKeggG_100_20210113143144'
    dfMetsInfer_fw_KG3 = 'yeast7_mappingKeggG_91_99_20210113143144_checked'

    dfMetsInfer_fw_C1 = 'yeast7_mappingChebi_empty_20210113143144'
    dfMetsInfer_fw_C2 = 'yeast7_mappingChebi_100_20210113143144'
    dfMetsInfer_fw_C3 = 'yeast7_mappingChebi_91_99_20210113143144_checked'
elif testModel == '3':
    ## Yeast 8
    outputFileName = 'yeast8'
    dfMetsInfer_classicMethod = 'y8_metabolites_wInferredIds_20201219102046'
    dfMetsInfer_fw_metacyc1 = 'yeast8_mappingMetaCyc_empty_20210111141255'
    dfMetsInfer_fw_metacyc2 = 'yeast8_mappingMetaCyc_100_20210111141255'
    dfMetsInfer_fw_metacyc3 = 'yeast8_mappingMetaCyc_91_99_20210111141255_checked'

    dfMetsInfer_fw_KC1 = 'yeast8_mappingKeggC_empty_20210112090407'
    dfMetsInfer_fw_KC2 = 'yeast8_mappingKeggC_100_20210112090407'
    dfMetsInfer_fw_KC3 = 'yeast8_mappingKeggC_91_99_20210112090407_checked'

    dfMetsInfer_fw_KG1 = 'yeast8_mappingKeggG_empty_20210112090407'
    dfMetsInfer_fw_KG2 = 'yeast8_mappingKeggG_100_20210112090407'
    dfMetsInfer_fw_KG3 = 'yeast8_mappingKeggG_91_99_20210112090407_checked'

    dfMetsInfer_fw_C1 = 'yeast8_mappingChebi_empty_20210112090407'
    dfMetsInfer_fw_C2 = 'yeast8_mappingChebi_100_20210112090407'
    dfMetsInfer_fw_C3 = 'yeast8_mappingChebi_91_99_20210112090407_checked'

elif testModel == '4':
    ## HMRcore
    outputFileName = 'hmrCore'
    dfMetsInfer_classicMethod = 'hmrCore_metabolites_wInferredIds_20210112092737'

    dfMetsInfer_fw_metacyc1 = 'hmrCore_mappingMetaCyc_empty_20210112162153'
    dfMetsInfer_fw_metacyc2 = 'hmrCore_mappingMetaCyc_100_20210112162153'
    dfMetsInfer_fw_metacyc3 = 'hmrCore_mappingMetaCyc_91_99_20210112162153_checked'

    dfMetsInfer_fw_KC1 = 'hmrCore_mappingKeggC_empty_20210112162153'
    dfMetsInfer_fw_KC2 = 'hmrCore_mappingKeggC_100_20210112162153'
    dfMetsInfer_fw_KC3 = 'hmrCore_mappingKeggC_91_99_20210112162153_checked'

    dfMetsInfer_fw_KG1 = 'hmrCore_mappingKeggG_empty_20210112162153'
    dfMetsInfer_fw_KG2 = 'hmrCore_mappingKeggG_100_20210112162153'
    dfMetsInfer_fw_KG3 = 'hmrCore_mappingKeggG_91_99_20210112162153_checked'

    dfMetsInfer_fw_C1 = 'hmrCore_mappingChebi_empty_20210112162153'
    dfMetsInfer_fw_C2 = 'hmrCore_mappingChebi_100_20210112162153'
    dfMetsInfer_fw_C3 = 'hmrCore_mappingChebi_91_99_20210112162153_checked'

start = time.time()

dMetMapping = {}

### METACYC
dfMetaCyc = pd.read_csv(os.path.join(OUTDIR,'metacyc_compounds_20201216152513.csv'), sep='\t', dtype = str)

df1 = pd.read_csv(os.path.join(OUTDIR, dfMetsInfer_fw_metacyc1 + '.tsv'), sep='\t', index_col=0)
df1['Matches'] = df1['Matches'].apply(ast.literal_eval)
df2 = pd.read_csv(os.path.join(OUTDIR, dfMetsInfer_fw_metacyc2 + '.tsv'), sep='\t', index_col=0)
df2['Matches'] = df2['Matches'].apply(ast.literal_eval)
df3 = pd.read_csv(os.path.join(OUTDIR, dfMetsInfer_fw_metacyc3 + '.tsv'), sep='\t', index_col=0)
df3['Matches'] = df3['Matches'].apply(ast.literal_eval)

dfAll = pd.concat([df1,df2,df3])

for row in dfAll.itertuples():
    # print('matches\t', row.Matches)
    dfIsolateMatches = dfMetaCyc[dfMetaCyc['Name'].isin(row.Matches)]
    if dfIsolateMatches.empty is True:
        lCorrespondences = []
    else:
        lCorrespondences = list(dfIsolateMatches['Id'].dropna()) + list(dfIsolateMatches['ChebiId'].dropna()) + list(dfIsolateMatches['KeggId'].dropna())
    if row.Name not in dMetMapping:
        dMetMapping[row.Name] = lCorrespondences
    else:
        dMetMapping[row.Name] += lCorrespondences

### KEGG
dfKeggC = pd.read_csv(os.path.join(OUTDIR,'kegg_compounds_20201216150802.csv'), sep='\t', dtype = str)
dfKeggG = pd.read_csv(os.path.join(OUTDIR,'kegg_glycans_20201216150802.csv'), sep='\t', dtype = str)
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
    # print('matches\t', row.Matches)
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
    # print('matches\t', row.Matches)
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


## CHEBI
dfChebiNames = pd.read_csv(os.path.join(RAWDIR, 'chebi_names_20201216.tsv'), sep = '\t', dtype=str)
dfChebiUniprot = pd.read_csv(os.path.join(RAWDIR, 'chebi_uniprot_20201216.tsv'), sep = '\t', header=None, dtype=str, names=['ID','NAME_wSymbols', 'NAME'])
dfChebiCompounds =  pd.read_csv(os.path.join(RAWDIR, 'chebi_compounds_20201216.tsv'), sep = '\t', dtype=str)

df1 = pd.read_csv(os.path.join(OUTDIR, dfMetsInfer_fw_C1 + '.tsv'), sep='\t', index_col=0)
df1['Matches'] = df1['Matches'].apply(ast.literal_eval)
df2 = pd.read_csv(os.path.join(OUTDIR, dfMetsInfer_fw_C2 + '.tsv'), sep='\t', index_col=0)
df2['Matches'] = df2['Matches'].apply(ast.literal_eval)
df3 = pd.read_csv(os.path.join(OUTDIR, dfMetsInfer_fw_C3 + '.tsv'), sep='\t', index_col=0)
df3['Matches'] = df3['Matches'].apply(ast.literal_eval)

dfAll = pd.concat([df1,df2,df3])

for row in dfAll.itertuples():
    # print('matches\t', row.Matches)
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
    # print('matches\t', row.Matches)
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
    # print('matches\t', row.Matches)
    dfIsolateMatches = dfChebiCompounds[dfChebiCompounds['NAME'].isin(row.Matches)]
    if dfIsolateMatches.empty is True:
        lCorrespondences = []
    else:
        lCorrespondences = list(dfIsolateMatches['ID'].dropna())
    if row.Name.strip() not in dMetMapping:
        dMetMapping[row.Name.strip()] = lCorrespondences
    else:
        dMetMapping[row.Name.strip()] += lCorrespondences

## elimino corrispondenze duplicate per ogni metabolita
dMetMapping_woDuplicates = {}
for met in dMetMapping:
    dMetMapping_woDuplicates[met.strip()] = gL.unique(dMetMapping[met])

dfMatches = pd.DataFrame(dMetMapping_woDuplicates.items(), columns=['Name', 'Identifiers'])
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingFuzzy_' + timeStamp + '.tsv'), sep='\t', index = False)

end = time.time()
print('elapsed time\t', end-start, '\tseconds')

## confronto col metodo classico
dfFuzzyMatches = pd.read_csv(os.path.join(OUTDIR, outputFileName + '_mappingFuzzy_' + timeStamp + '.tsv'), sep='\t')
dfFuzzyMatches['Identifiers_fuzzy'] = dfFuzzyMatches['Identifiers'].apply(ast.literal_eval)
print('dfFuzzyMatches\n', dfFuzzyMatches.columns, '\n')

dfClassicMatches = pd.read_csv(os.path.join(OUTDIR, dfMetsInfer_classicMethod + '.csv'), sep='\t')
dfClassicMatches['Identifiers_classic'] = dfClassicMatches['Identifiers'].apply(ast.literal_eval)
dfClassicMatches['Name'] = dfClassicMatches.Name.str.strip()

dfMerge = pd.merge(dfFuzzyMatches, dfClassicMatches, on = 'Name')
dfMerge['Identifiers'] = dfMerge.apply(lambda row: gL.unique(row['Identifiers_fuzzy'] + row['Identifiers_classic']), axis=1)

dfMerge_filter = dfMerge[['Name', 'Identifiers_fuzzy', 'Identifiers_classic', 'Identifiers']]

empty = 0
find = 0
for row in dfMerge_filter.itertuples():
    if row.Identifiers_fuzzy == [] and row.Identifiers_classic == []:
        empty += 1
    elif row.Identifiers_fuzzy != [] or row.Identifiers_classic != []:
        find += 1

print('empty\t', empty)
print('find\t', find)

dfMerge_filter.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingFuzzyAndClassic_' + timeStamp + '.tsv'), sep='\t', index = False)
