# -*- coding: utf-8 -*-
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
    prefix_modelName = 'recon3D'
elif testModel == 'y7':
    ## Yeast 7
    prefix_modelName = 'y7'
elif testModel == 'y8':
    ## Yeast 8
    prefix_modelName = 'y8'
elif testModel == 'hmr':
    ## HMRcore
    prefix_modelName = 'hmrCore'
elif testModel == 'ownData':
    ## specify your input data
    prefix_modelName = ''

dMetMapping = {}

# joining outputs from MetaCyc
dfMetaCyc = pd.read_csv(os.path.join(RAWDIR,'metacyc_compounds_20201216152513.csv'), sep='\t', dtype = str)

df1 = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingMetaCyc_empty.tsv'), sep='\t', index_col=0)
df1['Matches'] = df1['Matches'].apply(ast.literal_eval)
df2 = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingMetaCyc_100.tsv'), sep='\t', index_col=0)
df2['Matches'] = df2['Matches'].apply(ast.literal_eval)
df3 = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingMetaCyc_91_99.tsv'), sep='\t')
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

df1 = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingKeggC_empty.tsv'), sep='\t', index_col=0)
df1['Matches'] = df1['Matches'].apply(ast.literal_eval)
df2 = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingKeggC_100.tsv'), sep='\t', index_col=0)
df2['Matches'] = df2['Matches'].apply(ast.literal_eval)
df3 = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingKeggC_91_99.tsv'), sep='\t')
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


df1 = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingKeggG_empty.tsv'), sep='\t', index_col=0)
df1['Matches'] = df1['Matches'].apply(ast.literal_eval)
df2 = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingKeggG_100.tsv'), sep='\t', index_col=0)
df2['Matches'] = df2['Matches'].apply(ast.literal_eval)
df3 = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingKeggG_91_99.tsv'), sep='\t')
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

df1 = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingChebi_empty.tsv'), sep='\t', index_col=0)
df1['Matches'] = df1['Matches'].apply(ast.literal_eval)
df2 = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingChebi_100.tsv'), sep='\t', index_col=0)
df2['Matches'] = df2['Matches'].apply(ast.literal_eval)
df3 = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingChebi_91_99.tsv'), sep='\t')
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
dfMatches.to_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingFuzzy.tsv'), sep='\t', index = False)


## Joining FuzzyWuzzy and output from Step1
dfFuzzyMatches = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingFuzzy.tsv'), sep='\t')
dfFuzzyMatches['Identifiers_fuzzy'] = dfFuzzyMatches['Identifiers'].apply(ast.literal_eval)

dfClassicMatches = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_metabolites_wInferredIds.csv'), sep='\t')
dfClassicMatches['Identifiers_classic'] = dfClassicMatches['Identifiers'].apply(ast.literal_eval)
dfClassicMatches['Name'] = dfClassicMatches.Name.str.strip()

dfMerge = pd.merge(dfFuzzyMatches, dfClassicMatches, on = 'Name')
dfMerge['Identifiers'] = dfMerge.apply(lambda row: gL.unique(row['Identifiers_fuzzy'] + row['Identifiers_classic']), axis=1)

dfMerge_filter = dfMerge[['Name', 'Identifiers_fuzzy', 'Identifiers_classic', 'Identifiers']]
dfMerge_filter.to_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingFuzzyAndClassic.tsv'), sep='\t', index = False)
