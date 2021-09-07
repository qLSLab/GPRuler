import sys
import genericLib as gL
import os
import pandas as pd
import time
from ast import literal_eval

timeStamp = gL.getTimeStamp()

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]
FIGUREDIR = workingDirs[4]
DFDIR = workingDirs[5]

testModel = sys.argv[1]

recon if the pipeline starts from Recon 3D model,
y7 if the pipeline starts from Yeast 7 model, y8 if the pipeline starts from Yeast 8 model,
. if the pipeline starts from HMRcore model, or ownData if the pipeline starts from an existing SBML model added by the user.


## Input data
if testModel == 'hmr':
    ## HMRcore
    modelName = 'hmrCore'
    inputFuzzy = 'hmrCore_mappingFuzzyAndClassic_20210112165030.tsv'
elif testModel == 'recon':
    ## Recon3D
    modelName = 'recon3D'
    inputFuzzy = 'recon3_mappingFuzzyAndClassic_20210119085007.tsv'
elif testModel == 'y7':
    ## Yeast7
    modelName = 'yeast7'
    inputFuzzy = 'yeast7_mappingFuzzyAndClassic_20210113162642.tsv'
elif testModel == 'y8':
    ## Yeast7
    modelName = 'yeast8'
    inputFuzzy = 'yeast8_mappingFuzzyAndClassic_20210113162740.tsv'

## STEP 1. Get metabolites identifiers
dfchem_prop = pd.read_csv('chem_prop.tsv', sep = '\t', dtype=str, skiprows = 348, names = ['ID', 'name', 'reference', 'formula', 'charge', 'mass', 'InChI', 'InChIKey', 'SMILES'])
dfchem_xref = pd.read_csv('chem_xref.tsv', sep = '\t', dtype=str, skiprows = 348, names = ['source', 'ID', 'description'])

lDescription = []
for row in dfchem_xref.itertuples():
    if '||' not in row.description:
        lDescription.append([row.description])
    else:
        lDescription.append(row.description.split('||'))

dfchem_xref['description_splt'] = lDescription

dfchem_merge = pd.merge(dfchem_prop, dfchem_xref, on = 'ID')
dfchem_merge['name'] = dfchem_merge['name'].str.lower()

## Load list of metabolites from model
dfmets =  pd.read_csv(os.path.join(OUTDIR, inputFuzzy), sep = '\t', dtype=str)
dfmets['Name'] = dfmets['Name'].str.lower()

## get identifiers associated to common metabolites
inMetaNetX = gL.intersect(dfmets['Name'].tolist(), dfchem_merge['name'].tolist())
dfchem_merge_group = dfchem_merge[dfchem_merge['name'].isin(dfmets['Name'].tolist())]
dfFinal = dfchem_merge_group.groupby('ID')[['source', 'name']].agg(list).reset_index()

## test how many metabolites of the input model are not included in the column 'name' of dfchem_merge dataframe
notInMetaNetX = gL.difference(dfmets['Name'].tolist(), dfchem_merge['name'].tolist())

dfFinal_p2 = pd.DataFrame({'ID':[], 'source':[], 'description_splt':[]})
dfFilter = dfchem_merge[dfchem_merge.apply(lambda x: notInMetaNetX[0] in [el.lower() for el in x['description_splt']] , axis=1)].reset_index()
if dfFilter.empty == True:
    print('NOT FOUND:\t', notInMetaNetX[0])
else:
    dfFinal_p2 = dfFilter[['ID', 'source', 'description_splt']]
    dfFinal_p2 = dfFinal_p2.groupby('ID')[['source', 'description_splt']].agg(list).reset_index()
    print(dfFinal_p2)

for n in notInMetaNetX[1:]:
    dfFilter = dfchem_merge[dfchem_merge.apply(lambda x: n in [el.lower() for el in x['description_splt']] , axis=1)].reset_index()
    if dfFilter.empty == True:
        print('NOT FOUND:\t', n)
    else:
        dfFinal_tmp = dfFilter[['ID', 'source', 'description_splt']]
        dfFinal_tmp = dfFinal_tmp.groupby('ID')[['source', 'description_splt']].agg(list).reset_index()
        dfFinal_p2 = pd.concat([dfFinal_p2, dfFinal_tmp]).reset_index(drop = True)

dfFinal_p2.rename(columns={"description_splt": "name"}, inplace = True)
dfFinal_p2 = dfFinal_p2.reset_index(drop = True)

dfFinal_p2['name'] = dfFinal_p2.name.apply(lambda x: sum(x, []))
dfFinal_p2 = dfFinal_p2.explode('name')

dfFinal_p2['name'] = dfFinal_p2['name'].str.lower()
dfFinal_p2 = dfFinal_p2[dfFinal_p2['name'].isin(notInMetaNetX)]

dfFinal_p2 = dfFinal_p2.drop_duplicates(subset = ['ID', 'name'])
dfFinal_p2 = dfFinal_p2.reset_index(drop = True)

## join dfFinal and dfFinal_p2
dfAllMets = pd.concat([dfFinal, dfFinal_p2])
dfAllMets.to_csv(os.path.join(OUTDIR, modelName + '_mappingMetaNetX_20210901.csv'), sep = '\t')

## STEP 2. Comparison of MetaNetX and Fuzzy Wuzzy output
dfMetsFromMetaNetX = pd.read_csv(os.path.join(OUTDIR, modelName + '_mappingMetaNetX_20210901.csv'), sep = '\t', dtype=str, index_col = 0)

lName = []
for row in dfMetsFromMetaNetX.itertuples():
    if row.name.startswith("['"):
        name = literal_eval(row.name)
        name = gL.unique(name)
        lName.append(name)
    elif row.name.startswith('["'):
        name = literal_eval(row.name)
        name = gL.unique(name)
        lName.append(name)
    else:
        lName.append([row.name])

dfMetsFromMetaNetX['name'] = lName
dfMetsFromMetaNetX = dfMetsFromMetaNetX.explode('name')
dfMetsFromMetaNetX.to_csv(os.path.join(OUTDIR, modelName + '_mappingMetaNetX_20210901.csv'), sep = '\t')
dfmetsFuzzy =  pd.read_csv(os.path.join(OUTDIR, inputFuzzy), sep = '\t', dtype=str)
dfmetsFuzzy['Name'] = dfmetsFuzzy['Name'].str.lower()

lMetsModel = dfmetsFuzzy['Name'].tolist()

dfComparison = pd.merge(dfmetsFuzzy, dfMetsFromMetaNetX, left_on = 'Name', right_on = 'name', how = 'inner')
dfComparison['Identifiers'] = dfComparison['Identifiers'].apply(literal_eval)
dfComparison['source'] = dfComparison['source'].apply(literal_eval)

## fix identifiers to find the correct match, when present
lfuzzy = []
lmeta = []
for row in dfComparison.itertuples():
    lfuzzy_tmp = []
    for r1 in row.Identifiers:
        if r1.startswith('|'):
            lfuzzy_tmp.append(r1[1:-1])
        else:
            lfuzzy_tmp.append(r1)
    lfuzzy_tmp = gL.unique(lfuzzy_tmp)
    lfuzzy.append(lfuzzy_tmp)
    lmeta_tmp = []
    for r2 in row.source:
        if r2.startswith('CHEBI:') or r2.startswith('chebi:') or r2.startswith('kegg.compound:') or r2.startswith('metacyc.compound:'):
            lmeta_tmp.append(r2.split(':')[1])
    lmeta_tmp = gL.unique(lmeta_tmp)
    lmeta.append(lmeta_tmp)

dfComparison['fuzzy'] = lfuzzy
dfComparison['metaNet'] = lmeta
dfComparison.to_csv(os.path.join(OUTDIR, modelName + '_comparisonMetaNetX_FuzzyWuzzy.csv'), sep = '\t', index = False)

dfCommon = dfComparison[dfComparison.apply(lambda x: len(list(x['fuzzy'])) != 0 and len(list(x['metaNet'])) != 0, axis=1)]
dfCommon = dfCommon[dfCommon.apply(lambda x: len(list(set(x['fuzzy']) & set(x['metaNet']))) != 0 and len(list(set(x['fuzzy']) - set(x['metaNet']))) != 0 and len(list(set(x['metaNet']) - set(x['fuzzy']))) != 0, axis=1)]
lcommonMets = gL.unique(dfCommon['Name'])
print('Case 1: common identifiers for the same metabolite in addition to identifiers provided by our methodology not included in MetaNetX output and viceversa:\t', len(lcommonMets))

dfCommon = dfComparison[dfComparison.apply(lambda x: len(list(x['fuzzy'])) != 0 and len(list(x['metaNet'])) != 0, axis=1)]
dfCommon = dfCommon[dfCommon.apply(lambda x: len(list(set(x['fuzzy']) & set(x['metaNet']))) != 0 and len(list(set(x['fuzzy']) - set(x['metaNet']))) != 0 and len(list(set(x['metaNet']) - set(x['fuzzy']))) == 0, axis=1)]
lcommonMets = gL.unique(dfCommon['Name'])
print('Case 2: common identifiers for the same metabolite in addition to identifiers provided by only our methodology not included in MetaNetX output:\t', len(lcommonMets))

dfCommon = dfComparison[dfComparison.apply(lambda x: len(list(x['fuzzy'])) != 0 and len(list(x['metaNet'])) != 0, axis=1)]
dfCommon = dfCommon[dfCommon.apply(lambda x: len(list(set(x['fuzzy']) & set(x['metaNet']))) != 0 and len(list(set(x['fuzzy']) - set(x['metaNet']))) == 0 and len(list(set(x['metaNet']) - set(x['fuzzy']))) != 0, axis=1)]
lcommonMets = gL.unique(dfCommon['Name'])
print('Case 3: common identifiers for the same metabolite in addition to identifiers provided by only MetaNetX not included in our methodology output:\t', len(lcommonMets))

dfCommon = dfComparison[dfComparison.apply(lambda x: len(list(x['fuzzy'])) != 0 and len(list(x['metaNet'])) != 0, axis=1)]
dfCommon = dfCommon[dfCommon.apply(lambda x: len(list(set(x['fuzzy']) & set(x['metaNet']))) != 0 and len(list(set(x['fuzzy']) - set(x['metaNet']))) == 0 and len(list(set(x['metaNet']) - set(x['fuzzy']))) == 0, axis=1)]
lcommonMets = gL.unique(dfCommon['Name'])
print('Case 4: perfect match between identifiers for the same metabolite provided by MetaNetX and our methodology:\t', len(lcommonMets))

dfCommon = dfComparison[dfComparison.apply(lambda x: len(list(x['fuzzy'])) != 0 and len(list(x['metaNet'])) != 0, axis=1)]
dfCommon = dfCommon[dfCommon.apply(lambda x: len(list(set(x['fuzzy']) & set(x['metaNet']))) == 0, axis=1)]
lcommonMets = gL.unique(dfCommon['Name'])
print('Case 5: no intersection between identifiers for the same metabolite provided by MetaNetX and our methodology:\t', len(lcommonMets))

dfCommon = dfComparison[dfComparison.apply(lambda x: len(list(x['fuzzy'])) != 0 and len(list(x['metaNet'])) == 0, axis=1)]
lcommonMets = gL.unique(dfCommon['Name'])
print('Case 6: identifiers for the same metabolite provided by only our methodology:\t', len(lcommonMets))

dfCommon = dfComparison[dfComparison.apply(lambda x: len(list(x['fuzzy'])) == 0 and len(list(x['metaNet'])) != 0, axis=1)]
lcommonMets = gL.unique(dfCommon['Name'])
print('Case 7: identifiers for the same metabolite provided by only MetaNetX:\t', len(lcommonMets))
