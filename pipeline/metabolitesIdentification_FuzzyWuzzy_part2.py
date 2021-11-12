# -*- coding: utf-8 -*-
import os
import sys
import ast
import pandas as pd
import genericLib as gL
import metabolitesLib as metL

# setting working dirs
workingDirs = gL.setWorkingDirs()
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

### METACYC
dfMappingMetacyc_all = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingMetaCyc_allResults.tsv'), sep = '\t', index_col=0)
dfMatches100, dfMatches91_99, dfMatchesEmpty = metL.filterFWoutputs(dfMappingMetacyc_all)
dfMatches100.to_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingMetaCyc_100.tsv'), sep='\t')
dfMatches91_99.to_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingMetaCyc_91_99.tsv'), sep='\t')
dfMatchesEmpty.to_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingMetaCyc_empty.tsv'), sep='\t')

## Automatic interaction with the user to curate the list of putative matches with score between 91 and 99
dfMatches91_99 = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingMetaCyc_91_99.tsv'), sep="\t", index_col = 0)
dfMatches91_99['Matches'] = dfMatches91_99['Matches'].apply(ast.literal_eval)

curatedFWFile = open(os.path.join(OUTDIR, prefix_modelName + '_mappingMetaCyc_91_99.tsv'), mode='w')
gL.writeLineByLineToFile(curatedFWFile, ['Name', 'Matches'], '\t')

for row in dfMatches91_99.itertuples():
    lMatches = []
    for mat in row.Matches:
        print('Target metabolite:\t', row.Name)
        print('Proposed match:\t', mat, '\n')
        choice = input('Do you want to keep this proposed metabolite name inside the final list: yes (y) or no (n)? ')
        if choice == 'y':
            lMatches.append(mat)

    gL.writeLineByLineToFile(curatedFWFile, [row.Name, lMatches], '\t')

curatedFWFile.close()


## KEGG
dfMappingKeggC_all = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingKeggC_allResults.tsv'), sep = '\t', index_col=0)
dfMappingKeggG_all = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingKeggG_allResults.tsv'), sep = '\t', index_col=0)
dfMatches100, dfMatches91_99, dfMatchesEmpty = metL.filterFWoutputs(dfMappingKeggC_all)
dfMatches100.to_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingKeggC_100.tsv'), sep='\t')
dfMatches91_99.to_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingKeggC_91_99.tsv'), sep='\t')
dfMatchesEmpty.to_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingKeggC_empty.tsv'), sep='\t')

## Automatic interaction with the user to curate the list of putative matches with score between 91 and 99
dfMatches91_99 = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingKeggC_91_99.tsv'), sep="\t", index_col = 0)
dfMatches91_99['Matches'] = dfMatches91_99['Matches'].apply(ast.literal_eval)

curatedFWFile = open(os.path.join(OUTDIR, prefix_modelName + '_mappingKeggC_91_99.tsv'), mode='w')
gL.writeLineByLineToFile(curatedFWFile, ['Name', 'Matches'], '\t')

for row in dfMatches91_99.itertuples():
    lMatches = []
    for mat in row.Matches:
        print('Target metabolite:\t', row.Name)
        print('Proposed match:\t', mat, '\n')
        choice = input('Do you want to keep this proposed metabolite name inside the final list: yes (y) or no (n)? ')
        if choice == 'y':
            lMatches.append(mat)

    gL.writeLineByLineToFile(curatedFWFile, [row.Name, lMatches], '\t')

curatedFWFile.close()

dfMatches100, dfMatches91_99, dfMatchesEmpty = metL.filterFWoutputs(dfMappingKeggG_all)
dfMatches100.to_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingKeggG_100.tsv'), sep='\t')
dfMatches91_99.to_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingKeggG_91_99.tsv'), sep='\t')
dfMatchesEmpty.to_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingKeggG_empty.tsv'), sep='\t')

## Automatic interaction with the user to curate the list of putative matches with score between 91 and 99
dfMatches91_99 = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingKeggG_91_99.tsv'), sep="\t", index_col = 0)
dfMatches91_99['Matches'] = dfMatches91_99['Matches'].apply(ast.literal_eval)

curatedFWFile = open(os.path.join(OUTDIR, prefix_modelName + '_mappingKeggG_91_99.tsv'), mode='w')
gL.writeLineByLineToFile(curatedFWFile, ['Name', 'Matches'], '\t')

for row in dfMatches91_99.itertuples():
    lMatches = []
    for mat in row.Matches:
        print('Target metabolite:\t', row.Name)
        print('Proposed match:\t', mat, '\n')
        choice = input('Do you want to keep this proposed metabolite name inside the final list: yes (y) or no (n)? ')
        if choice == 'y':
            lMatches.append(mat)

    gL.writeLineByLineToFile(curatedFWFile, [row.Name, lMatches], '\t')

curatedFWFile.close()

## CHEBI
dfMappingChebi_all = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingChebi_allResults.tsv'), sep = '\t', index_col=0)
dfMatches100, dfMatches91_99, dfMatchesEmpty = metL.filterFWoutputs(dfMappingChebi_all)
dfMatches100.to_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingChebi_100.tsv'), sep='\t')
dfMatches91_99.to_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingChebi_91_99.tsv'), sep='\t')
dfMatchesEmpty.to_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingChebi_empty.tsv'), sep='\t')

## Automatic interaction with the user to curate the list of putative matches with score between 91 and 99
dfMatches91_99 = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingChebi_91_99.tsv'), sep="\t", index_col = 0)
dfMatches91_99['Matches'] = dfMatches91_99['Matches'].apply(ast.literal_eval)

curatedFWFile = open(os.path.join(OUTDIR, prefix_modelName + '_mappingChebi_91_99.tsv'), mode='w')
gL.writeLineByLineToFile(curatedFWFile, ['Name', 'Matches'], '\t')

for row in dfMatches91_99.itertuples():
    lMatches = []
    for mat in row.Matches:
        print('Target metabolite:\t', row.Name)
        print('Proposed match:\t', mat, '\n')
        choice = input('Do you want to keep this proposed metabolite name inside the final list: yes (y) or no (n)? ')
        if choice == 'y':
            lMatches.append(mat)

    gL.writeLineByLineToFile(curatedFWFile, [row.Name, lMatches], '\t')

curatedFWFile.close()
