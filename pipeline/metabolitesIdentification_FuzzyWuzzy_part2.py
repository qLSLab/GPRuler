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
    dfMappingMetaCyc = 'recon3_mappingMetaCyc_allResults_20210112093144'
    dfMappingKeggC = 'recon3_mappingKeggC_allResults_20210112093144'
    dfMappingKeggG = 'recon3_mappingKeggG_allResults_20210112093144'
    dfMappingChebi = 'recon3_mappingChebi_allResults_20210112093144'
    outputFileName = 'recon3'
elif testModel == '2':
    ## Yeast 7
    dfMappingMetaCyc = 'yeast7_mappingMetaCyc_allResults_20210112093252'
    dfMappingKeggC = 'yeast7_mappingKeggC_allResults_20210112093252'
    dfMappingKeggG = 'yeast7_mappingKeggG_allResults_20210112093252'
    dfMappingChebi = 'yeast7_mappingChebi_allResults_20210112093252'
    outputFileName = 'yeast7'
elif testModel == '3':
    ## Yeast 8
    dfMappingMetaCyc = 'yeast8_mappingMetaCyc_allResults_20210111095310'
    dfMappingKeggC = 'yeast8_mappingKeggC_allResults_20210111095310'
    dfMappingKeggG = 'yeast8_mappingKeggG_allResults_20210111095310'
    dfMappingChebi = 'yeast8_mappingChebi_allResults_20210111095310'
    outputFileName = 'yeast8'
elif testModel == '4':
    ## HMRcore
    dfMappingMetaCyc = 'hmrCore_mappingMetaCyc_allResults_20210112093341'
    dfMappingKeggC = 'hmrCore_mappingKeggC_allResults_20210112093341'
    dfMappingKeggG = 'hmrCore_mappingKeggG_allResults_20210112093341'
    dfMappingChebi = 'hmrCore_mappingChebi_allResults_20210112093341'
    outputFileName = 'hmrCore'

### METACYC
dfMappingMetacyc_all = pd.read_csv(os.path.join(OUTDIR, dfMappingMetaCyc + '.tsv'), sep = '\t', index_col=0)

dMappingMetaCyc_100 = {}
dMappingMetaCyc_91_99 = {}
dMappingMetaCyc_empty = {}

for index, row in dfMappingMetacyc_all.iterrows():
    col = 0
    threshold = 100
    print('met\t', index)
    while col < 10 and threshold > 90:
        try:
            tMatch = ast.literal_eval(row[col])
            match = tMatch[0]
            threshold = tMatch[1]
        except:
            tMatch = row[col]
            lMatch = row[col].split(',')
            match = lMatch[0].strip()[lMatch[0].find('(')+1:]
            threshold = int(lMatch[1].strip()[:lMatch[1].find(')')-1])

        if threshold == 100:
            if index not in dMappingMetaCyc_100:
                dMappingMetaCyc_100[index] = [match]
            else:
                dMappingMetaCyc_100[index] += [match]
            col = 10
        elif threshold < 100 and threshold > 90:
            if index not in dMappingMetaCyc_91_99:
                dMappingMetaCyc_91_99[index] = [match]
            else:
                dMappingMetaCyc_91_99[index] += [match]
            col += 1
        elif threshold <= 90:
            if index not in dMappingMetaCyc_empty:
                dMappingMetaCyc_empty[index] = []
            else:
                dMappingMetaCyc_empty[index] += []
            col = 10

dfMatches = pd.DataFrame(dMappingMetaCyc_100.items(), columns=['Name', 'Matches'])
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingMetaCyc_100_' + timeStamp + '.tsv'), sep='\t')

dfMatches = pd.DataFrame(dMappingMetaCyc_91_99.items(), columns=['Name', 'Matches'])
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingMetaCyc_91_99_' + timeStamp + '.tsv'), sep='\t')

dfMatches = pd.DataFrame(dMappingMetaCyc_empty.items(), columns=['Name', 'Matches'])
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingMetaCyc_empty_' + timeStamp + '.tsv'), sep='\t')

## KEGG
# dfKeggC = pd.read_csv(os.path.join(OUTDIR,'kegg_compounds_20201216150802.csv'), sep='\t', usecols=['Id', 'Name'])
# dfKeggG = pd.read_csv(os.path.join(OUTDIR,'kegg_glycans_20201216150802.csv'), sep='\t', usecols=['Id', 'Name'])

dfMappingKeggC_all = pd.read_csv(os.path.join(OUTDIR, dfMappingKeggC + '.tsv'), sep = '\t', index_col=0)
dfMappingKeggG_all = pd.read_csv(os.path.join(OUTDIR, dfMappingKeggG + '.tsv'), sep = '\t', index_col=0)

dMappingKeggC_100 = {}
dMappingKeggC_91_99 = {}
dMappingKeggC_empty = {}

for index, row in dfMappingKeggC_all.iterrows():
    col = 0
    threshold = 100
    print('met\t', index)
    while col < 10 and threshold > 90:
        try:
            tMatch = ast.literal_eval(row[col])
            match = tMatch[0]
            threshold = tMatch[1]
        except:
            tMatch = row[col]
            lMatch = row[col].split(',')
            match = lMatch[0].strip()[lMatch[0].find('(')+1:]
            threshold = int(lMatch[1].strip()[:lMatch[1].find(')')-1])

        if threshold == 100:
            if index not in dMappingKeggC_100:
                dMappingKeggC_100[index] = [match]
            else:
                dMappingKeggC_100[index] += [match]
            col = 10
        elif threshold < 100 and threshold > 90:
            if index not in dMappingKeggC_91_99:
                dMappingKeggC_91_99[index] = [match]
            else:
                dMappingKeggC_91_99[index] += [match]
            col += 1
        elif threshold <= 90:
            if index not in dMappingKeggC_empty:
                dMappingKeggC_empty[index] = []
            else:
                dMappingKeggC_empty[index] += []
            col = 10

dfMatches = pd.DataFrame(dMappingKeggC_100.items(), columns=['Name', 'Matches'])
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingKeggC_100_' + timeStamp + '.tsv'), sep='\t')

dfMatches = pd.DataFrame(dMappingKeggC_91_99.items(), columns=['Name', 'Matches'])
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingKeggC_91_99_' + timeStamp + '.tsv'), sep='\t')

dfMatches = pd.DataFrame(dMappingKeggC_empty.items(), columns=['Name', 'Matches'])
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingKeggC_empty_' + timeStamp + '.tsv'), sep='\t')

dMappingKeggG_100 = {}
dMappingKeggG_91_99 = {}
dMappingKeggG_empty = {}

for index, row in dfMappingKeggG_all.iterrows():
    col = 0
    threshold = 100
    print('met\t', index)
    while col < 10 and threshold > 90:
        try:
            tMatch = ast.literal_eval(row[col])
            match = tMatch[0]
            threshold = tMatch[1]
        except:
            tMatch = row[col]
            lMatch = row[col].split(',')
            match = lMatch[0].strip()[lMatch[0].find('(')+1:]
            threshold = int(lMatch[1].strip()[:lMatch[1].find(')')-1])

        if threshold == 100:
            if index not in dMappingKeggG_100:
                dMappingKeggG_100[index] = [match]
            else:
                dMappingKeggG_100[index] += [match]
            col = 10
        elif threshold < 100 and threshold > 90:
            if index not in dMappingKeggG_91_99:
                dMappingKeggG_91_99[index] = [match]
            else:
                dMappingKeggG_91_99[index] += [match]
            col += 1
        elif threshold <= 90:
            if index not in dMappingKeggG_empty:
                dMappingKeggG_empty[index] = []
            else:
                dMappingKeggG_empty[index] += []
            col = 10

dfMatches = pd.DataFrame(dMappingKeggG_100.items(), columns=['Name', 'Matches'])
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingKeggG_100_' + timeStamp + '.tsv'), sep='\t')

dfMatches = pd.DataFrame(dMappingKeggG_91_99.items(), columns=['Name', 'Matches'])
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingKeggG_91_99_' + timeStamp + '.tsv'), sep='\t')

dfMatches = pd.DataFrame(dMappingKeggG_empty.items(), columns=['Name', 'Matches'])
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingKeggG_empty_' + timeStamp + '.tsv'), sep='\t')

## CHEBI
dfMappingChebi_all = pd.read_csv(os.path.join(OUTDIR, dfMappingChebi + '.tsv'), sep = '\t', index_col=0)

dMappingChebi_100 = {}
dMappingChebi_91_99 = {}
dMappingChebi_empty = {}

for index, row in dfMappingChebi_all.iterrows():
    col = 0
    threshold = 100
    print('met\t', index)
    while col < 10 and threshold > 90:
        try:
            tMatch = ast.literal_eval(row[col])
            match = tMatch[0]
            threshold = tMatch[1]
        except:
            tMatch = row[col]
            lMatch = row[col].split(',')
            match = lMatch[0].strip()[lMatch[0].find('(')+1:]
            threshold = int(lMatch[1].strip()[:lMatch[1].find(')')-1])

        if threshold == 100:
            if index not in dMappingChebi_100:
                dMappingChebi_100[index] = [match]
            else:
                dMappingChebi_100[index] += [match]
            col = 10
        elif threshold < 100 and threshold > 90:
            if index not in dMappingChebi_91_99:
                dMappingChebi_91_99[index] = [match]
            else:
                dMappingChebi_91_99[index] += [match]
            col += 1
        elif threshold <= 90:
            if index not in dMappingChebi_empty:
                dMappingChebi_empty[index] = []
            else:
                dMappingChebi_empty[index] += []
            col = 10

dfMatches = pd.DataFrame(dMappingChebi_100.items(), columns=['Name', 'Matches'])
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingChebi_100_' + timeStamp + '.tsv'), sep='\t')

dfMatches = pd.DataFrame(dMappingChebi_91_99.items(), columns=['Name', 'Matches'])
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingChebi_91_99_' + timeStamp + '.tsv'), sep='\t')

dfMatches = pd.DataFrame(dMappingChebi_empty.items(), columns=['Name', 'Matches'])
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingChebi_empty_' + timeStamp + '.tsv'), sep='\t')
