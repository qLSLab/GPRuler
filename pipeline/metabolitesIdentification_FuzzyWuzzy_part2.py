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
    dfMappingMetaCyc = 'recon3_mappingMetaCyc_allResults'
    dfMappingKeggC = 'recon3_mappingKeggC_allResults'
    dfMappingKeggG = 'recon3_mappingKeggG_allResults'
    dfMappingChebi = 'recon3_mappingChebi_allResults'
    outputFileName = 'recon3'
elif testModel == 'y7':
    ## Yeast 7
    dfMappingMetaCyc = 'yeast7_mappingMetaCyc_allResults'
    dfMappingKeggC = 'yeast7_mappingKeggC_allResults'
    dfMappingKeggG = 'yeast7_mappingKeggG_allResults'
    dfMappingChebi = 'yeast7_mappingChebi_allResults'
    outputFileName = 'yeast7'
elif testModel == 'y8':
    ## Yeast 8
    dfMappingMetaCyc = 'yeast8_mappingMetaCyc_allResults'
    dfMappingKeggC = 'yeast8_mappingKeggC_allResults'
    dfMappingKeggG = 'yeast8_mappingKeggG_allResults'
    dfMappingChebi = 'yeast8_mappingChebi_allResults'
    outputFileName = 'yeast8'
elif testModel == 'hmr':
    ## HMRcore
    dfMappingMetaCyc = 'hmrCore_mappingMetaCyc_allResults'
    dfMappingKeggC = 'hmrCore_mappingKeggC_allResults'
    dfMappingKeggG = 'hmrCore_mappingKeggG_allResults'
    dfMappingChebi = 'hmrCore_mappingChebi_allResults'
    outputFileName = 'hmrCore'
elif testModel == 'ownData':
    ## specify your input data
    dfMappingMetaCyc = ''
    dfMappingKeggC = ''
    dfMappingKeggG = ''
    dfMappingChebi = ''
    outputFileName = ''

### METACYC
dfMappingMetacyc_all = pd.read_csv(os.path.join(OUTDIR, dfMappingMetaCyc + '.tsv'), sep = '\t', index_col=0)
dfMatches100, dfMatches91_99, dfMatchesEmpty = metL.filterFWoutputs(dfMappingMetacyc_all)
dfMatches100.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingMetaCyc_100.tsv'), sep='\t')
dfMatches91_99.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingMetaCyc_91_99.tsv'), sep='\t')
dfMatchesEmpty.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingMetaCyc_empty.tsv'), sep='\t')

## KEGG
dfMappingKeggC_all = pd.read_csv(os.path.join(OUTDIR, dfMappingKeggC + '.tsv'), sep = '\t', index_col=0)
dfMappingKeggG_all = pd.read_csv(os.path.join(OUTDIR, dfMappingKeggG + '.tsv'), sep = '\t', index_col=0)
dfMatches100, dfMatches91_99, dfMatchesEmpty = metL.filterFWoutputs(dfMappingKeggC_all)
dfMatches100.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingKeggC_100.tsv'), sep='\t')
dfMatches91_99.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingKeggC_91_99.tsv'), sep='\t')
dfMatchesEmpty.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingKeggC_empty.tsv'), sep='\t')

dfMatches100, dfMatches91_99, dfMatchesEmpty = metL.filterFWoutputs(dfMappingKeggG_all)
dfMatches100.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingKeggG_100.tsv'), sep='\t')
dfMatches91_99.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingKeggG_91_99.tsv'), sep='\t')
dfMatchesEmpty.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingKeggG_empty.tsv'), sep='\t')

## CHEBI
dfMappingChebi_all = pd.read_csv(os.path.join(OUTDIR, dfMappingChebi + '.tsv'), sep = '\t', index_col=0)
dfMatches100, dfMatches91_99, dfMatchesEmpty = metL.filterFWoutputs(dfMappingChebi_all)
dfMatches100.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingChebi_100.tsv'), sep='\t')
dfMatches91_99.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingChebi_91_99.tsv'), sep='\t')
dfMatchesEmpty.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingChebi_empty.tsv'), sep='\t')
