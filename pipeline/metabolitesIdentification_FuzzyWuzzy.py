import pandas as pd
import fuzzywuzzy as fw
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
import sys
from ast import literal_eval
import sys
import genericLib as gL
import os

def applyFW2Db(lMets2Search, lPutativeTargets):
    '''
    Apply Fuzzy Wuzzy to a list of compounds.
    Parameters:
    - lMets2Search: list of metabolites from the input model.
    - lPutativeTargets: list of compounds to look through.
    '''
    dMapping = {}
    for toMatch in lMets2Search:
        matches = process.extract(toMatch, lPutativeTargets, limit=10)
        if toMatch in dMapping:
            dMapping[toMatch] += matches
        else:
            dMapping[toMatch] = matches
    dfMatches = pd.DataFrame(dMapping).T
    return dfMatches

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]

# setting input data
testModel = sys.argv[1]
if testModel == 'recon':
    ## Recon 3
    modelXml = 'Recon3D_301_20200923'
    dfmetsInfo = 'recon3D_metabolites'
    outputFileName = 'recon3'
elif testModel == 'y7':
    ## Yeast 7
    modelXml = 'yeast_7.6_cobra'
    dfmetsInfo = 'y7_metabolites'
    outputFileName = 'yeast7'
elif testModel == 'y8':
    ## Yeast 8
    modelXml = 'yeast8'
    dfmetsInfo = 'y8_metabolites'
    outputFileName = 'yeast8'
elif testModel == 'hmr':
    ## HMRcore
    modelXml = 'HMRcore_20200328_wReconNames'
    dfmetsInfo = 'hmrCore_metabolites'
    outputFileName = 'hmrCore'
elif testModel == 'ownData':
    ## specify your input data
    modelXml = ''
    dfmetsInfo = ''
    outputFileName = ''

# extracting metabolites info from the input model
df = pd.read_csv(os.path.join(OUTDIR, dfmetsInfo + '_' + timeStamp + '.csv'), sep = '\t', dtype = str)

dfChebiDb = pd.read_csv(os.path.join(RAWDIR, 'chebi_database_accession_20201216.tsv.bz2'), sep = '\t', compression='bz2', dtype = {'ID': str, 'COMPOUND_ID': str, 'ACCESSION_NUMBER': str})
dfChebiRelations = pd.read_csv(os.path.join(RAWDIR, 'chebi_relation_20201216.tsv.bz2'), sep = '\t', compression='bz2', dtype = {'ID': str, 'TYPE': str, 'INIT_ID': str, 'FINAL_ID': str})

dfChebiInchi = pd.read_csv(os.path.join(RAWDIR, 'chebiId_inchi_20201216.tsv.bz2'), sep = '\t', compression='bz2', dtype=str)
dfChebiInchi['InChI_splitted'] = dfChebiInchi.InChI.str.split('/')
dfChebiInchi['InChI_formula'] = dfChebiInchi.InChI.str.split('/').str[1]

if testModel == 'y7' or testModel == 'y8':
    df['toMatch'] = df.Name.str.replace("\[(\w*\s*)+\]$", "")
    lMets2Search = list(df['toMatch'])
else:
    lMets2Search = list(df['Name'])
lMets2Search = gL.unique(lMets2Search)

# query MetaCyc compounds
dfMetaCyc = pd.read_csv(os.path.join(RAWDIR,'metacyc_compounds_20201216152513.csv'), sep='\t', usecols=['Id', 'Name'])
choicesMetaCyc = dfMetaCyc.Name.to_list()
dfMatches = applyFW2Db(lMets2Search, choicesMetaCyc)
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingMetaCyc_allResults.tsv'), sep='\t')

# query KEGG compounds
dfKeggC = pd.read_csv(os.path.join(RAWDIR,'kegg_compounds_20201216150802.csv'), sep='\t', usecols=['Id', 'Name'])
dfKeggG = pd.read_csv(os.path.join(RAWDIR,'kegg_glycans_20201216150802.csv'), sep='\t', usecols=['Id', 'Name'])
dfKeggC['Name'] = dfKeggC['Name'].apply(literal_eval)
dfKeggG['Name'] = dfKeggG['Name'].apply(literal_eval)
dfKeggC = dfKeggC.explode('Name')
dfKeggG = dfKeggG.explode('Name')
choicesKeggC = dfKeggC.Name.to_list() + dfKeggC.Id.to_list()
dfMatches = applyFW2Db(lMets2Search, choicesKeggC)
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingKeggC_allResults.tsv'), sep='\t')

choicesKeggG = dfKeggG.Name.to_list() + dfKeggG.Id.to_list()
dfMatches = applyFW2Db(lMets2Search, choicesKeggG)
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingKeggG_allResults.tsv'), sep='\t')

# query ChEBI compounds
dfChebiNames = pd.read_csv(os.path.join(RAWDIR, 'chebi_names_20201216.tsv.bz2'), sep = '\t', compression='bz2', dtype=str)
dfChebiUniprot = pd.read_csv(os.path.join(RAWDIR, 'chebi_uniprot_20201216.tsv.bz2'), sep = '\t',compression='bz2', header=None, dtype=str, names=['ID','NAME_wSymbols', 'NAME'])
dfChebiCompounds =  pd.read_csv(os.path.join(RAWDIR, 'chebi_compounds_20201216.tsv.bz2'), sep = '\t', compression='bz2', dtype=str)

choicesChebi = dfChebiNames.NAME.to_list() + dfChebiUniprot.NAME.to_list() + dfChebiCompounds.NAME.to_list()
dfMatches = applyFW2Db(lMets2Search, choicesChebi)
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingChebi_allResults_' + timeStamp + '.tsv'), sep='\t')
