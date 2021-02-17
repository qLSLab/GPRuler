import pandas as pd
import fuzzywuzzy as fw
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
import sys
from ast import literal_eval
import keggLib as kL
import sys
import genericLib as gL
import os
import time
import xmlLib2 as xL

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
    modelXml = 'Recon3D_301_20200923'
    dfmetsInfo = 'recon3D_metabolites'
    # dfmetsInfo = 'Recon3D_301_enriched'
    outputFileName = 'recon3'
elif testModel == '2':
    ## Yeast 7
    modelXml = 'yeast_7.6_cobra'
    # dfmetsInfo = 'y7_mets_20200814'
    dfmetsInfo = 'y7_metabolites'
    outputFileName = 'yeast7'
elif testModel == '3':
    ## Yeast 8
    modelXml = 'yeast8'
    dfmetsInfo = 'y8_metabolites'
    outputFileName = 'yeast8'
elif testModel == '4':
    ## HMRcore
    modelXml = 'HMRcore_20200328_wReconNames'
    dfmetsInfo = 'hmrCore_metabolites'
    outputFileName = 'hmrCore'


cId, cName, cKegg, cChebi, cPubchem, cBoundaryCondition, cChemicalFormula, cInchi = xL.getMetsInfoGEM(os.path.join(RAWDIR, modelXml + ".xml"))
df = pd.DataFrame({'Id': cId, 'Name': cName, 'KeggId': cKegg, 'ChebiId': cChebi, 'PubchemId': cPubchem, 'boundaryCondition': cBoundaryCondition, 'chemicalFormula': cChemicalFormula, 'Inchi': cInchi})
df.to_csv(os.path.join(OUTDIR, dfmetsInfo + '_' + timeStamp + '.csv'), sep = '\t', index = False)

dfChebiDb = pd.read_csv(os.path.join(RAWDIR, 'chebi_database_accession_20201216.tsv'), sep = '\t', dtype = {'ID': str, 'COMPOUND_ID': str, 'ACCESSION_NUMBER': str})
dfChebiRelations = pd.read_csv(os.path.join(RAWDIR, 'chebi_relation_20201216.tsv'), sep = '\t', dtype = {'ID': str, 'TYPE': str, 'INIT_ID': str, 'FINAL_ID': str})

dfChebiInchi = pd.read_csv(os.path.join(RAWDIR, 'chebiId_inchi_20201216.tsv'), sep = '\t', dtype=str)
dfChebiInchi['InChI_splitted'] = dfChebiInchi.InChI.str.split('/')
dfChebiInchi['InChI_formula'] = dfChebiInchi.InChI.str.split('/').str[1]

if testModel == '2' or testModel == '3':
    df['toMatch'] = df.Name.str.replace("\[(\w*\s*)+\]$", "")
    lMets2Search = list(df['toMatch'])
else:
    lMets2Search = list(df['Name'])
lMets2Search = gL.unique(lMets2Search)
print('len lMets2Search\t', len(lMets2Search))


### METACYC
dMapping2_AllResultsMetaCyc = {}
dMapping2metacyc = {}
dfMetaCyc = pd.read_csv(os.path.join(OUTDIR,'metacyc_compounds_20201216152513.csv'), sep='\t', usecols=['Id', 'Name'])

choicesMetaCyc = dfMetaCyc.Name.to_list()
for toMatch in lMets2Search:
    matches = process.extract(toMatch, choicesMetaCyc, limit=10)

    if toMatch in dMapping2_AllResultsMetaCyc:
        dMapping2_AllResultsMetaCyc[toMatch] += matches
    else:
        dMapping2_AllResultsMetaCyc[toMatch] = matches

dfMatches = pd.DataFrame(dMapping2_AllResultsMetaCyc).T
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingMetaCyc_allResults_' + timeStamp + '.tsv'), sep='\t')

### KEGG
dMapping2_AllResultskeggC = {}
dMapping2metacyc_keggC = {}

dMapping2_AllResultskeggG = {}
dMapping2metacyc_keggG = {}
dfKeggC = pd.read_csv(os.path.join(OUTDIR,'kegg_compounds_20201216150802.csv'), sep='\t', usecols=['Id', 'Name'])
dfKeggG = pd.read_csv(os.path.join(OUTDIR,'kegg_glycans_20201216150802.csv'), sep='\t', usecols=['Id', 'Name'])
dfKeggC['Name'] = dfKeggC['Name'].apply(literal_eval)
dfKeggG['Name'] = dfKeggG['Name'].apply(literal_eval)
dfKeggC = dfKeggC.explode('Name')
dfKeggG = dfKeggG.explode('Name')

choicesKeggC = dfKeggC.Name.to_list() + dfKeggC.Id.to_list()
for toMatch in lMets2Search:
    matches = process.extract(toMatch, choicesKeggC, limit=10)

    if toMatch in dMapping2_AllResultskeggC:
        dMapping2_AllResultskeggC[toMatch] += matches
    else:
        dMapping2_AllResultskeggC[toMatch] = matches

choicesKeggG = dfKeggG.Name.to_list() + dfKeggG.Id.to_list()
for toMatch in lMets2Search:
    matches = process.extract(toMatch, choicesKeggG, limit=10)

    if toMatch in dMapping2_AllResultskeggG:
        dMapping2_AllResultskeggG[toMatch] += matches
    else:
        dMapping2_AllResultskeggG[toMatch] = matches

dfMatches = pd.DataFrame(dMapping2_AllResultskeggC).T
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingKeggC_allResults_' + timeStamp + '.tsv'), sep='\t')

dfMatches = pd.DataFrame(dMapping2_AllResultskeggG).T
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingKeggG_allResults_' + timeStamp + '.tsv'), sep='\t')

## CHEBI
dMapping2_AllResultsChebi = {}
dMapping2metacyc_chebi = {}
dfChebiNames = pd.read_csv(os.path.join(RAWDIR, 'chebi_names_20201216.tsv'), sep = '\t', dtype=str)
dfChebiUniprot = pd.read_csv(os.path.join(RAWDIR, 'chebi_uniprot_20201216.tsv'), sep = '\t', header=None, dtype=str, names=['ID','NAME_wSymbols', 'NAME'])
dfChebiCompounds =  pd.read_csv(os.path.join(RAWDIR, 'chebi_compounds_20201216.tsv'), sep = '\t', dtype=str)

choicesChebi = dfChebiNames.NAME.to_list() + dfChebiUniprot.NAME.to_list() + dfChebiCompounds.NAME.to_list()
for toMatch in lMets2Search:
    matches = process.extract(toMatch, choicesChebi, limit=10)

    if toMatch in dMapping2_AllResultsChebi:
        dMapping2_AllResultsChebi[toMatch] += matches
    else:
        dMapping2_AllResultsChebi[toMatch] = matches

dfMatches = pd.DataFrame(dMapping2_AllResultsChebi).T
dfMatches.to_csv(os.path.join(OUTDIR, outputFileName + '_mappingChebi_allResults_' + timeStamp + '.tsv'), sep='\t')

end = time.time()
print('elapsed time\t', end-start, '\tseconds')
