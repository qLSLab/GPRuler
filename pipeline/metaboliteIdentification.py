import keggLib as kL
import sys
import genericLib as gL
import os
import pandas as pd
import gprLib as gprL
import xmlLib as xL
import time
import re

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
elif testModel == 'y7':
    ## Yeast 7
    modelXml = 'yeast_7.6_cobra'
    dfmetsInfo = 'y7_metabolites'
elif testModel == 'y8':
    ## Yeast 8
    modelXml = 'yeast8'
    dfmetsInfo = 'y8_metabolites'
elif testModel == 'hmr':
    ## HMRcore
    modelXml = 'HMRcore_20200328_wReconNames'
    dfmetsInfo = 'hmrCore_metabolites'
elif testModel == 'ownData':
    ## specify your input data
    modelXml = ''
    dfmetsInfo = ''

# extracting metabolites info from the input model
cId, cName, cKegg, cChebi, cPubchem, cBoundaryCondition, cChemicalFormula, cInchi = xL.getMetsInfoGEM(os.path.join(RAWDIR, modelXml + ".xml"))
df = pd.DataFrame({'Id': cId, 'Name': cName, 'KeggId': cKegg, 'ChebiId': cChebi, 'PubchemId': cPubchem, 'boundaryCondition': cBoundaryCondition, 'chemicalFormula': cChemicalFormula, 'Inchi': cInchi})
df.to_csv(os.path.join(OUTDIR, dfmetsInfo + '.csv'), sep = '\t', index = False)

# infer metabolites identifiers
dfChebiNames = pd.read_csv(os.path.join(RAWDIR, 'chebi_names_20201216.tsv'), sep = '\t', dtype=str)
dfChebiUniprot = pd.read_csv(os.path.join(RAWDIR, 'chebi_uniprot_20201216.tsv'), sep = '\t', header=None, dtype=str, names=['ID','NAME_wSymbols', 'NAME'])
dfChebiCompounds =  pd.read_csv(os.path.join(RAWDIR, 'chebi_compounds_20201216.tsv'), sep = '\t', dtype=str)

dfChebiDb = pd.read_csv(os.path.join(RAWDIR, 'chebi_database_accession_20201216.tsv'), sep = '\t', dtype = {'ID': str, 'COMPOUND_ID': str, 'ACCESSION_NUMBER': str})
dfChebiRelations = pd.read_csv(os.path.join(RAWDIR, 'chebi_relation_20201216.tsv'), sep = '\t', dtype = {'ID': str, 'TYPE': str, 'INIT_ID': str, 'FINAL_ID': str})

dfChebiInchi = pd.read_csv(os.path.join(RAWDIR, 'chebiId_inchi_20201216.tsv'), sep = '\t', dtype=str)
dfChebiInchi['InChI_splitted'] = dfChebiInchi.InChI.str.split('/')
dfChebiInchi['InChI_formula'] = dfChebiInchi.InChI.str.split('/').str[1]

if testModel == 'y7' or testModel == 'y8':
    lMets2Search_complete = list(df['Name'])
    lMets2Search = []
    for el in lMets2Search_complete:
        lPositions = [p.start() for p in re.finditer('\[', el)]
        pos = lPositions[-1]
        met = el[:pos].strip()
        lMets2Search.append(met)
elif testModel == 'recon' or testModel == 'hmr':
    lMets2Search = list(df['Name'])
lMets2Search = gL.unique(lMets2Search)
print('len lMets2Search\t', len(lMets2Search))

dfChebiNames['NAME'] = dfChebiNames['NAME'].str.lower()
dfChebiCompounds['NAME'] = dfChebiCompounds['NAME'].str.lower()
dfChebiUniprot['NAME'] = dfChebiUniprot['NAME'].str.lower()

dizMet2Ids = {}
for met in lMets2Search:
    keggId = ''
    if testModel == 'y7' or testModel == 'y8':
        dfMet = df[df['Name'].str.startswith(met + ' ' + '[')]
    else:
        dfMet = df[df['Name'] == met]
    dfMet = dfMet.reset_index(drop = True)
    inchiOriginal = dfMet.iloc[0]['Inchi']
    keggId = kL.extractKeggIdComp(met.lower(), dfChebiNames, dfChebiDb, dfChebiRelations, dfChebiUniprot, dfChebiCompounds, dfChebiInchi, inchiOriginal)
    dizMet2Ids[met] = keggId

dfMet2Ids = pd.DataFrame(dizMet2Ids.items(), columns=['Name', 'Identifiers'])
dfMet2Ids.to_csv(os.path.join(OUTDIR, dfmetsInfo + '_wInferredIds.csv'), sep = '\t', index = False)
