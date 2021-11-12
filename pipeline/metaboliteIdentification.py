# -*- coding: utf-8 -*-
import sys
import pandas as pd
import genericLib as gL
import metabolitesLib as metL

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]

# setting input data
testModel = sys.argv[1]
if testModel == 'recon':
    ## Recon 3
    modelXmlFile = 'Recon3D_301_20200923'
    prefix_modelName = 'recon3D'
    includeCompartment = False
elif testModel == 'y7':
    ## Yeast 7
    modelXmlFile = 'yeast_7.6_cobra'
    prefix_modelName = 'y7'
    includeCompartment = True
elif testModel == 'y8':
    ## Yeast 8
    modelXmlFile = 'yeast8'
    prefix_modelName = 'y8'
    includeCompartment = True
elif testModel == 'hmr':
    ## HMRcore
    modelXmlFile = 'HMRcore_20200328_wReconNames'
    prefix_modelName = 'hmrCore'
    includeCompartment = False
elif testModel == 'ownData':
    ## specify your input data
    modelXmlFile = 'toymodel'
    prefix_modelName = 'toy_metabolites'
    includeCompartment = False

# Extract metabolites info from the input model
fileName = gL.pathFilename(RAWDIR, modelXmlFile + ".xml")
cId, cName, cKegg, cChebi, cPubchem, cBoundaryCondition, cChemicalFormula, cInchi = metL.getMetsInfoGEM(
    fileName)
df = pd.DataFrame({
    'Id': cId,
    'Name': cName,
    'KeggId': cKegg,
    'ChebiId': cChebi,
    'PubchemId': cPubchem,
    'boundaryCondition': cBoundaryCondition,
    'chemicalFormula': cChemicalFormula,
    'Inchi': cInchi
})
fileName = gL.pathFilename(OUTDIR, prefix_modelName + '_metabolites.csv')
df.to_csv(fileName, sep='\t', index=False)

# Infer metabolites identifiers
fileName = gL.pathFilename(RAWDIR, 'chebi_names_20201216.tsv.bz2')
dfChebiNames = pd.read_csv(fileName, sep='\t', compression='bz2', dtype=str)

fileName = gL.pathFilename(RAWDIR, 'chebi_uniprot_20201216.tsv.bz2')
dfChebiUniprot = pd.read_csv(fileName,
                             sep='\t',
                             header=None,
                             dtype=str,
                             compression='bz2',
                             names=['ID', 'NAME_wSymbols', 'NAME'])

fileName = gL.pathFilename(RAWDIR, 'chebi_compounds_20201216.tsv.bz2')
dfChebiCompounds = pd.read_csv(fileName,
                               sep='\t',
                               compression='bz2',
                               dtype=str)

fileName = gL.pathFilename(RAWDIR, 'chebi_database_accession_20201216.tsv.bz2')
dfChebiDb = pd.read_csv(fileName,
                        sep='\t',
                        compression='bz2',
                        dtype={
                            'ID': str,
                            'COMPOUND_ID': str,
                            'ACCESSION_NUMBER': str
                        })
fileName = gL.pathFilename(RAWDIR, 'chebi_relation_20201216.tsv.bz2')
dfChebiRelations = pd.read_csv(fileName,
                               sep='\t',
                               compression='bz2',
                               dtype={
                                   'ID': str,
                                   'TYPE': str,
                                   'INIT_ID': str,
                                   'FINAL_ID': str
                               })

fileName = gL.pathFilename(RAWDIR, 'chebiId_inchi_20201216.tsv.bz2')
dfChebiInchi = pd.read_csv(fileName, sep='\t', compression='bz2', dtype=str)
dfChebiInchi['InChI_splitted'] = dfChebiInchi.InChI.str.split('/')
dfChebiInchi['InChI_formula'] = dfChebiInchi.InChI.str.split('/').str[1]

if includeCompartment is True:
    df['toMatch'] = df.Name.str.replace("\[(\w*\s*)+\]$", "", regex=True)
    df['toMatch'] = df.toMatch.str.strip()
    lMets2Search = list(df['toMatch'])
else:
    lMets2Search = list(df['Name'])
lMets2Search = gL.unique(lMets2Search)

dfChebiNames['NAME'] = dfChebiNames['NAME'].str.lower()
dfChebiCompounds['NAME'] = dfChebiCompounds['NAME'].str.lower()
dfChebiUniprot['NAME'] = dfChebiUniprot['NAME'].str.lower()

dizMet2Ids = {}
for met in lMets2Search:
    keggId = ''
    if includeCompartment is True:
        dfMet = df[df['Name'].str.startswith(met + ' ' + '[')]
    else:
        dfMet = df[df['Name'] == met]
    dfMet = dfMet.reset_index(drop=True)
    if dfMet.empty == False:
        inchiOriginal = dfMet.iloc[0]['Inchi']
    else:
        inchiOriginal = ''
    keggId = metL.extractKeggIdComp(met.lower(), dfChebiNames, dfChebiDb,
                                    dfChebiRelations, dfChebiUniprot,
                                    dfChebiCompounds, dfChebiInchi,
                                    inchiOriginal)
    dizMet2Ids[met] = keggId

dfMet2Ids = pd.DataFrame(dizMet2Ids.items(), columns=['Name', 'Identifiers'])
fileName = gL.pathFilename(OUTDIR, '_metabolites_wInferredIds.csv')
dfMet2Ids.to_csv(fileName, sep='\t', index=False)
