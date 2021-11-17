# -*- coding: utf-8 -*-
import sys
import pandas as pd
import genericLib as gL
import metabolitesLib as metL

# setting working dirs
dDirNames = {'raw': 'rawData', 'log': 'logs', 'out': 'outputs'}

(wrkDir, dDirs) = gL.setWorkingDirs(wrkDir=None, dizDirNames=dDirNames)
RAWDIR = dDirs['raw']
OUTDIR = dDirs['out']
LOGDIR = dDirs['log']

timeStamp = gL.getTimeStamp()

# setting input data
dModelPrms = gL.loadModelParams(sys.argv, timeStamp, dDirs)
print(dModelPrms)
logStream = gL.logFileOpen(logDIR=LOGDIR,
                           timeStamp=timeStamp,
                           basename=dModelPrms['basenameStr'])
sToLog = 'Input params\nmodel filename: ' + dModelPrms['filename'] + '\n'
sToLog += 'model prefix: ' + dModelPrms['prefix'] + '\n'
sToLog += 'include Compartement: ' + str(dModelPrms['incComp']) + '\n'
gL.toLog(logStream, sToLog)

# Extract metabolites info from the input model
cId, cName, cKegg, cChebi, cPubchem, cBoundaryCondition, cChemicalFormula, cInchi = metL.getMetsInfoGEM(
    dModelPrms['path'])
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
# write metabolites report file
fileName = gL.pathFilename(OUTDIR, dModelPrms['prefix'] + '_metabolites.csv')
df.to_csv(fileName, sep='\t', index=False)

# Infer metabolites identifiers
fileName = gL.pathJoinOrExit(RAWDIR, gL.dChebiFiles['names'])
dfChebiNames = pd.read_csv(fileName, sep='\t', compression='bz2', dtype=str)

fileName = gL.pathJoinOrExit(RAWDIR, gL.dChebiFiles['uniprot'])
dfChebiUniprot = pd.read_csv(fileName,
                             sep='\t',
                             header=None,
                             dtype=str,
                             compression='bz2',
                             names=['ID', 'NAME_wSymbols', 'NAME'])

fileName = gL.pathJoinOrExit(RAWDIR, gL.dChebiFiles['compounds'])
dfChebiCompounds = pd.read_csv(fileName,
                               sep='\t',
                               compression='bz2',
                               dtype=str)

fileName = gL.pathJoinOrExit(RAWDIR, gL.dChebiFiles['accession'])
dfChebiDb = pd.read_csv(fileName,
                        sep='\t',
                        compression='bz2',
                        dtype={
                            'ID': str,
                            'COMPOUND_ID': str,
                            'ACCESSION_NUMBER': str
                        })

fileName = gL.pathJoinOrExit(RAWDIR, gL.dChebiFiles['relation'])
dfChebiRelations = pd.read_csv(fileName,
                               sep='\t',
                               compression='bz2',
                               dtype={
                                   'ID': str,
                                   'TYPE': str,
                                   'INIT_ID': str,
                                   'FINAL_ID': str
                               })

fileName = gL.pathJoinOrExit(RAWDIR, gL.dChebiFiles['inchi'])
dfChebiInchi = pd.read_csv(fileName, sep='\t', compression='bz2', dtype=str)
dfChebiInchi['InChI_splitted'] = dfChebiInchi.InChI.str.split('/')
dfChebiInchi['InChI_formula'] = dfChebiInchi.InChI.str.split('/').str[1]

if dModelPrms['incComp'] is True:
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
    if dModelPrms['incComp'] is True:
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
fileName = gL.pathFilename(
    OUTDIR, dModelPrms['prefix'] + '_metabolites_wInferredIds.csv')
dfMet2Ids.to_csv(fileName, sep='\t', index=False)

logStream.close()
