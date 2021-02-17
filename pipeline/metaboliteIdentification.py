from Bio.KEGG.REST import *
import keggLib as kL
import sys
import genericLib as gL
import os
import pandas as pd
import gprLib as gprL
import xmlLib2 as xL
from ast import literal_eval
from nltk.tokenize import RegexpTokenizer
# from libchebipy._chebi_entity import ChebiEntity
import chebiLib as cL
import time
# import libchebipy
import re

def get_jaccard_sim(str1, str2):
    a = set(str1.split())
    b = set(str2.split())
    c = a.intersection(b)
    return float(len(c)) / (len(a) + len(b) - len(c))

timeStamp = gL.getTimeStamp()

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]
FIGUREDIR = workingDirs[4]
DFDIR = workingDirs[5]

## INPUTS:
# - metabolites list

start = time.time()

testModel = sys.argv[1]
if testModel == '1':
    ## Recon 3
    modelXml = 'Recon3D_301_20200923'
    dfmetsInfo = 'recon3D_metabolites'
    # dfmetsInfo = 'Recon3D_301_enriched'
elif testModel == '2':
    ## Yeast 7
    modelXml = 'yeast_7.6_cobra'
    # dfmetsInfo = 'y7_mets_20200814'
    dfmetsInfo = 'y7_metabolites'
elif testModel == '3':
    ## Yeast 8
    modelXml = 'yeast8'
    dfmetsInfo = 'y8_metabolites'
elif testModel == '4':
    ## HMRcore
    modelXml = 'HMRcore_20200328_wReconNames'
    dfmetsInfo = 'hmrCore_metabolites'
# elif testModel == '5':
#     ## Ecoli
#     modelXml = 'iML1515-ROS'
#     # dfmetsInfo = 'eColi_metabolites_enriched'
#     dfmetsInfo = 'eColi_metabolites'



## Step1. estrazione info mets da modello + inferimento dei kegg id
start = time.time()
cId, cName, cKegg, cChebi, cPubchem, cBoundaryCondition, cChemicalFormula, cInchi = xL.getMetsInfoGEM(os.path.join(RAWDIR, modelXml + ".xml"))
df = pd.DataFrame({'Id': cId, 'Name': cName, 'KeggId': cKegg, 'ChebiId': cChebi, 'PubchemId': cPubchem, 'boundaryCondition': cBoundaryCondition, 'chemicalFormula': cChemicalFormula, 'Inchi': cInchi})
df.to_csv(os.path.join(OUTDIR, dfmetsInfo + '_' + timeStamp + '.csv'), sep = '\t', index = False)


dfChebiNames = pd.read_csv(os.path.join(RAWDIR, 'chebi_names_20201216.tsv'), sep = '\t', dtype=str)
dfChebiUniprot = pd.read_csv(os.path.join(RAWDIR, 'chebi_uniprot_20201216.tsv'), sep = '\t', header=None, dtype=str, names=['ID','NAME_wSymbols', 'NAME'])
dfChebiCompounds =  pd.read_csv(os.path.join(RAWDIR, 'chebi_compounds_20201216.tsv'), sep = '\t', dtype=str)

dfChebiDb = pd.read_csv(os.path.join(RAWDIR, 'chebi_database_accession_20201216.tsv'), sep = '\t', dtype = {'ID': str, 'COMPOUND_ID': str, 'ACCESSION_NUMBER': str})
dfChebiRelations = pd.read_csv(os.path.join(RAWDIR, 'chebi_relation_20201216.tsv'), sep = '\t', dtype = {'ID': str, 'TYPE': str, 'INIT_ID': str, 'FINAL_ID': str})

dfChebiInchi = pd.read_csv(os.path.join(RAWDIR, 'chebiId_inchi_20201216.tsv'), sep = '\t', dtype=str)
dfChebiInchi['InChI_splitted'] = dfChebiInchi.InChI.str.split('/')
dfChebiInchi['InChI_formula'] = dfChebiInchi.InChI.str.split('/').str[1]

if testModel == '2' or testModel == '3':
    lMets2Search_complete = list(df['Name'])
    lMets2Search = []
    for el in lMets2Search_complete:
        lPositions = [p.start() for p in re.finditer('\[', el)]
        pos = lPositions[-1]
        met = el[:pos].strip()
        lMets2Search.append(met)
elif testModel == '1' or testModel == '4':
    lMets2Search = list(df['Name'])
lMets2Search = gL.unique(lMets2Search)
print('len lMets2Search\t', len(lMets2Search))

dfChebiNames['NAME'] = dfChebiNames['NAME'].str.lower()
dfChebiCompounds['NAME'] = dfChebiCompounds['NAME'].str.lower()
dfChebiUniprot['NAME'] = dfChebiUniprot['NAME'].str.lower()

dizMet2Ids = {}
for met in lMets2Search:
    keggId = ''
    #############################
    ## pulitura nome
    # met = "Decanoate (N-C10:0)"
    ###############################

    if testModel == '2' or testModel == '3':
        dfMet = df[df['Name'].str.startswith(met + ' ' + '[')]
    else:
        dfMet = df[df['Name'] == met]
    print('met:\t', met)
    dfMet = dfMet.reset_index(drop = True)
    print('Original keggId\t', dfMet.iloc[0]['KeggId'])
    inchiOriginal = dfMet.iloc[0]['Inchi']
    print('inchiOriginal\n', inchiOriginal, '\n')
    keggId = kL.extractKeggIdComp(met.lower(), dfChebiNames, dfChebiDb, dfChebiRelations, dfChebiUniprot, dfChebiCompounds, dfChebiInchi, inchiOriginal)
    ## ricerca
    print('----------Kegg ID inferito:\t', keggId, '\n')
    dizMet2Ids[met] = keggId
    print('\n\n')

dfMet2Ids = pd.DataFrame(dizMet2Ids.items(), columns=['Name', 'Identifiers'])
dfMet2Ids.to_csv(os.path.join(OUTDIR, dfmetsInfo + '_wInferredIds_' + timeStamp + '.csv'), sep = '\t', index = False)
end = time.time()
print('elapsed time\t', end-start, '\tseconds')
