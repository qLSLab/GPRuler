# -*- coding: utf-8 -*-
import os
import sys
import pandas as pd
import genericLib as gL
from ast import literal_eval
import metabolitesLib as metL

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
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

# extracting metabolites info from the input model
df = pd.read_csv(os.path.join(OUTDIR, prefix_modelName + '_metabolites.csv'), sep = '\t', dtype = str)

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
dfMatches = metL.applyFW2Db(lMets2Search, choicesMetaCyc)
dfMatches.to_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingMetaCyc_allResults.tsv'), sep='\t')

# query KEGG compounds
dfKeggC = pd.read_csv(os.path.join(RAWDIR,'kegg_compounds_20201216150802.csv'), sep='\t', usecols=['Id', 'Name'])
dfKeggG = pd.read_csv(os.path.join(RAWDIR,'kegg_glycans_20201216150802.csv'), sep='\t', usecols=['Id', 'Name'])
dfKeggC['Name'] = dfKeggC['Name'].apply(literal_eval)
dfKeggG['Name'] = dfKeggG['Name'].apply(literal_eval)
dfKeggC = dfKeggC.explode('Name')
dfKeggG = dfKeggG.explode('Name')
choicesKeggC = dfKeggC.Name.to_list() + dfKeggC.Id.to_list()
dfMatches = metL.applyFW2Db(lMets2Search, choicesKeggC)
dfMatches.to_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingKeggC_allResults.tsv'), sep='\t')

choicesKeggG = dfKeggG.Name.to_list() + dfKeggG.Id.to_list()
dfMatches = metL.applyFW2Db(lMets2Search, choicesKeggG)
dfMatches.to_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingKeggG_allResults.tsv'), sep='\t')

# query ChEBI compounds
dfChebiNames = pd.read_csv(os.path.join(RAWDIR, 'chebi_names_20201216.tsv.bz2'), sep = '\t', compression='bz2', dtype=str)
dfChebiUniprot = pd.read_csv(os.path.join(RAWDIR, 'chebi_uniprot_20201216.tsv.bz2'), sep = '\t',compression='bz2', header=None, dtype=str, names=['ID','NAME_wSymbols', 'NAME'])
dfChebiCompounds =  pd.read_csv(os.path.join(RAWDIR, 'chebi_compounds_20201216.tsv.bz2'), sep = '\t', compression='bz2', dtype=str)

choicesChebi = dfChebiNames.NAME.to_list() + dfChebiUniprot.NAME.to_list() + dfChebiCompounds.NAME.to_list()
dfMatches = metL.applyFW2Db(lMets2Search, choicesChebi)
dfMatches.to_csv(os.path.join(OUTDIR, prefix_modelName + '_mappingChebi_allResults.tsv'), sep='\t')
