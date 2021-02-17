
from Bio.KEGG.REST import *
import keggLib as kL
import sys
import genericLib as gL
import os
import pandas as pd
import gprLib as gprL
import xmlLib as xL
from ast import literal_eval
from nltk.tokenize import RegexpTokenizer
import chebiLib as cL
import time
import re
from ast import literal_eval
import cobra as cb
import numpy as np
import testKeggCompund as testKLib

def explodeFilterAndDrop(df, oldCol, lMets):
    df[oldCol + '_exploded'] = df[oldCol].values
    df = df.explode(oldCol + '_exploded')
    df = df[df[oldCol + '_exploded'].notna()]
    # df = df.reset_index(drop = True)
    newDf = df[df[oldCol + '_exploded'].isin(lMets)]
    # newDf = newDf.reset_index(drop = True)
    newDf[oldCol + '_t'] = newDf[oldCol].apply(tuple)
    # newDf = newDf.drop_duplicates(subset=[oldCol + '_t'])
    return newDf

def filterwExtractRows(df, dfR, dfP,  colNameR, colNameP):
    lDfs = []
    df['filter'] = df.apply(lambda row: extractRows(dfR, row[colNameR]), axis=1)
    if df['filter'].any() == True:
        result = df[df['filter'] == True]
        lDfs.append(result)

    df['filter'] = df.apply(lambda row: extractRows(dfR, row[colNameP]), axis=1)
    if df['filter'].any() == True:
        result = df[df['filter'] == True]
        lDfs.append(result)

    df['filter'] = df.apply(lambda row: extractRows(dfP, row[colNameR]), axis=1)
    if df['filter'].any() == True:
        result = df[df['filter'] == True]
        lDfs.append(result)

    df['filter'] = df.apply(lambda row: extractRows(dfP, row[colNameP]), axis=1)
    if df['filter'].any() == True:
        result = df[df['filter'] == True]
        lDfs.append(result)
    return lDfs


def findPutativeRxns(dfEqualMets, colL, colR, db):
    lPutativeRxns = []
    if db == 'metacyc' or db == 'kegg':
        dfEqualMets['fromRxn2Putative1'] = dfEqualMets.apply(lambda row: findReaction_fromRxn2Putative(lReactants_ids, lProducts_ids, row[colL], row[colR]), axis=1)
        dfEqualMets['fromPutative2Rxn1'] = dfEqualMets.apply(lambda row: findReaction_fromPutative2Rxn(lReactants_ids, lProducts_ids, row[colL], row[colR]), axis=1)
        dfEqualMets['fromRxn2Putative2'] = dfEqualMets.apply(lambda row: findReaction_fromRxn2Putative(lProducts_ids, lReactants_ids, row[colL], row[colR]), axis=1)
        dfEqualMets['fromPutative2Rxn2'] = dfEqualMets.apply(lambda row: findReaction_fromPutative2Rxn(lProducts_ids, lReactants_ids, row[colL], row[colR]), axis=1)

    elif db == 'rhea':
        dfEqualMets['fromRxn2Putative1'] = dfEqualMets.apply(lambda row: findReaction_rhea(lReactants_ids, lProducts_ids, row[colL], row[colR]), axis=1)
        dfEqualMets['fromPutative2Rxn1'] = dfEqualMets.apply(lambda row: findReaction_rhea(row[colL], row[colR], lReactants_ids, lProducts_ids), axis=1)
        dfEqualMets['fromRxn2Putative2'] = dfEqualMets.apply(lambda row: findReaction_rhea(lReactants_ids, lProducts_ids, row[colR], row[colL]), axis=1)
        dfEqualMets['fromPutative2Rxn2'] = dfEqualMets.apply(lambda row: findReaction_rhea(row[colR], row[colL], lReactants_ids, lProducts_ids), axis=1)

    dfMatchesEqual_L = dfEqualMets[((dfEqualMets['fromRxn2Putative1'] == True) & (dfEqualMets['fromPutative2Rxn1'] == True)) | ((dfEqualMets['fromRxn2Putative2'] == True) & (dfEqualMets['fromPutative2Rxn2'] == True))]

    if dfMatchesEqual_L.empty == False:
        lPutativeRxns = list(dfMatchesEqual_L['MetaCycId'].dropna()) + list(dfMatchesEqual_L['RheaId'].dropna()) + list(dfMatchesEqual_L['KeggId_fromKegg'].dropna()) + [el for l in list(dfMatchesEqual_L['OtherRheaId_fromKegg'].dropna()) for el in l] + [el for l in list(dfMatchesEqual_L['RheaId_master'].dropna()) for el in l] + [el for l in list(dfMatchesEqual_L['RheaId_lr'].dropna()) for el in l] + [el for l in list(dfMatchesEqual_L['RheaId_rl'].dropna()) for el in l] + [el for l in list(dfMatchesEqual_L['RheaId_bi'].dropna()) for el in l] + [el for l in list(dfMatchesEqual_L['KeggId_y'].dropna()) for sublist in l for el in sublist] + [el for l in list(dfMatchesEqual_L['MetaCycId_fromRhea'].dropna()) for el in l.split(',')]
        lPutativeRxns = gL.unique(lPutativeRxns)

    return lPutativeRxns


def findRxnsAfterFilter(lDfs2Concat, colNameS, colNameP, db):
    lPutativeRxns = []
    dfAllDBs_copy_All = pd.concat(lDfs2Concat)
    # print('IS iN dfAllDBs_copy_All R11680\t', 'R11680' in list(dfAllDBs_copy_All['KeggId_fromKegg']))

    dfAllDBs_copy_All[colNameS + '_t'] = dfAllDBs_copy_All[colNameS].apply(tuple)
    dfAllDBs_copy_All[colNameP + '_t'] = dfAllDBs_copy_All[colNameP].apply(tuple)
    # dfAllDBs_copy_All = dfAllDBs_copy_All.drop_duplicates(subset=[colNameS + '_t',colNameP + '_t'])
    # dfAllDBs_copy_All = dfAllDBs_copy_All.reset_index(drop = True)
    # print('dfFINALE kegg\t', dfAllDBs_copy_All.shape)

    lDfs = filterwExtractRows(dfAllDBs_copy_All, dfR, dfP, colNameS, colNameP)
    if len(lDfs) != 0:
        df = pd.concat(lDfs)
        # print('IS iN df R11680\t', 'R11680' in list(df['KeggId_fromKegg']))
        if db == 'metacyc':
            if df.empty is False:
                lPutativeRxns = findPutativeRxns(df, colNameS, colNameP, 'metacyc')
        elif db == 'kegg':
            if df.empty is False:
                lPutativeRxns = findPutativeRxns(df, colNameS, colNameP, 'kegg')
        elif db == 'rhea':
            if df.empty is False:
                lPutativeRxns = findPutativeRxns(df, colNameS[:-5], colNameP[:-5], 'rhea')

    lPutativeRxns = gL.unique(lPutativeRxns)
    return lPutativeRxns

def filterRows(dfIniziale, colNameS, colNameP, lReactants_ids_cofactors, lReactants_ids_internal, lProducts_ids_cofactors, lProducts_ids_internal):
    lDfs2Concat = []
    dfAllDBs_copy_filterL = dfIniziale.copy()
    dfAllDBs_copy_filterR = dfIniziale.copy()

    dfAllDBs_copy_filterL_filterCofactors = explodeFilterAndDrop(dfAllDBs_copy_filterL, colNameS, lReactants_ids_cofactors)
    # print('internal 1\t', dfAllDBs_copy_filterL_filterCofactors.shape)
    if dfAllDBs_copy_filterL_filterCofactors.empty is False:
        if len(lReactants_ids_internal) != 0:
            dfAllDBs_copy_filterL_filterInternal = explodeFilterAndDrop(dfAllDBs_copy_filterL_filterCofactors, colNameS, lReactants_ids_internal)
        else:
            dfAllDBs_copy_filterL_filterInternal = dfAllDBs_copy_filterL_filterCofactors.copy()
    else:
        if len(lReactants_ids_internal) != 0:
            dfAllDBs_copy_filterL_filterInternal = explodeFilterAndDrop(dfAllDBs_copy_filterL, colNameS, lReactants_ids_internal)
        else:
            dfAllDBs_copy_filterL_filterInternal = dfAllDBs_copy_filterL.copy()
    # print('internal 2\t', dfAllDBs_copy_filterL_filterInternal.shape)

    if dfAllDBs_copy_filterL_filterInternal.empty is False:
        lDfs2Concat.append(dfAllDBs_copy_filterL_filterInternal)

    dfAllDBs_copy_filterR_filterCofactors = explodeFilterAndDrop(dfAllDBs_copy_filterR, colNameP, lReactants_ids_cofactors)
    if dfAllDBs_copy_filterR_filterCofactors.empty is False:
        if len(lReactants_ids_internal) != 0:
            dfAllDBs_copy_filterR_filterInternal =  explodeFilterAndDrop(dfAllDBs_copy_filterR_filterCofactors, colNameP, lReactants_ids_internal)
        else:
            dfAllDBs_copy_filterR_filterInternal = dfAllDBs_copy_filterR_filterCofactors.copy()
    else:
        if len(lReactants_ids_internal) != 0:
            dfAllDBs_copy_filterR_filterInternal = explodeFilterAndDrop(dfAllDBs_copy_filterR, colNameP, lReactants_ids_internal)
        else:
            dfAllDBs_copy_filterR_filterInternal = dfAllDBs_copy_filterR.copy()
    # print('internal R\t', dfAllDBs_copy_filterR_filterInternal.shape)

    if dfAllDBs_copy_filterR_filterInternal.empty is False:
        lDfs2Concat.append(dfAllDBs_copy_filterR_filterInternal)


    dfAllDBs_copy_filterL = dfIniziale.copy()
    dfAllDBs_copy_filterR = dfIniziale.copy()

    dfAllDBs_copy_filterL_filterCofactors = explodeFilterAndDrop(dfAllDBs_copy_filterL, colNameS, lProducts_ids_cofactors)
    if dfAllDBs_copy_filterL_filterCofactors.empty is False:
        if len(lProducts_ids_internal) != 0:
            dfAllDBs_copy_filterL_filterInternal = explodeFilterAndDrop(dfAllDBs_copy_filterL_filterCofactors, colNameS, lProducts_ids_internal)
        else:
            dfAllDBs_copy_filterL_filterInternal = dfAllDBs_copy_filterL_filterCofactors.copy()
    else:
        if len(lProducts_ids_internal) != 0:
            dfAllDBs_copy_filterL_filterInternal = explodeFilterAndDrop(dfAllDBs_copy_filterL, colNameS, lProducts_ids_internal)
        else:
            dfAllDBs_copy_filterL_filterInternal = dfAllDBs_copy_filterL.copy()
    # print('internal L\t', dfAllDBs_copy_filterL_filterInternal.shape)

    if dfAllDBs_copy_filterL_filterInternal.empty is False:
        lDfs2Concat.append(dfAllDBs_copy_filterL_filterInternal)

    dfAllDBs_copy_filterR_filterCofactors = explodeFilterAndDrop(dfAllDBs_copy_filterR, colNameP, lProducts_ids_cofactors)
    if dfAllDBs_copy_filterR_filterCofactors.empty is False:
        if len(lProducts_ids_internal) != 0:
            dfAllDBs_copy_filterR_filterInternal = explodeFilterAndDrop(dfAllDBs_copy_filterR_filterCofactors, colNameP, lProducts_ids_internal)
        else:
            dfAllDBs_copy_filterR_filterInternal = dfAllDBs_copy_filterR_filterCofactors.copy()
    else:
        if len(lProducts_ids_internal) != 0:
            dfAllDBs_copy_filterR_filterInternal = explodeFilterAndDrop(dfAllDBs_copy_filterR, colNameP, lProducts_ids_internal)
        else:
            dfAllDBs_copy_filterR_filterInternal = dfAllDBs_copy_filterR.copy()
    # print('internal R\t', dfAllDBs_copy_filterR_filterInternal.shape)

    if dfAllDBs_copy_filterR_filterInternal.empty is False:
        lDfs2Concat.append(dfAllDBs_copy_filterR_filterInternal)

    return lDfs2Concat

def containsDesideredMets(rLeft, rRight, lista):
    rL_set = set(rLeft)
    rR_set = set(rRight)
    lista_set = set(lista)
    if (rL_set & lista_set) | (rR_set & lista_set):
        return True
    else:
        return False

def flatList(lMets):
    lMets = [el for l in lMets for el in l]
    return lMets

def removeProtons(lMets, rxnType):
    if rxnType != '|TRANSPORT|' and '|PROTON|' in lMets:
        lMets.remove('|PROTON|')
    return lMets

def removeProtonsKegg(dMets, rxnType):
    lMets = []
    for k in list(dMets.keys()):
        # print('k\t', k)
        # print('estratto\t', testKLib.KEGGCompundMatch(k), '\n')
        lMets.append(testKLib.KEGGCompundMatch(k))
    # lMets = []
    # for k in list(dMets.keys()):
    #     if '(' in k:
    #         lMets.append(k[:k.find('(')])
    #     else:
    #         lMets.append(k)
    if rxnType != '|TRANSPORT|' and 'C00080' in lMets:
        lMets.remove('C00080')
    return lMets

def removeProtonsRhea(ldizMets, rxnType):
    lMets = []
    for diz in ldizMets:
        if rxnType != '|TRANSPORT|':
            lMets.append([lcomps for lcomps in list(diz.values()) if '24636' not in lcomps and '15378' not in lcomps])
        else:
            lMets.append(list(diz.values()))
    return lMets


def findReaction_fromRxn2Putative(lReactants_ids, lProducts_ids, lLeft, lRight):
    ## vengono confrontati il primo argomento col terzo, e il secondo col quarto
    lSideSx = []
    # print('lReactants_ids\t', lReactants_ids, '\n')
    for item in lReactants_ids:
        # print('item\t', item)
        # print('lLeft\t', lLeft)
        lSideSx.append(len(gL.intersect(item, list(lLeft))) != 0)
    # print('lSideSx\t', lSideSx)

    lSideDx = []
    for item in lProducts_ids:
        lSideDx.append(len(gL.intersect(item, list(lRight))) != 0)
    # print('lSideDx\t', lSideDx)

    if all(lSideSx) == True and all(lSideDx) == True:
        return True
    else:
        return False

def findReaction_fromPutative2Rxn(lReactants_ids, lProducts_ids, lLeft, lRight):
    ## vengono confrontati il primo argomento col terzo, e il secondo col quarto
    lSideSx = []
    for item in list(lLeft):
        lSideSx_internal = []
        for lComps in lReactants_ids:
            lSideSx_internal.append(item in lComps)
        lSideSx.append(any(lSideSx_internal))
    # print('lSideSx\t', lSideSx)

    lSideDx = []
    for item in list(lRight):
        lSideDx_internal = []
        for lComps in lProducts_ids:
            lSideDx_internal.append(item in lComps)
        lSideDx.append(any(lSideDx_internal))
    # print('lSideDx\t', lSideDx)

    if all(lSideSx) == True and all(lSideDx) == True:
        return True
    else:
        return False

def applyFindReaction_rhea(lS, lP, lReactants_ids, lProducts_ids):
    idx = 0 ## cosi mi muovo su ogni reazione putativa (in rhea trovo associate a volte piu' di 1 rxn)
    foundRhea = False
    while idx < len(lS) and foundRhea == False:
        if all(r != [] for r in lS[idx]) == True and all(p != [] for p in lP[idx]) == True:
            if (len(lReactants_ids) == len(lS[idx]) and len(lProducts_ids) == len(lP[idx])) or (len(lProducts_ids) == len(lS[idx]) and len(lReactants_ids) == len(lP[idx])):
                found_fromRxn2Putative1 = findReaction_rhea(lReactants_ids, lProducts_ids, lS[idx], lP[idx])
                found_fromPutative2Rxn1 = findReaction_rhea(lS[idx], lP[idx], lReactants_ids, lProducts_ids)
                found_fromRxn2Putative2 = findReaction_rhea(lProducts_ids, lReactants_ids, lS[idx], lP[idx])
                found_fromPutative2Rxn2 = findReaction_rhea(lS[idx], lP[idx], lProducts_ids, lReactants_ids)
                if (found_fromRxn2Putative1 == True and found_fromPutative2Rxn1 == True) or (found_fromRxn2Putative2 == True and found_fromPutative2Rxn2 == True):
                    foundRhea = True
        idx += 1
    return foundRhea


def findReaction_rhea(lReactants_ids, lProducts_ids, lLeft, lRight):
    ## vengono confrontati il primo argomento col terzo, e il secondo col quarto
    lSideSx = []
    for item in lReactants_ids:
        lSideSx_internal = []
        for possibleComp in lLeft:
            lSideSx_internal.append(len(gL.intersect(item, possibleComp)) != 0)
        lSideSx.append(any(lSideSx_internal))
    # print('lSideSx\t', lSideSx)

    lSideDx = []
    for item in lProducts_ids:
        lSideDx_internal = []
        for possibleComp in lRight:
            lSideDx_internal.append(len(gL.intersect(item, possibleComp)) != 0)
        lSideDx.append(any(lSideDx_internal))
    # print('lSideDx\t', lSideDx)

    if all(lSideSx) == True and all(lSideDx) == True:
        return True
    else:
        return False

def checkEqualSubsProds(lLeft, lRight):
    # print('lLeft\t',lLeft, '\t----\t', type(lLeft))
    # print('lRight\t',lRight, '\t----\t', type(lRight))
    if len(gL.intersect(lLeft, lRight)) != 0:
        return True
    else:
        return False

def detectEqualSubsProds(lLeft, lRight):
    if len(gL.intersect(lLeft, lRight)) != 0:
        return gL.intersect(lLeft, lRight)
    else:
        return []

def removeEqualCompound(l2Remove, lMets):
    # print('cont\t', l2Remove, lMets)
    # print('diff\t', gL.difference(lMets, l2Remove))
    # print('tipo\t', type(gL.difference(lMets, l2Remove)), '\n')
    return gL.difference(lMets, l2Remove)

# def removeEqualCompoundFromRight(l2Remove, lRight):
#     return gL.difference(lRight, l2Remove)


def extractRows(df, column2Search):
    out = df.isin(column2Search)
    return out.any().all()


timeStamp = gL.getTimeStamp()

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]
FIGUREDIR = workingDirs[4]
DFDIR = workingDirs[5]
LOGDIR = workingDirs[7]

## INPUTS:
# - metabolites list

start = time.time()

testModel = sys.argv[1]
if testModel == 'recon3':
    ## Recon 3
    modelXml = 'Recon3D_301_20200923'
    dfmetsInfo = 'recon3D_metabolites_20201218172131'
    dfmetsIds = 'recon3_mappingFuzzyAndClassic_20210119085007'
    dfrxnsInfo = 'recon3D_reactions'
elif testModel == 'y7':
    ## Yeast 7
    modelXml = 'yeast_7.6_cobra'
    dfmetsInfo = 'y7_metabolites_20201219101924'
    dfmetsIds = 'yeast7_mappingFuzzyAndClassic_20210113162642'
    dfrxnsInfo = 'y7_reactions'
elif testModel == 'y8':
    ## Yeast 8
    modelXml = 'yeast8'
    dfmetsInfo = 'y8_metabolites_20201219102046'
    dfmetsIds = 'yeast8_mappingFuzzyAndClassic_20210113162740'
    dfrxnsInfo = 'y8_reactions'
elif testModel == 'hmr':
    ## HMRcore
    modelXml = 'HMRcore_20200328_wReconNames'
    dfmetsInfo = 'hmrCore_metabolites_20210112092737'
    dfmetsIds = 'hmrCore_mappingFuzzyAndClassic_20210112165030'
    dfrxnsInfo = 'hmrCore_reactions'
# elif testModel == '4':
#     ## Ecoli
#     modelXml = 'iML1515-ROS'
#     # dfmetsInfo = 'eColi_metabolites_enriched'
#     dfmetsInfo = 'eColi_metabolites'


## Step 1. arricchire dfMetsFromModel con le info da dfMetsIdentifiers
## df metacyc
dfmetacyc = pd.read_csv(os.path.join(OUTDIR, 'metacyc_compounds_20201216152513.csv'), sep = '\t', dtype=str)

## df kegg compound - glycan
dfkeggC = pd.read_csv(os.path.join(OUTDIR, 'kegg_compounds_20201216150802.csv'), sep = '\t', dtype=str)
dfkeggC['ChebiId'] = dfkeggC['ChebiId'].apply(literal_eval)
dfkeggG = pd.read_csv(os.path.join(OUTDIR, 'kegg_glycans_20201216150802.csv'), sep = '\t', dtype=str)
dfkeggG['ChebiId'] = dfkeggG['ChebiId'].apply(literal_eval)

dfMetsFromModel = pd.read_csv(os.path.join(OUTDIR, dfmetsInfo + '.csv'), sep = '\t', dtype=str)

dfMetsIdentifiers = pd.read_csv(os.path.join(OUTDIR, dfmetsIds + '.tsv'), sep = '\t')
dfMetsIdentifiers['Name'] = dfMetsIdentifiers['Name'].str.strip()
dfMetsIdentifiers['Identifiers'] = dfMetsIdentifiers['Identifiers'].apply(literal_eval)
dizMetsIdentifiers = dfMetsIdentifiers.set_index('Name')['Identifiers'].to_dict()
# print(dizMetsIdentifiers)

for k in dizMetsIdentifiers:
    lCompleteIds = []
    # print(k, '\t', dizMetsIdentifiers[k])
    lCompleteIds += dizMetsIdentifiers[k]
    metacDb = dfmetacyc[dfmetacyc['Id'].isin(dizMetsIdentifiers[k])]
    if metacDb.empty is False:
        lCompleteIds += list(metacDb['ChebiId'].dropna())
        lCompleteIds += list(metacDb['KeggId'].dropna())
    chebiDb = dfmetacyc[dfmetacyc['ChebiId'].isin(dizMetsIdentifiers[k])]
    if chebiDb.empty is False:
        lCompleteIds += list(chebiDb['Id'].dropna())
        lCompleteIds += list(chebiDb['KeggId'].dropna())
    keggDb = dfmetacyc[dfmetacyc['KeggId'].isin(dizMetsIdentifiers[k])]
    if keggDb.empty is False:
        lCompleteIds += list(keggDb['Id'].dropna())
        lCompleteIds += list(keggDb['ChebiId'].dropna())

    keggDbC = dfkeggC[dfkeggC['Id'].isin(dizMetsIdentifiers[k])]
    if keggDbC.empty is False:
        lCompleteIds += [el for l in list(keggDbC['ChebiId'].dropna()) for el in l]

    dfkeggC_filter =dfkeggC[pd.DataFrame(dfkeggC.ChebiId.tolist()).isin(dizMetsIdentifiers[k]).any(1).values]
    if dfkeggC_filter.empty is False:
        lCompleteIds += list(dfkeggC_filter['Id'].dropna())

    keggDbG = dfkeggG[dfkeggG['Id'].isin(dizMetsIdentifiers[k])]
    if keggDbG.empty is False:
        lCompleteIds += [el for l in list(keggDbG['ChebiId'].dropna()) for el in l]

    keggDbG_filter =keggDbG[pd.DataFrame(keggDbG.ChebiId.tolist()).isin(dizMetsIdentifiers[k]).any(1).values]
    if keggDbG_filter.empty is False:
        lCompleteIds += list(keggDbG_filter['Id'].dropna())
    lCompleteIds = gL.unique(lCompleteIds)
    dizMetsIdentifiers[k] = lCompleteIds

lOriginalEmpty_newEmpty = 0
lOriginalEmpty_newFound = 0
lOriginalYes_newFound = 0
lOriginalYes_newNotFound = 0
lOriginalYes_newEmpty = 0
total = 0


if testModel == 'y7' or testModel == 'y8':
    dfMetsFromModel['Name'] = dfMetsFromModel.Name.str.replace("\[(\w*\s*)+\]$", "")

lIdentifiers = []
for row in dfMetsFromModel.itertuples():
    # print('--' + row.Name.strip() + '---')
    lInferredIds = dizMetsIdentifiers[row.Name.strip()]
    lIdentifiers.append(lInferredIds)
    # print('lInferredIds\t', lInferredIds)

    lIdsWithinModel = []
    if pd.isna(row.KeggId) is False:
        lIdsWithinModel += [row.KeggId]
    if pd.isna(row.ChebiId) is False:
        lIdsWithinModel += [row.ChebiId]
    # print('lIdsWithinModel\t', lIdsWithinModel)

    total += 1
    if lIdsWithinModel == [] and lInferredIds == []:
        lOriginalEmpty_newEmpty += 1
    elif lIdsWithinModel == [] and lInferredIds != []:
        lOriginalEmpty_newFound += 1
    elif lIdsWithinModel != [] and lInferredIds != []:
        if len(gL.intersect(lInferredIds, lIdsWithinModel)) != 0:
            lOriginalYes_newFound += 1
        else:
            lOriginalYes_newNotFound += 1
    elif lIdsWithinModel != [] and lInferredIds == []:
        lOriginalYes_newEmpty += 1

print('total\t', total)
print('lOriginalEmpty_newEmpty\t', lOriginalEmpty_newEmpty)
print('lOriginalEmpty_newFound\t', lOriginalEmpty_newFound)
print('lOriginalYes_newFound\t', lOriginalYes_newFound)
print('lOriginalYes_newNotFound\t', lOriginalYes_newNotFound)
print('lOriginalYes_newEmpty\t', lOriginalYes_newEmpty)

dfMetsFromModel['lIdentifiers'] = lIdentifiers
dfMetsFromModel.to_csv(os.path.join(OUTDIR, dfmetsInfo + '_enriched.csv'), sep = '\t', index = False)


## db chebi con tutti i parentali
dfChebiFormula = pd.read_csv(os.path.join(RAWDIR, 'chebi_chemical_data_20201216.tsv'), sep = '\t', dtype=str)

dfChebiCompounds_exploded = pd.read_csv(os.path.join(OUTDIR, 'chebi_compounds_20201216153117_exploded.csv'), sep = '\t', dtype=str) ## 1
dfChebiCompounds_exploded['ParentalChebiIds'] = dfChebiCompounds_exploded['ParentalChebiIds'].apply(literal_eval)
dfChebiCompounds_exploded['AllChebiIds'] = dfChebiCompounds_exploded['AllChebiIds'].apply(literal_eval)
dfChebiCompounds_exploded['ParentalKeggC'] = dfChebiCompounds_exploded['ParentalKeggC'].apply(literal_eval)
dfChebiCompounds_exploded['ParentalKeggG'] = dfChebiCompounds_exploded['ParentalKeggG'].apply(literal_eval)
dfChebiCompounds_exploded['ParentalMetacyc'] = dfChebiCompounds_exploded['ParentalMetacyc'].apply(literal_eval)
dfChebiCompounds_exploded['AllKeggC'] = dfChebiCompounds_exploded['AllKeggC'].apply(literal_eval)
dfChebiCompounds_exploded['AllKeggG'] = dfChebiCompounds_exploded['AllKeggG'].apply(literal_eval)
dfChebiCompounds_exploded['AllMetacyc'] = dfChebiCompounds_exploded['AllMetacyc'].apply(literal_eval)

dfChebiCompounds = pd.read_csv(os.path.join(OUTDIR, 'chebi_compounds_20201216153117.csv'), sep = '\t', dtype=str) ## 1
dfChebiCompounds['ParentalChebiIds'] = dfChebiCompounds['ParentalChebiIds'].apply(literal_eval)
dfChebiCompounds['AllChebiIds'] = dfChebiCompounds['AllChebiIds'].apply(literal_eval)
dfChebiCompounds['ParentalKeggC'] = dfChebiCompounds['ParentalKeggC'].apply(literal_eval)
dfChebiCompounds['ParentalKeggG'] = dfChebiCompounds['ParentalKeggG'].apply(literal_eval)
dfChebiCompounds['ParentalMetacyc'] = dfChebiCompounds['ParentalMetacyc'].apply(literal_eval)
dfChebiCompounds['AllKeggC'] = dfChebiCompounds['AllKeggC'].apply(literal_eval)
dfChebiCompounds['AllKeggG'] = dfChebiCompounds['AllKeggG'].apply(literal_eval)
dfChebiCompounds['AllMetacyc'] = dfChebiCompounds['AllMetacyc'].apply(literal_eval)

dfMetsFromModel = pd.read_csv(os.path.join(OUTDIR, dfmetsInfo + '_enriched.csv'), sep = '\t', dtype=str) ## 1

# dfMetsFromModel['boundaryCondition'] = dfMetsFromModel['boundaryCondition'].apply(literal_eval)
dfMetsFromModel['lIdentifiers'] = dfMetsFromModel['lIdentifiers'].apply(literal_eval)

if testModel == 'recon3':
    dfMetsFromModel['Id_converted'] = dfMetsFromModel['Id'].str.strip('M_')
    dfMetsFromModel['Id_converted'] = dfMetsFromModel['Id_converted'].str.replace('__91__', '[')
    dfMetsFromModel['Id_converted'] = dfMetsFromModel['Id_converted'].str.replace('__93', ']')
    dfMetsFromModel['Id_converted2'] = dfMetsFromModel['Id'].str.strip('M_')
    dfMetsFromModel['Id_converted2'] = dfMetsFromModel['Id_converted2'] + '__'
    # print(dfMetsFromModel.head())

dfAllDBs = pd.read_csv(os.path.join(OUTDIR, 'dfJoin_metacyc_kegg_rhea_20201218164915.csv'), sep = '\t', dtype=str) ## 3
dfAllDBs['ec_number'] = dfAllDBs['ec_number'].apply(literal_eval)
dfAllDBs['enzymatic_reaction'] = dfAllDBs['enzymatic_reaction'].fillna({i: '[]' for i in dfAllDBs.index})
dfAllDBs['enzymatic_reaction'] = dfAllDBs['enzymatic_reaction'].apply(literal_eval)
dfAllDBs['left'] = dfAllDBs['left'].apply(literal_eval)
dfAllDBs['right'] = dfAllDBs['right'].apply(literal_eval)
dfAllDBs['rxn_locations'] = dfAllDBs['rxn_locations'].apply(literal_eval)
dfAllDBs['subreactions'] = dfAllDBs['subreactions'].apply(literal_eval)
# dfAllDBs['names'] = dfAllDBs['names'].apply(literal_eval)
dfAllDBs['synonyms'] = dfAllDBs['synonyms'].apply(literal_eval)
dfAllDBs['enzymes_not_used'] = dfAllDBs['enzymes_not_used'].apply(literal_eval)
# dfAllDBs['systematic_name'] = dfAllDBs['systematic_name'].apply(literal_eval)
dfAllDBs['compartments_of_reaction'] = dfAllDBs['compartments_of_reaction'].apply(literal_eval)
dfAllDBs['enzymes_of_reaction'] = dfAllDBs['enzymes_of_reaction'].apply(literal_eval)
dfAllDBs['genes_of_reaction'] = dfAllDBs['genes_of_reaction'].apply(literal_eval)
dfAllDBs['specific_forms_of_rxn'] = dfAllDBs['specific_forms_of_rxn'].apply(literal_eval)
dfAllDBs['Subs_fromKegg'] = dfAllDBs['Subs_fromKegg'].apply(literal_eval)
dfAllDBs['Prods_fromKegg'] = dfAllDBs['Prods_fromKegg'].apply(literal_eval)

dfAllDBs['left_metacyc'] = dfAllDBs.apply(lambda row: removeProtons(row['left'], row['reaction_type']), axis=1)
dfAllDBs['right_metacyc'] = dfAllDBs.apply(lambda row: removeProtons(row['right'], row['reaction_type']), axis=1)
dfAllDBs['lSubs_fromKegg'] = dfAllDBs.apply(lambda row: removeProtonsKegg(row['Subs_fromKegg'], row['reaction_type']), axis=1)
dfAllDBs['lProds_fromKegg'] = dfAllDBs.apply(lambda row: removeProtonsKegg(row['Prods_fromKegg'], row['reaction_type']), axis=1)

dfAllDBs['OtherRheaId_fromKegg'] = dfAllDBs['OtherRheaId_fromKegg'].apply(literal_eval)
lotherrhea = []
for el in dfAllDBs['OtherRheaId_fromKegg']:
    if type(el) == list:
        lotherrhea.append(el)
    elif type(el) == int:
        lotherrhea.append([el])
dfAllDBs['OtherRheaId_fromKegg'] = lotherrhea

dfAllDBs['RheaId_master'] = dfAllDBs['RheaId_master'].fillna({i: '[]' for i in dfAllDBs.index})
dfAllDBs['RheaId_master'] = dfAllDBs['RheaId_master'].apply(literal_eval)
dfAllDBs['RheaId_lr'] = dfAllDBs['RheaId_lr'].fillna({i: '[]' for i in dfAllDBs.index})
dfAllDBs['RheaId_lr'] = dfAllDBs['RheaId_lr'].apply(literal_eval)
dfAllDBs['RheaId_rl'] = dfAllDBs['RheaId_rl'].fillna({i: '[]' for i in dfAllDBs.index})
dfAllDBs['RheaId_rl'] = dfAllDBs['RheaId_rl'].apply(literal_eval)
dfAllDBs['RheaId_bi'] = dfAllDBs['RheaId_bi'].fillna({i: '[]' for i in dfAllDBs.index})
dfAllDBs['RheaId_bi'] = dfAllDBs['RheaId_bi'].apply(literal_eval)
dfAllDBs['KeggId_y'] = dfAllDBs['KeggId_y'].fillna({i: '[]' for i in dfAllDBs.index})
dfAllDBs['KeggId_y'] = dfAllDBs['KeggId_y'].apply(literal_eval)
dfAllDBs['Outgoings'] = dfAllDBs['Outgoings'].fillna({i: '[]' for i in dfAllDBs.index})
dfAllDBs['Outgoings'] = dfAllDBs['Outgoings'].apply(literal_eval)
dfAllDBs['Reactants_name_id'] = dfAllDBs['Reactants_name_id'].fillna({i: '[]' for i in dfAllDBs.index})
dfAllDBs['Reactants_name_id'] = dfAllDBs['Reactants_name_id'].apply(literal_eval)
dfAllDBs['Products_name_id'] = dfAllDBs['Products_name_id'].fillna({i: '[]' for i in dfAllDBs.index})
dfAllDBs['Products_name_id'] = dfAllDBs['Products_name_id'].apply(literal_eval)

dfAllDBs['lReactants_name_id'] = dfAllDBs.apply(lambda row: removeProtonsRhea(row['Reactants_name_id'], row['reaction_type']), axis=1)
dfAllDBs['lProducts_name_id'] = dfAllDBs.apply(lambda row: removeProtonsRhea(row['Products_name_id'], row['reaction_type']), axis=1)
## colonna lReactants_name_id e lProducts_name_id: faccio si che ogni id sia su una riga diversa duplicando la corrispondente riga quando necessario (uso funzione explode di pandas)

list_cols = {'lReactants_name_id','lProducts_name_id'}
other_cols = list(set(dfAllDBs.columns) - set(list_cols))
exploded = [dfAllDBs[col].explode() for col in list_cols]
dfAllDBs_explode = pd.DataFrame(dict(zip(list_cols, exploded)))
dfAllDBs_explode = dfAllDBs[other_cols].merge(dfAllDBs_explode, how="right", left_index=True, right_index=True)
dfAllDBs_explode = dfAllDBs_explode.reset_index(drop = True)

dfAllDBs_explode['lReactants_name_id'] = dfAllDBs_explode['lReactants_name_id'].fillna({i: [] for i in dfAllDBs_explode.index})
dfAllDBs_explode['lProducts_name_id'] = dfAllDBs_explode['lProducts_name_id'].fillna({i: [] for i in dfAllDBs_explode.index})

dfAllDBs_explode['lReactants_name_id_flat'] = dfAllDBs_explode.apply(lambda row: flatList(row['lReactants_name_id']), axis=1)
dfAllDBs_explode['lProducts_name_id_flat'] = dfAllDBs_explode.apply(lambda row: flatList(row['lProducts_name_id']), axis=1)

## aggiungo colonna con lunghezza di alcune colonne
dfAllDBs_explode['left_metacyc_l'] = dfAllDBs_explode.left_metacyc.str.len()
dfAllDBs_explode['right_metacyc_l'] = dfAllDBs_explode.right_metacyc.str.len()

dfAllDBs_explode['lSubs_fromKegg_l'] = dfAllDBs_explode.lSubs_fromKegg.str.len()
dfAllDBs_explode['lProds_fromKegg_l'] = dfAllDBs_explode.lProds_fromKegg.str.len()

dfAllDBs_explode['lReactants_name_id_flat_l'] = dfAllDBs_explode.lReactants_name_id_flat.str.len()
dfAllDBs_explode['lProducts_name_id_flat_l'] = dfAllDBs_explode.lProducts_name_id_flat.str.len()

# dfAllDBs_explode_sort_left_metacyc = dfAllDBs_explode.sort_values(by=['left_metacyc_l'], ascending = False)
# dfAllDBs_explode_sort_right_metacyc = dfAllDBs_explode.sort_values(by=['right_metacyc_l'], ascending = False)
# dfAllDBs_explode_sort_lSubs_fromKegg = dfAllDBs_explode.sort_values(by=['lSubs_fromKegg_l'], ascending = False)
# dfAllDBs_explode_sort_lProds_fromKegg = dfAllDBs_explode.sort_values(by=['lProds_fromKegg_l'], ascending = False)
# dfAllDBs_explode_sort_lReactants_name_id_flat = dfAllDBs_explode.sort_values(by=['lReactants_name_id_flat_l'], ascending = False)
# dfAllDBs_explode_sort_lProducts_name_id_flat = dfAllDBs_explode.sort_values(by=['lProducts_name_id_flat_l'], ascending = False)

# print(dfAllDBs_explode.loc[[0]])


dgb_left_metacyc = dfAllDBs_explode.groupby(by=['left_metacyc_l']).indices
# print(dgb_left_metacyc[11].tolist(), '\n')
# print(dfAllDBs_explode.loc[dgb_left_metacyc[11].tolist()]['MetaCycId'], '\n')

dgb_right_metacyc = dfAllDBs_explode.groupby(by=['right_metacyc_l']).indices
# print(dgb_right_metacyc[11].tolist(), '\n')
# print(dfAllDBs_explode.loc[dgb_right_metacyc[11].tolist()]['MetaCycId'], '\n')

dgb_lSubs_fromKegg = dfAllDBs_explode.groupby(by=['lSubs_fromKegg_l']).indices
# print(dgb_lSubs_fromKegg[11].tolist(), '\n')
# print(dfAllDBs_explode.loc[dgb_lSubs_fromKegg[11].tolist()]['MetaCycId'], '\n')

dgb_lProds_fromKegg = dfAllDBs_explode.groupby(by=['lProds_fromKegg_l']).indices
# print(dgb_lProds_fromKegg[11].tolist(), '\n')
# print(dfAllDBs_explode.loc[dgb_lProds_fromKegg[11].tolist()]['MetaCycId'], '\n')

dgb_lReactants_name_id_flat = dfAllDBs_explode.groupby(by=['lReactants_name_id_flat_l']).indices
# print(dgb_lReactants_name_id_flat[11].tolist(), '\n')
# print(dfAllDBs_explode.loc[dgb_lReactants_name_id_flat[11].tolist()]['MetaCycId'], '\n')

dgb_lProducts_name_id_flat = dfAllDBs_explode.groupby(by=['lProducts_name_id_flat_l']).indices
# print(dgb_lProducts_name_id_flat[11].tolist(), '\n')
# print(dfAllDBs_explode.loc[dgb_lProducts_name_id_flat[11].tolist()]['MetaCycId'], '\n')


## creo df per le righe che hanno in comune S e P
dfAllDBs_equalSP_M = dfAllDBs_explode.copy()
dfAllDBs_equalSP_M = dfAllDBs_equalSP_M[dfAllDBs_equalSP_M['reaction_type'] != '|TRANSPORT|']
dfAllDBs_equalSP_M['equalSP'] = dfAllDBs_equalSP_M.apply(lambda row: checkEqualSubsProds(row['left'], row['right']), axis=1)
dfAllDBs_equalSP_filter_M = dfAllDBs_equalSP_M[dfAllDBs_equalSP_M['equalSP'] == True]
dfAllDBs_equalSP_filter_M['equalSP_mets'] = dfAllDBs_equalSP_filter_M.apply(lambda row: detectEqualSubsProds(row['left'], row['right']), axis=1)
dfAllDBs_equalSP_filter_M['left_woEqualMet'] = dfAllDBs_equalSP_filter_M.apply(lambda row: removeEqualCompound(row['equalSP_mets'], row['left']), axis=1)
dfAllDBs_equalSP_filter_M['right_woEqualMet'] = dfAllDBs_equalSP_filter_M.apply(lambda row: removeEqualCompound(row['equalSP_mets'], row['right']), axis=1)
print('dfAllDBs_equalSP_filter_M\t', dfAllDBs_equalSP_filter_M.shape)

dfAllDBs_equalSP_K = dfAllDBs_explode.copy()
dfAllDBs_equalSP_K = dfAllDBs_equalSP_K[dfAllDBs_equalSP_K['reaction_type'] != '|TRANSPORT|']
dfAllDBs_equalSP_K['equalSP'] = dfAllDBs_equalSP_K.apply(lambda row: checkEqualSubsProds(row['lSubs_fromKegg'], row['lProds_fromKegg']), axis=1)
dfAllDBs_equalSP_filter_K = dfAllDBs_equalSP_K[dfAllDBs_equalSP_K['equalSP'] == True]
dfAllDBs_equalSP_filter_K['equalSP_mets'] = dfAllDBs_equalSP_filter_K.apply(lambda row: detectEqualSubsProds(row['lSubs_fromKegg'], row['lProds_fromKegg']), axis=1)

dfAllDBs_equalSP_filter_K['left_woEqualMet'] = dfAllDBs_equalSP_filter_K.apply(lambda row: removeEqualCompound(row['equalSP_mets'], row['lSubs_fromKegg']), axis=1)
dfAllDBs_equalSP_filter_K['right_woEqualMet'] = dfAllDBs_equalSP_filter_K.apply(lambda row: removeEqualCompound(row['equalSP_mets'], row['lProds_fromKegg']), axis=1)

# for row in dfAllDBs_equalSP_filter_K.itertuples():
#     print(row.KeggId_fromKegg, '\t', row.equalSP_mets,'\t', row.lSubs_fromKegg,'\t', row.lProds_fromKegg, row.left_woEqualMet,'\n')
#
# sys.exit()
# print('dfAllDBs_equalSP_filter_K\t', dfAllDBs_equalSP_filter_K['left_woEqualMet'].tolist())


dfAllDBs_equalSP_R = dfAllDBs_explode.copy()
dfAllDBs_equalSP_R = dfAllDBs_equalSP_R[dfAllDBs_equalSP_R['reaction_type'] != '|TRANSPORT|']
dfAllDBs_equalSP_R['equalSP'] = dfAllDBs_equalSP_R.apply(lambda row: checkEqualSubsProds(row['lReactants_name_id_flat'], row['lProducts_name_id_flat']), axis=1)
dfAllDBs_equalSP_filter_R = dfAllDBs_equalSP_R[dfAllDBs_equalSP_R['equalSP'] == True]
dfAllDBs_equalSP_filter_R['equalSP_mets'] = dfAllDBs_equalSP_filter_R.apply(lambda row: detectEqualSubsProds(row['lReactants_name_id_flat'], row['lProducts_name_id_flat']), axis=1)
dfAllDBs_equalSP_filter_R['left_woEqualMet'] = dfAllDBs_equalSP_filter_R.apply(lambda row: removeEqualCompound(row['equalSP_mets'], row['lReactants_name_id_flat']), axis=1)
dfAllDBs_equalSP_filter_R['right_woEqualMet'] = dfAllDBs_equalSP_filter_R.apply(lambda row: removeEqualCompound(row['equalSP_mets'], row['lProducts_name_id_flat']), axis=1)
print('dfAllDBs_equalSP_filter_R\t', dfAllDBs_equalSP_filter_R.shape)

###############################################################
## Step 2. estrarre da modello info su rxns (eseguire solo la prima volta)
dfRxns = xL.getRxnsInfoGEM(os.path.join(RAWDIR, modelXml + ".xml"))

model = cb.io.read_sbml_model(os.path.join(RAWDIR, modelXml + '.xml'))

# verifico quali reazioni sono trasporti
lNames = []
lTransport = []
lMetsTrasportati = []
for row in dfRxns.itertuples():

    if row.Rxn.startswith('R_'):
        rxnId = row.Rxn[2:]
    else:
        rxnId = row.Rxn
    # print('rxnId\t', rxnId)
    rxn = model.reactions.get_by_id(rxnId)
    lNames.append(rxn.name)
    # lReactants = []
    dReactants = {}
    # print('rxn.reactants\t', rxn.reactants)
    for reactant in rxn.reactants:
        # print(reactant.id, '\t', reactant.compartment)
        if testModel == 'y7' or testModel == 'y8':
            # print('prima\t', reactant.name)
            dReactants[reactant.id] = reactant.name[:reactant.name.rfind('[')].strip()
            # print('Dopo\t', reactant.name[:reactant.name.rfind('[')].strip())
        else:
            dReactants[reactant.id] = reactant.name
            # print('reactant.name\t', reactant.name)
        # if '__91__' in reactant.id and '__93__' in reactant.id:
        #     lReactants.append(reactant.id.split('__91')[0])
    # print('dReactants\t', dReactants)
    lReactants = list(dReactants.values())
    lReactants.sort()
    # print('lReactants\t', lReactants)

    # lProducts = []
    dProducts = {}
    # for product in rxn.products:
    #     if '__91__' in product.id and '__93__' in product.id:
    #         lProducts.append(product.id.split('__91')[0])
    for product in rxn.products:
        # print(product.id, '\t', product.compartment)
        if testModel == 'y7' or testModel == 'y8':
            # print('prima\t', product.name)
            dProducts[product.id] = product.name[:product.name.rfind('[')].strip()
            # print('Dopo\t', product.name[:product.name.rfind('[')].strip())
        else:
            dProducts[product.id] = product.name
            # print('product.name\t', product.name)
    # print('dProducts\t', dProducts)
    lProducts = list(dProducts.values())
    lProducts.sort()
    # print('lProducts\t', lProducts)

    # print(lReactants, '\t', lProducts)
    # print('Intersect:\t', gL.intersect(lReactants, lProducts))
    if len(gL.intersect(lReactants, lProducts)) != 0:
    # if len(lReactants) == len(lProducts) and lReactants == lProducts:
        lTransport.append(True)
        metTrasportati = gL.intersect(lReactants, lProducts)
        lMetsTrasportati.append(metTrasportati)
        # print('lMetsTrasportati\t', metTrasportati)
        # print('Is Transport')
    else:
        lTransport.append(False)
        lMetsTrasportati.append([])
    # print('\n')

dfRxns['Name'] = lNames
dfRxns['IsTransport'] = lTransport
dfRxns['trasportedMets'] = lMetsTrasportati

initialLetter = all(el.startswith('M_') for el in list(dfMetsFromModel['Id']))
# print('initialLetter\t', initialLetter)

# verifico quali reazioni sono exchange
lExchange = []
for row in dfRxns.itertuples():
    if row.Rxn.startswith('R_'):
        rxnId = row.Rxn[2:]
    else:
        rxnId = row.Rxn
    # print('rxnId\t', rxnId)
    rxn = model.reactions.get_by_id(rxnId)
    dReactants = {}
    # lReactants = []
    for reactant in rxn.reactants:
        # print('reactant\t', reactant.id)
        if initialLetter is True:
        # if modelName == 'Recon3':
            outMet = dfMetsFromModel.loc[dfMetsFromModel['Id'] == 'M_' + reactant.id]
        else:
            outMet = dfMetsFromModel.loc[dfMetsFromModel['Id'] == reactant.id]
        # outMet = dfMetsFromModel.loc[dfMetsFromModel['Id_converted2'] == reactant.id]
        outMet = outMet.reset_index(drop=True)
        # print('outMetR\n', outMet)
        # print('Is empty\t', outMet.empty)
        dReactants[reactant.id] = outMet.iloc[0]['boundaryCondition']
        # lReactants.append(reactant.id)
    # print('dReactants\t', dReactants)
    dProducts = {}
    for product in rxn.products:
        if initialLetter is True:
        # if modelName == 'Recon3':
            outMet = dfMetsFromModel.loc[dfMetsFromModel['Id'] == 'M_' + product.id]
        else:
            outMet = dfMetsFromModel.loc[dfMetsFromModel['Id'] == product.id]
        # outMet = dfMetsFromModel.loc[dfMetsFromModel['Id_converted2'] == product.id]
        outMet = outMet.reset_index(drop=True)
        # print('outMetP\n', outMet)
        dProducts[product.id] = outMet.iloc[0]['boundaryCondition']
        # lProducts.append(product.id)
    # print('dProducts\t', dProducts)
    if len(dReactants) == 0 or len(dProducts) == 0:
        lExchange.append(True)
        # print('Is Exchange')
    elif any(v == True for v in list(dReactants.values())) is True or any(v == True for v in list(dReactants.values())) is 'True' or any(v == True for v in list(dReactants.values())) is 'true' or any(v == True for v in list(dProducts.values())) is True or any(v == True for v in list(dProducts.values())) is 'True' or any(v == True for v in list(dProducts.values())) is 'true':
        lExchange.append(True)
        # print('Is Exchange')
    else:
        lExchange.append(False)
    # print('\n')

dfRxns['IsExchange'] = lExchange

# estraggo GPR da modello
lRules = []
for row in dfRxns.itertuples():
    if row.Rxn.startswith('R_'):
        rxnId = row.Rxn[2:]
    else:
        rxnId = row.Rxn
    rxn = model.reactions.get_by_id(rxnId)
    rule = rxn.gene_reaction_rule
    # print('rule\t', rule)
    if '__46__' in rule:
        s = '__46__'
        pos = rule.find(s)
        # print('pos\t', pos)
        while pos != -1:
            rule = rule.replace(s, '.')
            pos = rule.find(s)
            # print('pos\t', pos)
        # print('Rule:\t', rule)
    lRules.append(rule)

dfRxns['GPRrule'] = lRules

dfRxns.to_csv(os.path.join(OUTDIR, dfrxnsInfo + '.csv'), sep = '\t', index = False)

############################


### Step 3. Associare putativo identificativo alla reazione
model = cb.io.read_sbml_model(os.path.join(RAWDIR, modelXml + '.xml'))
dfRxns = pd.read_csv(os.path.join(OUTDIR, dfrxnsInfo + '.csv'), sep = '\t', dtype = {'Rxn': str, 'KeggId': str, 'GPR': str, 'Name': str, 'IsTransport': bool, 'IsExchange': bool, 'GPRrule': str})
dfRxns['EC number'] = dfRxns['EC number'].apply(literal_eval)
dfRxns['trasportedMets'] = dfRxns['trasportedMets'].apply(literal_eval)

if testModel == 'recon3':
    dfRxns['Rxn_converted'] = dfRxns['Rxn'].str.replace('__91__', '[')
    dfRxns['Rxn_converted'] = dfRxns['Rxn_converted'].str.replace('__93__', ']')

lcofactors = ['44215', '|NAD|', 'C00003', '15846', '57540', 'C00001', '|WATER|', '15377', '16908', 'C00004', '|NADH|', '57945', '25805', '|OXYGEN-MOLECULE|', '15379', 'C00007', '|NADP|', '58349', '|B-HEP-1:5|', 'C00667', '18009', 'C00006',
            '16474', '77312', '|NADPH|', '|CPD-16005|', 'C20745', '77177', '57783', 'C00005', 'CPD-16005', '15996', 'C00044', '37565', '|GTP|', 'C00010', '57287', '|CO-A|', '15346', 'CARBON-DIOXIDE', 'C00011', '|CARBON-DIOXIDE|', '16526',
            'AMMONIUM', 'C01342', '|AMMONIUM|', '|AMMONIA|', 'C00014', 'AMMONIA', '16134', '28938', 'C00002', '15422', '22258', '|ATP|', '30616', '456216', '73342', '|ADP|', '22251', 'C00008', 'G11113', '22252', '16761']

# lcofactors = ['44215', '|NAD|', 'C00003', '15846', '57540', '16908', 'C00004', '|NADH|', '57945', '25805', '|OXYGEN-MOLECULE|', '15379', 'C00007',
#             '|NADP|', '58349', '|B-HEP-1:5|', 'C00667', '18009', 'C00006', '16474', '77312', '|NADPH|', '|CPD-16005|', 'C20745', '77177', '57783', 'C00005', 'CPD-16005',
#             'C00010', '57287', '|CO-A|', '15346', 'CARBON-DIOXIDE', 'C00011', '|CARBON-DIOXIDE|', '16526', '22258', '73342', '22252']

# lcofactors = ['44215', '|NAD|', 'C00003', '15846', '57540', '16908', 'C00004', '|NADH|', '57945', '25805', '|OXYGEN-MOLECULE|', '15379', 'C00007',
#             '|NADP|', '58349', '|B-HEP-1:5|', 'C00667', '18009', 'C00006', '16474', '77312', '|NADPH|', '|CPD-16005|', 'C20745', '77177', '57783', 'C00005', 'CPD-16005',
#             '22258', '73342', '22252']


if testModel == 'y7' or testModel == 'y8':
    dfMetsFromModel['Name'] = dfMetsFromModel['Name'].str.strip()

outFile = open(os.path.join(OUTDIR, dfrxnsInfo + '_wIds_' + timeStamp + '.csv'), mode='w')
gL.writeLineByLineToFile(outFile, ['RxnId', 'PutativeIdentifiers'], '\t')

lcofactors = gL.unique(lcofactors)
lIds_all = []
lIdentifiersRxn_all = []
for rowRxn in dfRxns.itertuples():
    isTransport = rowRxn.IsTransport
    isExchange = rowRxn.IsExchange
    if isTransport == True:
        if rowRxn.Rxn.startswith('R_'): ####### Rxn_converted
            rxnId = rowRxn.Rxn[2:] ####### Rxn_converted
        else:
            rxnId = rowRxn.Rxn ####### Rxn_converted

        print('rxnId\t', rxnId)
        rxn = model.reactions.get_by_id(rxnId)
        print('ID\t', rxn.id)
        print('isTransport\t', isTransport)
        print('isExchange\t', isExchange)
        lIds_all.append(rxn.id)
        lIdentifiersRxn = [] ## lista dove salvo tutti gli identificativi che trovo associati alla rxn corrente
        lReactants = []
        if isTransport == True:
            lReactants = rowRxn.trasportedMets
        else:
            for r in rxn.reactants:
                # print('react\t', r.id)
                lReactants.append(r.id)
        lProducts = []
        if isTransport == True:
            lProducts = rowRxn.trasportedMets
        else:
            for r in rxn.products:
                lProducts.append(r.id)
        print('lReactants\t', lReactants)
        print('lProducts\t', lProducts)


        if isTransport == True:
            dfMetsFromModel_reactants = dfMetsFromModel[dfMetsFromModel['Name'].isin(lReactants)] ### Id_converted
        else:
            if testModel == 'recon3' or testModel == 'hmr':
                dfMetsFromModel_reactants = dfMetsFromModel[dfMetsFromModel['Id'].isin(['M_' + m for m in lReactants])] ### Id_converted
            else:
                dfMetsFromModel_reactants = dfMetsFromModel[dfMetsFromModel['Id'].isin(lReactants)] ### Id_converted
        # dfMetsFromModel_reactants = dfMetsFromModel[dfMetsFromModel['Id_converted'].isin(lReactants)] ### Id
        print('dfMetsFromModel_reactants\n', dfMetsFromModel_reactants, '\n')
        # for row in dfMetsFromModel_reactants.itertuples():
        #     print('dfMetsFromModel_reactants\t', row.Id, '\t', row.Name, '\t', row.lIdentifiers, '\n')
        lReactants_ids_original = list(dfMetsFromModel_reactants['lIdentifiers'].dropna())

        if isTransport == False and isExchange == False:
            lReactants_ids_original = [l for l in lReactants_ids_original if 'C00080' not in l] ## tolgo protone perchp non sempre incluso a seconda di come espressa la rxn nei vari db
        print('lReactants_ids_original\t', lReactants_ids_original)
        lReactants_ids = []
        for lIdReactant in lReactants_ids_original:
            # print('prima lIdReactant:\t', len(lIdReactant))
            lCompsFound = []
            lCompsFound += lIdReactant
            # l = [el for el in l if el.startswith('|') == False and el.startswith('C') == False and el.startswith('G') == False]
            chebiIdentifiers = gL.intersect(lIdReactant,list(dfChebiCompounds['Id']))

            dfChebiExplodedFilter = dfChebiCompounds_exploded[dfChebiCompounds_exploded['ParentalChebiIds_exploded'].isin(chebiIdentifiers)]
            dfChebiFilter = dfChebiCompounds[dfChebiCompounds['Id'].isin(chebiIdentifiers)]

            lAllChebi2Search = []
            if dfChebiExplodedFilter.empty is False:
                for riga in list(dfChebiExplodedFilter['ParentalChebiIds']):
                    lAllChebi2Search += riga
                for riga in list(dfChebiExplodedFilter['AllChebiIds']):
                    lAllChebi2Search += riga

            if dfChebiFilter.empty is False:
                for riga in list(dfChebiFilter['ParentalChebiIds']):
                    lAllChebi2Search += riga
                for riga in list(dfChebiFilter['AllChebiIds']):
                    lAllChebi2Search += riga
            lAllChebi2Search = gL.unique(lAllChebi2Search)
            # print('lAllChebi2Search\t', lAllChebi2Search, '\n')
            dfChebiFormula_filtered = dfChebiFormula[dfChebiFormula['COMPOUND_ID'].isin(lAllChebi2Search)]
            if dfChebiFormula_filtered.empty is False:
                dfChebiFormula_filtered = dfChebiFormula_filtered.reset_index(drop = True)
                dfChebiFormula_filtered_toKeep = dfChebiFormula_filtered[dfChebiFormula_filtered['TYPE'] == 'FORMULA']
                if dfChebiFormula_filtered_toKeep.empty is False:
                    # print('Con formula\t', list(dfChebiFormula_filtered_toKeep['COMPOUND_ID']), '\n')
                    lCompsFound += list(dfChebiFormula_filtered_toKeep['COMPOUND_ID'])

                    for riga in list(dfChebiFormula_filtered_toKeep['COMPOUND_ID']):
                        dfChebiExplodedFilter2 = dfChebiCompounds_exploded[dfChebiCompounds_exploded['ParentalChebiIds_exploded'].isin(list(dfChebiFormula_filtered_toKeep['COMPOUND_ID']))]
                        for riga2 in list(dfChebiExplodedFilter2['ParentalChebiIds']):
                            # print('len riga2\t', len(riga2))
                            lCompsFound += riga2
                            # print('Eventuali alternativi \t', riga2, '\n')

            lCompsFound += [element for el in list(dfChebiFilter['ParentalKeggC']) for element in el] + [element for el in list(dfChebiFilter['ParentalKeggG']) for element in el] + [element for el in list(dfChebiFilter['ParentalMetacyc']) for element in el] + [element for el in list(dfChebiFilter['AllKeggC']) for element in el] + [element for el in list(dfChebiFilter['AllKeggG']) for element in el] + [element for el in list(dfChebiFilter['AllMetacyc']) for element in el]
            # print('lCompsFound LEN\t', len(lCompsFound))
            lCompsFound = gL.unique(lCompsFound)
            # print('lCompsFound FINE\t', len(lCompsFound))
            lReactants_ids.append(lCompsFound)
        print('lReactants_ids\n', len(lReactants_ids), '\n')

        if isTransport == True:
            dfMetsFromModel_products = dfMetsFromModel[dfMetsFromModel['Name'].isin(lProducts)] ### Id_converted
        else:
            if testModel == 'recon3' or testModel == 'hmr':
                dfMetsFromModel_products = dfMetsFromModel[dfMetsFromModel['Id'].isin(['M_' + m for m in lProducts])]### Id_converted
            else:
                dfMetsFromModel_products = dfMetsFromModel[dfMetsFromModel['Id'].isin(lProducts)]### Id_converted
        # dfMetsFromModel_products = dfMetsFromModel[dfMetsFromModel['Id_converted'].isin(lProducts)]### Id_converted
        # print('dfMetsFromModel_products\n', dfMetsFromModel_products, '\n')
        # for row in dfMetsFromModel_products.itertuples():
        #     print('dfMetsFromModel_products\t', row.Id, '\t', row.Name, '\t', row.lIdentifiers, '\n')
        lProducts_ids_original = list(dfMetsFromModel_products['lIdentifiers'].dropna())
        if isTransport == False and isExchange == False:
            lProducts_ids_original = [l for l in lProducts_ids_original if 'C00080' not in l] ## tolgo protone perchp non sempre incluso a seconda di come espressa la rxn nei vari db

        print('lProducts_ids_original\t', lProducts_ids_original)
        lProducts_ids = []
        for lIdProduct in lProducts_ids_original:
            lCompsFound = []
            lCompsFound += lIdProduct
            # l = [el for el in l if el.startswith('|') == False and el.startswith('C') == False and el.startswith('G') == False]
            chebiIdentifiers = gL.intersect(lIdProduct,list(dfChebiCompounds['Id']))
            dfChebiFilter = dfChebiCompounds[dfChebiCompounds['Id'].isin(chebiIdentifiers)]
            dfChebiExplodedFilter = dfChebiCompounds_exploded[dfChebiCompounds_exploded['ParentalChebiIds_exploded'].isin(chebiIdentifiers)]

            lAllChebi2Search = []
            if dfChebiExplodedFilter.empty is False:
                for riga in list(dfChebiExplodedFilter['ParentalChebiIds']):
                    lAllChebi2Search += riga
                for riga in list(dfChebiExplodedFilter['AllChebiIds']):
                    lAllChebi2Search += riga

            if dfChebiFilter.empty is False:
                for riga in list(dfChebiFilter['ParentalChebiIds']):
                    lAllChebi2Search += riga
                for riga in list(dfChebiFilter['AllChebiIds']):
                    lAllChebi2Search += riga
            lAllChebi2Search = gL.unique(lAllChebi2Search)
            dfChebiFormula_filtered = dfChebiFormula[dfChebiFormula['COMPOUND_ID'].isin(lAllChebi2Search)]
            if dfChebiFormula_filtered.empty is False:
                dfChebiFormula_filtered = dfChebiFormula_filtered.reset_index(drop = True)
                dfChebiFormula_filtered_toKeep = dfChebiFormula_filtered[dfChebiFormula_filtered['TYPE'] == 'FORMULA']
                if dfChebiFormula_filtered_toKeep.empty is False:
                    # print('Con formula\t', list(dfChebiFormula_filtered_toKeep['COMPOUND_ID']), '\n')
                    lCompsFound += list(dfChebiFormula_filtered_toKeep['COMPOUND_ID'])
                    for riga in list(dfChebiFormula_filtered_toKeep['COMPOUND_ID']):
                        dfChebiExplodedFilter2 = dfChebiCompounds_exploded[dfChebiCompounds_exploded['ParentalChebiIds_exploded'].isin(list(dfChebiFormula_filtered_toKeep['COMPOUND_ID']))]
                        for riga2 in list(dfChebiExplodedFilter2['ParentalChebiIds']):
                            lCompsFound += riga2
                            # print('Eventuali alternativi \t', riga2, '\n')

            lCompsFound += [element for el in list(dfChebiFilter['ParentalKeggC']) for element in el] + [element for el in list(dfChebiFilter['ParentalKeggG']) for element in el] + [element for el in list(dfChebiFilter['ParentalMetacyc']) for element in el] + [element for el in list(dfChebiFilter['AllKeggC']) for element in el] + [element for el in list(dfChebiFilter['AllKeggG']) for element in el] + [element for el in list(dfChebiFilter['AllMetacyc']) for element in el]
            lCompsFound = gL.unique(lCompsFound)
            lProducts_ids.append(lCompsFound)
        print('lProducts_ids\n', len(lProducts_ids), '\n')

        dfR = pd.DataFrame(lReactants_ids).transpose()
        dfP = pd.DataFrame(lProducts_ids).transpose()

        if all(r != [] for r in lReactants_ids) == True and all(p != [] for p in lProducts_ids) == True:

            lReactants_ids_cofactors = []
            lReactants_ids_internal = []
            for lReagenti in lReactants_ids:
                if len(gL.intersect(lcofactors, lReagenti)) != 0:
                    lReactants_ids_cofactors += lReagenti
                else:
                    lReactants_ids_internal += lReagenti

            lReactants_ids_cofactors = gL.unique(lReactants_ids_cofactors)
            lReactants_ids_internal = gL.unique(lReactants_ids_internal)
            print('lReactants_ids_cofactors\t', len(lReactants_ids_cofactors))
            print('lReactants_ids_internal\t', len(lReactants_ids_internal))

            lProducts_ids_cofactors = []
            lProducts_ids_internal = []
            for lProdotti in lProducts_ids:
                if len(gL.intersect(lcofactors, lProdotti)) != 0:
                    lProducts_ids_cofactors += lProdotti
                else:
                    lProducts_ids_internal += lProdotti
            lProducts_ids_cofactors = gL.unique(lProducts_ids_cofactors)
            lProducts_ids_internal = gL.unique(lProducts_ids_internal)
            print('lProducts_ids_cofactors\t', len(lProducts_ids_cofactors))
            print('lProducts_ids_internal\t', len(lProducts_ids_internal))
            dfAllDBs_copy = dfAllDBs_explode.copy()
            if isTransport == True:
                dfAllDBs_copy = dfAllDBs_copy[dfAllDBs_copy['reaction_type'] == '|TRANSPORT|']
            else:
                dfAllDBs_copy = dfAllDBs_copy[dfAllDBs_copy['reaction_type'] != '|TRANSPORT|']

            # l2Search = []
            # for jj in lReactants_ids:
            #     l2Search += gL.unique(jj)
            # for jj in lProducts_ids:
            #     l2Search+= gL.unique(jj)
            # l2Search = gL.unique(l2Search)
            # print('l2Search\n', l2Search, '\n')

            # print('START\t', dfAllDBs_copy.shape)
            lIdxs = []
            if len(lReactants_ids) in dgb_left_metacyc and len(lProducts_ids) in dgb_right_metacyc:
                lCommonIdxs = gL.intersect(dgb_left_metacyc[len(lReactants_ids)].tolist(), dgb_right_metacyc[len(lProducts_ids)].tolist())
                lIdxs += lCommonIdxs

            if len(lReactants_ids) in dgb_right_metacyc and len(lProducts_ids) in dgb_left_metacyc:
                lCommonIdxs = gL.intersect(dgb_right_metacyc[len(lReactants_ids)].tolist(), dgb_left_metacyc[len(lProducts_ids)].tolist())
                lIdxs += lCommonIdxs

            if len(lReactants_ids) in dgb_lSubs_fromKegg and len(lProducts_ids) in dgb_lProds_fromKegg:
                lCommonIdxs = gL.intersect(dgb_lSubs_fromKegg[len(lReactants_ids)].tolist(), dgb_lProds_fromKegg[len(lProducts_ids)].tolist())
                lIdxs += lCommonIdxs

            if len(lReactants_ids) in dgb_lProds_fromKegg and len(lProducts_ids) in dgb_lSubs_fromKegg:
                lCommonIdxs = gL.intersect(dgb_lProds_fromKegg[len(lReactants_ids)].tolist(), dgb_lSubs_fromKegg[len(lProducts_ids)].tolist())
                lIdxs += lCommonIdxs

            if len(lReactants_ids) in dgb_lReactants_name_id_flat and len(lProducts_ids) in dgb_lProducts_name_id_flat:
                lCommonIdxs = gL.intersect(dgb_lReactants_name_id_flat[len(lReactants_ids)].tolist(), dgb_lProducts_name_id_flat[len(lProducts_ids)].tolist())
                lIdxs += lCommonIdxs

            if len(lReactants_ids) in dgb_lProducts_name_id_flat and len(lProducts_ids) in dgb_lReactants_name_id_flat:
                lCommonIdxs = gL.intersect(dgb_lProducts_name_id_flat[len(lReactants_ids)].tolist(), dgb_lReactants_name_id_flat[len(lProducts_ids)].tolist())
                lIdxs += lCommonIdxs

            lIdxs = gL.unique(lIdxs)
            # print('lIdxs\t', len(lIdxs))

            lDfIdx = dfAllDBs_copy.index.tolist()
            # print('lDfIdx\t', len(gL.intersect(lIdxs, lDfIdx)))

            dfAllDBs_copy_filter = dfAllDBs_copy.loc[gL.intersect(lIdxs, lDfIdx)]
            # print('IS IN\t', 'R00105' in list(dfAllDBs_copy_filter['KeggId_fromKegg']))
            # sys.exit()

            ## Metacyc
            print('\nMetacyc')
            lDfs2Concat = filterRows(dfAllDBs_copy_filter, 'left_metacyc', 'right_metacyc', lReactants_ids_cofactors, lReactants_ids_internal, lProducts_ids_cofactors, lProducts_ids_internal)
            if len(lDfs2Concat) != 0:
                lPutativeRxns = findRxnsAfterFilter(lDfs2Concat, 'left_metacyc', 'right_metacyc', 'metacyc')
                lIdentifiersRxn += lPutativeRxns

            dfAllDBs_equalSP_filter_M_copy = dfAllDBs_equalSP_filter_M.copy()
            lPutativeRxns = findPutativeRxns(dfAllDBs_equalSP_filter_M_copy, 'left_woEqualMet', 'right', 'metacyc')
            lIdentifiersRxn += lPutativeRxns
            dfAllDBs_equalSP_filter_M_copy = dfAllDBs_equalSP_filter_M.copy()
            lPutativeRxns = findPutativeRxns(dfAllDBs_equalSP_filter_M_copy, 'left', 'right_woEqualMet', 'metacyc')
            lIdentifiersRxn += lPutativeRxns

            dfAllDBs_equalSP_filter_K_copy = dfAllDBs_equalSP_filter_K.copy()
            lPutativeRxns = findPutativeRxns(dfAllDBs_equalSP_filter_K_copy, 'left_woEqualMet', 'lProds_fromKegg', 'metacyc')
            lIdentifiersRxn += lPutativeRxns
            dfAllDBs_equalSP_filter_K_copy = dfAllDBs_equalSP_filter_K.copy()
            lPutativeRxns = findPutativeRxns(dfAllDBs_equalSP_filter_K_copy, 'lSubs_fromKegg', 'right_woEqualMet', 'metacyc')
            lIdentifiersRxn += lPutativeRxns

            dfAllDBs_equalSP_filter_R_copy = dfAllDBs_equalSP_filter_R.copy()
            lPutativeRxns = findPutativeRxns(dfAllDBs_equalSP_filter_R_copy, 'left_woEqualMet', 'lProducts_name_id_flat', 'metacyc')
            lIdentifiersRxn += lPutativeRxns
            dfAllDBs_equalSP_filter_R_copy = dfAllDBs_equalSP_filter_R.copy()
            lPutativeRxns = findPutativeRxns(dfAllDBs_equalSP_filter_R_copy, 'lReactants_name_id_flat', 'right_woEqualMet', 'metacyc')
            lIdentifiersRxn += lPutativeRxns

            ## kegg
            print('\nkegg')
            print('IS iN R11680\t', 'R11680' in list(dfAllDBs_copy_filter['KeggId_fromKegg']))
            lDfs2Concat = filterRows(dfAllDBs_copy_filter, 'lSubs_fromKegg', 'lProds_fromKegg', lReactants_ids_cofactors, lReactants_ids_internal, lProducts_ids_cofactors, lProducts_ids_internal)
            if len(lDfs2Concat) != 0:
                lPutativeRxns = findRxnsAfterFilter(lDfs2Concat, 'lSubs_fromKegg', 'lProds_fromKegg', 'kegg')
                lIdentifiersRxn += lPutativeRxns

            dfAllDBs_equalSP_filter_M_copy = dfAllDBs_equalSP_filter_M.copy()
            lPutativeRxns = findPutativeRxns(dfAllDBs_equalSP_filter_M_copy, 'left_woEqualMet', 'right', 'kegg')
            lIdentifiersRxn += lPutativeRxns
            dfAllDBs_equalSP_filter_M_copy = dfAllDBs_equalSP_filter_M.copy()
            lPutativeRxns = findPutativeRxns(dfAllDBs_equalSP_filter_M_copy, 'left', 'right_woEqualMet', 'kegg')
            lIdentifiersRxn += lPutativeRxns

            dfAllDBs_equalSP_filter_K_copy = dfAllDBs_equalSP_filter_K.copy()
            lPutativeRxns = findPutativeRxns(dfAllDBs_equalSP_filter_K_copy, 'left_woEqualMet', 'lProds_fromKegg', 'kegg')
            lIdentifiersRxn += lPutativeRxns
            dfAllDBs_equalSP_filter_K_copy = dfAllDBs_equalSP_filter_K.copy()
            lPutativeRxns = findPutativeRxns(dfAllDBs_equalSP_filter_K_copy, 'lSubs_fromKegg', 'right_woEqualMet', 'kegg')
            lIdentifiersRxn += lPutativeRxns

            dfAllDBs_equalSP_filter_R_copy = dfAllDBs_equalSP_filter_R.copy()
            lPutativeRxns = findPutativeRxns(dfAllDBs_equalSP_filter_R_copy, 'left_woEqualMet', 'lProducts_name_id_flat', 'kegg')
            lIdentifiersRxn += lPutativeRxns
            dfAllDBs_equalSP_filter_R_copy = dfAllDBs_equalSP_filter_R.copy()
            lPutativeRxns = findPutativeRxns(dfAllDBs_equalSP_filter_R_copy, 'lReactants_name_id_flat', 'right_woEqualMet', 'kegg')
            lIdentifiersRxn += lPutativeRxns

            ## rhea
            print('\nrhea')
            lDfs2Concat = filterRows(dfAllDBs_copy_filter, 'lReactants_name_id_flat', 'lProducts_name_id_flat', lReactants_ids_cofactors, lReactants_ids_internal, lProducts_ids_cofactors, lProducts_ids_internal)

            if len(lDfs2Concat) != 0:
                lPutativeRxns = findRxnsAfterFilter(lDfs2Concat, 'lReactants_name_id_flat', 'lProducts_name_id_flat', 'rhea')
                lIdentifiersRxn += lPutativeRxns

            dfAllDBs_equalSP_filter_M_copy = dfAllDBs_equalSP_filter_M.copy()
            lPutativeRxns = findPutativeRxns(dfAllDBs_equalSP_filter_M_copy, 'left_woEqualMet', 'right', 'rhea')
            lIdentifiersRxn += lPutativeRxns
            dfAllDBs_equalSP_filter_M_copy = dfAllDBs_equalSP_filter_M.copy()
            lPutativeRxns = findPutativeRxns(dfAllDBs_equalSP_filter_M_copy, 'left', 'right_woEqualMet', 'rhea')
            lIdentifiersRxn += lPutativeRxns

            dfAllDBs_equalSP_filter_K_copy = dfAllDBs_equalSP_filter_K.copy()
            lPutativeRxns = findPutativeRxns(dfAllDBs_equalSP_filter_K_copy, 'left_woEqualMet', 'lProds_fromKegg', 'rhea')
            lIdentifiersRxn += lPutativeRxns
            dfAllDBs_equalSP_filter_K_copy = dfAllDBs_equalSP_filter_K.copy()
            lPutativeRxns = findPutativeRxns(dfAllDBs_equalSP_filter_K_copy, 'lSubs_fromKegg', 'right_woEqualMet', 'rhea')
            lIdentifiersRxn += lPutativeRxns

            dfAllDBs_equalSP_filter_R_copy = dfAllDBs_equalSP_filter_R.copy()
            lPutativeRxns = findPutativeRxns(dfAllDBs_equalSP_filter_R_copy, 'left_woEqualMet', 'lProducts_name_id_flat', 'rhea')
            lIdentifiersRxn += lPutativeRxns
            dfAllDBs_equalSP_filter_R_copy = dfAllDBs_equalSP_filter_R.copy()
            lPutativeRxns = findPutativeRxns(dfAllDBs_equalSP_filter_R_copy, 'lReactants_name_id_flat', 'right_woEqualMet', 'rhea')
            lIdentifiersRxn += lPutativeRxns


        lIdentifiersRxn = gL.unique(lIdentifiersRxn)
        print('lIdentifiersRxn\n', lIdentifiersRxn, '\n\n')

        ## non posso rimuovere reazioni gia' identificate per via di reazioni duplicate in multipli compartimenti
        # df1 = dfAllDBs_explode[dfAllDBs_explode['MetaCycId'].isin(lIdentifiersRxn)]
        # df2 = dfAllDBs_explode[dfAllDBs_explode['KeggId_fromKegg'].isin(lIdentifiersRxn)]
        # df3 = explodeFilterAndDrop(dfAllDBs_explode, 'RheaId_master', lIdentifiersRxn)
        #
        # lMatchIdxs = df1.index.tolist() + df2.index.tolist() + df3.index.tolist()
        # lMatchIdxs = gL.unique(lMatchIdxs)
        # dfAllDBs_explode = dfAllDBs_explode.drop(dfAllDBs_explode.index[lMatchIdxs])
        lIdentifiersRxn_all.append(lIdentifiersRxn)
        gL.writeLineByLineToFile(outFile, [rxn.id, lIdentifiersRxn], '\t')

outFile.close()


dfRxns['RxnId'] = lIds_all


dfRxnWId = pd.read_csv(os.path.join(OUTDIR, dfrxnsInfo + '_wIds_' + timeStamp + '.csv'), sep = '\t', dtype=str)
print('dfRxnWId\n', dfRxnWId.head())

dfMerge = pd.merge(dfRxns, dfRxnWId, on = 'RxnId')
dfMerge.to_csv(os.path.join(OUTDIR, dfrxnsInfo + '_enriched_' + timeStamp + '.csv'), sep = '\t', index = False)


# dfRxns['PutativeIdentifiers'] = lIdentifiersRxn_all
# dfRxns.to_csv(os.path.join(OUTDIR, dfrxnsInfo + '_enriched_' + timeStamp + '.csv'), sep = '\t', index = False)

end = time.time()
print('elapsed time\t', end-start, '\tseconds')
