# -*- coding: utf-8 -*-
import pandas as pd
import genericLib as gL
import re
from lxml import etree as ET

def explodeFilterAndDrop(df, oldCol, lMets):
    df[oldCol + '_exploded'] = df[oldCol].values
    df = df.explode(oldCol + '_exploded')
    df = df[df[oldCol + '_exploded'].notna()]
    newDf = df[df[oldCol + '_exploded'].isin(lMets)]
    newDf[oldCol + '_t'] = newDf[oldCol].apply(tuple)
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


def findRxnsAfterFilter(lDfs2Concat, colNameS, colNameP, db, dfR, dfP):
    lPutativeRxns = []
    dfAllDBs_copy_All = pd.concat(lDfs2Concat)
    dfAllDBs_copy_All[colNameS + '_t'] = dfAllDBs_copy_All[colNameS].apply(tuple)
    dfAllDBs_copy_All[colNameP + '_t'] = dfAllDBs_copy_All[colNameP].apply(tuple)

    lDfs = filterwExtractRows(dfAllDBs_copy_All, dfR, dfP, colNameS, colNameP)
    if len(lDfs) != 0:
        df = pd.concat(lDfs)
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
        lMets.append(KEGGCompundMatch(k))
    if rxnType != '|TRANSPORT|' and 'C00080' in lMets:
        lMets.remove('C00080')
    return lMets

def KEGGCompundMatch(strKEGG):
    """
    Split the KEGG string of the equation into a strKEGGcompound and the other
    """
    if 'C' in strKEGG:
        tokens = re.search(r'(C\d{5})', strKEGG)
        match = tokens.group(1)
    else:
        tokens = re.search(r'(G\d{5})', strKEGG)
        match = tokens.group(1)
    return match


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
    for item in lReactants_ids:
        lSideSx.append(len(gL.intersect(item, list(lLeft))) != 0)

    lSideDx = []
    for item in lProducts_ids:
        lSideDx.append(len(gL.intersect(item, list(lRight))) != 0)

    if all(lSideSx) == True and all(lSideDx) == True:
        return True
    else:
        return False

def findReaction_fromPutative2Rxn(lReactants_ids, lProducts_ids, lLeft, lRight):
    ## The first vs. the third and the second vs. the fourth arguments of the function are compared
    lSideSx = []
    for item in list(lLeft):
        lSideSx_internal = []
        for lComps in lReactants_ids:
            lSideSx_internal.append(item in lComps)
        lSideSx.append(any(lSideSx_internal))

    lSideDx = []
    for item in list(lRight):
        lSideDx_internal = []
        for lComps in lProducts_ids:
            lSideDx_internal.append(item in lComps)
        lSideDx.append(any(lSideDx_internal))

    if all(lSideSx) == True and all(lSideDx) == True:
        return True
    else:
        return False

def applyFindReaction_rhea(lS, lP, lReactants_ids, lProducts_ids):
    idx = 0
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
    ## The first vs. the third and the second vs. the fourth arguments of the function are compared
    lSideSx = []
    for item in lReactants_ids:
        lSideSx_internal = []
        for possibleComp in lLeft:
            lSideSx_internal.append(len(gL.intersect(item, possibleComp)) != 0)
        lSideSx.append(any(lSideSx_internal))

    lSideDx = []
    for item in lProducts_ids:
        lSideDx_internal = []
        for possibleComp in lRight:
            lSideDx_internal.append(len(gL.intersect(item, possibleComp)) != 0)
        lSideDx.append(any(lSideDx_internal))

    if all(lSideSx) == True and all(lSideDx) == True:
        return True
    else:
        return False

def checkEqualSubsProds(lLeft, lRight):
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
    return gL.difference(lMets, l2Remove)

def extractRows(df, column2Search):
    out = df.isin(column2Search)
    return out.any().all()

def getRxnsInfoGEM(model):
    rId=[]
    rKegg=[]
    rGpr = []
    rEc = []
    tree = ET.parse(model)
    root = tree.getroot()
    ## Definition of useful tags
    listRxnsTag = ''
    annotTag = ''
    liTag = ''
    pTag = ''
    for element in root.getiterator():
        if element.tag.endswith('listOfReactions'):
            listRxnsTag = element.tag
        elif element.tag.endswith('annotation'):
            annotTag = element.tag
        elif element.tag.endswith('li'):
            liTag = element.tag
        elif element.tag.endswith('p'):
            pTag = element.tag

    for rxns in root.getiterator(listRxnsTag):
        for rxn in rxns:
            rId.append(rxn.attrib['id'])
            keggId = ''
            ecNumber = []
            if annotTag != '':
                for annotation in rxn.getiterator(annotTag):
                    if annotation.tag != '' and liTag != '':
                        for element in annotation.getiterator(liTag):
                            for k,v in element.attrib.items():
                                if 'kegg.reaction' in v:
                                    keggId = v[-6:]
                                if 'ec-code' in v:
                                    ecNumber.append(v.split('/')[-1])
            rKegg.append(keggId)
            rEc.append(ecNumber)
            rule = ''
            if pTag != '':
                for gpr in rxn.getiterator(pTag):
                    if "GENE_ASSOCIATION" in gpr.text:
                        rule = gpr.text.split(':')[1].strip()
            rGpr.append(rule)
    dfRxns = pd.DataFrame({'Rxn': rId, 'KeggId': rKegg, 'EC number': rEc, 'GPR':rGpr})
    return dfRxns

def extractChebiIds(l):
    lChebi = []
    for el in l:
        lChebi.append(el.split(';')[0][6:])
    lChebi = gL.unique(lChebi)
    return lChebi


def findTC_transports(lR, r):
    lFindR = [len(gL.intersect(l,r)) != 0 for l in lR]

    if any(lFindR) == True:
        return True
    else:
        return False

def checkNadNadpDependencies_and(nadp, nad, id2Search, id2Search_complete, lAnd, dGenesFromKegg):
    if nadp == True:
        if id2Search in dGenesFromKegg:
            dGeneK = dGenesFromKegg[id2Search]
        else:
            dGeneK = getKeggInfo(id2Search_complete)
            dGenesFromKegg[id2Search] = dGeneK
        if dGeneK != 400 and dGeneK != 404 and 'DEFINITION' in dGeneK:
            if '(NAD(+))' not in dGeneK['DEFINITION'] and '(NAD+)' not in dGeneK['DEFINITION'] and '(NADP(+))' not in dGeneK['DEFINITION'] and '(NADP+)' not in dGeneK['DEFINITION']:
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            elif ('(NAD(+))' in dGeneK['DEFINITION'] or '(NAD+)' in dGeneK['DEFINITION']) and ('(NADP(+))' in dGeneK['DEFINITION'] or '(NADP+)' in dGeneK['DEFINITION']):
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            elif ('(NAD(+))' not in dGeneK['DEFINITION'] and '(NAD+)' not in dGeneK['DEFINITION']) and ('(NADP(+))' in dGeneK['DEFINITION'] or '(NADP+)' in dGeneK['DEFINITION']):
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            else:
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
        if dGeneK != 400 and dGeneK != 404 and 'ORTHOLOGY' in dGeneK:
            if '(NAD(+))' not in dGeneK['ORTHOLOGY'] and '(NAD+)' not in dGeneK['ORTHOLOGY'] and '(NADP(+))' not in dGeneK['ORTHOLOGY'] and '(NADP+)' not in dGeneK['ORTHOLOGY']:
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            elif ('(NAD(+))' in dGeneK['ORTHOLOGY'] or '(NAD+)' in dGeneK['ORTHOLOGY']) and ('(NADP(+))' in dGeneK['ORTHOLOGY'] or '(NADP+)' in dGeneK['ORTHOLOGY']):
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            elif ('(NAD(+))' not in dGeneK['ORTHOLOGY'] and '(NAD+)' not in dGeneK['ORTHOLOGY']) and ('(NADP(+))' in dGeneK['ORTHOLOGY'] or '(NADP+)' in dGeneK['ORTHOLOGY']):
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            else:
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
    if nad == True:
        if id2Search in dGenesFromKegg:
            dGeneK = dGenesFromKegg[id2Search]
        else:
            dGeneK = getKeggInfo(id2Search_complete)
            dGenesFromKegg[id2Search] = dGeneK
        if dGeneK != 400 and dGeneK != 404 and 'DEFINITION' in dGeneK:
            if '(NAD(+))' not in dGeneK['DEFINITION'] and '(NAD+)' not in dGeneK['DEFINITION'] and '(NADP(+))' not in dGeneK['DEFINITION'] and '(NADP+)' not in dGeneK['DEFINITION']:
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            elif ('(NAD(+))' in dGeneK['DEFINITION'] or '(NAD+)' in dGeneK['DEFINITION']) and ('(NADP(+))' in dGeneK['DEFINITION'] or '(NADP+)' in dGeneK['DEFINITION']):
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            elif ('(NAD(+))' in dGeneK['DEFINITION'] or '(NAD+)' in dGeneK['DEFINITION']) and ('(NADP(+))' not in dGeneK['DEFINITION'] and '(NADP+)' not in dGeneK['DEFINITION']):
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            else:
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
        if dGeneK != 400 and dGeneK != 404 and 'ORTHOLOGY' in dGeneK:
            if '(NAD(+))' not in dGeneK['ORTHOLOGY'] and '(NAD+)' not in dGeneK['ORTHOLOGY'] and '(NADP(+))' not in dGeneK['ORTHOLOGY'] and '(NADP+)' not in dGeneK['ORTHOLOGY']:
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            elif ('(NAD(+))' in dGeneK['ORTHOLOGY'] or '(NAD+)' in dGeneK['ORTHOLOGY']) and ('(NADP(+))' in dGeneK['ORTHOLOGY'] or '(NADP+)' in dGeneK['ORTHOLOGY']):
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            elif ('(NAD(+))' in dGeneK['ORTHOLOGY'] or '(NAD+)' in dGeneK['ORTHOLOGY']) and ('(NADP(+))' not in dGeneK['ORTHOLOGY'] and '(NADP+)' not in dGeneK['ORTHOLOGY']):
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            else:
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
    if nadp == False and nad == False:
        if id2Search not in lAnd:
            lAnd.append(id2Search)
    return lAnd, dGenesFromKegg

def checkNadNadpDependencies_or(nadp, nad, id2Search, id2Search_complete, lMetaEnzOR, dGenesFromKegg):
    if nadp == True:
        if id2Search in dGenesFromKegg:
            dGeneK = dGenesFromKegg[id2Search]
        else:
            dGeneK = getKeggInfo(id2Search_complete)
            dGenesFromKegg[id2Search] = dGeneK
        if dGeneK != 400 and dGeneK != 404 and 'DEFINITION' in dGeneK:
            if '(NAD(+))' not in dGeneK['DEFINITION'] and '(NAD+)' not in dGeneK['DEFINITION'] and '(NADP(+))' not in dGeneK['DEFINITION'] and '(NADP+)' not in dGeneK['DEFINITION']:
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            elif ('(NAD(+))' in dGeneK['DEFINITION'] or '(NAD+)' in dGeneK['DEFINITION']) and ('(NADP(+))' in dGeneK['DEFINITION'] or '(NADP+)' in dGeneK['DEFINITION']):
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            elif ('(NAD(+))' not in dGeneK['DEFINITION'] and '(NAD+)' not in dGeneK['DEFINITION']) and ('(NADP(+))' in dGeneK['DEFINITION'] or '(NADP+)' in dGeneK['DEFINITION']):
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            else:
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
        if dGeneK != 400 and dGeneK != 404 and 'ORTHOLOGY' in dGeneK:
            if '(NAD(+))' not in dGeneK['ORTHOLOGY'] and '(NAD+)' not in dGeneK['ORTHOLOGY'] and '(NADP(+))' not in dGeneK['ORTHOLOGY'] and '(NADP+)' not in dGeneK['ORTHOLOGY']:
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            elif ('(NAD(+))' in dGeneK['ORTHOLOGY'] or '(NAD+)' in dGeneK['ORTHOLOGY']) and ('(NADP(+))' in dGeneK['ORTHOLOGY'] or '(NADP+)' in dGeneK['ORTHOLOGY']):
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            elif ('(NAD(+))' not in dGeneK['ORTHOLOGY'] and '(NAD+)' not in dGeneK['ORTHOLOGY']) and ('(NADP(+))' in dGeneK['ORTHOLOGY'] or '(NADP+)' in dGeneK['ORTHOLOGY']):
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            else:
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
    if nad == True:
        if id2Search in dGenesFromKegg:
            dGeneK = dGenesFromKegg[id2Search]
        else:
            dGeneK = getKeggInfo(id2Search_complete)
            dGenesFromKegg[id2Search] = dGeneK
        if dGeneK != 400 and dGeneK != 404 and 'DEFINITION' in dGeneK:
            if '(NAD(+))' not in dGeneK['DEFINITION'] and '(NAD+)' not in dGeneK['DEFINITION'] and '(NADP(+))' not in dGeneK['DEFINITION'] and '(NADP+)' not in dGeneK['DEFINITION']:
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            elif ('(NAD(+))' in dGeneK['DEFINITION'] or '(NAD+)' in dGeneK['DEFINITION']) and ('(NADP(+))' in dGeneK['DEFINITION'] or '(NADP+)' in dGeneK['DEFINITION']):
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            elif ('(NAD(+))' in dGeneK['DEFINITION'] or '(NAD+)' in dGeneK['DEFINITION']) and ('(NADP(+))' not in dGeneK['DEFINITION'] and '(NADP+)' not in dGeneK['DEFINITION']):
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            else:
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
        if dGeneK != 400 and dGeneK != 404 and 'ORTHOLOGY' in dGeneK:
            if '(NAD(+))' not in dGeneK['ORTHOLOGY'] and '(NAD+)' not in dGeneK['ORTHOLOGY'] and '(NADP(+))' not in dGeneK['ORTHOLOGY'] and '(NADP+)' not in dGeneK['ORTHOLOGY']:
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            elif ('(NAD(+))' in dGeneK['ORTHOLOGY'] or '(NAD+)' in dGeneK['ORTHOLOGY']) and ('(NADP(+))' in dGeneK['ORTHOLOGY'] or '(NADP+)' in dGeneK['ORTHOLOGY']):
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            elif ('(NAD(+))' in dGeneK['ORTHOLOGY'] or '(NAD+)' in dGeneK['ORTHOLOGY']) and ('(NADP(+))' not in dGeneK['ORTHOLOGY'] and '(NADP+)' not in dGeneK['ORTHOLOGY']):
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            else:
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
    if nadp == False and nad == False:
        if id2Search not in lMetaEnzOR:
            lMetaEnzOR.append([id2Search])
    return lMetaEnzOR, dGenesFromKegg

def getKeggInfo(rxnString):
    try:
        k = kegg.KEGG()
        info = k.get(rxnString)
        dizInfo = k.parse(info)
    except:
        dizInfo = {}
    return dizInfo
