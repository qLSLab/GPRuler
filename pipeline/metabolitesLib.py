# -*- coding: utf-8 -*-
import re
import pandas as pd
import genericLib as gL
from Bio.KEGG.REST import *
from fuzzywuzzy import process

def queryKeggCompound(db, met, metOriginale, lSymbols):
    request = kegg_find(database = db, query = met)
    lAllResults = request.readlines()
    lkeggId_perfectMatch = []
    dRequest = {}
    for r in lAllResults:
        lsComps = r.strip('\n').split("\t")
        if lsComps != [''] and lsComps[0].startswith('cpd'):
            alternativeNames = lsComps[1].split(';')
            dRequest[lsComps[0][4:]] = alternativeNames
        elif lsComps != [''] and lsComps[0].startswith('gl'):
            alternativeNames = lsComps[1].split(';')
            dRequest[lsComps[0][3:]] = alternativeNames
    lkeggId_perfectMatch += [kId for kId, comps in dRequest.items() for comp in comps if metOriginale == comp.strip().lower()]
    if len(lkeggId_perfectMatch) == 0:
        for kId, comps in dRequest.items():
            for comp in comps:
                dfl1 = gL.extractRegexFromItem(metOriginale, r"([0-9a-z]+)")
                l1 = list(dfl1[0])
                l1.sort()
                dfl2 = gL.extractRegexFromItem(comp.strip().lower(), r"([0-9a-z]+)")
                l2 = list(dfl2[0])
                l2.sort()
                if l1 == l2:
                    lkeggId_perfectMatch.append(kId)
    if len(lkeggId_perfectMatch) == 0:
        lkeggId_perfectMatch = []
        for symbol in lSymbols:
            dRequest_copy = dRequest.copy()
            for kId in dRequest_copy:
                for comp in dRequest_copy[kId]:
                    comp = comp.strip().lower()
                    lPositions = [p.start() for p in re.finditer(symbol, comp)]
                    for posFound in lPositions:
                        try:
                            comp_subst = comp[:posFound] + ' ' + comp[posFound + 1:]
                        except:
                            comp_subst = comp[:posFound] + '+' + comp[posFound + 1:]
                        if comp_subst == met:
                            lkeggId_perfectMatch += [kId]
                        else:
                            comp_subst_woSpaces = replaceSpacesWithPlusAndSearch(comp_subst)
                            if comp_subst_woSpaces == met:
                                lkeggId_perfectMatch += [kId]
                            else:
                                comp_subst_woSpaces = replaceSpacesWithNoneAndSearch(comp_subst)
                                if comp_subst_woSpaces == met:
                                    lkeggId_perfectMatch += [kId]
                    try:
                        comp_Allsubst = re.sub(symbol, ' ', comp)
                    except:
                        comp_Allsubst = re.sub(symbol, '+', comp)
                    if comp_Allsubst == met:
                        lkeggId_perfectMatch += [kId]
                    else:
                        comp_Allsubst_woSpaces = replaceSpacesWithPlusAndSearch(comp_Allsubst)
                        if comp_Allsubst_woSpaces == met:
                            lkeggId_perfectMatch += [kId]
                        else:
                            comp_Allsubst_woSpaces = replaceSpacesWithNoneAndSearch(comp_Allsubst)
                            if comp_Allsubst_woSpaces == met:
                                lkeggId_perfectMatch += [kId]
    return lkeggId_perfectMatch

def replaceSpacesWithPlusAndSearch(metOriginale):
    metNew = re.sub(r'\s+', '+', metOriginale)
    return metNew

def substituteSymbolsAndSearch(met, lSymbols, db):
    lOutputs = []
    lSymbols2Change = [symbol for symbol in lSymbols if symbol in met]
    for k in lSymbols2Change:
        lPositions = [p.start() for p in re.finditer(k, met)]
        for posFound in lPositions:
            metSubstituted = met[:posFound] + '+' + met[posFound + 1:]
            try:
                lkeggId_perfectMatch = queryKeggCompound(db, metSubstituted, metSubstituted, lSymbols)
                lOutputs += lkeggId_perfectMatch
            except:
                metSubstituted = replaceSpacesWithPlusAndSearch(metSubstituted)
                try:
                    lkeggId_perfectMatch = queryKeggCompound(db, metSubstituted, metSubstituted, lSymbols)
                    lOutputs += lkeggId_perfectMatch
                except:
                    print('No results')

        met_AllsubstAll = re.sub(k, '+', met)
        try:
            lkeggId_perfectMatch = queryKeggCompound(db, met_AllsubstAll, met_AllsubstAll, lSymbols)
        except:
            met_AllsubstAll = replaceSpacesWithPlusAndSearch(met_AllsubstAll)
            try:
                lkeggId_perfectMatch = queryKeggCompound(db, met_AllsubstAll, met_AllsubstAll, lSymbols)
                lOutputs += lkeggId_perfectMatch
            except:
                print('No results')
        else:
            lOutputs += lkeggId_perfectMatch
    lOutputs = gL.unique(lOutputs)
    return lOutputs

def substituteWordsAndSearch(met, dWords, db, lSymbols):
    lOutputs = []
    lWords2Change = [k for k,v in dWords.items() if k in met]
    met_substAll = (met + '.')[:-1]
    for k in lWords2Change:
        for v in dWords[k]:
            metSubstituted = met.replace(k,v)
            met_substAll = met_substAll.replace(k,v)
            try:
                lkeggId_perfectMatch = queryKeggCompound(db, metSubstituted, metSubstituted, lSymbols)
            except:
                metSubstituted = replaceSpacesWithPlusAndSearch(metSubstituted)
                try:
                    lkeggId_perfectMatch = queryKeggCompound(db, metSubstituted, metSubstituted, lSymbols)
                    lOutputs += lkeggId_perfectMatch
                except:
                    print('No results')
            else:
                lOutputs += lkeggId_perfectMatch
            try:
                lOutputs_changeSymbol = substituteSymbolsAndSearch(metSubstituted, lSymbols, db)
                lOutputs += lOutputs_changeSymbol
            except:
                print('No results')
    try:
        lkeggId_perfectMatch = queryKeggCompound(db, met_substAll, met_substAll, lSymbols)
    except:
        met_substAll = replaceSpacesWithPlusAndSearch(met_substAll)
        try:
            lkeggId_perfectMatch = queryKeggCompound(db, met_substAll, met_substAll, lSymbols)
            lOutputs += lkeggId_perfectMatch
        except:
            print('No results')
    else:
        lOutputs += lkeggId_perfectMatch
    lOutputs = gL.unique(lOutputs)
    return lOutputs

def replaceSpacesWithNoneAndSearch(metOriginale):
    metNew = re.sub(r'\s+', '', metOriginale)
    return metNew

def name2ChebiIds_perfectMatch(met, dfChebiNames, dfChebiUniprot, dfChebiCompounds):
    lchebiId = []
    dfSearch_comps = dfChebiCompounds[dfChebiCompounds['NAME'] == met]
    if dfSearch_comps.empty is False:
        lchebiId += list(dfSearch_comps['ID'])

    dfSearch_names = dfChebiNames[dfChebiNames['NAME'] == met]
    if dfSearch_names.empty is False:
        lchebiId += list(dfSearch_names['COMPOUND_ID'])

    dfSearch_uniprot = dfChebiUniprot[dfChebiUniprot['NAME'] == met]
    if dfSearch_uniprot.empty is False:
        lchebiId += list(dfSearch_uniprot['ID'])
    return lchebiId

def fromChebi2KeggAndMetaCyc(chebiId, dfChebiDb, dfChebiRelations):
    dbRef = dfChebiDb[dfChebiDb['COMPOUND_ID'].isin(chebiId)]
    dbRef = dbRef.reset_index(drop = True)
    keggIdentifier = []
    possibleStrings = ['KEGG COMPOUND accession', 'KEGG GLYCAN accession', 'MetaCyc accession']
    for possibleString in possibleStrings:
        if dbRef.empty is False and any(ref == possibleString for ref in dbRef['TYPE']) is True:
            crossDb = dbRef[dbRef['TYPE'] == possibleString]
            keggIdentifier += list(crossDb['ACCESSION_NUMBER'])
        else:
            dfRel = dfChebiRelations[(dfChebiRelations['INIT_ID'].isin(chebiId)) | (dfChebiRelations['FINAL_ID'].isin(chebiId))]
            if dfRel.empty is False and any(rel.startswith('is_conjugate_') for rel in dfRel['TYPE']) is True:
                dfconj = dfRel[dfRel['TYPE'].str.startswith('is_conjugate_')]
                for row in dfconj.itertuples():
                    if row.INIT_ID in chebiId:
                        conj = row.FINAL_ID
                        dbRef = dfChebiDb[dfChebiDb['COMPOUND_ID'] == conj]
                        if dbRef.empty is False and any(ref == possibleString for ref in dbRef['TYPE']) is True:
                            crossDb = dbRef[dbRef['TYPE'] == possibleString]
                            keggIdentifier += list(crossDb['ACCESSION_NUMBER'])
                    elif row.FINAL_ID in chebiId:
                        conj = row.INIT_ID
                        dbRef = dfChebiDb[dfChebiDb['COMPOUND_ID'] == conj]
                        if dbRef.empty is False and any(ref == possibleString for ref in dbRef['TYPE']) is True:
                            crossDb = dbRef[dbRef['TYPE'] == possibleString]
                            keggIdentifier += list(crossDb['ACCESSION_NUMBER'])
    return keggIdentifier

def name2ChebiIds_subsSymbolsInChebiOutput(met, lSymbols, dfChebiNames, dfChebiUniprot, dfChebiCompounds):
    lchebiId = []
    s = 0
    dfChebiCompounds_copy_all = dfChebiCompounds.copy()
    dfChebiNames_copy_all = dfChebiNames.copy()
    dfChebiUniprot_copy_all = dfChebiUniprot.copy()
    met_substAll = (met + '.')[:-1]
    while s < len(lSymbols):
        met_substSingle = re.sub(lSymbols[s], '', met)
        met_substAll = re.sub(lSymbols[s], '', met_substAll)
        dfChebiCompounds_copy = dfChebiCompounds.copy()
        dfChebiCompounds_copy['NAME'] = dfChebiCompounds_copy.NAME.str.replace(lSymbols[s], '')
        dfChebiCompounds_copy['NAME'] = dfChebiCompounds_copy.NAME.str.strip()
        dfChebiCompounds_copy_all['NAME'] = dfChebiCompounds_copy_all.NAME.str.replace(lSymbols[s], '')
        dfChebiCompounds_copy_all['NAME'] = dfChebiCompounds_copy_all.NAME.str.strip()
        dfSearch_comps = dfChebiCompounds_copy[dfChebiCompounds_copy['NAME'] == met_substSingle.strip()]
        if dfSearch_comps.empty is False:
            lchebiId += list(dfSearch_comps['ID'])
        dfChebiNames_copy = dfChebiNames.copy()
        dfChebiNames_copy['NAME'] = dfChebiNames_copy.NAME.str.replace(lSymbols[s], '')
        dfChebiNames_copy['NAME'] = dfChebiNames_copy.NAME.str.strip()
        dfChebiNames_copy_all['NAME'] = dfChebiNames_copy_all.NAME.str.replace(lSymbols[s], '')
        dfChebiNames_copy_all['NAME'] = dfChebiNames_copy_all.NAME.str.strip()
        dfSearch_names = dfChebiNames_copy[dfChebiNames_copy['NAME'] == met_substSingle.strip()]
        if dfSearch_names.empty is False:
            lchebiId += list(dfSearch_names['COMPOUND_ID'])
        dfChebiUniprot_copy = dfChebiUniprot.copy()
        dfChebiUniprot_copy['NAME'] = dfChebiUniprot_copy.NAME.str.replace(lSymbols[s], '')
        dfChebiUniprot_copy['NAME'] = dfChebiUniprot_copy.NAME.str.strip()
        dfChebiUniprot_copy_all['NAME'] = dfChebiUniprot_copy_all.NAME.str.replace(lSymbols[s], '')
        dfChebiUniprot_copy_all['NAME'] = dfChebiUniprot_copy_all.NAME.str.strip()
        dfSearch_unip = dfChebiUniprot_copy[dfChebiUniprot_copy['NAME'] == met_substSingle.strip()]
        if dfSearch_unip.empty is False:
            lchebiId += list(dfSearch_unip['ID'])
        s += 1

    met_substAll = re.sub(r'\s+', '', met_substAll.strip())
    dfSearch_comps_all = dfChebiCompounds_copy_all[dfChebiCompounds_copy_all['NAME'] == met_substAll.strip()]
    if dfSearch_comps_all.empty is False:
        lchebiId += list(dfSearch_comps_all['ID'])
    dfSearch_names_all = dfChebiNames_copy_all[dfChebiNames_copy_all['NAME'] == met_substAll.strip()]
    if dfSearch_names_all.empty is False:
        lchebiId += list(dfSearch_names_all['COMPOUND_ID'])
    dfSearch_unip_all = dfChebiUniprot_copy_all[dfChebiUniprot_copy_all['NAME'] == met_substAll.strip()]
    if dfSearch_unip_all.empty is False:
        lchebiId += list(dfSearch_unip_all['ID'])

    return lchebiId

def name2ChebiIds_subsWords(met, dWords, lSymbols, dfChebiNames, dfChebiUniprot, dfChebiCompounds, dfChebiDb, dfChebiRelations):
    lAllId = []
    lWords2Change = [k for k,v in dWords.items() if k in met]
    met_substAll = (met + '.')[:-1]
    for k in lWords2Change:
        for v in dWords[k]:
            metSubstituted = met.replace(k,v)
            met_substAll = met_substAll.replace(k,v)
            lchebiId = name2ChebiIds_perfectMatch(metSubstituted, dfChebiNames, dfChebiUniprot, dfChebiCompounds)
            lAllId += lchebiId
            ## From the list of ChEBI identifiers, search the associated KEGG and MetaCyc identifiers
            keggIdentifiers = fromChebi2KeggAndMetaCyc(lchebiId, dfChebiDb, dfChebiRelations)
            lAllId += keggIdentifiers
            if len(keggIdentifiers) == 0:
                lchebiId = name2ChebiIds_subsSymbolsInChebiOutput(metSubstituted, lSymbols, dfChebiNames, dfChebiUniprot, dfChebiCompounds)
                lAllId += lchebiId
                keggIdentifiers = fromChebi2KeggAndMetaCyc(lchebiId, dfChebiDb, dfChebiRelations)
                lAllId += keggIdentifiers

    lchebiId = name2ChebiIds_perfectMatch(met_substAll, dfChebiNames, dfChebiUniprot, dfChebiCompounds)
    lAllId += lchebiId
    ## From the list of ChEBI identifiers, search the associated KEGG and MetaCyc identifiers
    keggIdentifiers = fromChebi2KeggAndMetaCyc(lchebiId, dfChebiDb, dfChebiRelations)
    lAllId += keggIdentifiers
    if len(keggIdentifiers) == 0:
        lchebiId = name2ChebiIds_subsSymbolsInChebiOutput(met_substAll, lSymbols, dfChebiNames, dfChebiUniprot, dfChebiCompounds)
        lAllId += lchebiId
        keggIdentifiers = fromChebi2KeggAndMetaCyc(lchebiId, dfChebiDb, dfChebiRelations)
        lAllId += keggIdentifiers

    return lAllId

def splitChebiResultsOnSymbol(df, colName, colId, symbol):
    lchebiId = []
    df[colName + '_par'] = df[colName].str.split(symbol)
    for row in df.itertuples():
        if pd.isna(getattr(row, colName + '_par')) is False:
            putativeChebi = symbol.join(getattr(row, colName + '_par')[:-1])
            if putativeChebi == met:
                lchebiId.append(getattr(row, colId))
    return lchebiId

def extractKeggIdComp(met, dfChebiNames, dfChebiDb, dfChebiRelations, dfChebiUniprot, dfChebiCompounds, dfChebiInchi, inchiOriginal):
    dWords = {'ic acid':['ate'], 'ate':['ic acid'], 'bisphosphate':['diphosphate'], 'diphosphate':['bisphosphate'],
                'aminium': ['amine'], 'amine': ['aminium'], 'ammonia': ['nh4+', 'nh3'], 'ammonium': ['nh4+', 'nh3'], 'proton': ['H+'],
                'adenosine triphosphate': ['atp'], 'adenosine diphosphate': ['adp'], 'coenzyme a': ['coa'], 'coa': ['coenzyme a'], 'apotransferin': ['apotransferrin'],
                "adenosine-5'-diphosphate": ['adp'], "adenosine 5'-diphosphate": ['adp'], "uridine-5'-diphosphate": ['udp'], "uridine 5'-diphosphate": ['udp'],
                "deoxyuridine-5'-diphosphate": ['dudp'], "deoxyuridine-5'-triphosphate": ['dutp'],
                'acp': ['[acp]', 'acyl-carrier protein', '[acyl-carrier protein]', 'acyl-carrier-protein', '[acyl-carrier-protein]'], '[acp]': ['acp'], "cytidine-5'-monophosphate": ['cmp'], "cytidine 5'-monophosphate": ['cmp'],
                '-ld-pe-pool': [''], '-ld-ps-pool': [''], '-ld-pc-pool': [''], '-ld-pe-pool': [''], '-ld-tg1-pool': [''], '-ld-tg2-pool': [''], '-ld-tg3-pool': [''], '-ld-pi-pool': [''], '-vldl-pool': [''], '-bile-pc-pool': [''],
                '-uptake-pool': [''], '-pool': [''], 'ide': ['ic acid']}
    lSymbols = ['/', '-', "'", '\[', '\]', '\(', '\)']
    keggId =  []

    lDbs = ['compound', 'glycan']
    for database in lDbs:
        try:
            lkeggId_perfectMatch = queryKeggCompound(database, met, met, lSymbols)
        except:
            try:
                metNew = replaceSpacesWithPlusAndSearch(met)
                lkeggId_perfectMatch = queryKeggCompound(database, metNew, met, lSymbols)
            except:
                if any(symbol in met for symbol in lSymbols) is True:
                    lOutputs = substituteSymbolsAndSearch(met, lSymbols, database)
                    keggId += lOutputs
                if any(k in met for k,v in dWords.items()) is True:
                    lOutputs = substituteWordsAndSearch(met, dWords, database, lSymbols)
                    keggId += lOutputs
            else:
                keggId += lkeggId_perfectMatch
        else:
            if len(lkeggId_perfectMatch) != 0:
                keggId += lkeggId_perfectMatch
            if any(symbol in met for symbol in lSymbols) is True:
                lOutputs = substituteSymbolsAndSearch(met, lSymbols, database)
                keggId += lOutputs
            if any(k in met for k,v in dWords.items()) is True:
                lOutputs = substituteWordsAndSearch(met, dWords, database, lSymbols)
                keggId += lOutputs

        try:
            metNew_woSpaces = replaceSpacesWithNoneAndSearch(met)
            lkeggId_perfectMatch = queryKeggCompound(database, metNew_woSpaces, met, lSymbols)
        except:
            print('not found')
        else:
            keggId += lkeggId_perfectMatch

        lPositions = [p.start() for p in re.finditer(' ', met)]
        for posFound in lPositions:
            met_subst = met[:posFound] + met[posFound + 1:]
            try:
                lkeggId_perfectMatch = queryKeggCompound(database, met_subst, met_subst, lSymbols)
            except:
                print('not found')
            else:
                keggId += lkeggId_perfectMatch

        met_subst = (met + '.')[:-1]
        for posFound in lPositions:
            met_subst = met_subst[:posFound] + met_subst[posFound + 1:]
        try:
            lkeggId_perfectMatch = queryKeggCompound(database, met_subst, met_subst, lSymbols)
        except:
            print('not found')
        else:
            keggId += lkeggId_perfectMatch

    ## search the metabolite in ChEBI
    lchebiId = name2ChebiIds_perfectMatch(met, dfChebiNames, dfChebiUniprot, dfChebiCompounds)
    keggId += lchebiId
    ## And from the list of ChEBI identifiers, search information in ChEBI about KEGG and MetaCyc identifiers
    keggIdentifiers = fromChebi2KeggAndMetaCyc(lchebiId, dfChebiDb, dfChebiRelations)
    keggId += keggIdentifiers
    ## If the perfect match is not found, try to manipulate the name of the metabolite
    lchebiId = name2ChebiIds_subsSymbolsInChebiOutput(met, lSymbols, dfChebiNames, dfChebiUniprot, dfChebiCompounds)
    keggId += lchebiId
    keggIdentifier = fromChebi2KeggAndMetaCyc(lchebiId, dfChebiDb, dfChebiRelations)
    keggId += keggIdentifier
    ## try to substitute words within the metabolite name to find the match in ChEBI or if not found try again to manipulate the string
    lAllId = name2ChebiIds_subsWords(met, dWords, lSymbols, dfChebiNames, dfChebiUniprot, dfChebiCompounds, dfChebiDb, dfChebiRelations)
    keggId += lAllId

    ## try to search the string before the parenthesis or the comma symbol
    if len(keggId) == 0:
        lchebiId_fromChebiCompounds = splitChebiResultsOnSymbol(dfChebiCompounds, 'NAME', 'ID', '(')
        lchebiId_fromChebiNames = splitChebiResultsOnSymbol(dfChebiNames, 'NAME', 'COMPOUND_ID', '(')
        lchebiId_fromChebiUniprot = splitChebiResultsOnSymbol(dfChebiUniprot, 'NAME', 'ID', '(')
        lchebiId = lchebiId_fromChebiCompounds + lchebiId_fromChebiNames + lchebiId_fromChebiUniprot
        keggId += lchebiId
        keggIdentifier = fromChebi2KeggAndMetaCyc(lchebiId, dfChebiDb, dfChebiRelations)
        keggId += keggIdentifier

    if len(keggId) == 0:
        lchebiId_fromChebiCompounds = splitChebiResultsOnSymbol(dfChebiCompounds, 'NAME', 'ID', ',')
        lchebiId_fromChebiNames = splitChebiResultsOnSymbol(dfChebiNames, 'NAME', 'COMPOUND_ID', ',')
        lchebiId_fromChebiUniprot = splitChebiResultsOnSymbol(dfChebiUniprot, 'NAME', 'ID', ',')
        lchebiId = lchebiId_fromChebiCompounds + lchebiId_fromChebiNames + lchebiId_fromChebiUniprot
        keggId += lchebiId
        keggIdentifier = fromChebi2KeggAndMetaCyc(lchebiId, dfChebiDb, dfChebiRelations)
        keggId += keggIdentifier

    if len(keggId) == 0 and pd.isna(inchiOriginal) is False and inchiOriginal != '':
        lParts = inchiOriginal.split('/')
        originalInchi_formula = lParts[1]
        i = 2
        originalInchi_atomConnection = ''
        originalInchi_hydrogen = ''
        while i < len(lParts):
            if lParts[i].startswith('c'):
                originalInchi_atomConnection = lParts[i]
            elif lParts[i].startswith('h'):
                originalInchi_hydrogen = lParts[i]
            i+=1

        lchebiId = []
        for row in dfChebiInchi.itertuples():
            currentChebiInchi = row.InChI_splitted
            currentInchi_formula = currentChebiInchi[1]
            i = 2
            currentInchi_atomConnection = ''
            currentInchi_hydrogen = ''
            while i < len(currentChebiInchi):
                if currentChebiInchi[i].startswith('c'):
                    currentInchi_atomConnection = currentChebiInchi[i]
                elif currentChebiInchi[i].startswith('h'):
                    currentInchi_hydrogen = currentChebiInchi[i]
                i+=1
            if originalInchi_formula != '' and currentInchi_formula != '' and originalInchi_atomConnection != '' and currentInchi_atomConnection != '' and originalInchi_hydrogen != '' and currentInchi_hydrogen != '':
                if originalInchi_formula == currentInchi_formula and originalInchi_atomConnection == currentInchi_atomConnection and originalInchi_hydrogen == currentInchi_hydrogen:
                    lchebiId.append(row.CHEBI_ID)

        keggId += lchebiId
        ## From the list of ChEBI identifiers search the associated KEGG and MetaCyc identifiers
        keggIdentifiers = fromChebi2KeggAndMetaCyc(lchebiId, dfChebiDb, dfChebiRelations)
        keggId += keggIdentifiers

    if len(keggId) == 0:
        for database in lDbs:
            if '(' in met:
                lPositions = [p.start() for p in re.finditer('\(', met)]
                met_untilParenthesis = met[:lPositions[-1]].strip()
                try:
                    lkeggId_perfectMatch = queryKeggCompound(database, met_untilParenthesis, met_untilParenthesis, lSymbols)
                except:
                    if any(symbol in met_untilParenthesis for symbol in lSymbols) is True:
                        lOutputs = substituteSymbolsAndSearch(met_untilParenthesis, lSymbols, database)
                        keggId += lOutputs
                    if any(k in met_untilParenthesis for k,v in dWords.items()) is True:
                        lOutputs = substituteWordsAndSearch(met_untilParenthesis, dWords, database, lSymbols)
                        keggId += lOutputs
                else:
                    keggId += lkeggId_perfectMatch

            if ',' in met:
                lPositions = [p.start() for p in re.finditer(',', met)]
                met_untilParenthesis = met[:lPositions[-1]].strip()
                try:
                    lkeggId_perfectMatch = queryKeggCompound(database, met_untilParenthesis, met_untilParenthesis, lSymbols)
                except:
                    if any(symbol in met_untilParenthesis for symbol in lSymbols) is True:
                        lOutputs = substituteSymbolsAndSearch(met_untilParenthesis, lSymbols, database)
                        keggId += lOutputs
                    if any(k in met_untilParenthesis for k,v in dWords.items()) is True:
                        lOutputs = substituteWordsAndSearch(met_untilParenthesis, dWords, database, lSymbols)
                        keggId += lOutputs
                else:
                    keggId += lkeggId_perfectMatch

    keggId = gL.unique(keggId)
    keggId = [k for k in keggId if k != '']
    return keggId

def getMetsInfoGEM(model):
    """ get all info about metabolites in GEMS"""
    cId=[]
    cName=[]
    cKegg=[]
    cChebi = []
    cPubchem = []
    cBoundaryCondition = []
    cChemicalFormula = []
    cInchi = []
    tree = ET.parse(model)
    root = tree.getroot()
    ## Definition of useful tags
    annotTag = ''
    liTag = ''
    notesTag = ''
    bodyTag = ''
    for  element in root.getiterator():
        if element.tag.endswith('listOfSpecies'):
            listSpeciesTag = element.tag
        elif element.tag.endswith('annotation'):
            annotTag = element.tag
        elif element.tag.endswith('li'):
            liTag = element.tag
        elif element.tag.endswith('notes'):
            notesTag = element.tag
        elif element.tag.endswith('body'):
            bodyTag = element.tag

    for species in root.getiterator(listSpeciesTag):
        for el in species:
            keggCompound = ''
            chebiCompound = ''
            pubchemCompound = ''
            chemicalFormula = ''
            inchi = ''
            cId.append(el.attrib['id'])
            cBoundaryCondition.append(el.attrib['boundaryCondition'])
            lChemicalFormula = [k for k, v in el.attrib.items() if k.endswith('chemicalFormula')]
            if len(lChemicalFormula) != 0:
                chemicalFormula = el.attrib[lChemicalFormula[0]]
            cChemicalFormula.append(chemicalFormula)
            name = ''
            if 'name' in el.attrib:
                name = el.attrib['name']
            cName.append(name)
            if annotTag != '':
                for ell in el.getiterator(annotTag):
                    if ell.tag != '' and liTag != '':
                        for annot in ell.getiterator(liTag):
                            for k,v in annot.attrib.items():
                                if 'kegg.compound' in v:
                                    keggCompound = v[-6:]
                                if '/chebi/' in v:
                                    posStart = v.find('CHEBI:')
                                    chebiCompound = v[posStart+len('CHEBI:'):]
                                    print('chebiCompound\t', chebiCompound)
                                if 'pubchem.compound' in v:
                                    pubchemCompound = v.split('/')[-1]
                                if '/inchi/' in v:
                                    posStart = v.find('/inchi/')
                                    inchi = v[posStart+len('/inchi/'):]
            if (inchi == '' or keggCompound == '') and notesTag != '':
                for ell in el.getiterator(notesTag):
                    if ell.tag != '' and bodyTag != '':
                        for annot in ell.getiterator(bodyTag):
                            for p in annot:
                                if p.text.startswith('metKEGGID:'):
                                    keggCompound = p.text[10:]
                                elif p.text.startswith('metInChIString:'):
                                    inchi = p.text[15:]
            cInchi.append(inchi)
            cKegg.append(keggCompound)
            cChebi.append(chebiCompound)
            cPubchem.append(pubchemCompound)
    return cId, cName, cKegg, cChebi, cPubchem, cBoundaryCondition, cChemicalFormula, cInchi

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

def filterFWoutputs(dfMappingInput):
    dMapping_100 = {}
    dMapping_91_99 = {}
    dMapping_empty = {}
    for index, row in dfMappingInput.iterrows():
        col = 0
        threshold = 100
        while col < 10 and threshold > 90:
            try:
                tMatch = ast.literal_eval(row[col])
                match = tMatch[0]
                threshold = tMatch[1]
            except:
                tMatch = row[col]
                lMatch = row[col].split(',')
                match = lMatch[0].strip()[lMatch[0].find('(')+1:]
                threshold = int(lMatch[1].strip()[:lMatch[1].find(')')-1])

            if threshold == 100:
                if index not in dMapping_100:
                    dMapping_100[index] = [match]
                else:
                    dMapping_100[index] += [match]
                col = 10
            elif threshold < 100 and threshold > 90:
                if index not in dMapping_91_99:
                    dMapping_91_99[index] = [match]
                else:
                    dMapping_91_99[index] += [match]
                col += 1
            elif threshold <= 90:
                if index not in dMapping_empty:
                    dMapping_empty[index] = []
                else:
                    dMapping_empty[index] += []
                col = 10

    dfMatches100 = pd.DataFrame(dMapping_100.items(), columns=['Name', 'Matches'])
    dfMatches91_99 = pd.DataFrame(dMapping_91_99.items(), columns=['Name', 'Matches'])
    dfMatchesEmpty = pd.DataFrame(dMapping_empty.items(), columns=['Name', 'Matches'])
    return dfMatches100, dfMatches91_99, dfMatchesEmpty
