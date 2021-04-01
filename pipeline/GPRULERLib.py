# -*- coding: utf-8 -*-
import os
import sys
import json
import requests
import xmltodict
# from sympy import *
import pandas as pd
import re as reModule
import genericLib as gL
import bioservices.kegg as kegg
from bioservices import UniProt
from nltk.tokenize import word_tokenize
from nltk.tokenize import RegexpTokenizer
import nltk
nltk.download('punkt')

def getUniprotAndComplexPortalData(df, gene='uniprotId'):
    '''
    This function retrieves from Uniprot and Complex Portal databases protein quaternary structures,
    protein macromolecular complexes and protein-protein interactions established by a given
    metabolic gene.
    Input:
    - df: df kegg2Uniprot identifiers conversion.
    - gene: column name of df variable containing Uniprot identifiers
    Output:
    - Initial dataframe enriched with information retrieved from Uniprot and Complex Portal.
    '''
    u = UniProt(verbose=False)
    lAllProteinNames = []
    lAllGeneNames =  []
    lUniprotIds = []
    lsubunits = []
    lComplexPortal = []
    lAllComplexPortal_uniprotId = []
    lAllComplexPortal_protName = []
    lFunction = []
    for el in df[gene]:
        res = u.search("%s" % el, frmt="xml")
        proteinNames = []
        geneNames = []
        complexPortal = ''
        subunitDescription = ''
        uniprotIds = []
        complexPortal_uniprotId = []
        complexPortal_protName = []
        if res != '':
            dRes = xmltodict.parse(res)
            try:
                for ref in dRes['uniprot']['entry']['dbReference']:
                    if ref['@type'] == 'ComplexPortal':
                        complexPortal = ref['@id']
            except:
                complexPortal = ''

            if complexPortal != '':
                url = 'https://www.ebi.ac.uk/intact/complex-ws/complex/' + complexPortal
                response = requests.get(url, headers={"Accept": "application/json"})
                jsonData = json.loads(response.text)
                try:
                    for participant in jsonData['participants']:
                        if participant['interactorType'] == 'protein':
                            complexPortal_uniprotId.append(participant['identifier'])
                            complexPortal_protName.append(participant['name'])
                except:
                    print('no complex portal element')

            try:
                for k, val in dRes['uniprot']['entry']['protein']['recommendedName'].items():
                    if k == 'fullName' or k == 'shortName':
                        if type(val) == str:
                            proteinNames.append(val)
                        else:
                            proteinNames.append(val['#text'])
            except:
                print('no protein recommended name')
            try:
                for diz in dRes['uniprot']['entry']['protein']['alternativeName']:
                    for k, val in diz.items():
                        if k == 'fullName' or k == 'shortName':
                            if type(val) == str:
                                proteinNames.append(val)
                            else:
                                proteinNames.append(val['#text'])
            except:
                print('no protein alternative name')
            try:
                for diz in dRes['uniprot']['entry']['gene']['name']:
                    for k, val in diz.items():
                        if k == '#text':
                            geneNames.append(val)
            except:
                print('no gene name')

            function = ''
            try:
                for el in dRes['uniprot']['entry']['comment']:
                    if el['@type'] == 'function' and type(el['text']) == str:
                        function = el['text']
                    elif el['@type'] == 'function' and type(el['text']) != str:
                        function = el['text']['#text']
            except:
                print('no function')

            subunitDescription = ''
            if 'comment' in dRes['uniprot']['entry'].keys() and type(dRes['uniprot']['entry']['comment']) == list:
                try:
                    for diz in dRes['uniprot']['entry']['comment']:
                        if diz['@type'] == 'subunit':
                            for k, val in diz.items():
                                if k == 'text':
                                    subunitDescription = val['#text']
                except:
                    for el in dRes['uniprot']['entry']['comment']:
                        if el['@type'] == 'subunit':
                            subunitDescription = el['text']
            try:
                if type(dRes['uniprot']['entry']['accession']) != str:
                    uniprotIds = dRes['uniprot']['entry']['accession']
                else:
                    uniprotIds = [dRes['uniprot']['entry']['accession']] ## questa è una lista
            except:
                uniprotIds = []
        lAllProteinNames.append(proteinNames)
        lAllGeneNames.append(geneNames)
        lUniprotIds.append(uniprotIds)
        lFunction.append(function)
        lsubunits.append(subunitDescription)
        lComplexPortal.append(complexPortal)
        lAllComplexPortal_uniprotId.append(complexPortal_uniprotId)
        lAllComplexPortal_protName.append(complexPortal_protName)
    df['proteinNames'] = lAllProteinNames
    df['geneNames'] = lAllGeneNames
    df['txt_subunit'] = lsubunits
    df['function'] = lFunction
    df['id_uniprot'] = lUniprotIds
    df['complexPortal'] = lComplexPortal
    df['complexPortal_uniprotId'] = lAllComplexPortal_uniprotId
    df['complexPortal_protName'] = lAllComplexPortal_protName
    u1mer = []
    for i in range(0, len(df)):
        struttura = [word for word in word_tokenize(df['txt_subunit'][i], "english", False) if 'mer' in el]
        if struttura != []:
            u1mer.append(','.join(struttura))
        else:
            u1mer.append('')
    df['structure'] = u1mer
    return(df)

def textMiningFromUniprot(df):
    '''
    This function performs text mining from "Function" and "Interaction" sections of Uniprot database
    retrieved by getUniprotAndComplexPortalData function.
    Input:
    - df: dataframe generated by getUniprotAndComplexPortalData function
    Output:
    - Enriched input dataframe where for each metabolic gene product genes participating in the same
    protein complexes, together with genes coding for isoforms are identified.
    '''
    tokenizer = RegexpTokenizer(r'\w+')
    lOtherIsoforms = []
    lAllSubunitsAllGenes = []
    lSubunitsFromName = []
    lInteract = []
    lby_similarity = []
    lgene_from_structure = []
    lIsoformIndication = []
    lSubunitIndication = []
    lSameEnzymeMembership = []
    lredundancies = []
    for row in df.itertuples():
        funzione = row.function
        testo = row.txt_subunit
        lPossibleSubunitsFromName = []

        namesWithIsoformIndication = False
        lAllNames = row.proteinNames + row.geneNames
        for el in lAllNames:
            if ('isoform' in el) or ('isozyme' in el) or ('isoenzyme' in el):
                namesWithIsoformIndication = True

        namesWithProteinIndication = []
        lnameEnzyme = []
        possibleWords = ['subunit', 'chain']
        subunitIndication = False
        for p in possibleWords:
            namesWithSubunitsIndication = [el for el in lAllNames if p in el]
            if len(namesWithSubunitsIndication) != 0:
                subunitIndication = True
            for n in namesWithSubunitsIndication:
                if n.endswith(p) or p + ',' in n:
                    lPossibleSubunitsFromName.append(n.split(p)[0].split()[-1].strip(','))
                    lnameEnzyme.append(' '.join(n.split(p)[0].split()[:-1]).lower().strip())
                else:
                    lPossibleSubunitsFromName.append(n.split(p)[1].split()[0].strip(','))
                    lnameEnzyme.append(n.split(p)[0].lower().strip())

        sameEnzymeMembership = False # Returns True or False depending on whether indications of at least one of enzyme names within Function section of Uniprot is found
        pos = [s.start() for s in reModule.finditer('is one of', funzione.lower())]
        for p in pos:
            funzione_sub = funzione[p:]
            if any(el in funzione_sub for el in lnameEnzyme) is True:
                sameEnzymeMembership = True
        redundancy = []
        pos = [s.start() for s in reModule.finditer('redundant with', funzione.lower())]
        for p in pos:
            funzione_sub = funzione[p:].split('.')[0]
            redundancy += [word for word in word_tokenize(funzione_sub, 'english', False) if word.isupper() and len(word) > 2]

        pos = [s.start() for s in reModule.finditer('mer', testo.lower())]
        genefromstructure = []
        for p in pos:
            testo_sub = testo[p:]
            try:
                genefromstructure += [word for word in word_tokenize(testo_sub.split('mer')[1].split(".")[0], 'english', False) if word.isupper() and len(word) > 2]
            except:
                print('no gene structure')

        pos = [s.start() for s in reModule.finditer('interact', testo.lower())]
        interactPart = []
        for p in pos:
            testo_sub = testo[p:]
            try:
                interactPart += [word for word in word_tokenize(testo_sub.split("interact")[1].split(".")[0], 'english', False) if word.isupper() and len(word) > 2]
            except:
                interactPart += [word for word in word_tokenize(testo_sub.split("Interact")[1].split(".")[0], 'english', False) if word.isupper() and len(word) > 2]

        bysimilarity = []
        if "by similarity" in testo.lower():
            try:
                bysimilarity += [word for word in word_tokenize(row.txt_subunit.split("By similarity")[0].split(".")[1], 'english', False) if word.isupper() and len(word) > 2]
            except:
                bysimilarity += [word for word in word_tokenize(row.txt_subunit.split("By similarity")[0], 'english', False) if word.isupper() and len(word) > 2]

        lInteract.append(interactPart)
        lby_similarity.append(bysimilarity)
        lgene_from_structure.append(genefromstructure)
        lIsoformIndication.append(namesWithIsoformIndication)

        lOtherSubunits = []
        lIsoforms = []

        lChainWords = ['chain', 'chains']
        for chain in lChainWords:
            pos = [s.start() for s in reModule.finditer(chain, testo.lower())]
            for p in reversed(pos):
                testo_sub = testo[:p]
                lOtherSubunits.append(testo_sub.split()[-1])
                lOtherSubunits.append(testo_sub.split()[-2])

        lOtherSubunitsWords = ['complex', 'component', 'composed of', 'bound to', 'mer of', 'mer with',
                               'mers of', 'mers with', 'interacts with','interact with','binds','bind',
                               'together with', 'comprises', 'comprise', 'comprising', 'containing',
                               'part of a complex with', 'found in a complex with', 'form', 'forms',
                               'binds to', 'bind to', 'consisting of', 'consists of','consist of',
                               'associates with', 'merization with', 'identified in complexes with',
                               'contains', 'types of subunits', 'part of', 'interacts','interact',
                               'associate with', 'formed of', 'merizes with', 'merize with',
                               'mers including', 'mer including', 'complex', 'component']

        # Split entire text by individual sentences and analyse each one
        lTestoSplt = testo.strip('.').split('.')
        for frase in lTestoSplt:
            ## Select from lOtherSubunitsWords list only the words that are present in the sentence under consideration
            wordsSelection = [el for el in lOtherSubunitsWords if el in frase.lower()]
            columns = ['proteinNames', 'geneNames']
            for col in columns:
                orPosition = frase.find(' or ')
                frase_before = frase[:orPosition]
                frase_after = frase[orPosition-1+len(' or '):]
                while (orPosition != -1) and ((any(frase_before.endswith(el) for el in getattr(row, col)) is True) or (any(frase_after.startswith(el) for el in getattr(row, col)) is True)):
                    tmpIsoformList = []
                    foundWords = [el for el in getattr(row, col) if el in frase]
                    wordsToRemove = []
                    for f in foundWords:
                        lfraseStrip = frase.strip(f).strip()
                        lfraseSplt = lfraseStrip.split(f)
                        if len(lfraseSplt) == 2:
                            el0 = lfraseSplt[0].strip()
                            el1 = lfraseSplt[1].strip()
                            if el0.endswith(' or') is True and col == 'proteinNames':
                                selectedIsoform = el0.strip('or').strip().split(' ')[-1]
                                tmpIsoformList.append(selectedIsoform)
                                lIsoforms.append(selectedIsoform)
                                wordsToRemove.append(selectedIsoform + ' or ')
                            if el1.startswith('or ') is True and col == 'proteinNames':
                                selectedIsoform = el1.strip('or').strip().split(' ')[0]
                                tmpIsoformList.append(selectedIsoform)
                                lIsoforms.append(selectedIsoform)
                                wordsToRemove.append(' or ' + selectedIsoform)
                            if el0.endswith(' or') is True and col == 'geneNames':
                                out0 = el0.strip('or').strip().split()[-1]
                                dfout0 = gL.extractRegexFromItem(out0, r"([A-Za-z0-9-]+)")
                                lout0 = list(dfout0[0])
                                tmpIsoformList += lout0
                                lIsoforms += lout0
                                wordsToRemove.append(lout0[0] + ' or ')
                            if el1.startswith('or ') is True and col == 'geneNames':
                                out1 = el1.strip('or').strip().split()[0]
                                dfout1 = gL.extractRegexFromItem(out1, r"([A-Za-z0-9-]+)")
                                lout1 = list(dfout1[0])
                                tmpIsoformList += lout1
                                lIsoforms += lout1
                                wordsToRemove.append(' or ' + lout1[0])
                    if len(wordsToRemove) != 0:
                        for w in wordsToRemove:
                            frase = reModule.sub(w, "", frase)
                        orPosition = frase.find(' or ')
                        if orPosition != -1:
                            frase_before = frase[:orPosition]
                            frase_after = frase[orPosition-1+len(' or '):]
                        else:
                            frase_before = frase[:]
                            frase_after = frase[:]
                    else:
                        orPosition = -1

            ## Once updated the sentence deprived of genes identified as isoforms, look for AND relationships
            for w in wordsSelection:
                pos = [s.start() for s in reModule.finditer(w, frase.lower())]
                for p in pos:
                    frase_sub = frase[p:]
                    if 'isoforms of' not in frase_sub:
                        try:
                            frase_sub_splitted = reModule.split(w, frase_sub.lower())[1]
                        except:
                            frase_sub_splitted = reModule.split(w[0].upper() + w[1:], frase_sub)[1]
                        frase_sub_splitted_tokenized = tokenizer.tokenize(frase_sub_splitted)
                        lOtherSubunits += frase_sub_splitted_tokenized

        lOtherIsoforms.append(lIsoforms)
        lSubunitsFromName.append(lPossibleSubunitsFromName + namesWithProteinIndication)
        lAllSubunitsAllGenes.append(lOtherSubunits)
        lSubunitIndication.append(subunitIndication)
        lSameEnzymeMembership.append(sameEnzymeMembership)
        lredundancies.append(redundancy)

    df['interact'] = lInteract
    df['by_similarity'] = lby_similarity
    df['gene_from_structure'] = lgene_from_structure
    df['otherSubunits'] = lAllSubunitsAllGenes
    df['subunitsFromName'] = lSubunitsFromName
    df['otherIsoforms'] = lOtherIsoforms
    df['isoformIndication'] = lIsoformIndication
    df['subunitIndication'] = lSubunitIndication
    df['sameEnzymeMembership'] = lSameEnzymeMembership
    df['redundancy'] = lredundancies
    return(df)

def getStringData(df, gene='id_uniprot'):
    '''
    This function retrieves from STRING database known and predicted protein-protein interactions
    established by each queried metabolic gene.
    Input:
    - df: dataframe generated by the previous step;
    - id_prot: column name of uniprot identifiers of genes. By default it is set equal to 'id_uniprot'.
    Output:
    - df: enriched input dataframe with information retrieved from STRING database.
    '''
    dizUniprotString = {}
    dizlInteractors = {}
    stringSubunits = []
    for row in df.itertuples():
        lStringInteractors = []
        for uniprotId in getattr(row, gene):
            if uniprotId not in dizUniprotString:
                originalUniprotNames = row.proteinNames + row.geneNames
                url = "https://string-db.org/api/json/network?identifiers=" + uniprotId
                response = requests.get(url, verify=False)
                while response.status_code == 524:
                    response = requests.get(url, verify=False)
                net = response.json()
                lInteractors = []
                original = ''
                if type(net) == list:
                    for i in range(0, len(net)):
                        if (net[i]['preferredName_A'] == uniprotId) or (any(net[i]['preferredName_A'] == name for name in originalUniprotNames) is True):
                            original = net[i]['preferredName_A']
                            lInteractors.append(net[i]['preferredName_B'])
                        elif (net[i]['preferredName_B'] == uniprotId) or (any(net[i]['preferredName_B'] == name for name in originalUniprotNames) is True):
                            lInteractors.append(net[i]['preferredName_A'])
                dizUniprotString[uniprotId] = (original, lInteractors)
            else:
                original = dizUniprotString[uniprotId][0]
                lInteractors = dizUniprotString[uniprotId][1]
            if original != '':
                for interactor in lInteractors:
                    if interactor not in dizlInteractors:
                        url = 'https://string-db.org//api/json/enrichment?identifiers=' + original + '%0d' + interactor
                        response = requests.get(url, verify=False)
                        while response.status_code == 524:
                            response = requests.get(url, verify=False)
                        diz = response.json()
                        dizlInteractors[interactor] = diz
                    else:
                        diz = dizlInteractors[interactor]
                    for stringElement in diz:
                        if stringElement['category'] == 'Component' and ("complex" in str(stringElement['description']) or "chain" in str(stringElement['description'])) and original in stringElement['inputGenes']:
                            complesso = stringElement['description']
                            lStringInteractors.append(interactor)
        lStringInteractors = gL.unique(lStringInteractors)
        stringSubunits.append(lStringInteractors)

    df['stringSubunits'] = stringSubunits
    return(df)

def getKeggData(df, organism, gene='keggId'):
    '''
    This function retrieves from KEGG database protein isoforms established by each queried metabolic gene.
    Input:
    - df: dataframe generated by getStringData function;
    - gene: column name of KEGG identifiers of genes;
    - organism: KEGG code of the organism under investigation.
    Output:
    - df: enriched input dataframe with information retrieved from KEGG database.
    '''
    k = kegg.KEGG()
    list_iso = []
    lNames = []
    for el in df[gene]:
        iso = ""
        try:
            find = str(k.parse(k.get(organism + ":" + str(el)))['ORTHOLOGY'].keys()).split("'", 1)[1].split("'", 1)[0]
        except:
            find = "Nessuna isoforma trovata"
            name = ''

        try:
            name = k.parse(k.get(organism + ":" + str(el)))['NAME'][0].split(',')
        except:
            name = ''

        isoform = ""
        try:
            isoform = k.parse(k.get(find))['GENES'][organism.upper()]
            for i in range(0, len(isoform.split(" "))):
                iso = isoform.split(" ")[i].split("(", 1)[1].split(")", 1)[0] + " " + iso
        except:
            iso = ""
        list_iso.append(iso)
        lNames.append(name)

    lFinalIsoforms = []
    for item in list_iso:
        lFinalIsoforms.append(list(word_tokenize(item)))

    df['isoform'] = lFinalIsoforms
    df['geneName_fromKEGG'] = lNames
    return (df)

def mergeData(df):
    '''
    This function joins all retrieved information from explored databases and assembles them to generate the final GPR rule.
    Input:
    - df: dataframe generated by getKeggData function.
    Output:
    - dfFinal: dataframe where for each input gene the list of its AND and OR relationships is returned.
    '''
    list_and = []
    list_or = []
    lTotalGenes = []
    lTotalGenes = [item for elem in list(df['geneName_fromKEGG']) for item in elem] + [item for elem in list(df['proteinNames']) for item in elem] \
                   + [item for elem in list(df['geneNames']) for item in elem] + [item for elem in list(df['id_uniprot']) for item in elem] \
                   + [item for elem in list(df['subunitsFromName']) for item in elem]
    lTotalGenes = gL.unique(lTotalGenes)

    lTotalAndSubs = lTotalGenes + list(df['subunitsFromName'])

    dErroneousNames = {'SLCA7A7': 'SLC7A7', 'SLCA7A11': 'SLC7A11'}
    for r in df.itertuples():
        ## clean all information retrieved from the explored databases by removing those elements that are not included into the lTotalGenes list
        finallDipComplex = []
        finallDipBinary = []

        dfinallDipBinary_names = {}
        for f in finallDipBinary:
            out = df.loc[df.uniprotId == f]
            for o in out.itertuples():
                dfinallDipBinary_names[f] = gL.unique([o.uniprotId] + o.proteinNames + o.geneNames + o.id_uniprot + o.geneName_fromKEGG)

        lComplexPortal_unipId = [x for x in list(r.complexPortal_uniprotId) if x not in list(r.id_uniprot)]
        finallComplexPortal_unipId = [x for x in lComplexPortal_unipId if x in lTotalAndSubs] # Select only the isoforms falling within input genes list

        lComplexPortal_protName = [x for x in list(r.complexPortal_protName) if x not in list(r.id_uniprot)]
        finallComplexPortal_protName = [x for x in lComplexPortal_protName if x in lTotalAndSubs]

        try:
            ldfStructure = [x for x in list(r.gene_from_structure) if x in lTotalAndSubs]
        except:
            ldfStructure = []

        try:
            ldfInteract = [x for x in list(r.interact) if x in lTotalAndSubs]
        except:
            ldfInteract = []

        try:
            ldfSimilarity = [x for x in list(r.by_similarity) if x in lTotalAndSubs]
        except:
            ldfSimilarity = []

        lAllIsoforms = list(r.otherIsoforms) + list(r.isoform)
        lAllIsoforms = gL.unique(lAllIsoforms)

        finalStringSubs = [x for x in list(r.stringSubunits) if x in lTotalAndSubs]

        lOtherSubunits = [x for x in list(r.otherSubunits) if x not in list(r.id_uniprot) + list(r.proteinNames) + list(r.geneNames) +
                          list(r.subunitsFromName) + list(r.geneName_fromKEGG)]
        finalOtherSubunits = [x for x in lOtherSubunits if x in lTotalAndSubs]

        finallDipBinary_woIsoforms = []
        for k in dfinallDipBinary_names:
            if all(el not in lAllIsoforms for el in dfinallDipBinary_names[k]) is True:
                finallDipBinary_woIsoforms.append(k)

        subunit = []
        subunit += finallDipComplex + finallDipBinary_woIsoforms + finallComplexPortal_unipId + finallComplexPortal_protName
        subunit += ldfInteract + ldfSimilarity + gL.difference(finalStringSubs, lAllIsoforms) + finalOtherSubunits
        subunit = gL.unique(subunit)

        # uso il dizionario di nomi che sono sbagliati e che è necessario correggere
        finalSubunitSet = []
        for s in subunit:
            if s in dErroneousNames:
                finalSubunitSet.append(dErroneousNames[s])
            else:
                finalSubunitSet.append(s)

        lIsoforms = [x for x in list(r.otherIsoforms) + list(r.isoform) if
                     x not in list(r.id_uniprot) + list(r.proteinNames) + list(r.geneNames) +
                     list(r.subunitsFromName) + list(r.geneName_fromKEGG)]
        finalIsoforms = [x for x in lIsoforms if x in lTotalAndSubs]

        list_and.append(finalSubunitSet)
        list_or.append(finalIsoforms)

    dfFinal = pd.DataFrame({'gene': df['geneNames'], 'uniprotId': df['id_uniprot']})
    dfFinal['AND'] = list_and
    dfFinal['OR'] = list_or
    return(dfFinal)

def putativeOrganisms(organism):
    df = getOrgs()
    out = df[df['orgName'].str.contains(organism)]
    dOut = {}
    for i, row in out.iterrows():
        dOut[row['orgName']] = row['orgCode']
    return dOut

def getOrgs():
    k = kegg.KEGG()
    Org = k.list('organism')
    spltOrg = Org.split("\n")
    dfAllOrgs = pd.DataFrame(columns=['Tnumber', 'orgCode', 'orgName', 'phylogeny'])
    for el in range(0,len(spltOrg)-1):
        sngItem = pd.DataFrame([spltOrg[el].split("\t")], columns=['Tnumber', 'orgCode', 'orgName', 'phylogeny'])
        dfAllOrgs = dfAllOrgs.append(sngItem, ignore_index=True)
    return dfAllOrgs

def generateOrganismSpecificRegex(keggGeneList):
    # Compose the regex according to characters present in KEGG gene identifiers.
    ## Scorrere lista di geni e per ognuno stabilire se ci sono solo caratteri maiuscoli, solo minuscoli o misti; verificare poi se ci sono numeri
    lower = [gene for gene in keggGeneList if gene.islower()]
    upper = [gene for gene in keggGeneList if gene.isupper()]
    mixed = [gene for gene in keggGeneList if gene.islower() and not gene.isupper()]

    if len(mixed) > 0:
        p1 = 'A-Za-z'
    elif len(mixed) == 0 and len(lower) == 0 and len(upper) > 0:
        p1 = 'A-Z'
    elif len(mixed) == 0 and len(lower) > 0 and len(upper) == 0:
        p1 = 'a-z'
    else:
        p1 = ''

    numeric = False
    for gene in keggGeneList:
        if any(c.isdigit() for c in gene) is True: #True if contains at least one numeric character
            numeric = True

    if numeric is True:
        p2 = '0-9'
    else:
        p2 = ''

    # Return all non-alphanumeric characters
    noAlphaNumeric = []
    for gene in keggGeneList:
        noAlphaNumeric += reModule.findall(r'\W+', gene)

    noAlphaNumeric = gL.unique(noAlphaNumeric)
    p3 = ''.join(noAlphaNumeric)

    # Construct the regex
    regexOrgSpecific = '([' + p1 + p2 + p3 + ']+)'
    return regexOrgSpecific
