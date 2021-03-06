# -*- coding: utf-8 -*-
import time
import pandas as pd
import bioservices.kegg as kegg
from Bio.KEGG.REST import *
import os
import sys
from bioservices import UniProt


from nltk.tokenize import word_tokenize
import xmltodict
import re as reModule
import requests
import itertools as it
from sympy import *
import json
from nltk.tokenize import RegexpTokenizer
from ast import literal_eval
import os
import sys
import nltk
nltk.download('punkt')
import time

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
        print('el:\t', el)
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
    print('lsubunits\n', lsubunits, '\n', len(lsubunits), '\n')
    print('df txt_subunit \n', df['txt_subunit'])
    u1mer = []
    for i in range(0, len(df)):
        print('word:\n', df['txt_subunit'][i])
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
        print('checkpoint1')

        sameEnzymeMembership = False # Returns True or False depending on whether indications of at least one of enzyme names within Function section of Uniprot is found
        pos = [s.start() for s in reModule.finditer('is one of', funzione.lower())]
        for p in pos:
            funzione_sub = funzione[p:]
            if any(el in funzione_sub for el in lnameEnzyme) is True:
                sameEnzymeMembership = True
        print('checkpoint2')

        redundancy = []
        pos = [s.start() for s in reModule.finditer('redundant with', funzione.lower())]
        for p in pos:
            funzione_sub = funzione[p:].split('.')[0]
            redundancy += [word for word in word_tokenize(funzione_sub, 'english', False) if word.isupper() and len(word) > 2]

        print('checkpoint3')

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

        print('checkpoint4')
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

        print('checkpoint5')

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
            print('frase\n', frase, '\n')
            for col in columns:
                orPosition = frase.find(' or ')
                print('orPosition\t', orPosition)
                frase_before = frase[:orPosition]
                frase_after = frase[orPosition-1+len(' or '):]
                print('frase_before\t', frase_before)
                print('frase_after\t', frase_after)
                print('find\t', getattr(row, col), '\n')
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
                                dfout0 = extractRegexFromItem(out0, r"([A-Za-z0-9-]+)")
                                lout0 = list(dfout0[0])
                                tmpIsoformList += lout0
                                lIsoforms += lout0
                                wordsToRemove.append(lout0[0] + ' or ')
                            if el1.startswith('or ') is True and col == 'geneNames':
                                out1 = el1.strip('or').strip().split()[0]
                                dfout1 = extractRegexFromItem(out1, r"([A-Za-z0-9-]+)")
                                lout1 = list(dfout1[0])
                                tmpIsoformList += lout1
                                lIsoforms += lout1
                                wordsToRemove.append(' or ' + lout1[0])
                    if len(wordsToRemove) != 0:
                        for w in wordsToRemove:
                            frase = reModule.sub(w, "", frase)
                        print('frase\t', frase)
                        orPosition = frase.find(' or ')
                        if orPosition != -1:
                            frase_before = frase[:orPosition]
                            frase_after = frase[orPosition-1+len(' or '):]
                        else:
                            frase_before = frase[:]
                            frase_after = frase[:]
                        print('orPosition\t', orPosition)
                        print('frase_before\t', frase_before)
                        print('frase_after\t', frase_after, '\n')
                    else:
                        orPosition = -1

            ## Once updated the sentence deprived of genes identified as isoforms, look for AND relationships
            print('Frase\t', frase)
            for w in wordsSelection:
                print('w\t', w)
                pos = [s.start() for s in reModule.finditer(w, frase.lower())]
                # print('pos\t', pos)
                for p in pos:
                    frase_sub = frase[p:]
                    print('frase_sub\t', frase_sub)
                    if 'isoforms of' not in frase_sub:
                        try:
                            frase_sub_splitted = reModule.split(w, frase_sub.lower())[1]
                            print('frase_sub_splitted\t', frase_sub_splitted)
                        except:
                            frase_sub_splitted = reModule.split(w[0].upper() + w[1:], frase_sub)[1]
                            print('frase_sub_splitted2\t', frase_sub_splitted)
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

def getDipData(df, id_prot='id_uniprot'):
    '''
    This function retrieves from DIP database protein-protein interactions and protein macromolecular complexes
    established by each queried metabolic gene.
    Input:
    - df: dataframe generated by textMiningFromUniprot function;
    - id_prot: column name of uniprot identifiers of genes. By default it is set equal to 'id_uniprot'.
    Output:
    - df: enriched input dataframe with information retrieved from DIP database.
    '''
    u = UniProt(verbose=False)
    lAllDip = []
    lAllDipComplex = []
    lAllDipBinary = []
    dizUniprotDip = {}
    diztmpDipId = {}
    for uniprotIds in df[id_prot]:
        currentDipId = []
        currentDipComplex = []
        currentDipBinary = []
        for unip in uniprotIds:
            tmpDipId = []
            if unip not in dizUniprotDip:
                res = u.search("%s" % unip, frmt="xml")
                dizUniprotDip[unip] = res
            else:
                res = dizUniprotDip[unip]
            if res != '':
                dRes = xmltodict.parse(res)
                try:
                    for ref in dRes['uniprot']['entry']['dbReference']:
                        if ref['@type'] == 'DIP':
                            foundEl = "".join(reModule.findall("[0-9]",ref['@id']))
                            currentDipId.append(foundEl)
                            tmpDipId.append(foundEl)
                except:
                    dipId = ''

            for singleDipId in tmpDipId:
                if singleDipId not in diztmpDipId:
                    url1 = 'https://dip.doe-mbi.ucla.edu/dip/Browse.cgi?PK=' + str(singleDipId) + '&MD=1'
                    response = requests.get(url1, verify=False)
                    test1 = response.content
                    dfs1 = pd.read_html(test1)

                    url0 = 'https://dip.doe-mbi.ucla.edu/dip/Browse.cgi?PK=' + str(singleDipId) + '&MD=0'
                    response = requests.get(url0, verify=False)
                    test0 = response.content
                    dfs0 = pd.read_html(test0)

                    diztmpDipId[singleDipId] = (dfs1, dfs0)
                else:
                    dfs1 = diztmpDipId[singleDipId][0]
                    dfs0 = diztmpDipId[singleDipId][1]

                for i in range(0, len(dfs1[8]['Cross Reference']['SWISSPROT'])):
                    if str(dfs1[8]['Cross Reference']['SWISSPROT'][i]) != "nan":
                        currentDipComplex.append(dfs1[8]['Cross Reference']['SWISSPROT'][i])


                for i in range(0, len(dfs0[8]['Cross Reference']['SWISSPROT'])):
                    if str(dfs0[8]['Cross Reference']['SWISSPROT'][i]) != "nan":
                        currentDipBinary.append(dfs0[8]['Cross Reference']['SWISSPROT'][i])

        lAllDip.append(currentDipId)
        lAllDipComplex.append(currentDipComplex)
        lAllDipBinary.append(currentDipBinary)
    df['DIPids'] = lAllDip
    df['complex_fromDIP'] = lAllDipComplex
    df['binary_fromDIP'] = lAllDipBinary
    return(df)

def getStringData(df, gene='id_uniprot'):
    '''
    This function retrieves from STRING database known and predicted protein-protein interactions
    established by each queried metabolic gene.
    Input:
    - df: dataframe generated by getDipData function;
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
            print('uniprotId\t', uniprotId)
            if uniprotId not in dizUniprotString:
                originalUniprotNames = row.proteinNames + row.geneNames
                url = "https://string-db.org/api/json/network?identifiers=" + uniprotId
                response = requests.get(url, verify=False)
                while response.status_code == 524:
                    response = requests.get(url, verify=False)
                    print('response.status_code\t', response.status_code)
                net = response.json()
                print('type net\t', type(net))
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
                        print('interactor\t', interactor)
                        url = 'https://string-db.org//api/json/enrichment?identifiers=' + original + '%0d' + interactor
                        response = requests.get(url, verify=False)
                        while response.status_code == 524:
                            response = requests.get(url, verify=False)
                            print('response.status_code\t', response.status_code)
                        diz = response.json()
                        dizlInteractors[interactor] = diz
                    else:
                        diz = dizlInteractors[interactor]
                    for stringElement in diz:
                        print('stringElement\t', stringElement)
                        if stringElement['category'] == 'Component' and ("complex" in str(stringElement['description']) or "chain" in str(stringElement['description'])) and original in stringElement['inputGenes']:
                            complesso = stringElement['description']
                            lStringInteractors.append(interactor)
        lStringInteractors = unique(lStringInteractors)
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

def mergeData(df, dip):
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
    lTotalGenes = unique(lTotalGenes)

    lTotalAndSubs = lTotalGenes + list(df['subunitsFromName'])

    dErroneousNames = {'SLCA7A7': 'SLC7A7', 'SLCA7A11': 'SLC7A11'}
    for r in df.itertuples():
        ## clean all information retrieved from the explored databases by removing those elements that are not included into the lTotalGenes list

        if dip == 'yes':
            lDipComplex = [x for x in list(r.complex_fromDIP) if x not in list(r.id_uniprot)]
            finallDipComplex = intersect(lDipComplex, lTotalGenes)
            lDipBinary = [x for x in list(r.binary_fromDIP) if x not in list(r.id_uniprot)]
            finallDipBinary = intersect(lDipBinary, lTotalGenes)
            finallDipBinary = unique(finallDipBinary)
        else:
            finallDipComplex = []
            finallDipBinary = []


        dfinallDipBinary_names = {}
        for f in finallDipBinary:
            out = df.loc[df.uniprotId == f]
            for o in out.itertuples():
                dfinallDipBinary_names[f] = unique([o.uniprotId] + o.proteinNames + o.geneNames + o.id_uniprot + o.geneName_fromKEGG)

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
        lAllIsoforms = unique(lAllIsoforms)

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
        subunit += ldfInteract + ldfSimilarity + difference(finalStringSubs, lAllIsoforms) + finalOtherSubunits
        subunit = unique(subunit)

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


def pathJoinCheck(dir2add, rootPath='.'):
    """Check the existence of a path and build it if necessary."""
    path = os.path.join(rootPath, dir2add)
    if not os.path.exists(path):
        os.makedirs(path)
    return path

def setWorkingDirs(mainDir=None):
    """Set the working directories inputData, outputData.

    If arguments are omitted current working directory replace them.

    Keyword arguments:
     dataDir: The directory to look for the input files.
     outDir:  The directory to write the output files.

    Return:
    (dataDir, outDir, mapDir, reportDir, modelDir, figureDir, dataframeDir, logDir, blastDir)
    """
    cwDir = os.getcwd()
    if mainDir is None:
        mainDir = cwDir
    else:
        mainDir = os.path.abspath(mainDir) + os.sep
    inputDir = pathJoinCheck('inputData', mainDir)
    outputDir = pathJoinCheck('outputData', mainDir)
    return mainDir, inputDir, outputDir

def getTimeStamp():
    return time.strftime('%Y%m%d%H%M%S', time.localtime())

def unique(a):
	""" return the list with duplicate elements removed """
	return list(set(a))

def intersect(a, b):
	""" return the intersection of two lists """
	return list(set(a) & set(b))

def difference(a, b):
    """ return the difference of two lists """
    return list(set(a) - set(b))

def writeLineByLineToFile(stream, dataToWrite):
    stream.write('\t'.join(str(s) for s in dataToWrite) + '\n')

def extractRegexFromItem(item, regex):
    sItem = pd.Series(data=item)
    if sItem.empty == False:
        dfItem = sItem.str.extractall(regex)
    else:
        dfItem = []
    return dfItem

def getOrgs():
    k = kegg.KEGG()
    Org = k.list('organism')
    spltOrg = Org.split("\n")
    dfAllOrgs = pd.DataFrame(columns=['Tnumber', 'orgCode', 'orgName', 'phylogeny'])
    for el in range(0,len(spltOrg)-1):
        sngItem = pd.DataFrame([spltOrg[el].split("\t")], columns=['Tnumber', 'orgCode', 'orgName', 'phylogeny'])
        dfAllOrgs = dfAllOrgs.append(sngItem, ignore_index=True)
    return dfAllOrgs

def putativeOrganisms(organism):
    df = getOrgs()
    out = df[df['orgName'].str.contains(organism)]
    dOut = {}
    for i, row in out.iterrows():
        dOut[row['orgName']] = row['orgCode']
    return dOut

def getKeggInfo(query):
    dizInfo = {}
    try:
        k = kegg.KEGG()
        info = k.get(query)
        dizInfo = k.parse(info)
    except dizInfo == 400 or dizInfo == 404:
            dizInfo = {}
    return dizInfo, info

def kegg2UniprotGenesId(organismCode, model, dirPath):
    ## Conversion of kegg identifiers to uniprot identifiers
    kegg2uniprot = kegg_conv('uniprot', organismCode)
    dOrgKegg2Uniprot = {}
    for el in kegg2uniprot.readlines():
        elSplt = el.strip().split('\t')
        if elSplt[0].split(':')[1] not in dOrgKegg2Uniprot.keys():
            dOrgKegg2Uniprot[elSplt[0].split(':')[1]] = [elSplt[1].split(':')[1]]
        else:
            dOrgKegg2Uniprot[elSplt[0].split(':')[1]] += [elSplt[1].split(':')[1]]

    outFile = open(os.path.join(dirPath, model + '_Kegg2UniprotGenes.csv'), mode='w')
    writeLineByLineToFile(outFile, ['keggId', 'uniprotId'])

    for k, v in dOrgKegg2Uniprot.items():
        for vv in v:
            writeLineByLineToFile(outFile, [k, vv])
    outFile.close()

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

    noAlphaNumeric = unique(noAlphaNumeric)
    p3 = ''.join(noAlphaNumeric)

    # Construct the regex
    regexOrgSpecific = '([' + p1 + p2 + p3 + ']+)'
    return regexOrgSpecific


def isCoeff(s):
    """Determine if a string splitted on the spaces the first element is the
    stoichiometric coefficient or not.
    Example: if string is "2 A" return True; if string is "A" it returns False"""
    answer = reModule.match('((\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)$', s)
    if answer is not None:
        # print('\t match ', answer.group(0)
        return True
    else:
        return False

def dizReaProd(s):
    termini = s.strip().split(' + ')
    # print('termini', termini
    diz = {}
    for termine in termini:
        coeffMetabo = termine.split()
        # print('coeffMetabo', coeffMetabo
        coeff = coeffMetabo[0]
        if isCoeff(coeff) is True:
            coeff = float(coeffMetabo[0])
            metabolita = ' '.join(coeffMetabo[1:])
        else:
            metabolita = ' '.join(coeffMetabo)
            coeff = 1.0
        # print('coeff', coeff
        # print('metabolita', metabolita
        diz[metabolita] = coeff
    return diz
