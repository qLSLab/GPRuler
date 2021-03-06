# -*- coding: utf-8 -*-
import requests
import json

def getDefinitionFromGO(evidenceCode):
    requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/eco/terms/ECO%3A" + evidenceCode.split(':')[1]
    r = requests.get(requestURL, headers={ "Accept" : "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()
    responseBody = json.loads(r.text)
    if 'results' in responseBody and 'name' in responseBody['results'][0]:
        defEco = responseBody['results'][0]['name']
        if 'automatic' not in defEco:
            return 'manual'
        else:
            return 'automatic'

def getGOsearch(query):
    requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query=" + query.replace(' ', '%20') + "&limit=25&page=1"
    r = requests.get(requestURL, headers={ "Accept" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    responseBody = json.loads(r.text)
    id = ''
    if 'results' in responseBody:
        find = False
        i = 0
        while i < len(responseBody['results']) and find == False:
            if 'name' in responseBody['results'][i] and responseBody['results'][i]['name'].lower() == query:
                find = True
                id = responseBody['results'][i]['id']
            i += 1
    return id

def getGOAncestors(goTerm):
    requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO%3A" + goTerm.split(':')[1] + "/ancestors?relations=is_a%2Cpart_of%2Coccurs_in%2Cregulates"
    r = requests.get(requestURL, headers={ "Accept" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    responseBody = json.loads(r.text)
    if 'ancestors' in responseBody['results'][0]:
        lAncestors = responseBody['results'][0]['ancestors']
    else:
        lAncestors = []

    return lAncestors


def getGOName(goTerm):
    requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO%3A" + goTerm.split(':')[1]
    r = requests.get(requestURL, headers={ "Accept" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    nome = ''
    responseBody = json.loads(r.text)
    if 'name' in responseBody['results'][0]:
        nome = responseBody['results'][0]['name']
    return nome

def getCompartmentInfo(model):
    dCompartments = {}
    tree = ET.parse(model)
    root = tree.getroot()
    ## Definition of useful tags
    for element in root.getiterator():
        if element.tag.endswith('compartment'):
            compartmentTag = element.tag
            # print('compartmentTag\t', compartmentTag)
    for compartment in root.getiterator(compartmentTag):
        if 'id' in compartment.attrib:
            k = compartment.attrib['id']
        if 'name' in compartment.attrib:
            v = compartment.attrib['name']
        dCompartments[k] = v
    return dCompartments


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
