# -*- coding: utf-8 -*-
import cobra as cb
import os
import sys
import time
import genericLib as gL
import pandas as pd
from ast import literal_eval
from bioservices import UniProt
import xmltodict
# import etreeLib as eL
from lxml import etree as ET
import requests
import json
import xmlLib as xL
import RESTmoduleModified as RESTmod

def getDefinitionFromGO(evidenceCode):
    requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/eco/terms/ECO%3A" + evidenceCode.split(':')[1]
    r = requests.get(requestURL, headers={ "Accept" : "application/json"})

    if not r.ok:
        r.raise_for_status()
        print('errore ECO')
        sys.exit()
    responseBody = json.loads(r.text)
    print(responseBody.keys())
    if 'results' in responseBody and 'name' in responseBody['results'][0]:
        defEco = responseBody['results'][0]['name']
        print('Def ECO:\t', defEco)
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
        # print('results\n', responseBody['results'])
        find = False
        i = 0
        while i < len(responseBody['results']) and find == False:
            if 'name' in responseBody['results'][i] and responseBody['results'][i]['name'].lower() == query:
                find = True
                id = responseBody['results'][i]['id']
            i += 1
    print('ID\t', id)
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




timeStamp = gL.getTimeStamp()

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]
FIGUREDIR = workingDirs[4]
DFDIR = workingDirs[5]
LOGDIR = workingDirs[7]

start = time.time()

testModel = sys.argv[1]
if testModel == 'recon3':
    ## Recon 3
    modelXml = 'Recon3D_301_20200923'
    dfrxnsInfo = 'recon3D_genes2Compartments'
    dfrxns2Genes = 'recon3D_reactions_wGenes_20210208210444'
    orgCode = 'hsa'
    ## questa lista è da settare in base all'organizzazione del modello
    lCompartmentsOrganization = {'cytoplasm': ['cytoplasm','preribosome', 'cytosol', 'endosome', 'cytoskeleton', 'mating projection'],
                                'extracellular': ['cell periphery', 'boundary', 'extracellular','plasma membrane', 'cell envelope'],
                                'golgi': ['golgi apparatus', 'golgi', 'golgi apparatus membrane', 'golgi membrane'],
                                'lysosome': ['lysosome', 'lysosomal membrane'],
                                'mitochondrion': ['mitochondrion', 'mitochondrial matrix', 'mitochondrial membrane','mitochondrion inner membrane', 'mitochondrion inner membrane', 'mitochondrial inner membrane', 'mitochondrial intermembrane space', 'mitochondrial outer membrane'],
                                'nucleus': ['nucleus', 'nuclear lumen', 'nucleoplasm', 'nuclear envelope', 'replication compartment'],
                                'endoplasmic reticulum': ['endoplasmic reticulum', 'endoplasmic reticulum membrane'],
                                'peroxisome': ['peroxisome', 'peroxisomal membrane']}
elif testModel == 'y7':
    ## Yeast 7
    modelXml = 'yeast_7.6_cobra'
    dfrxnsInfo = 'y7_genes2Compartments'
    dfrxns2Genes = 'y7_reactions_wGenes_20210208143051'
    orgCode = 'sce'
    ## questa lista è da settare in base all'organizzazione del modello
    lCompartmentsOrganization = {'extracellular': ['cell periphery', 'boundary', 'extracellular','plasma membrane', 'cell envelope'],
                                'cytoplasm': ['cytoplasm','preribosome', 'cytosol', 'endosome', 'cytoskeleton', 'mating projection'],
                                'endoplasmic reticulum': ['endoplasmic reticulum'], 'endoplasmic reticulum membrane': ['endoplasmic reticulum membrane'],
                                'golgi': ['golgi apparatus', 'golgi'], 'golgi membrane': ['golgi apparatus membrane', 'golgi membrane'],
                                'lipid droplet': ['lipid particle', 'lipid droplet'],
                                'mitochondrial membrane': ['mitochondrial membrane','mitochondrion inner membrane', 'mitochondrion inner membrane', 'mitochondrial inner membrane', 'mitochondrial intermembrane space', 'mitochondrial outer membrane'],
                                'mitochondrion': ['mitochondrion', 'mitochondrial matrix'],
                                'nucleus': ['nucleus', 'nuclear lumen', 'nucleoplasm', 'nuclear envelope', 'replication compartment'],
                                'peroxisome': ['peroxisome'], 'peroxisomal membrane': ['peroxisomal membrane'],
                                 'vacuole': ['vacuole'], 'vacuolar membrane': ['vacuolar membrane']}
elif testModel == 'y8':
    ## Yeast 8
    modelXml = 'yeast8'
    dfrxnsInfo = 'y8_genes2Compartments'
    dfrxns2Genes = 'y8_reactions_wGenes_20210208143053'
    orgCode = 'sce'
    ## questa lista è da settare in base all'organizzazione del modello
    lCompartmentsOrganization = {'cytoplasm': ['cytoplasm','preribosome', 'cytosol', 'endosome', 'cytoskeleton', 'mating projection'],
                                'cell envelope': ['plasma membrane', 'cell envelope'],
                                'extracellular': ['cell periphery', 'boundary', 'extracellular'],
                                'endoplasmic reticulum': ['endoplasmic reticulum'], 'endoplasmic reticulum membrane': ['endoplasmic reticulum membrane'],
                                'golgi': ['golgi apparatus', 'golgi'], 'golgi membrane': ['golgi apparatus membrane', 'golgi membrane'],
                                'lipid droplet': ['lipid particle', 'lipid droplet'],
                                'mitochondrial membrane': ['mitochondrial membrane','mitochondrion inner membrane', 'mitochondrion inner membrane', 'mitochondrial inner membrane', 'mitochondrial intermembrane space', 'mitochondrial outer membrane'],
                                'mitochondrion': ['mitochondrion', 'mitochondrial matrix'],
                                'nucleus': ['nucleus', 'nuclear lumen', 'nucleoplasm', 'nuclear envelope', 'replication compartment'],
                                'peroxisome': ['peroxisome', 'peroxisomal membrane'],
                                 'vacuole': ['vacuole'], 'vacuolar membrane': ['vacuolar membrane']}
elif testModel == 'hmr':
    ## HMRcore
    modelXml = 'HMRcore_20200328_wReconNames'
    dfrxnsInfo = 'hmrCore_genes2Compartments'
    dfrxns2Genes = 'hmrCore_reactions_wGenes_20210208094226'
    orgCode = 'hsa'
    ## questa lista è da settare in base all'organizzazione del modello
    lCompartmentsOrganization = {'cytoplasm': ['cytoplasm','preribosome', 'cytosol', 'endosome', 'cytoskeleton', 'mating projection'],
                                'extracellular': ['cell periphery', 'boundary', 'extracellular', 'plasma membrane', 'cell envelope'],
                                'mitochondrion': ['mitochondrion', 'mitochondrial matrix', 'mitochondrial membrane','mitochondrion inner membrane', 'mitochondrion inner membrane', 'mitochondrial inner membrane', 'mitochondrial intermembrane space', 'mitochondrial outer membrane']}


# model = cb.io.read_sbml_model(os.path.join(RAWDIR, modelXml + '.xml'))
#
# dModelCompartments = xL.getCompartmentInfo(os.path.join(RAWDIR, modelXml + '.xml'))
# print('dModelCompartments\n', dModelCompartments)
#
#
# dRxn2Compartments = {}
# for r in model.reactions:
#     # print(r.id, '\tCompartments: ', r.get_compartments())
#     lRxnComps = []
#     for comp in r.get_compartments():
#         lRxnComps.append(dModelCompartments[comp])
#     # print('lRxnComps\t', lRxnComps)
#     dRxn2Compartments[r.id] = lRxnComps

## extraggo per geni trovati le info sui compartimenti in cui sono annotati
dfModelRxns2Genes = pd.read_csv(os.path.join(OUTDIR, dfrxns2Genes + '.csv'), sep = '\t', dtype = {'Rxn': str, 'KeggId': str, 'GPR': str, 'Name': str, 'IsTransport': bool, 'IsExchange': bool, 'GPRrule': str})
dfModelRxns2Genes['lGenes'] = dfModelRxns2Genes['lGenes'].apply(literal_eval)
print('dfModelRxns2Genes\t', dfModelRxns2Genes.shape)

lAllGenes = []
for row in dfModelRxns2Genes.itertuples():
    if row.lGenes != []:
        for l in row.lGenes:
            lAllGenes += l

lAllGenes = gL.unique(lAllGenes)
print('Len lAllGenes\t', len(lAllGenes))

## conversione kegg a uniprot
dGene2Uniprot = {}
uniprot2Org = RESTmod.kegg_conv(orgCode, "uniprot").readlines()
for line in uniprot2Org:
    separo = line.strip().split('\t')
    ncbi = separo[1].split(':')[1]
    unip = separo[0].split(':')[1]
    if ncbi not in dGene2Uniprot:
        dGene2Uniprot[ncbi] = [unip]
    else:
        dGene2Uniprot[ncbi] += [unip]


dCompartments = {'cytosolic small ribosomal subunit': 'ribosome', 'Slx1-Slx4 complex': 'nucleus', 'cytosol': 'cytoplasm',
                'membrane': 'plasma membrane', 'U4/U6 x U5 tri-snRNP complex': 'nucleus', 'TRAPP complex': 'plasma membrane',
                'Cdc73/Paf1 complex': 'nucleoplasm', 'cell': 'cytoplasm', 'spliceosomal complex': 'nucleus',
                'mismatch repair complex': 'nucleus', 'small nucleolar ribonucleoprotein complex': 'cytoplasm',
                'MLL3/4 complex': 'nucleoplasm', 'Set1C/COMPASS complex': 'nucleoplasm', 'motile cilium': 'cilium',
                'TIM23 mitochondrial import inner membrane translocase complex': 'mitochondrial membrane',
                'dynein complex': 'cytoskeleton', 'proteasome complex': 'cytoplasm', 'mitochondrial membrane': 'mitochondrial membrane',
                'prefoldin complex': 'cytoplasm', 'DNA replication factor C complex': 'nucleus', 'chromosome': 'nucleus',
                'TORC1 complex': 'cytoplasm', 'TORC2 complex': 'cytoplasm', 'eukaryotic translation initiation factor 3 complex': 'cytoplasm',
                'vacuolar proton-transporting V-type ATPase, V1 domain': 'vacuolar membrane', 'proton-transporting V-type ATPase, V0 domain': 'vacuolar membrane',
                'proton-transporting V-type ATPase, V0 domain': 'plasma membrane', 'mitochondrial respiratory chain complex III': 'mitochondrial membrane',
                'vacuolar proton-transporting V-type ATPase, V0 domain': 'vacuolar membrane', 'proton-transporting V-type ATPase, V1 domain': 'vacuolar membrane',
                'proton-transporting V-type ATPase, V1 domain': 'plasma membrane', 'mitochondrial inner membrane': 'mitochondrial membrane',
                'cell membrane': 'plasma membrane', 'vacuole membrane': 'vacuolar membrane','mitochondrion inner membrane': 'mitochondrial membrane',
                'mitochondrion membrane': 'mitochondrial membrane', 'mitochondrion matrix': 'mitochondrion',
                'preribosome, small subunit precursor': 'ribosome', 'eukaryotic 43S preinitiation complex': 'cytoplasm', 'eukaryotic 48S preinitiation complex': 'cytoplasm',
                'small ribosomal subunit': 'ribosome', 'nuclear pore': 'nuclear envelope', 'SAGA complex': 'nucleoplasm', 'dna replication preinitiation complex': 'nucleoplasm',
                'proteasome regulatory particle, base subcomplex': 'cytoplasm', 'chromosome, telomeric region': 'nucleus',
                'U2 snRNP': 'nucleus', 'U5 snRNP': 'nucleus', 'U1 snRNP': 'nucleus', 'U4/U6 x U5 tri-snRNP complex, U4 snRNP': 'nucleus',
                'cohesin complex': 'nucleus', 'COP9 signalosome': 'nucleus', 'Arp2/3 protein complex': 'cytoskeleton',
                'nuclear pore outer ring': 'nuclear envelope', 'small-subunit processome': 'ribosome', 'COPII vesicle coat': 'Golgi apparatus',
                'protein phosphatase type 2A complex': 'cytoplasm', 'Mre11 complex': 'nucleus', 'GMP reductase complex': 'cytoplasm',
                'proton-transporting two-sector ATPase complex, catalytic domain': 'mitochondrion', 'pyruvate dehydrogenase complex': 'mitochondrial matrix',
                'extrinsic component of mitochondrial inner membrane': 'mitochondrial membrane', 'phosphopyruvate hydratase complex': 'cytoplasm',
                'vacuolar proton-transporting V-type ATPase complex': 'vacuolar membrane',
                'proton-transporting ATP synthase complex, catalytic core F(1)': 'plasma membrane',
                'proton-transporting ATP synthase complex, catalytic core F(1)': 'mitochondrial membrane',
                'glycerol-3-phosphate dehydrogenase complex': 'cytoplasm', 'glycine cleavage complex': 'cytoplasm',
                'proton-transporting ATP synthase complex, coupling factor F(o)': 'plasma membrane',
                'proton-transporting ATP synthase complex, coupling factor F(o)': 'mitochondrial membrane',
                'oxoglutarate dehydrogenase complex': 'mitochondrion', 'GPI-anchor transamidase complex': 'endoplasmic reticulum membrane',
                'respirasome': 'mitochondrial membrane',
                'glycosylphosphatidylinositol-N-acetylglucosaminyltransferase (GPI-GnT) complex': 'endoplasmic reticulum membrane',
                'integral component of endoplasmic reticulum membrane': 'endoplasmic reticulum membrane',
                'mcm complex':'cytoplasm', 'mcm core complex':'cytoplasm'}

dEvidenceCodes = {'manual': ['ECO:0000269', 'ECO:0000303', 'ECO:0000305', 'ECO:0000250', 'ECO:0000255','ECO:0000256', 'ECO:0000259',
                'ECO:0000312', 'ECO:0007005', 'ECO:0000315','ECO:0000244', 'ECO:0000314', 'ECO:0000304', 'ECO:0000318', 'ECO:0000501'],
                'automatic': ['ECO:0000313', 'ECO:0000213']}

# ECO:0000501 dichiarato in uniprot (https://www.uniprot.org/help/gene_ontology) che molte delle info qui presenti sono state curate manualmente.

lPossibleCompartments = ['nucleus', 'nuclear lumen', 'nucleoplasm', 'nuclear envelope', 'mitochondrion inner membrane', 'mitochondrial inner membrane', 'mitochondrial membrane',
                        'mitochondrion', 'cytoplasm', 'vacuole', 'vacuolar membrane', 'endoplasmic reticulum', 'lipid droplet', 'plasma membrane', 'mitochondrial matrix',
                        'preribosome', 'cell periphery', 'mitochondrial intermembrane space', 'Golgi apparatus', 'peroxisome', 'endoplasmic reticulum membrane',
                        'cytoskeleton', 'mating projection', 'replication compartment', 'peroxisomal membrane', 'mitochondrial outer membrane', 'golgi membrane',
                        'golgi apparatus membrane', 'Golgi membrane', 'lysosomal membrane', 'lysosome']



dname2GO = {}
dAnc2Name = {}
outFile = open(os.path.join(OUTDIR, dfrxnsInfo + '_' + timeStamp + '.csv'), mode='w')
gL.writeLineByLineToFile(outFile, ['Gene', 'lCompartments'], '\t')

u = UniProt(verbose=False)

for gene in lAllGenes:
    print('gene\t', gene)
    lDefinitiveCompartment = []
    if gene in dGene2Uniprot:
        lCompartments = []
        lUniprotIds = dGene2Uniprot[gene]
        for unip in lUniprotIds:
            print('unip\t', unip)
            uniprotSearch = u.search("%s" % unip, frmt="xml")
            if uniprotSearch != '':
                dUniprotSearch = xmltodict.parse(uniprotSearch)
                dEcoEvidence = {}
                ## costruisco dizionario degli evidence codes
                if 'evidence' in dUniprotSearch['uniprot']['entry']:
                    if type(dUniprotSearch['uniprot']['entry']['evidence']) == list:
                        for evidence in dUniprotSearch['uniprot']['entry']['evidence']:
                            # print('evidence\n', evidence)
                            for k, v in evidence.items():
                                if k == '@key':
                                    code = v
                                if k == '@type':
                                    ecoCode = v
                            dEcoEvidence[str(code)] = ecoCode
                    else:
                        # print('Evidence\n', dUniprotSearch['uniprot']['entry']['evidence'])
                        for k, v in dUniprotSearch['uniprot']['entry']['evidence'].items():
                            if k == '@key':
                                code = v
                            if k == '@type':
                                ecoCode = v
                        dEcoEvidence[str(code)] = ecoCode
                print('dEcoEvidence\t', dEcoEvidence, '\n')
                ## locations da GO terms
                if 'dbReference' in dUniprotSearch['uniprot']['entry']:
                    for el in dUniprotSearch['uniprot']['entry']['dbReference']:
                        if el['@type'] == 'GO':
                            location = ''
                            for member in el['property']:
                                if member['@type'] == 'term' and member['@value'].startswith('C:'):
                                    location = member['@value'][2:]
                                if member['@type'] == 'evidence':
                                    evidence = member['@value'].split()
                                    # print()

                            print('location1\t', location, '\tevidence GO term\t', evidence)
                            lAllEco = dEvidenceCodes['manual'] + dEvidenceCodes['automatic']
                            evidence2Search = gL.difference(evidence, lAllEco)
                            for eco in evidence2Search:
                                # print('eco\t', eco, '\t', eco.split(':')[1])

                                goDefinition = getDefinitionFromGO(eco)
                                if goDefinition == 'manual':
                                    # print('Manuale')
                                    dEvidenceCodes['manual'] += [eco]

                            if location != '' and any(eco in dEvidenceCodes['manual'] for eco in evidence) is True:
                                lCompartments.append(location.lower())
                                print('Dentro\t', '\n')
                            elif location != '' and all(eco not in dEvidenceCodes['manual'] for eco in evidence) is True and all(eco not in dEvidenceCodes['automatic'] for eco in evidence) is True:
                                print('Fuori\t', location, '\t', evidence, '\n')
                ## locations da Uniprot
                print('uniprot')
                if 'comment' in dUniprotSearch['uniprot']['entry'] and type(dUniprotSearch['uniprot']['entry']['comment']) == list:
                    for el in dUniprotSearch['uniprot']['entry']['comment']:
                        if el['@type'] == 'subcellular location':
                            if 'subcellularLocation' in el.keys() and type(el['subcellularLocation']) == list:
                                print('uniprotLoc\t',  el['subcellularLocation'], '\tTIPO\t', type(el['subcellularLocation']))
                                for member in el['subcellularLocation']:
                                    if member == 'location':
                                        for localiz in el['subcellularLocation'][member]:
                                            print('localizU\t', localiz)
                                            if localiz == '@evidence':
                                                lEvidences = el['subcellularLocation'][member][localiz].split()
                                            elif localiz == '#text':
                                                location = el['subcellularLocation'][member][localiz]

                                        print('Location2:\t', location, '\t', 'Evidence:\t', lEvidences, '\n')
                                        lEcos = gL.unique([dEcoEvidence[code] for code in lEvidences])

                                        lAllEco = dEvidenceCodes['manual'] + dEvidenceCodes['automatic']
                                        evidence2Search = gL.difference(lEcos, lAllEco)
                                        for eco in evidence2Search:
                                            # print('eco\t', eco, '\t', eco.split(':')[1])
                                            goDefinition = getDefinitionFromGO(eco)
                                            if goDefinition == 'manual':
                                                # print('Manuale')
                                                dEvidenceCodes['manual'] += [eco]

                                        if any(eco in dEvidenceCodes['manual'] for eco in lEcos) is True:
                                            lCompartments.append(location.lower())
                                            print('dentro2')
                                        elif all(eco not in dEvidenceCodes['manual'] for eco in lEcos) is True and all(eco not in dEvidenceCodes['automatic'] for eco in lEcos) is True:
                                            print('Fuori1\t', lEcos, '\n')
                                    elif 'location' in member.keys():
                                        print('MEMBER\t', member, '\nLOC\t', member['location'], '\t',  type(member['location']), '\n')
                                        if type(member['location']) == list:
                                            for el in member['location']:
                                                print('EL\t', el)
                                                if '@evidence' in el and '#text' in el:
                                                    lEvidences =  el['@evidence'].split()
                                                    print('EV\t',lEvidences)

                                                    lEcos = gL.unique([dEcoEvidence[code] for code in lEvidences])
                                                    if any(eco in dEvidenceCodes['manual'] for eco in lEcos) is True:
                                                        lCompartments.append(el['#text'].lower())
                                                        print('TEXT\t', el['#text'])
                                                elif '@evidence' not in el and '#text' in el:
                                                    lCompartments.append(el['#text'].lower())
                                                    print('TEXTNOEVIDENCE\t', el['#text'])
                                        elif type(member['location']) == str:
                                            print('AGGIUNTO')
                                            lCompartments.append(member['location'].lower())
                                        else:
                                            print('KEYS\t', member['location'].keys(), '\t', member['location'].values())

                                            if '@evidence' in member['location'] and '#text' in member['location']:
                                                lEvidences =  member['location']['@evidence'].split()
                                                print('EV\t',lEvidences)

                                                lEcos = gL.unique([dEcoEvidence[code] for code in lEvidences])
                                                if any(eco in dEvidenceCodes['manual'] for eco in lEcos) is True:
                                                    lCompartments.append(member['location']['#text'].lower())
                                                    print('TEXT\t', member['location']['#text'])
                                            elif '@evidence' not in member['location'] and '#text' in member['location']:
                                                lCompartments.append(member['location']['#text'].lower())
                                                print('TEXTNOEVIDENCE\t', member['location']['#text'])


                            elif 'subcellularLocation' in el.keys() and type(el['subcellularLocation']) != list:
                                print('caso2')
                                # print('DIVERSO\n', el['subcellularLocation'], '\tTIPO\t', type(el['subcellularLocation']), '\n')
                                unipLoc = el['subcellularLocation']['location']
                                if type(unipLoc) == str:
                                    lCompartments.append(unipLoc.lower())
                                    # print('unipLoc1\t',unipLoc, '\n') ## no info su evidence
                                elif type(unipLoc) == list:
                                    # print('unipLoc2\t', unipLoc, '\n')
                                    for component in unipLoc:
                                        # print('Type component\t', type(component))
                                        if type(component) == str:
                                            # print('component\t', component)
                                            lCompartments.append(component.lower())
                                        else:
                                            if '#text' in component and '@evidence' in component:
                                                # print('text--\t', component['#text'])
                                                lEco = gL.unique([dEcoEvidence[code] for code in component['@evidence'].split()])
                                                # print('lEco3\t', lEco, '\n')

                                                lAllEco = dEvidenceCodes['manual'] + dEvidenceCodes['automatic']
                                                evidence2Search = gL.difference(lEco, lAllEco)
                                                for eco in evidence2Search:
                                                    print('eco\t', eco, '\t', eco.split(':')[1])
                                                    goDefinition = getDefinitionFromGO(eco)
                                                    if goDefinition == 'manual':
                                                        # print('Manuale')
                                                        dEvidenceCodes['manual'] += [eco]

                                                if any(eco in dEvidenceCodes['manual'] for eco in lEco) is True:
                                                    lCompartments.append(component['#text'].lower())
                                                    print('dentro2')
                                                elif all(eco not in dEvidenceCodes['manual'] for eco in lEco) is True and all(eco not in dEvidenceCodes['automatic'] for eco in lEco) is True:
                                                    print('Fuori2\t', lEco, '\n')
                                                # print('componentcK\t',component.keys(), '\n')
                                                # print('componentV\t',component.values(), '\n')
                                            elif '#text' in component and '@evidence' not in component:
                                                lCompartments.append(component['#text'].lower())
                                else:
                                    # print('unipLoc3\t', unipLoc, '\n')
                                    if '#text' in unipLoc and '@evidence' in unipLoc:
                                        # print('text\t', unipLoc['#text'])
                                        lEco = gL.unique([dEcoEvidence[code] for code in unipLoc['@evidence'].split()])
                                        # print('lEco\t', lEco, '\n')
                                        lAllEco = dEvidenceCodes['manual'] + dEvidenceCodes['automatic']
                                        evidence2Search = gL.difference(lEco, lAllEco)
                                        for eco in evidence2Search:
                                            # print('eco\t', eco, '\t', eco.split(':')[1])
                                            goDefinition = getDefinitionFromGO(eco)
                                            if goDefinition == 'manual':
                                                # print('Manuale')
                                                dEvidenceCodes['manual'] += [eco]

                                        if any(eco in dEvidenceCodes['manual'] for eco in lEco) is True:
                                            lCompartments.append(unipLoc['#text'].lower())
                                            print('dentro3')
                                        elif all(eco not in dEvidenceCodes['manual'] for eco in lEco) is True and all(eco not in dEvidenceCodes['automatic'] for eco in lEco) is True:
                                            print('Fuori3\t', lEco, '\n')

                                    elif '#text' in unipLoc and '@evidence' not in unipLoc:
                                        lCompartments.append(unipLoc['#text'].lower())
                            else:
                                print('caso3')
        # print('lCompartmentsFINALE\t', lCompartments)
        lCompartments_sistemati = []
        for comp in lCompartments:
            if comp in dCompartments:
                comp = dCompartments[comp]
            elif 'integral component of' in comp:
                comp = comp.replace('integral component of', '').strip()
                if comp in dCompartments:
                    comp = dCompartments[comp]
            else:
                print('New compartment:\t', comp)
            lCompartments_sistemati.append(comp)
        lCompartments_sistemati = gL.unique(lCompartments_sistemati)
        print('Lista locations:\t', lCompartments_sistemati, '\n')

        knownComps = gL.intersect(lCompartments_sistemati, lPossibleCompartments)
        print('knownComps\t', knownComps)
        unknown = gL.difference(lCompartments_sistemati, lPossibleCompartments)
        print('unknown\t', unknown)

        if len(unknown) != 0:
            for unk in unknown:
                # print('unk\t', unk)
                if unk not in dname2GO:
                    idGo = getGOsearch(unk)
                    # print('idGo\t', idGo)
                    if idGo != '':
                        lAncestors = getGOAncestors(idGo)
                        # print('lAncestors\t', lAncestors)
                        lnewComps = []
                        for anc in lAncestors:
                            # print('anc\t', anc)
                            if anc in dAnc2Name and dAnc2Name[anc] in lPossibleCompartments:
                                lnewComps.append(dAnc2Name[anc])
                                # print('Nome anc:\t', dAnc2Name[anc])
                            else:
                                name = getGOName(anc)
                                # print('name\t', name)
                                if name != '' and name in lPossibleCompartments:
                                    lnewComps.append(name)
                                    dAnc2Name[anc] = name
                                    # print('Nome anc:\t', name)
                        # print('lnewComps\t', lnewComps)
                        if len(lnewComps) == []:
                            dname2GO[unk] = ['cytoplasm']
                        else:
                            isCitopl = 'cytoplasm' in lnewComps or 'cytosol' in lnewComps
                            # print('isCitopl\t', isCitopl)
                            isMembrane = any('membrane' in c for c in lnewComps)
                            # print('isMembrane\t', isMembrane)
                            otherComp = [c for c in lnewComps if 'membrane' not in c and isCitopl is False]
                            # print('otherComp\t', otherComp)
                            if isMembrane is True:
                                # print('isM tryu\t', [c for c in lnewComps if 'membrane' in c])
                                dname2GO[unk] = [c for c in lnewComps if 'membrane' in c]
                            elif len(otherComp) != 0:
                                # print('otherComp non 0\t', otherComp)
                                dname2GO[unk] = otherComp
                            else:
                                # print('solo cit')
                                dname2GO[unk] = ['cytoplasm']
                        knownComps += dname2GO[unk]
                else:
                    # print('known comp\t', dname2GO[unk])
                    knownComps += dname2GO[unk]
        knownComps = gL.unique(knownComps)
        # print('knownComps\t', knownComps, '\n')
        lDefinitiveCompartment = []
        for defCompartment in knownComps:
            # print('foundComp\t', defCompartment)
            for k in lCompartmentsOrganization:
                if defCompartment in lCompartmentsOrganization[k]:
                    # print('def compartment\t', k)
                    lDefinitiveCompartment.append(k)
        lDefinitiveCompartment = gL.unique(lDefinitiveCompartment)
        print('lDefinitiveCompartment\t', lDefinitiveCompartment)
    gL.writeLineByLineToFile(outFile, [gene, lDefinitiveCompartment], '\t')

outFile.close()

end = time.time()
print('Elapsed time:\t', end-start, '\tseconds')
