# -*- coding: utf-8 -*-
import os
import sys
import xmltodict
import pandas as pd
import genericLib as gL
import genesLib as genesL
from ast import literal_eval
from bioservices import UniProt
import RESTmoduleModified as RESTmod

# setting working dirs
workingDirs = gL.setWorkingDirs()
OUTDIR = workingDirs[2]

# setting input data
testModel = sys.argv[1]
if testModel == 'recon3':
    ## Recon 3
    dfGenesLoc = 'recon3D_genes2Compartments'
    dfrxns2Genes = 'recon3D_reactions_wGenes'
    orgCode = 'hsa'
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
    dfGenesLoc = 'y7_genes2Compartments'
    dfrxns2Genes = 'y7_reactions_wGenes'
    orgCode = 'sce'
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
    dfGenesLoc = 'y8_genes2Compartments'
    dfrxns2Genes = 'y8_reactions_wGenes'
    orgCode = 'sce'
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
    dfGenesLoc = 'hmrCore_genes2Compartments'
    dfrxns2Genes = 'hmrCore_reactions_wGenes'
    orgCode = 'hsa'
    lCompartmentsOrganization = {'cytoplasm': ['cytoplasm','preribosome', 'cytosol', 'endosome', 'cytoskeleton', 'mating projection'],
                                'extracellular': ['cell periphery', 'boundary', 'extracellular', 'plasma membrane', 'cell envelope'],
                                'mitochondrion': ['mitochondrion', 'mitochondrial matrix', 'mitochondrial membrane','mitochondrion inner membrane', 'mitochondrion inner membrane', 'mitochondrial inner membrane', 'mitochondrial intermembrane space', 'mitochondrial outer membrane']}

elif testModel == 'ownData':
    ## specify your input data
    dfGenesLoc = ''
    dfrxns2Genes = ''
    orgCode = ''
    lCompartmentsOrganization = {}


# Extract for the identified genes all the corresponding compartment information
dfModelRxns2Genes = pd.read_csv(os.path.join(OUTDIR, dfrxns2Genes + '.csv'), sep = '\t', dtype = {'Rxn': str, 'KeggId': str, 'GPR': str, 'Name': str, 'IsTransport': bool, 'IsExchange': bool, 'GPRrule': str})
dfModelRxns2Genes['lGenes'] = dfModelRxns2Genes['lGenes'].apply(literal_eval)

lAllGenes = []
for row in dfModelRxns2Genes.itertuples():
    if row.lGenes != []:
        for l in row.lGenes:
            lAllGenes += l

lAllGenes = gL.unique(lAllGenes)

# Convert KEGG gene identifiers to Uniprot identifiers
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

lPossibleCompartments = ['nucleus', 'nuclear lumen', 'nucleoplasm', 'nuclear envelope', 'mitochondrion inner membrane', 'mitochondrial inner membrane', 'mitochondrial membrane',
                        'mitochondrion', 'cytoplasm', 'vacuole', 'vacuolar membrane', 'endoplasmic reticulum', 'lipid droplet', 'plasma membrane', 'mitochondrial matrix',
                        'preribosome', 'cell periphery', 'mitochondrial intermembrane space', 'Golgi apparatus', 'peroxisome', 'endoplasmic reticulum membrane',
                        'cytoskeleton', 'mating projection', 'replication compartment', 'peroxisomal membrane', 'mitochondrial outer membrane', 'golgi membrane',
                        'golgi apparatus membrane', 'Golgi membrane', 'lysosomal membrane', 'lysosome']

dname2GO = {}
dAnc2Name = {}
outFile = open(os.path.join(OUTDIR, dfGenesLoc + '.csv'), mode='w')
gL.writeLineByLineToFile(outFile, ['Gene', 'lCompartments'], '\t')

u = UniProt(verbose=False)

for gene in lAllGenes:
    lDefinitiveCompartment = []
    if gene in dGene2Uniprot:
        lCompartments = []
        lUniprotIds = dGene2Uniprot[gene]
        for unip in lUniprotIds:
            uniprotSearch = u.search("%s" % unip, frmt="xml")
            if uniprotSearch != '':
                dUniprotSearch = xmltodict.parse(uniprotSearch)
                dEcoEvidence = {}
                # Construct dictionary of the Evidence codes
                if 'evidence' in dUniprotSearch['uniprot']['entry']:
                    if type(dUniprotSearch['uniprot']['entry']['evidence']) == list:
                        for evidence in dUniprotSearch['uniprot']['entry']['evidence']:
                            for k, v in evidence.items():
                                if k == '@key':
                                    code = v
                                if k == '@type':
                                    ecoCode = v
                            dEcoEvidence[str(code)] = ecoCode
                    else:
                        for k, v in dUniprotSearch['uniprot']['entry']['evidence'].items():
                            if k == '@key':
                                code = v
                            if k == '@type':
                                ecoCode = v
                        dEcoEvidence[str(code)] = ecoCode
                # Retrieve localization information from GO terms
                if 'dbReference' in dUniprotSearch['uniprot']['entry']:
                    for el in dUniprotSearch['uniprot']['entry']['dbReference']:
                        if el['@type'] == 'GO':
                            location = ''
                            for member in el['property']:
                                if member['@type'] == 'term' and member['@value'].startswith('C:'):
                                    location = member['@value'][2:]
                                if member['@type'] == 'evidence':
                                    evidence = member['@value'].split()

                            lAllEco = dEvidenceCodes['manual'] + dEvidenceCodes['automatic']
                            evidence2Search = gL.difference(evidence, lAllEco)
                            for eco in evidence2Search:
                                goDefinition = genesL.getDefinitionFromGO(eco)
                                if goDefinition == 'manual':
                                    dEvidenceCodes['manual'] += [eco]

                            if location != '' and any(eco in dEvidenceCodes['manual'] for eco in evidence) is True:
                                lCompartments.append(location.lower())

                # Retrieve localization information from Uniprot
                if 'comment' in dUniprotSearch['uniprot']['entry'] and type(dUniprotSearch['uniprot']['entry']['comment']) == list:
                    for el in dUniprotSearch['uniprot']['entry']['comment']:
                        if el['@type'] == 'subcellular location':
                            if 'subcellularLocation' in el.keys() and type(el['subcellularLocation']) == list:
                                for member in el['subcellularLocation']:
                                    if member == 'location':
                                        for localiz in el['subcellularLocation'][member]:
                                            if localiz == '@evidence':
                                                lEvidences = el['subcellularLocation'][member][localiz].split()
                                            elif localiz == '#text':
                                                location = el['subcellularLocation'][member][localiz]

                                        lEcos = gL.unique([dEcoEvidence[code] for code in lEvidences])
                                        lAllEco = dEvidenceCodes['manual'] + dEvidenceCodes['automatic']
                                        evidence2Search = gL.difference(lEcos, lAllEco)
                                        for eco in evidence2Search:
                                            goDefinition = genesL.getDefinitionFromGO(eco)
                                            if goDefinition == 'manual':
                                                dEvidenceCodes['manual'] += [eco]

                                        if any(eco in dEvidenceCodes['manual'] for eco in lEcos) is True:
                                            lCompartments.append(location.lower())
                                    elif 'location' in member.keys():
                                        if type(member['location']) == list:
                                            for el in member['location']:
                                                if '@evidence' in el and '#text' in el:
                                                    lEvidences =  el['@evidence'].split()
                                                    lEcos = gL.unique([dEcoEvidence[code] for code in lEvidences])
                                                    if any(eco in dEvidenceCodes['manual'] for eco in lEcos) is True:
                                                        lCompartments.append(el['#text'].lower())
                                                elif '@evidence' not in el and '#text' in el:
                                                    lCompartments.append(el['#text'].lower())
                                        elif type(member['location']) == str:
                                            lCompartments.append(member['location'].lower())
                                        else:
                                            if '@evidence' in member['location'] and '#text' in member['location']:
                                                lEvidences =  member['location']['@evidence'].split()
                                                lEcos = gL.unique([dEcoEvidence[code] for code in lEvidences])
                                                if any(eco in dEvidenceCodes['manual'] for eco in lEcos) is True:
                                                    lCompartments.append(member['location']['#text'].lower())
                                            elif '@evidence' not in member['location'] and '#text' in member['location']:
                                                lCompartments.append(member['location']['#text'].lower())

                            elif 'subcellularLocation' in el.keys() and type(el['subcellularLocation']) != list:
                                unipLoc = el['subcellularLocation']['location']
                                if type(unipLoc) == str:
                                    lCompartments.append(unipLoc.lower())
                                elif type(unipLoc) == list:
                                    for component in unipLoc:
                                        if type(component) == str:
                                            lCompartments.append(component.lower())
                                        else:
                                            if '#text' in component and '@evidence' in component:
                                                lEco = gL.unique([dEcoEvidence[code] for code in component['@evidence'].split()])
                                                lAllEco = dEvidenceCodes['manual'] + dEvidenceCodes['automatic']
                                                evidence2Search = gL.difference(lEco, lAllEco)
                                                for eco in evidence2Search:
                                                    goDefinition = genesL.getDefinitionFromGO(eco)
                                                    if goDefinition == 'manual':
                                                        dEvidenceCodes['manual'] += [eco]

                                                if any(eco in dEvidenceCodes['manual'] for eco in lEco) is True:
                                                    lCompartments.append(component['#text'].lower())
                                            elif '#text' in component and '@evidence' not in component:
                                                lCompartments.append(component['#text'].lower())
                                else:
                                    if '#text' in unipLoc and '@evidence' in unipLoc:
                                        lEco = gL.unique([dEcoEvidence[code] for code in unipLoc['@evidence'].split()])
                                        lAllEco = dEvidenceCodes['manual'] + dEvidenceCodes['automatic']
                                        evidence2Search = gL.difference(lEco, lAllEco)
                                        for eco in evidence2Search:
                                            goDefinition = genesL.getDefinitionFromGO(eco)
                                            if goDefinition == 'manual':
                                                dEvidenceCodes['manual'] += [eco]

                                        if any(eco in dEvidenceCodes['manual'] for eco in lEco) is True:
                                            lCompartments.append(unipLoc['#text'].lower())

                                    elif '#text' in unipLoc and '@evidence' not in unipLoc:
                                        lCompartments.append(unipLoc['#text'].lower())
        lCompartments_sistemati = []
        for comp in lCompartments:
            if comp in dCompartments:
                comp = dCompartments[comp]
            elif 'integral component of' in comp:
                comp = comp.replace('integral component of', '').strip()
                if comp in dCompartments:
                    comp = dCompartments[comp]
            lCompartments_sistemati.append(comp)
        lCompartments_sistemati = gL.unique(lCompartments_sistemati)

        knownComps = gL.intersect(lCompartments_sistemati, lPossibleCompartments)
        unknown = gL.difference(lCompartments_sistemati, lPossibleCompartments)

        if len(unknown) != 0:
            for unk in unknown:
                if unk not in dname2GO:
                    idGo = genesL.getGOsearch(unk)
                    if idGo != '':
                        lAncestors = genesL.getGOAncestors(idGo)
                        lnewComps = []
                        for anc in lAncestors:
                            if anc in dAnc2Name and dAnc2Name[anc] in lPossibleCompartments:
                                lnewComps.append(dAnc2Name[anc])
                            else:
                                name = genesL.getGOName(anc)
                                if name != '' and name in lPossibleCompartments:
                                    lnewComps.append(name)
                                    dAnc2Name[anc] = name
                        if len(lnewComps) == []:
                            dname2GO[unk] = ['cytoplasm']
                        else:
                            isCitopl = 'cytoplasm' in lnewComps or 'cytosol' in lnewComps
                            isMembrane = any('membrane' in c for c in lnewComps)
                            otherComp = [c for c in lnewComps if 'membrane' not in c and isCitopl is False]
                            if isMembrane is True:
                                dname2GO[unk] = [c for c in lnewComps if 'membrane' in c]
                            elif len(otherComp) != 0:
                                dname2GO[unk] = otherComp
                            else:
                                dname2GO[unk] = ['cytoplasm']
                        knownComps += dname2GO[unk]
                else:
                    knownComps += dname2GO[unk]
        knownComps = gL.unique(knownComps)
        lDefinitiveCompartment = []
        for defCompartment in knownComps:
            for k in lCompartmentsOrganization:
                if defCompartment in lCompartmentsOrganization[k]:
                    lDefinitiveCompartment.append(k)
        lDefinitiveCompartment = gL.unique(lDefinitiveCompartment)
    gL.writeLineByLineToFile(outFile, [gene, lDefinitiveCompartment], '\t')

outFile.close()
