# -*- coding: utf-8 -*-
import os
import sys
import cobra as cb
import pandas as pd
import genericLib as gL
import genesLib as genesL
from ast import literal_eval

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]
FIGUREDIR = workingDirs[4]
DFDIR = workingDirs[5]
LOGDIR = workingDirs[7]

# setting input data
testModel = sys.argv[1]
if testModel == 'recon3':
    ## Recon 3
    modelXml = 'Recon3D_301_20200923'
    dfGenes2Comp = 'recon3D_genes2Compartments'
    dfrxns2Genes = 'recon3D_reactions_wGenes'
    outputFileName = 'recon3D_genes2Compartments_wFilter'
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
    dfGenes2Comp = 'y7_genes2Compartments'
    dfrxns2Genes = 'y7_reactions_wGenes'
    outputFileName = 'y7_genes2Compartments_wFilter'
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
    dfGenes2Comp = 'y8_genes2Compartments'
    dfrxns2Genes = 'y8_reactions_wGenes'
    outputFileName = 'y8_genes2Compartments_wFilter'
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
    dfGenes2Comp = 'hmrCore_genes2Compartments'
    dfrxns2Genes = 'hmrCore_reactions_wGenes'
    outputFileName = 'hmrCore_genes2Compartments_wFilter'
    ## questa lista è da settare in base all'organizzazione del modello
    lCompartmentsOrganization = {'cytoplasm': ['cytoplasm','preribosome', 'cytosol', 'endosome', 'cytoskeleton', 'mating projection'],
                                'extracellular': ['cell periphery', 'boundary', 'extracellular', 'plasma membrane', 'cell envelope'],
                                'mitochondrion': ['mitochondrion', 'mitochondrial matrix', 'mitochondrial membrane','mitochondrion inner membrane', 'mitochondrion inner membrane', 'mitochondrial inner membrane', 'mitochondrial intermembrane space', 'mitochondrial outer membrane']}
elif testModel == 'ownData':
    ## specify your input data
    modelXml = ''
    dfGenes2Comp = ''
    dfrxns2Genes = ''
    outputFileName = ''
    lCompartmentsOrganization = {}

model = cb.io.read_sbml_model(os.path.join(RAWDIR, modelXml + '.xml'))

# Extract compartment information from the model
dModelCompartments = genesL.getCompartmentInfo(os.path.join(RAWDIR, modelXml + '.xml'))

dRxn2Compartments = {}
for r in model.reactions:
    lRxnComps = []
    for comp in r.get_compartments():
        for k in lCompartmentsOrganization:
            if dModelCompartments[comp].lower() in lCompartmentsOrganization[k]:
                lRxnComps.append(k)
    dRxn2Compartments[r.id] = lRxnComps

# Filter genes associated to each reaction according to the associated compartment
dfGenes2Compartment =  pd.read_csv(os.path.join(OUTDIR, dfGenes2Comp + '.csv'), sep = '\t', dtype = {'Gene': str})
dfGenes2Compartment['lCompartments'] = dfGenes2Compartment['lCompartments'].apply(literal_eval)
dGenes2Compartment = dfGenes2Compartment.set_index('Gene')['lCompartments'].to_dict()

dfModelRxns2Genes = pd.read_csv(os.path.join(OUTDIR, dfrxns2Genes + '.csv'), sep = '\t', dtype = {'Rxn': str, 'KeggId': str, 'GPR': str, 'Name': str, 'IsTransport': bool, 'IsExchange': bool, 'GPRrule': str})
dfModelRxns2Genes['lGenes'] = dfModelRxns2Genes['lGenes'].apply(literal_eval)

lGenesFiltered_all = []
for row in dfModelRxns2Genes.itertuples():
    lCompRxnModel = gL.unique(dRxn2Compartments[row.Rxn_conv])
    if len(row.lGenes) != 0 and len(lCompRxnModel) != 0:
        lAllGene_currentRxn = []
        for l in row.lGenes:
            lAllGene_currentRxn += l
        lAllGene_currentRxn = gL.unique(lAllGene_currentRxn)
        lGenes2Remove = []
        for g in lAllGene_currentRxn:
            if len(dGenes2Compartment[g]) != 0 and len(gL.intersect(gL.unique(dGenes2Compartment[g]), lCompRxnModel)) == 0:
                lGenes2Remove.append(g)
        lGenesFiltered = [[el for el in l if el not in lGenes2Remove] for l in row.lGenes]
        lGenesFiltered = [subL for subL in lGenesFiltered if subL != []]
        lGenesFiltered_all.append(lGenesFiltered)
    else:
        lGenesFiltered_all.append(row.lGenes)

dfModelRxns2Genes['lGenes_filtered'] = lGenesFiltered_all
dfModelRxns2Genes.to_csv(os.path.join(OUTDIR, outputFileName + '.csv'), sep = '\t', index = False)
