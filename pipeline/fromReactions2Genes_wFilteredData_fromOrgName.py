# -*- coding: utf-8 -*-
import os
import sys
import pandas as pd
import genericLib as gL
from ast import literal_eval

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]

# setting input data
dfGenes2Comp = ''
dfrxns2Genes = ''
modelName = ''
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


# Filter genes associated to each reaction according to the associated compartment
dfGenes2Compartment =  pd.read_csv(os.path.join(OUTDIR, dfGenes2Comp + '.csv'), sep = '\t', dtype = {'Gene': str})
dfGenes2Compartment['lCompartments'] = dfGenes2Compartment['lCompartments'].apply(literal_eval)
dGenes2Compartment = dfGenes2Compartment.set_index('Gene')['lCompartments'].to_dict()

dfModelRxns2Genes = pd.read_csv(os.path.join(OUTDIR, dfrxns2Genes + '.csv'), sep = '\t', dtype = {'RxnId': str})
dfModelRxns2Genes['Genes'] = dfModelRxns2Genes['Genes'].apply(literal_eval)

dRxn2GeneswLocation = {}
for row in dfModelRxns2Genes.itertuples():

    ## get all the genes of the current reaction
    lAllGene_currentRxn = []
    for l in row.Genes:
        lAllGene_currentRxn += l
    lAllGene_currentRxn = gL.unique(lAllGene_currentRxn)

    # reduce dGenes2Compartment to only consider genes of the current reaction
    dGenesCurrentRxn = {g: dGenes2Compartment[g] for g in lAllGene_currentRxn}

    # retrieved all the possible compartments of this reaction and create the same number of duplicates of the current reaction
    lRetrievedCompartments = []
    for g in dGenesCurrentRxn:
        lRetrievedCompartments += dGenesCurrentRxn[g]
    lRetrievedCompartments = gL.unique(lRetrievedCompartments)

    rxnSuffix = 1
    dRxn2AnnotatedCompartments = {}
    for putativeComp in lRetrievedCompartments:
        dRxn2AnnotatedCompartments[row.RxnId + '_' + str(rxnSuffix)] = putativeComp
        rxnSuffix += 1

    for duplicateRxn in dRxn2AnnotatedCompartments:
        lCompRxnModel = [dRxn2AnnotatedCompartments[duplicateRxn]]
        if len(row.Genes) != 0:
            lGenes2Remove = []
            for g in lAllGene_currentRxn:
                if len(dGenes2Compartment[g]) != 0 and len(gL.intersect(gL.unique(dGenes2Compartment[g]), lCompRxnModel)) == 0:
                    lGenes2Remove.append(g)
            lGenesFiltered = [[el for el in l if el not in lGenes2Remove] for l in row.Genes]
            lGenesFiltered = [subL for subL in lGenesFiltered if subL != []]
            dRxn2GeneswLocation[duplicateRxn] = lGenesFiltered
        else:
            dRxn2GeneswLocation[duplicateRxn] = row.Genes

    if len(lRetrievedCompartments) == 0:
        dRxn2GeneswLocation[row.RxnId] = row.Genes

dfRxn2GeneswLocation = pd.DataFrame(dRxn2GeneswLocation.items(), columns=['RxnId', 'Genes'])
dfRxn2GeneswLocation.to_csv(os.path.join(OUTDIR, modelName + '_Rxns2Genes.csv'), sep = '\t', index = False)
