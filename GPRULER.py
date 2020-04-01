#!/usr/bin/env python
# coding: utf-8

# ###### Useful library
import GPRULERLib as gprL
import pandas as pd
from nltk.tokenize import word_tokenize
from bioservices import UniProt
import bioservices.kegg as kegg
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

timeStamp = gprL.getTimeStamp()

#############################
## Step 1. Load input data ##
#############################
workingDirs = gprL.setWorkingDirs()
MAINDIR = workingDirs[0]
INPDIR = workingDirs[1]
OUTDIR = workingDirs[2]

start = time.time()

analysis = input('Do you want to use test models (1) or your own data (2)? ')
if analysis == '1':
    testModel = input('Which test model you wanto to use? HMRcore (1), Recon 3 (2), Yeast 7 (3), Yeast 8 (4) ')
    if testModel == '1':
        ## HMRcore
        dfrxns2Genes = 'HMRcore_Rxns2Genes'
        dfKegg2UniprotGenesFileName = 'HMRcore_Kegg2UniprotGenes'
        regexOrgSpecific = r"([0-9]+)"
        model = 'HMRcore'
        organismCode = 'hsa'
    elif testModel == '2':
        ## Recon 3
        dfrxns2Genes = 'Recon3_Rxns2Genes'
        dfKegg2UniprotGenesFileName = 'Recon3_Kegg2UniprotGenes'
        regexOrgSpecific = r"([0-9]+)"
        model = 'Recon3'
        organismCode = 'hsa'
    elif testModel == '3':
        ## Yeast 7
        dfrxns2Genes = 'Yeast7_Rxns2Genes'
        dfKegg2UniprotGenesFileName = 'Yeast7_Kegg2UniprotGenes'
        regexOrgSpecific = r"([A-Z0-9-]+)"
        model = 'Yeast7'
        organismCode = 'sce'
    elif testModel == '4':
        ## Yeast 8
        dfrxns2Genes = 'Yeast8_Rxns2Genes'
        dfKegg2UniprotGenesFileName = 'Yeast8_Kegg2UniprotGenes'
        regexOrgSpecific = r"([A-Z0-9-]+)"
        model = 'Yeast8'
        organismCode = 'sce'
    dfRxnToGenes = pd.read_csv(os.path.join(INPDIR, dfrxns2Genes + '.csv'), sep='\t')
    dfRxnToGenes['Genes'] = dfRxnToGenes['Genes'].apply(literal_eval)
    dfKegg2UniprotId = pd.read_csv(os.path.join(INPDIR, dfKegg2UniprotGenesFileName + '.csv'), sep="\t", dtype = {'keggId': str})
elif analysis == '2':
    inputType = input('Which type of input data you have? Organism Name (1) or Reactions-Genes associations (2)')
    if inputType == '1':
        model = input('Which is the model name? ')
        organismChoice = input('Do you have the organism name (1) or the KEGG code (2) of the organism under investigation? ')
        if organismChoice == '1':
            organism = input('Insert the organism name: ')
            dfPutativeOrgs = gprL.putativeOrganisms(organism)
            while len(dfPutativeOrgs) == 0:
                organism = input('Organism not found! Insert the organism name: ')
                dfPutativeOrgs = gprL.putativeOrganisms(organism)
            print('\nPutative organisms are:\n')
            for key, value in dfPutativeOrgs.items():
                print(key, '\t', value)
            print('\n')
            organismCode = input('Type the correct KEGG code among the returned ones: ')
        elif organismChoice == '2':
            organismCode = input('Insert the KEGG organism code: ')
        ## 'list' returns the entire list of organism genes
        k = kegg.KEGG()
        keggGenes = k.list(organismCode)
        keggGenesSplt = keggGenes.strip().split("\n")
        lOrganismGenes = []
        for gene in keggGenesSplt:
            lOrganismGenes.append(gene.split('\t')[0].split(organismCode + ':')[1])

        dGene2RxnsList = {}
        cont = 1
        geneFile = open(os.path.join(INPDIR, model + '_GeneId2Rxns.csv'), mode='w')
        gprL.writeLineByLineToFile(geneFile, ['GeneId', 'Rxns'])
        rxnFile = open(os.path.join(INPDIR, model + '_RxnId2Equation.csv'), mode='w')
        gprL.writeLineByLineToFile(rxnFile, ['RxnId', 'Equation', 'Definition'])
        rxn2EcFile = open(os.path.join(INPDIR, model + '_RxnId2ECs.csv'), mode='w')
        gprL.writeLineByLineToFile(rxn2EcFile, ['RxnId', 'EC number'])
        dRxn2EcNumber = {}
        dCompleteRxns_equation = {} # to save all the reactions candidate for the model avoiding duplicated elements
        dCompleteRxns_definition = {} # to save all the reactions candidate for the model avoiding duplicated elements
        for gene in lOrganismGenes:
            ## get for each gene its BRITE information in order to select only metabolic genes
            dizInfo, info = gprL.getKeggInfo(organismCode + ':' + gene)
            lAssociatedRxns = []
            if 'BRITE' in dizInfo.keys() and 'Metabolism' in dizInfo['BRITE']:
                lOrths = list(dizInfo['ORTHOLOGY'].keys())
                for orth in list(dizInfo['ORTHOLOGY'].values()):
                    posOpen = orth.find('[EC:')
                    posClose = orth[posOpen:].find(']')
                    ecNumber = orth[posOpen:posClose+posOpen+1][1:-1]
                    lEcNumbers = list(gprL.extractRegexFromItem(ecNumber, r"([0-9.-]+)")[0])
                    if len(lEcNumbers) == 0:
                        for orth in lOrths:
                            dizInfoOrth, infoOrth = gprL.getKeggInfo('ko:' + orth)
                            if 'DBLINKS' in dizInfoOrth and 'RN' in dizInfoOrth['DBLINKS']:
                                lSelectedRxns = list(gprL.extractRegexFromItem(dizInfoOrth['DBLINKS']['RN'], r"([A-Z0-9]+)")[0])
                                for rxnName in lSelectedRxns:
                                    if rxnName not in dRxn2EcNumber:
                                        dRxn2EcNumber[rxnName] = [orth]
                                    else:
                                        dRxn2EcNumber[rxnName] += [orth]

                                    if rxnName.startswith('R') is True and rxnName not in dCompleteRxns_equation:
                                        dizInfoRxn, infoRxn = gprL.getKeggInfo('rn:' + rxnName)
                                        if 'EQUATION' in dizInfoRxn.keys():
                                            rxnEquation = dizInfoRxn['EQUATION']
                                            rxnDefinition = dizInfoRxn['DEFINITION']
                                        else:
                                            rxnEquation = ''
                                            rxnDefinition = ''

                                        dCompleteRxns_equation[rxnName] = rxnEquation
                                        dCompleteRxns_definition[rxnName] = rxnDefinition
                                        gprL.writeLineByLineToFile(rxnFile, [rxnName, rxnEquation, rxnDefinition])

                                    if rxnName not in lAssociatedRxns:
                                        lAssociatedRxns.append(rxnName)
                    else:
                        for ec in lEcNumbers:
                            dizInfoEc, infoEc = gprL.getKeggInfo('ec:' + ec)
                            try:
                                if 'ALL_REAC' not in dizInfoEc.keys() and 'REACTION' in dizInfoEc.keys():
                                    equazione = ' '. join(dizInfoEc['REACTION'])
                                    eqSplt = equazione.split('=')
                                    eqSplt_subs = eqSplt[0]
                                    eqSplt_prods = eqSplt[1]
                                    dSubs = gprL.dizReaProd(eqSplt_subs)
                                    dProds = gprL.dizReaProd(eqSplt_prods)

                                    lSubs = []
                                    for s in dizInfoEc['SUBSTRATE']:
                                        if 'CPD:' in s:
                                            pos = s.find('CPD:')
                                            substringa = s[pos + len('CPD:'):]
                                            lSubs += list(gprL.extractRegexFromItem(substringa, r"([A-Z0-9]+)")[0])
                                        else:
                                            lSubs += [s]

                                    i = 0
                                    dSubs_kegg = {}
                                    for k, v in dSubs.items():
                                        dSubs_kegg[lSubs[i]] = v
                                        i += 1

                                    lsubstratesParts = []
                                    for k, v in dSubs.items():
                                        lsubstratesParts.append(str(v) + ' ' + str(k))
                                    substrates = ' + '.join(lsubstratesParts)

                                    lsubstratesParts_kegg = []
                                    for k, v in dSubs_kegg.items():
                                        lsubstratesParts_kegg.append(str(v) + ' ' + str(k))
                                    substrates_kegg = ' + '.join(lsubstratesParts_kegg)

                                    lProds = []
                                    for p in dizInfoEc['PRODUCT']:
                                        if 'CPD:' in p:
                                            pos = p.find('CPD:')
                                            substringa = p[pos + len('CPD:'):]
                                            lProds += list(gprL.extractRegexFromItem(substringa, r"([A-Z0-9]+)")[0])
                                        else:
                                            lProds += [p]

                                    i = 0
                                    dProds_kegg = {}
                                    for k, v in dProds.items():
                                        dProds_kegg[lProds[i]] = v
                                        i += 1

                                    lproductsParts = []
                                    for k, v in dProds.items():
                                        lproductsParts.append(str(v) + ' ' + str(k))
                                    products = ' + '.join(lproductsParts)

                                    lproductsParts_kegg = []
                                    for k, v in dProds_kegg.items():
                                        lproductsParts_kegg.append(str(v) + ' ' + str(k))
                                    products_kegg = ' + '.join(lproductsParts_kegg)

                                    rxnName = 'R' + str(cont)
                                    rxnEquation = substrates + ' <=> ' + products
                                    rxnEquation_kegg = substrates_kegg + ' <=> ' + products_kegg

                                    if rxnName not in dRxn2EcNumber:
                                        dRxn2EcNumber[rxnName] = [ec]
                                    else:
                                        dRxn2EcNumber[rxnName] += [ec]

                                    if rxnEquation not in dCompleteRxns_equation.values():
                                        dCompleteRxns[rxnName] = rxnEquation_kegg
                                        gprL.writeLineByLineToFile(rxnFile, [rxnName, rxnEquation_kegg, rxnEquation])
                                        cont += 1
                                        lAssociatedRxns.append(rxnName)
                                    else:
                                        for k, v in dCompleteRxns_equation.items():
                                            if v == rxnEquation_kegg and k not in lAssociatedRxns:
                                                lAssociatedRxns.append(k)
                                elif 'ALL_REAC' in dizInfoEc.keys():
                                    for rxn in dizInfoEc['ALL_REAC']:
                                        if rxnName not in dRxn2EcNumber:
                                            dRxn2EcNumber[rxnName] = [ec]
                                        else:
                                            dRxn2EcNumber[rxnName] += [ec]

                                        lSelectedRxns = list(gprL.extractRegexFromItem(rxn, r"([A-Z0-9]+)")[0])
                                        for rxnName in lSelectedRxns:
                                            if rxnName.startswith('R') is True and rxnName not in dCompleteRxns:
                                                dizInfoRxn, infoRxn = gprL.getKeggInfo('rn:' + rxnName)
                                                if 'EQUATION' in dizInfoRxn.keys():
                                                    rxnEquation = dizInfoRxn['EQUATION']
                                                    rxnDefinition = dizInfoRxn['DEFINITION']
                                                else:
                                                    rxnEquation = ''
                                                    rxnDefinition = ''

                                                dCompleteRxns_equation[rxnName] = rxnEquation
                                                dCompleteRxns_definition[rxnName] = rxnDefinition
                                                gprL.writeLineByLineToFile(rxnFile, [rxnName, rxnEquation, rxnDefinition])

                                            if rxnName not in lAssociatedRxns:
                                                lAssociatedRxns.append(rxnName)
                            except:
                                for orth in lOrths:
                                    dizInfoOrth, infoOrth = gprL.getKeggInfo('ko:' + orth)
                                    if 'DBLINKS' in dizInfoOrth and 'RN' in dizInfoOrth['DBLINKS']:
                                        lSelectedRxns = list(gprL.extractRegexFromItem(dizInfoOrth['DBLINKS']['RN'], r"([A-Z0-9]+)")[0])
                                        for rxnName in lSelectedRxns:
                                            if rxnName not in dRxn2EcNumber:
                                                dRxn2EcNumber[rxnName] = [orth]
                                            else:
                                                dRxn2EcNumber[rxnName] += [orth]

                                            if rxnName.startswith('R') is True and rxnName not in dCompleteRxns_equation:
                                                dizInfoRxn, infoRxn = gprL.getKeggInfo('rn:' + rxnName)
                                                if 'EQUATION' in dizInfoRxn.keys():
                                                    rxnEquation = dizInfoRxn['EQUATION']
                                                    rxnDefinition = dizInfoRxn['DEFINITION']
                                                else:
                                                    rxnEquation = ''
                                                    rxnDefinition = ''

                                                dCompleteRxns_equation[rxnName] = rxnEquation
                                                dCompleteRxns_definition[rxnName] = rxnDefinition
                                                gprL.writeLineByLineToFile(rxnFile, [rxnName, rxnEquation, rxnDefinition])

                                            if rxnName not in lAssociatedRxns:
                                                lAssociatedRxns.append(rxnName)
            lAssociatedRxns = gprL.unique(lAssociatedRxns)
            dGene2RxnsList[gene] = lAssociatedRxns
            gprL.writeLineByLineToFile(geneFile, [gene, lAssociatedRxns])
        rxnFile.close()
        geneFile.close()

        for k, v in dRxn2EcNumber.items():
            gprL.writeLineByLineToFile(rxn2EcFile, [k, v])
        rxn2EcFile.close()

        flatdGene2RxnsList_values = [rxn for associatedRxns in dGene2RxnsList.values() for rxn in associatedRxns]
        flatdGene2RxnsList_values = gprL.unique(flatdGene2RxnsList_values)
        dfrxns2Genes = model + '_Rxns2Genes'
        outFile = open(os.path.join(INPDIR, dfrxns2Genes + '.csv'), mode='w')
        gprL.writeLineByLineToFile(outFile, ['Rxn', 'Genes'])
        for r in flatdGene2RxnsList_values:
            lAssociatedGenes = []
            for k,v in dGene2RxnsList.items():
                if r in v:
                    lAssociatedGenes.append(k)
            lAssociatedGenes = gprL.unique(lAssociatedGenes)
            gprL.writeLineByLineToFile(outFile, [r, lAssociatedGenes])
        outFile.close()

        # load the file reporting for each reaction the associated metabolic genes
        dfRxnToGenes = pd.read_csv(os.path.join(INPDIR, dfrxns2Genes + '.csv'), sep='\t')
        dfRxnToGenes['Genes'] = dfRxnToGenes['Genes'].apply(literal_eval)

        lOrganismGenesSet = []
        for gene in dfRxnToGenes['Genes']:
            lOrganismGenesSet += gene
        lOrganismGenesSet = gprL.unique(lOrganismGenesSet)

        ## Convert the list of all KEGG genes identifiers of the target organism to their corresponding UNIPROT identifiers
        gprL.kegg2UniprotGenesId(organismCode, model, INPDIR)

        ## Load the list of all KEGG genes identifiers of the target organism with their corresponding UNIPROT identifiers
        dfKegg2UniprotGenesFileName = model + '_Kegg2UniprotGenes'
        dfCompleteKegg2UniprotId = pd.read_csv(os.path.join(INPDIR, dfKegg2UniprotGenesFileName + '.csv'), sep="\t", dtype = {'keggId': str})

        ## Select from the dfCompleteKegg2UniprotId dataframe only the metabolic genes included in 'lOrganismGenesSet' list
        dfKegg2UniprotId = dfCompleteKegg2UniprotId[dfCompleteKegg2UniprotId.keggId.isin(lOrganismGenesSet)]
        dfKegg2UniprotId = dfKegg2UniprotId.reset_index(drop=True)

        regex = input('Do you want to manually insert your regex (1) or infer it from your data (2)? ')
        if regex == '1':
            regexOrgSpecific = input('Insert your regex: ')
        elif regex == '2':
            regexOrgSpecific = gprL.generateOrganismSpecificRegex(dfKegg2UniprotId['keggId'])

    elif inputType == '2':
        model = input('Which is the model name? ')
        # 'dfrxns2Genes' is the file name of the dataframe expliciting for each reaction the catalyst genes
        dfrxns2Genes = model + '_Rxns2Genes'
        dfRxnToGenes = pd.read_csv(os.path.join(INPDIR, dfrxns2Genes + '.csv'), sep='\t') # df in cui per ogni reazione viene esplicitata la lista di geni responsabili della sua catalisi
        dfRxnToGenes['Genes'] = dfRxnToGenes['Genes'].apply(literal_eval)
        organismChoice = input('Do you have the organism name (1) or the KEGG code (2) of the organism under investigation? ')
        if organismChoice == '1':
            organism = input('Insert the organism name: ')
            dfPutativeOrgs = gprL.putativeOrganisms(organism)
            while len(dfPutativeOrgs) == 0:
                organism = input('Organism not found! Insert the organism name: ')
                dfPutativeOrgs = gprL.putativeOrganisms(organism)
            print('\nPutative organisms are:\n')
            for key, value in dfPutativeOrgs.items():
                print(key, '\t', value)
            print('\n')
            organismCode = input('Type the correct KEGG code among the returned ones: ')
        elif organismChoice == '2':
            organismCode = input('Insert the KEGG organism code: ')

        lOrganismGenesSet = []
        for rxn in dfRxnToGenes['Genes']:
            lOrganismGenesSet += rxn
        lOrganismGenesSet = gprL.unique(lOrganismGenesSet)

        ## Convert the list of all KEGG genes identifiers of the target organism to their corresponding UNIPROT identifiers
        gprL.kegg2UniprotGenesId(organismCode, model, INPDIR)

        ## Load the list of all KEGG genes identifiers of the target organism with their corresponding UNIPROT identifiers
        dfKegg2UniprotGenesFileName = model + '_Kegg2UniprotGenes'
        dfCompleteKegg2UniprotId = pd.read_csv(os.path.join(INPDIR, dfKegg2UniprotGenesFileName + '.csv'), sep="\t", dtype = {'keggId': str})

        ## Select from the dfCompleteKegg2UniprotId dataframe only the metabolic genes included in 'lOrganismGenesSet' list
        dfKegg2UniprotId = dfCompleteKegg2UniprotId[dfCompleteKegg2UniprotId.keggId.isin(lOrganismGenesSet)]
        dfKegg2UniprotId = dfKegg2UniprotId.reset_index(drop=True)

        regex = input('Do you want to manually insert your regex (1) or infer it from your data (2)? ')
        if regex == '1':
            regexOrgSpecific = input('Insert your regex: ')
        elif regex == '2':
            regexOrgSpecific = gprL.generateOrganismSpecificRegex(dfKegg2UniprotId['keggId'])

#############################################################
## Step 2. Execute getUniprotAndComplexPortalData function ##
#############################################################
print('Get data from Uniprot and Complex Portal: DOING')
dfData = gprL.getUniprotAndComplexPortalData(dfKegg2UniprotId)
print('Get data from Uniprot and Complex Portal: DONE\n')

#############################################################
## Step 3. Execute textMiningFromUniprot function ##
#############################################################
print('Do text mining: DOING')
dfData = gprL.textMiningFromUniprot(dfData)
print('Do text mining: DONE\n')

#############################################################
## Step 4. Execute getDipData function ##
#############################################################
print('Get data from DIP: DOING')
dfData = gprL.getDipData(dfData)
print('Get data from DIP: DONE\n')

#############################################################
## Step 5. Execute getStringData function ##
#############################################################
print('Get data from STRING: DOING')
dfData = gprL.getStringData(dfData)
print('Get data from STRING: DONE\n')

#############################################################
## Step 6. Execute getKeggData function ##
#############################################################
print('Get data from KEGG: DOING')
dfData = gprL.getKeggData(dfData,organismCode)
print('Get data from KEGG: DONE\n')
dfData.to_csv(os.path.join(OUTDIR, model + '_GenesData_' + timeStamp + '.csv'), sep="\t", index = False)

#############################################################
## Step 7. Execute mergeData function ##
#############################################################
print('Merge data: DOING')
dfData = gprL.mergeData(dfData)
print('Merge data: DONE\n')
dfData.to_csv(os.path.join(OUTDIR, model + '_GenesRelationships_' + timeStamp + '.csv'), sep="\t", index = False)

#############################################################
## Step 8. Generate final GPR rules ##
#############################################################
print('Generation of rules: DOING')
dfGenesRelationships = pd.read_csv(os.path.join(OUTDIR, model + '_GenesRelationships_' + timeStamp + '.csv'), sep="\t")
dfCompleteGenesInfo = pd.read_csv(os.path.join(OUTDIR, model + '_GenesData_' + timeStamp + '.csv'), sep="\t")

dfGenesRelationships['gene'] = dfGenesRelationships['gene'].apply(literal_eval)
dfGenesRelationships['uniprotId'] = dfGenesRelationships['uniprotId'].apply(literal_eval)
dfGenesRelationships['AND'] = dfGenesRelationships['AND'].apply(literal_eval)
dfGenesRelationships['OR'] = dfGenesRelationships['OR'].apply(literal_eval)

dfCompleteGenesInfo['proteinNames'] = dfCompleteGenesInfo['proteinNames'].apply(literal_eval)
dfCompleteGenesInfo['geneNames'] = dfCompleteGenesInfo['geneNames'].apply(literal_eval)
dfCompleteGenesInfo['id_uniprot'] = dfCompleteGenesInfo['id_uniprot'].apply(literal_eval)
dfCompleteGenesInfo['subunitsFromName'] = dfCompleteGenesInfo['subunitsFromName'].apply(literal_eval)
dfCompleteGenesInfo['geneName_fromKEGG'] = dfCompleteGenesInfo['geneName_fromKEGG'].fillna('[]')
dfCompleteGenesInfo['geneName_fromKEGG'] = dfCompleteGenesInfo['geneName_fromKEGG'].apply(literal_eval)
dfCompleteGenesInfo['otherIsoforms'] = dfCompleteGenesInfo['otherIsoforms'].apply(literal_eval)
dfCompleteGenesInfo['isoform'] = dfCompleteGenesInfo['isoform'].apply(literal_eval)
dfCompleteGenesInfo['redundancy'] = dfCompleteGenesInfo['redundancy'].apply(literal_eval)


lRules = []
for i, row in dfRxnToGenes.iterrows():
    lCatalystGenes  = row['Genes']
    if len(lCatalystGenes) == 0:
        ruleLastVersion = ''
    elif len(lCatalystGenes) == 1:
        ruleLastVersion = str(lCatalystGenes[0])
    elif len(lCatalystGenes) > 1:
        combin = list(it.combinations(lCatalystGenes, 2))
        lCombin = []
        andFlag = False
        orFlag = False
        for c in combin:
            val0 = dfKegg2UniprotId.loc[dfKegg2UniprotId['keggId'] == c[0]]
            luniprot0 = []
            for i0, row0 in val0.iterrows():
                luniprot0.append(row0['uniprotId'])
            val1 = dfKegg2UniprotId.loc[dfKegg2UniprotId['keggId'] == c[1]]
            luniprot1 = []
            for i1, row1 in val1.iterrows():
                luniprot1.append(row1['uniprotId'])

            lcomplex0 = []
            lnames0 = []
            luniprot0_append = []
            lisoformsFromUniprot0 = []
            lIsoformsIndications0 = []
            lIsoformsFromKegg0 = []
            lRedundancy0 = []

            lcomplex1 = []
            lnames1 = []
            luniprot1_append = []
            lisoformsFromUniprot1 = []
            lIsoformsIndications1 = []
            lIsoformsFromKegg1 = []
            lRedundancy1 = []

            for u0 in luniprot0:
                resAllInfo = dfCompleteGenesInfo.loc[dfCompleteGenesInfo['uniprotId'] == u0]
                for i0, row0 in resAllInfo.iterrows():
                    lisoformsFromUniprot0 += row0['otherIsoforms'] # colonna "otherIsoforms" che contiene le isoforme da UniProt
                    lIsoformsFromKegg0 += row0['isoform'] # colonna "isoforms" che contiene le isoforme da KEGG
                    luniprot0_append += row0['id_uniprot'] # attacco anche la colonna con gli altri id_uniprot trovati
                    nomi0 = row0['proteinNames'] + row0['geneNames'] + row0['subunitsFromName'] + row0['geneName_fromKEGG'] # includo tutti i nomi trovati associati alla proteina
                    lnames0 += nomi0
                    lIsoformsIndications0.append(row0['isoformIndication'])
                    lRedundancy0 += row0['redundancy']

                comparison0 = [rule == luniprot0_append for rule in dfGenesRelationships['uniprotId']]
                res = dfGenesRelationships.loc[comparison0]
                for i0, r0 in res.iterrows():
                    lnames0 += r0['gene']
                    lcomplex0 += r0['AND']
            luniprot0 = gprL.unique(luniprot0 + luniprot0_append)
            lcomplex0 = gprL.unique(lcomplex0)
            lnames0 = gprL.unique(lnames0)
            lisoformsFromUniprot0 = gprL.unique(lisoformsFromUniprot0)
            lIsoformsIndications0 = gprL.unique(lIsoformsIndications0)
            lIsoformsFromKegg0 = gprL.unique(lIsoformsFromKegg0)

            for u1 in luniprot1:
                resAllInfo = dfCompleteGenesInfo[dfCompleteGenesInfo['uniprotId'] == u1]
                for i1, row1 in resAllInfo.iterrows():
                    lisoformsFromUniprot1 += row1['otherIsoforms'] # colonna "otherIsoforms" che contiene le isoforme da UniProt
                    lIsoformsFromKegg1 += row1['isoform'] # colonna "isoforms" che contiene le isoforme da KEGG
                    luniprot1_append += row1['id_uniprot'] # attacco anche la colonna con gli altri id_uniprot trovati
                    nomi1 = row1['proteinNames'] + row1['geneNames'] + row1['subunitsFromName'] + row1['geneName_fromKEGG'] # includo tutti i nomi trovati associati alla proteina
                    lnames1 += nomi1
                    lIsoformsIndications1.append(row1['isoformIndication'])
                    lRedundancy1 += row1['redundancy']


                comparison1 = [rule == luniprot1_append for rule in dfGenesRelationships['uniprotId']]
                res = dfGenesRelationships.loc[comparison1]
                for i1, r1 in res.iterrows():
                    lnames1 += r1['gene']
                    lcomplex1 += r1['AND']

            luniprot1 = gprL.unique(luniprot1 + luniprot1_append)
            lcomplex1 = gprL.unique(lcomplex1)
            lnames1 = gprL.unique(lnames1)
            lisoformsFromUniprot1 = gprL.unique(lisoformsFromUniprot1)
            lIsoformsIndications1 = gprL.unique(lIsoformsIndications1)
            lIsoformsFromKegg1 = gprL.unique(lIsoformsFromKegg1)

            if len(gprL.intersect(lisoformsFromUniprot0,luniprot1)) != 0 or len(gprL.intersect(lisoformsFromUniprot0, lnames1)) != 0 or len(gprL.intersect(lisoformsFromUniprot1, luniprot0)) != 0 or len(gprL.intersect(lisoformsFromUniprot1, lnames0)) != 0:
                lCombin.append('(' + str(c[0]) + ' or ' + str(c[1]) + ')')
                orFlag = True
            elif len(gprL.intersect(lRedundancy0,luniprot1)) != 0 or len(gprL.intersect(lRedundancy0, lnames1)) != 0 or len(gprL.intersect(lRedundancy1, luniprot0)) != 0 or len(gprL.intersect(lRedundancy1, lnames0)) != 0:
                lCombin.append('(' + str(c[0]) + ' or ' + str(c[1]) + ')')
                orFlag = True
            elif (any(el == True for el in lIsoformsIndications0) is True) and (any(el == True for el in lIsoformsIndications1) is True):
                orFlag = True
                lCombin.append('(' + str(c[0]) + ' or ' + str(c[1]) + ')')
            elif len(gprL.intersect(lcomplex0,luniprot1)) != 0 or len(gprL.intersect(lcomplex0, lnames1)) != 0 or len(gprL.intersect(lcomplex1, luniprot0)) != 0 or len(gprL.intersect(lcomplex1, lnames0)) != 0:
                lCombin.append('(' + str(c[0]) + ' and ' + str(c[1]) + ')')
                andFlag = True
            elif len(gprL.intersect(lIsoformsFromKegg0,luniprot1)) != 0 or len(gprL.intersect(lIsoformsFromKegg0, lnames1)) != 0 or len(gprL.intersect(lIsoformsFromKegg1, luniprot0)) != 0 or len(gprL.intersect(lIsoformsFromKegg1, lnames0)) != 0:
                lCombin.append('(' + str(c[0]) + ' or ' + str(c[1]) + ')')
                orFlag = True
            else:
                lCombin.append('(' + str(c[0]) + ' or ' + str(c[1]) + ')')
                orFlag = True

        if andFlag is True and orFlag is False:
            ruleLastVersion = ' and '.join(lCatalystGenes)
        elif andFlag is False and orFlag is True:
            ruleLastVersion = ' or '.join(lCatalystGenes)
        elif andFlag is True and orFlag is True:
            lParentheses_and = []
            lParentheses_or = []
            allGenes = []
            for l in lCombin:
                dfgenes = gprL.extractRegexFromItem(l, regexOrgSpecific)
                if ' and ' in l:
                    allGenes += list(dfgenes[0])
                    lParentheses_and.append(list(dfgenes[0]))
                elif ' or ' in l:
                    allGenes += list(dfgenes[0])
                    lParentheses_or.append(list(dfgenes[0]))

            lParentheses_and.sort() # lista di liste dove ogni lista interna contiene geni legati da stesso operatore nella stessa parentesi
            lParentheses_or.sort()
            allGenes = gprL.unique(allGenes)

            ruleParts = [lParentheses_and[0]] # lista di liste dove ogni lista interna contiene geni legati da stesso operatore nella stessa parentesi

            for i in range(1, len(lParentheses_and)):
                tmpI = lParentheses_and[i][:]
                coexist = False
                for p in ruleParts:  # calcoo inmtersezione con ogni elemento ddella lista e se positivo allora verifico che ci sia congruenza anche con altro elemento della lista
                    tmpP = p[:]
                    intersezione = gprL.intersect(tmpP, tmpI)
                    if intersezione != []:
                        diffPI = gprL.difference(tmpP, tmpI)
                        diffIP = gprL.difference(tmpI, tmpP)
                        prodotto = list(it.product(diffIP, diffPI))
                        present = True
                        for prod in prodotto:
                            prod_forward = [prod[0], prod[1]]
                            prod_backward = [prod[1], prod[0]]
                            if prod_forward not in lParentheses_and and prod_backward not in lParentheses_and:
                                present = False
                        lParentheses_and_reverse = [lParentheses_and[i][1], lParentheses_and[i][0]]
                        if present == True and lParentheses_and[i] not in ruleParts and lParentheses_and_reverse not in ruleParts:
                            p += diffIP
                            coexist = True
                if coexist == False:
                    ruleParts.append(lParentheses_and[i])

            if len(ruleParts) > 1:
                # Convert them to sets and use the set.intersection method, reducing over the list of sets
                # semplifico espressione logica isolando la parte comune
                intersected = list(set(ruleParts[0]).intersection(*ruleParts))

                ## raccolgo la parter di string in comune e le restanti parti le metto in OR
                ruleParts_remainingGenes = []
                for c in ruleParts:
                    ruleParts_remainingGenes.append(gprL.difference(c, intersected))

                ## poi verifico che le coppie in or dentro 'lParentheses_or' sono rispettate dentro commonReplaced e poi scrivo la regola finale
                ruleParts_remainingGenes_combs = list(it.combinations(range(0, len(ruleParts_remainingGenes)), 2))
                orCombinations = []
                for el in ruleParts_remainingGenes_combs:
                    if el[0] != el[1]:
                        orCombinations += list(it.product(ruleParts_remainingGenes[el[0]], ruleParts_remainingGenes[el[1]]))
                toRemove = []
                for o in orCombinations:
                    if o[0] == o[1]:
                        toRemove.append(o)
                orCombinations = [x for x in orCombinations if x not in toRemove]
                present = True
                for prod in orCombinations:
                    prod_forward = [prod[0], prod[1]]
                    prod_backward = [prod[1], prod[0]]
                    inside = False
                    for r in ruleParts:
                        if prod[0] in r and prod[1] in r:
                            inside = True
                    if prod_forward not in lParentheses_or and prod_backward not in lParentheses_or and inside == False: #and non coesistono in nessuna sottolista di ruleParts
                        present = False
                        ruleParts.append(prod_forward)

                if present == True:
                    lpartAnd = []
                    for r in ruleParts_remainingGenes:
                        if len(r) == 1:
                            lpartAnd.append(r[0])
                        else:
                            lpartAnd.append('(' + ' and '.join(r) + ')')
                    if intersected == []:
                        ruleFinale =  ' or '.join(lpartAnd)
                    else:
                        ruleFinale = ' and '.join(intersected) + ' and (' + ' or '.join(lpartAnd) + ')'
                andPairs = [item for elem in ruleParts for item in elem]
                andPairs = gprL.unique(andPairs)
                allMinusAnd = gprL.difference(allGenes, andPairs)
                combinToEvaluate = list(it.product(andPairs, allMinusAnd))
                isAllOr = True
                lRemainingOr = []
                for prod in combinToEvaluate:
                    prod_forward = [prod[0], prod[1]]
                    prod_backward = [prod[1], prod[0]]
                    inside2 = False
                    for r in ruleParts:
                        if prod[0] in r and prod[1] in r:
                            inside2 = True
                    if prod_forward not in lParentheses_or and prod_backward not in lParentheses_or and inside2 == False:
                        ruleParts.append(prod_forward)
                        isAllOr = False
                    else:
                        lRemainingOr.append(prod_forward)

                if isAllOr == True and allMinusAnd != []:
                    ruleFinale = '(' + ruleFinale + ') or ' + ' or '.join(allMinusAnd)

                ruleLastVersion = ruleFinale[:]
            elif len(ruleParts) == 1:
                andPairs = [item for elem in ruleParts for item in elem]
                andPairs = gprL.unique(andPairs)
                allMinusAnd = gprL.difference(allGenes, andPairs)
                combinToEvaluate = list(it.product(andPairs, allMinusAnd))
                present2 = True
                for prod in combinToEvaluate:
                    prod_forward = [prod[0], prod[1]]
                    prod_backward = [prod[1], prod[0]]
                    inside2 = False
                    for r in ruleParts:
                        if prod[0] in r and prod[1] in r:
                            inside2 = True
                    if prod_forward not in lParentheses_or and prod_backward not in lParentheses_or and inside2 == False:
                        present2 = False
                        ruleParts.append(prod_forward)

                tmpAndList = []
                for el in ruleParts:
                    tmpAndList.append('(' + ' and '.join(el) + ')')

                ruleFinale = ' or '.join(tmpAndList) + ' or ' + ' or '.join(allMinusAnd)
                ruleLastVersion = ruleFinale[:]

    lRules.append(ruleLastVersion)

dfRxnToGenes['GPR rule'] = lRules
dfRxnToGenes.to_csv(os.path.join(OUTDIR, model + '_gprRules_' + timeStamp + '.csv'), sep="\t", index = False)
print('Generation of rules: DONE\n')

end = time.time()
print('Elapsed time:\t', end-start, '\tseconds')
