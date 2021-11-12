#!/usr/bin/env python
# coding: utf-8
import os
import sys
# from sympy import *
import pandas as pd
import itertools as it
import GPRULERLib as gprL
import genesLib as genesL
import genericLib as gL
from ast import literal_eval

# setting working dirs
workingDirs = gL.setWorkingDirs()
OUTDIR = workingDirs[2]

# setting input data
testModel = sys.argv[1]
if testModel == 'hmr':
    ## HMRcore
    regexOrgSpecific = r"([0-9]+)"
    model = 'HMRcore'
    organismCode = 'hsa'
elif testModel == 'recon':
    ## Recon 3
    regexOrgSpecific = r"([0-9]+)"
    model = 'Recon3'
    organismCode = 'hsa'
elif testModel == 'y7':
    ## Yeast 7
    regexOrgSpecific = r"([A-Z0-9-]+)"
    model = 'Yeast7'
    organismCode = 'sce'
elif testModel == 'y8':
    ## Yeast 8
    regexOrgSpecific = r"([A-Z0-9-]+)"
    model = 'Yeast8'
    organismCode = 'sce'
elif testModel == 'ownData':
    ## specify your input data
    model = ''
    regex = input('Do you want to manually insert your regex (1) or infer it from your data (2)? ')
    if regex == '1':
        regexOrgSpecific = input('Insert your regex: ')
    elif regex == '2':
        regexOrgSpecific = gprL.generateOrganismSpecificRegex(dfKegg2UniprotId['keggId'])

    organismChoice = input('Do you have the organism name (1) or the KEGG code (2) of the organism under investigation? ')
    if organismChoice == '1':
        organism = input('Insert the organism name: ')
        dfPutativeOrgs = genesL.putativeOrganisms(organism)
        while len(dfPutativeOrgs) == 0:
            organism = input('Organism not found! Insert the organism name: ')
            dfPutativeOrgs = genesL.putativeOrganisms(organism)
        print('\nPutative organisms are:\n')
        for key, value in dfPutativeOrgs.items():
            print(key, '\t', value)
        print('\n')
        organismCode = input('Type the correct KEGG code among the returned ones: ')
    elif organismChoice == '2':
        organismCode = input('Insert the KEGG organism code: ')

dfRxnToGenes = pd.read_csv(os.path.join(OUTDIR, model + '_Rxns2Genes.csv'), sep='\t')
dfRxnToGenes['Genes'] = dfRxnToGenes['Genes'].apply(literal_eval)
dfKegg2UniprotId = pd.read_csv(os.path.join(OUTDIR, model + '_Kegg2UniprotGenes.csv'), sep="\t", dtype = {'keggId': str})

lOrganismGenesSet = []
for gene in dfRxnToGenes['Genes']:
    lRxnGenes = [el for g in gene for el in g]
    lOrganismGenesSet += lRxnGenes
lOrganismGenesSet = gL.unique(lOrganismGenesSet)

## Select from the dfKegg2UniprotId dataframe only the metabolic genes included in 'lOrganismGenesSet' list
dfKegg2UniprotId = dfKegg2UniprotId[dfKegg2UniprotId.keggId.isin(lOrganismGenesSet)]
dfKegg2UniprotId = dfKegg2UniprotId.reset_index(drop=True)

#############################################################
# Execute getUniprotAndComplexPortalData function
#############################################################
print('Get data from Uniprot and Complex Portal: DOING')
dfData = gprL.getUniprotAndComplexPortalData(dfKegg2UniprotId)
print('Get data from Uniprot and Complex Portal: DONE\n')

#############################################################
# Execute textMiningFromUniprot function
#############################################################
print('Do text mining: DOING')
dfData = gprL.textMiningFromUniprot(dfData)
print('Do text mining: DONE\n')

#############################################################
# Execute getStringData function
#############################################################
print('Get data from STRING: DOING')
dfData = gprL.getStringData(dfData)
print('Get data from STRING: DONE\n')

#############################################################
# Execute getKeggData function
#############################################################
print('Get data from KEGG: DOING')
dfData = gprL.getKeggData(dfData,organismCode)
print('Get data from KEGG: DONE\n')
dfData.to_csv(os.path.join(OUTDIR, model + '_GenesData.csv'), sep="\t", index = False)

#############################################################
# Execute mergeData function
#############################################################
dfData = pd.read_csv(os.path.join(OUTDIR, model + '_GenesData.csv'), sep = '\t', dtype = str)
dfData['geneName_fromKEGG'] = dfData['geneName_fromKEGG'].fillna('[]')
dfData['geneName_fromKEGG'] = dfData['geneName_fromKEGG'].apply(literal_eval)
dfData['proteinNames'] = dfData['proteinNames'].fillna('[]')
dfData['proteinNames'] = dfData['proteinNames'].apply(literal_eval)
dfData['geneNames'] = dfData['geneNames'].fillna('[]')
dfData['geneNames'] = dfData['geneNames'].apply(literal_eval)
dfData['id_uniprot'] = dfData['id_uniprot'].fillna('[]')
dfData['id_uniprot'] = dfData['id_uniprot'].apply(literal_eval)
dfData['subunitsFromName'] = dfData['subunitsFromName'].fillna('[]')
dfData['subunitsFromName'] = dfData['subunitsFromName'].apply(literal_eval)
dfData['complexPortal_uniprotId'] = dfData['complexPortal_uniprotId'].fillna('[]')
dfData['complexPortal_uniprotId'] = dfData['complexPortal_uniprotId'].apply(literal_eval)
dfData['complexPortal_protName'] = dfData['complexPortal_protName'].fillna('[]')
dfData['complexPortal_protName'] = dfData['complexPortal_protName'].apply(literal_eval)
dfData['gene_from_structure'] = dfData['gene_from_structure'].fillna('[]')
dfData['gene_from_structure'] = dfData['gene_from_structure'].apply(literal_eval)
dfData['interact'] = dfData['interact'].fillna('[]')
dfData['interact'] = dfData['interact'].apply(literal_eval)
dfData['by_similarity'] = dfData['by_similarity'].fillna('[]')
dfData['by_similarity'] = dfData['by_similarity'].apply(literal_eval)
dfData['isoform'] = dfData['isoform'].fillna('[]')
dfData['isoform'] = dfData['isoform'].apply(literal_eval)
dfData['otherIsoforms'] = dfData['otherIsoforms'].fillna('[]')
dfData['otherIsoforms'] = dfData['otherIsoforms'].apply(literal_eval)
dfData['stringSubunits'] = dfData['stringSubunits'].fillna('[]')
dfData['stringSubunits'] = dfData['stringSubunits'].apply(literal_eval)
dfData['otherSubunits'] = dfData['otherSubunits'].fillna('[]')
dfData['otherSubunits'] = dfData['otherSubunits'].apply(literal_eval)

print('Merge data: DOING')
dfData = gprL.mergeData(dfData)
print('Merge data: DONE\n')
dfData.to_csv(os.path.join(OUTDIR, model + '_GenesRelationships.csv'), sep="\t", index = False)

#############################################################
## Step 8. Generate final GPR rules ##
#############################################################
print('Generation of rules: DOING')
dfGenesRelationships = pd.read_csv(os.path.join(OUTDIR, model + '_GenesRelationships.csv'), sep="\t")
dfCompleteGenesInfo = pd.read_csv(os.path.join(OUTDIR, model + '_GenesData.csv'), sep="\t")

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
for row in dfRxnToGenes.itertuples():
    lCatalystGenes = row.Genes
    if len(lCatalystGenes) == 0:
        ruleLastVersion = ''
    elif len(lCatalystGenes) == 1:
        if lCatalystGenes[0] != []:
            ruleLastVersion = '(' + ' and '.join(lCatalystGenes[0])+ ')'
    elif len(lCatalystGenes) > 1:
        lCatalystGenes_listing = []
        lCatalystGenes_original = []
        for internalList in lCatalystGenes.copy():
            if internalList != [] and internalList not in lCatalystGenes_original:
                lCatalystGenes_original.append(internalList)
                for el in internalList:
                    lCatalystGenes_listing.append(el)
        lCatalystGenes_listing = gL.unique(lCatalystGenes_listing)
        lCatalystGenes = lCatalystGenes_listing.copy()
        lCatalystGenes = gL.unique(lCatalystGenes)
        if len(lCatalystGenes) == 1:
            ruleLastVersion = str(lCatalystGenes[0])
        else:
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
                        lisoformsFromUniprot0 += row0['otherIsoforms'] # the "otherIsoforms" column includes the isoforms retrieved from Uniprot
                        lIsoformsFromKegg0 += row0['isoform'] # the "isoforms" column includes the isoforms retrieved from KEGG
                        luniprot0_append += row0['id_uniprot'] # also used the Uniprot identifiers included in id_uniprot column
                        nomi0 = row0['proteinNames'] + row0['geneNames'] + row0['subunitsFromName'] + row0['geneName_fromKEGG'] # use all the names retrieved for the protein
                        lnames0 += nomi0
                        lIsoformsIndications0.append(row0['isoformIndication'])
                        lRedundancy0 += row0['redundancy']

                    comparison0 = [len(gL.intersect(rule,luniprot0_append)) != 0 for rule in dfGenesRelationships['uniprotId']]
                    res = dfGenesRelationships.loc[comparison0]
                    for i0, r0 in res.iterrows():
                        lnames0 += r0['gene']
                        lcomplex0 += r0['AND']
                luniprot0 = gL.unique(luniprot0 + luniprot0_append)
                lcomplex0 = gL.unique(lcomplex0)
                lnames0 = gL.unique(lnames0)
                lisoformsFromUniprot0 = gL.unique(lisoformsFromUniprot0)
                lIsoformsIndications0 = gL.unique(lIsoformsIndications0)
                lIsoformsFromKegg0 = gL.unique(lIsoformsFromKegg0)

                for u1 in luniprot1:
                    resAllInfo = dfCompleteGenesInfo[dfCompleteGenesInfo['uniprotId'] == u1]
                    for i1, row1 in resAllInfo.iterrows():
                        lisoformsFromUniprot1 += row1['otherIsoforms']  # the "otherIsoforms" column includes the isoforms retrieved from Uniprot
                        lIsoformsFromKegg1 += row1['isoform'] # the "isoforms" column includes the isoforms retrieved from KEGG
                        luniprot1_append += row1['id_uniprot']  # also used the Uniprot identifiers included in id_uniprot column
                        nomi1 = row1['proteinNames'] + row1['geneNames'] + row1['subunitsFromName'] + row1['geneName_fromKEGG'] # use all the names retrieved for the protein
                        lnames1 += nomi1
                        lIsoformsIndications1.append(row1['isoformIndication'])
                        lRedundancy1 += row1['redundancy']

                    comparison1 = [len(gL.intersect(rule,luniprot1_append)) != 0 for rule in dfGenesRelationships['uniprotId']]
                    res = dfGenesRelationships.loc[comparison1]
                    for i1, r1 in res.iterrows():
                        lnames1 += r1['gene']
                        lcomplex1 += r1['AND']

                luniprot1 = gL.unique(luniprot1 + luniprot1_append)
                lcomplex1 = gL.unique(lcomplex1)
                lnames1 = gL.unique(lnames1)
                lisoformsFromUniprot1 = gL.unique(lisoformsFromUniprot1)
                lIsoformsIndications1 = gL.unique(lIsoformsIndications1)
                lIsoformsFromKegg1 = gL.unique(lIsoformsFromKegg1)

                if len(gL.intersect(lisoformsFromUniprot0,luniprot1)) != 0 or len(gL.intersect(lisoformsFromUniprot0, lnames1)) != 0 or len(gL.intersect(lisoformsFromUniprot1, luniprot0)) != 0 or len(gL.intersect(lisoformsFromUniprot1, lnames0)) != 0:
                    lCombin.append('(' + str(c[0]) + ' or ' + str(c[1]) + ')')
                    orFlag = True
                elif len(gL.intersect(lRedundancy0,luniprot1)) != 0 or len(gL.intersect(lRedundancy0, lnames1)) != 0 or len(gL.intersect(lRedundancy1, luniprot0)) != 0 or len(gL.intersect(lRedundancy1, lnames0)) != 0:
                    lCombin.append('(' + str(c[0]) + ' or ' + str(c[1]) + ')')
                    orFlag = True
                elif (any(el == True for el in lIsoformsIndications0) is True) and (any(el == True for el in lIsoformsIndications1) is True):
                    orFlag = True
                    lCombin.append('(' + str(c[0]) + ' or ' + str(c[1]) + ')')
                elif len(gL.intersect(lcomplex0,luniprot1)) != 0 or len(gL.intersect(lcomplex0, lnames1)) != 0 or len(gL.intersect(lcomplex1, luniprot0)) != 0 or len(gL.intersect(lcomplex1, lnames0)) != 0:
                    lCombin.append('(' + str(c[0]) + ' and ' + str(c[1]) + ')')
                    andFlag = True
                elif len(gL.intersect(lIsoformsFromKegg0,luniprot1)) != 0 or len(gL.intersect(lIsoformsFromKegg0, lnames1)) != 0 or len(gL.intersect(lIsoformsFromKegg1, luniprot0)) != 0 or len(gL.intersect(lIsoformsFromKegg1, lnames0)) != 0:
                    lCombin.append('(' + str(c[0]) + ' or ' + str(c[1]) + ')')
                    orFlag = True
                else:
                    lCombin.append('(' + str(c[0]) + ' or ' + str(c[1]) + ')')
                    orFlag = True

            for internalL in lCatalystGenes_original:
                if len(internalL) > 1:
                    lCombin.append('(' + ' and '.join(internalL)+ ')')
                    andFlag = True
            if andFlag is True and orFlag is False:
                ruleLastVersion = ' and '.join(lCatalystGenes)
            elif andFlag is False and orFlag is True:
                ruleLastVersion = ' or '.join(lCatalystGenes)
            elif andFlag is True and orFlag is True:
                lParentheses_and = []
                lParentheses_or = []
                allGenes = []
                for l in lCombin:
                    dfgenes = gL.extractRegexFromItem(l, regexOrgSpecific)
                    if ' and ' in l:
                        allGenes += list(dfgenes[0])
                        lParentheses_and.append(list(dfgenes[0]))
                    elif ' or ' in l:
                        allGenes += list(dfgenes[0])
                        lParentheses_or.append(list(dfgenes[0]))

                lParentheses_and.sort() # list of lists where each internal list includes genes linked by the same operator in the same parenthesis
                lParentheses_or.sort()
                allGenes = gL.unique(allGenes)
                ruleParts = [lParentheses_and[0]]

                for i in range(1, len(lParentheses_and)):
                    tmpI = lParentheses_and[i].copy()
                    coexist = False
                    # Compute the intersection with each element of the list and, if found, check that consistency exists also with the other element of the list
                    for p in ruleParts:
                        tmpP = p.copy()
                        intersezione = gL.intersect(tmpP, tmpI)
                        if intersezione != []:
                            diffPI = gL.difference(tmpP, tmpI)
                            diffIP = gL.difference(tmpI, tmpP)
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
                    # Simplify the logic expression isolating the common part and joining the remnant parts with OR operator
                    intersected = list(set(ruleParts[0]).intersection(*ruleParts))
                    ruleParts_remainingGenes = []
                    for c in ruleParts:
                        ruleParts_remainingGenes.append(gL.difference(c, intersected))

                    ruleParts_remainingGenes_combs = list(it.combinations(range(0, len(ruleParts_remainingGenes)), 2))
                    orCombinations = []
                    ruleParts_remainingGenes_combs = gL.unique(ruleParts_remainingGenes_combs)
                    for el in ruleParts_remainingGenes_combs:
                        if el[0] != el[1]:
                            orCombinations += list(it.product(ruleParts_remainingGenes[el[0]], ruleParts_remainingGenes[el[1]]))
                    toRemove = []
                    orCombinations = gL.unique(orCombinations)
                    for o in orCombinations:
                        if o[0] == o[1]:
                            toRemove.append(o)
                    orCombinations = [x for x in orCombinations if x not in toRemove]
                    present = True
                    orCombinations = gL.unique(orCombinations)
                    for prod in orCombinations:
                        prod_forward = [prod[0], prod[1]]
                        prod_backward = [prod[1], prod[0]]
                        inside = False
                        for r in ruleParts:
                            if prod[0] in r and prod[1] in r:
                                inside = True
                        if prod_forward not in lParentheses_or and prod_backward not in lParentheses_or and inside == False:
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
                    andPairs = gL.unique(andPairs)
                    allMinusAnd = gL.difference(allGenes, andPairs)
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
                    andPairs = gL.unique(andPairs)
                    allMinusAnd = gL.difference(allGenes, andPairs)
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
dfRxnToGenes.to_csv(os.path.join(OUTDIR, model + '_gprRules.csv'), sep="\t", index = False)
print('Generation of rules: DONE\n')
