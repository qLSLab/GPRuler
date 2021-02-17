from Bio.KEGG.REST import *
import keggLib as kL
import sys
import genericLib as gL
import os
import pandas as pd
import gprLib as gprL
import xmlLib as xL
from ast import literal_eval
from nltk.tokenize import RegexpTokenizer
import chebiLib as cL
import time
import re
from ast import literal_eval
import cobra as cb
import numpy as np

def extractChebiIds(l):
    lChebi = []
    for el in l:
        lChebi.append(el.split(';')[0][6:])
    lChebi = gL.unique(lChebi)
    return lChebi

def findTC(lR, lP, r):
    lFindR = [len(gL.intersect(l,r)) != 0 for l in lR]
    lFindP = [len(gL.intersect(l,r)) != 0 for l in lP]

    if any(lFindR) == True and any(lFindP) == True:
        return True
    else:
        return False

def findTC_transports(lR, r):
    lFindR = [len(gL.intersect(l,r)) != 0 for l in lR]

    if any(lFindR) == True:
        return True
    else:
        return False

timeStamp = gL.getTimeStamp()

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]
FIGUREDIR = workingDirs[4]
DFDIR = workingDirs[5]

start = time.time()

testModel = sys.argv[1]
if testModel == 'recon3':
    ## Recon 3
    modelXml = 'Recon3D_301_20200923'
    dfmetsInfo = 'recon3D_metabolites_20201218172131'
    dfmetsIds = 'recon3_mappingFuzzyAndClassic_20210119085007'
    dfrxnsInfo = 'recon3D_reactions'
elif testModel == 'y7':
    ## Yeast 7
    modelXml = 'yeast_7.6_cobra'
    dfmetsInfo = 'y7_metabolites_20201219101924'
    dfmetsIds = 'yeast7_mappingFuzzyAndClassic_20210113162642'
    dfrxnsInfo = 'y7_reactions'
elif testModel == 'y8':
    ## Yeast 8
    modelXml = 'yeast8'
    dfmetsInfo = 'y8_metabolites_20201219102046'
    dfmetsIds = 'yeast8_mappingFuzzyAndClassic_20210113162740'
    dfrxnsInfo = 'y8_reactions'
elif testModel == 'hmr':
    ## HMRcore
    modelXml = 'HMRcore_20200328_wReconNames'
    dfmetsInfo = 'hmrCore_metabolites_20210112092737'
    dfmetsIds = 'hmrCore_mappingFuzzyAndClassic_20210112165030'
    dfrxnsInfo = 'hmrCore_reactions'
# elif testModel == '4':
#     ## Ecoli
#     modelXml = 'iML1515-ROS'
#     # dfmetsInfo = 'eColi_metabolites_enriched'
#     dfmetsInfo = 'eColi_metabolites'


model = cb.io.read_sbml_model(os.path.join(RAWDIR, modelXml + '.xml'))
dfRxns = pd.read_csv(os.path.join(OUTDIR, dfrxnsInfo + '.csv'), sep = '\t', dtype = {'Rxn': str, 'KeggId': str, 'GPR': str, 'Name': str, 'IsTransport': bool, 'IsExchange': bool, 'GPRrule': str})
dfRxns['EC number'] = dfRxns['EC number'].apply(literal_eval)
dfRxns['trasportedMets'] = dfRxns['trasportedMets'].apply(literal_eval)

dfTc2Subs_original = pd.read_csv(os.path.join(RAWDIR, 'tcdb_tc2substrates.csv'), sep = '\t', dtype=str, names= ['TCnumber', 'Substrates'])
dfTc2Uniprot = pd.read_csv(os.path.join(RAWDIR, 'tcdb_tc2Uniprot.csv'), sep = '\t', dtype=str, names= ['UniprotId', 'TCnumber'])

dfTc2Subs_original['Substrates_splt'] = dfTc2Subs_original['Substrates'].str.split('|')
dfTc2Subs_original['Substrates_converted'] = dfTc2Subs_original.apply(lambda row: extractChebiIds(row['Substrates_splt']), axis=1)
print(dfTc2Subs_original.iloc[0:5]['Substrates_converted'], '\n')

dfMetsFromModel = pd.read_csv(os.path.join(OUTDIR, dfmetsInfo + '_enriched.csv'), sep = '\t', dtype=str) ## 1
dfMetsFromModel['lIdentifiers'] = dfMetsFromModel['lIdentifiers'].apply(literal_eval)

if testModel == 'recon3':
    dfMetsFromModel['Id_converted'] = dfMetsFromModel['Id'].str.strip('M_')
    dfMetsFromModel['Id_converted'] = dfMetsFromModel['Id_converted'].str.replace('__91__', '[')
    dfMetsFromModel['Id_converted'] = dfMetsFromModel['Id_converted'].str.replace('__93', ']')
    dfMetsFromModel['Id_converted2'] = dfMetsFromModel['Id'].str.strip('M_')
    dfMetsFromModel['Id_converted2'] = dfMetsFromModel['Id_converted2'] + '__'
    # print(dfMetsFromModel.head())

## Step 1. arricchire dfMetsFromModel con le info da dfMetsIdentifiers
## df metacyc
dfmetacyc = pd.read_csv(os.path.join(OUTDIR, 'metacyc_compounds_20201216152513.csv'), sep = '\t', dtype=str)

## df kegg compound - glycan
dfkeggC = pd.read_csv(os.path.join(OUTDIR, 'kegg_compounds_20201216150802.csv'), sep = '\t', dtype=str)
dfkeggC['ChebiId'] = dfkeggC['ChebiId'].apply(literal_eval)
dfkeggG = pd.read_csv(os.path.join(OUTDIR, 'kegg_glycans_20201216150802.csv'), sep = '\t', dtype=str)
dfkeggG['ChebiId'] = dfkeggG['ChebiId'].apply(literal_eval)

dfMetsFromModel = pd.read_csv(os.path.join(OUTDIR, dfmetsInfo + '.csv'), sep = '\t', dtype=str)
dfMetsIdentifiers = pd.read_csv(os.path.join(OUTDIR, dfmetsIds + '.tsv'), sep = '\t')
dfMetsIdentifiers['Name'] = dfMetsIdentifiers['Name'].str.strip()
dfMetsIdentifiers['Identifiers'] = dfMetsIdentifiers['Identifiers'].apply(literal_eval)
dizMetsIdentifiers = dfMetsIdentifiers.set_index('Name')['Identifiers'].to_dict()
# print(dizMetsIdentifiers)

for k in dizMetsIdentifiers:
    lCompleteIds = []
    # print(k, '\t', dizMetsIdentifiers[k])
    lCompleteIds += dizMetsIdentifiers[k]
    metacDb = dfmetacyc[dfmetacyc['Id'].isin(dizMetsIdentifiers[k])]
    if metacDb.empty is False:
        lCompleteIds += list(metacDb['ChebiId'].dropna())
        lCompleteIds += list(metacDb['KeggId'].dropna())
    chebiDb = dfmetacyc[dfmetacyc['ChebiId'].isin(dizMetsIdentifiers[k])]
    if chebiDb.empty is False:
        lCompleteIds += list(chebiDb['Id'].dropna())
        lCompleteIds += list(chebiDb['KeggId'].dropna())
    keggDb = dfmetacyc[dfmetacyc['KeggId'].isin(dizMetsIdentifiers[k])]
    if keggDb.empty is False:
        lCompleteIds += list(keggDb['Id'].dropna())
        lCompleteIds += list(keggDb['ChebiId'].dropna())

    keggDbC = dfkeggC[dfkeggC['Id'].isin(dizMetsIdentifiers[k])]
    if keggDbC.empty is False:
        lCompleteIds += [el for l in list(keggDbC['ChebiId'].dropna()) for el in l]

    dfkeggC_filter =dfkeggC[pd.DataFrame(dfkeggC.ChebiId.tolist()).isin(dizMetsIdentifiers[k]).any(1).values]
    if dfkeggC_filter.empty is False:
        lCompleteIds += list(dfkeggC_filter['Id'].dropna())

    keggDbG = dfkeggG[dfkeggG['Id'].isin(dizMetsIdentifiers[k])]
    if keggDbG.empty is False:
        lCompleteIds += [el for l in list(keggDbG['ChebiId'].dropna()) for el in l]

    keggDbG_filter =keggDbG[pd.DataFrame(keggDbG.ChebiId.tolist()).isin(dizMetsIdentifiers[k]).any(1).values]
    if keggDbG_filter.empty is False:
        lCompleteIds += list(keggDbG_filter['Id'].dropna())
    lCompleteIds = gL.unique(lCompleteIds)
    dizMetsIdentifiers[k] = lCompleteIds

lOriginalEmpty_newEmpty = 0
lOriginalEmpty_newFound = 0
lOriginalYes_newFound = 0
lOriginalYes_newNotFound = 0
lOriginalYes_newEmpty = 0
total = 0

if testModel == 'y7' or testModel == 'y8':
    dfMetsFromModel['Name'] = dfMetsFromModel.Name.str.replace("\[(\w*\s*)+\]$", "")

lIdentifiers = []
for row in dfMetsFromModel.itertuples():

    lInferredIds = dizMetsIdentifiers[row.Name.strip()]
    lIdentifiers.append(lInferredIds)

    # # print(row.Name)
    # if testModel == '2' or testModel == '3':
    #     lPositions = [p.start() for p in re.finditer('\[', row.Name)]
    #     pos = lPositions[-1]
    #     met = row.Name[:pos].strip()
    #     lInferredIds = dizMetsIdentifiers[met]
    # else:
    #     lInferredIds = dizMetsIdentifiers[row.Name]
    # print('lInferredIds\t', lInferredIds)

    lIdsWithinModel = []
    if pd.isna(row.KeggId) is False:
        lIdsWithinModel += [row.KeggId]
    if pd.isna(row.ChebiId) is False:
        lIdsWithinModel += [row.ChebiId]
    # print('lIdsWithinModel\t', lIdsWithinModel)

    total += 1
    if lIdsWithinModel == [] and lInferredIds == []:
        lOriginalEmpty_newEmpty += 1
    elif lIdsWithinModel == [] and lInferredIds != []:
        lOriginalEmpty_newFound += 1
    elif lIdsWithinModel != [] and lInferredIds != []:
        if len(gL.intersect(lInferredIds, lIdsWithinModel)) != 0:
            lOriginalYes_newFound += 1
        else:
            lOriginalYes_newNotFound += 1
    elif lIdsWithinModel != [] and lInferredIds == []:
        lOriginalYes_newEmpty += 1

print('total\t', total)
print('lOriginalEmpty_newEmpty\t', lOriginalEmpty_newEmpty)
print('lOriginalEmpty_newFound\t', lOriginalEmpty_newFound)
print('lOriginalYes_newFound\t', lOriginalYes_newFound)
print('lOriginalYes_newNotFound\t', lOriginalYes_newNotFound)
print('lOriginalYes_newEmpty\t', lOriginalYes_newEmpty)


dfMetsFromModel['lIdentifiers'] = lIdentifiers
# dfMetsFromModel.to_csv(os.path.join(OUTDIR, dfmetsInfo + '_enriched.csv'), sep = '\t', index = False)

dfChebiFormula = pd.read_csv(os.path.join(RAWDIR, 'chebi_chemical_data_20201216.tsv'), sep = '\t', dtype=str)

dfChebiCompounds_exploded = pd.read_csv(os.path.join(OUTDIR, 'chebi_compounds_20201216153117_exploded.csv'), sep = '\t', dtype=str) ## 1
dfChebiCompounds_exploded['ParentalChebiIds'] = dfChebiCompounds_exploded['ParentalChebiIds'].apply(literal_eval)
dfChebiCompounds_exploded['AllChebiIds'] = dfChebiCompounds_exploded['AllChebiIds'].apply(literal_eval)
dfChebiCompounds_exploded['ParentalKeggC'] = dfChebiCompounds_exploded['ParentalKeggC'].apply(literal_eval)
dfChebiCompounds_exploded['ParentalKeggG'] = dfChebiCompounds_exploded['ParentalKeggG'].apply(literal_eval)
dfChebiCompounds_exploded['ParentalMetacyc'] = dfChebiCompounds_exploded['ParentalMetacyc'].apply(literal_eval)
dfChebiCompounds_exploded['AllKeggC'] = dfChebiCompounds_exploded['AllKeggC'].apply(literal_eval)
dfChebiCompounds_exploded['AllKeggG'] = dfChebiCompounds_exploded['AllKeggG'].apply(literal_eval)
dfChebiCompounds_exploded['AllMetacyc'] = dfChebiCompounds_exploded['AllMetacyc'].apply(literal_eval)

dfChebiCompounds = pd.read_csv(os.path.join(OUTDIR, 'chebi_compounds_20201216153117.csv'), sep = '\t', dtype=str) ## 1
dfChebiCompounds['ParentalChebiIds'] = dfChebiCompounds['ParentalChebiIds'].apply(literal_eval)
dfChebiCompounds['AllChebiIds'] = dfChebiCompounds['AllChebiIds'].apply(literal_eval)
dfChebiCompounds['ParentalKeggC'] = dfChebiCompounds['ParentalKeggC'].apply(literal_eval)
dfChebiCompounds['ParentalKeggG'] = dfChebiCompounds['ParentalKeggG'].apply(literal_eval)
dfChebiCompounds['ParentalMetacyc'] = dfChebiCompounds['ParentalMetacyc'].apply(literal_eval)
dfChebiCompounds['AllKeggC'] = dfChebiCompounds['AllKeggC'].apply(literal_eval)
dfChebiCompounds['AllKeggG'] = dfChebiCompounds['AllKeggG'].apply(literal_eval)
dfChebiCompounds['AllMetacyc'] = dfChebiCompounds['AllMetacyc'].apply(literal_eval)


dfMetsFromModel['Name'] = dfMetsFromModel['Name'].str.strip()


dUniprotId_grouped_all = []
lIds_all = []
for rowRxn in dfRxns.itertuples():
    dUniprotId_grouped = {}
    isTransport = rowRxn.IsTransport
    isExchange = rowRxn.IsExchange
    if isTransport == True or isExchange == True:
        if rowRxn.Rxn.startswith('R_'):
            rxnId = rowRxn.Rxn[2:]
        else:
            rxnId = rowRxn.Rxn
        rxn = model.reactions.get_by_id(rxnId)
        print('ID\t', rxn.id)
        print('isTransport\t', isTransport)
        print('isExchange\t', isExchange)
        lIds_all.append(rxn.id)
        lIdentifiersRxn = [] ## lista dove salvo tutti gli identificativi che trovo associati alla rxn corrente

        lReactants = []
        if len(rowRxn.trasportedMets) != 0:
            lReactants = rowRxn.trasportedMets
            # print('lReactants1\t', lReactants, '\t', len(lReactants))
            dfMetsFromModel_reactants = dfMetsFromModel[dfMetsFromModel['Name'].isin(lReactants)]
            # print('dfMetsFromModel_reactants1\t', dfMetsFromModel_reactants)

            dfMetsFromModel_reactants = dfMetsFromModel_reactants.drop_duplicates(subset = ['Name'])
            lReactants_ids_original = list(dfMetsFromModel_reactants['lIdentifiers'].dropna())
            # print('lReactants_ids_original1\t', lReactants_ids_original)

        else:
            for r in rxn.reactants:
                lReactants.append(r.id)
            # print('lReactants2\t', lReactants)
            if testModel == 'recon3' or testModel == 'hmr':
                dfMetsFromModel_reactants = dfMetsFromModel[dfMetsFromModel['Id'].isin(['M_' + m for m in lReactants])] ### Id_converted
            else:
                dfMetsFromModel_reactants = dfMetsFromModel[dfMetsFromModel['Id'].isin(lReactants)]
            # print('dfMetsFromModel_reactants2\t', dfMetsFromModel_reactants)

            dfMetsFromModel_reactants = dfMetsFromModel_reactants.drop_duplicates(subset = ['Id'])
            lReactants_ids_original = list(dfMetsFromModel_reactants['lIdentifiers'].dropna())
            # print('lReactants_ids_original2\t', lReactants_ids_original)


        lProducts = []
        if len(rowRxn.trasportedMets) != 0:
            lProducts = rowRxn.trasportedMets
            dfMetsFromModel_products = dfMetsFromModel[dfMetsFromModel['Name'].isin(lProducts)]### Id_converted

            # print('dfMetsFromModel_products\n', dfMetsFromModel_products, '\n')
            dfMetsFromModel_products = dfMetsFromModel_products.drop_duplicates(subset = ['Name'])
            lProducts_ids_original = list(dfMetsFromModel_products['lIdentifiers'].dropna())
        else:
            lProducts = []
            for r in rxn.products:
                lProducts.append(r.id)

            if testModel == 'recon3' or testModel == 'hmr':
                dfMetsFromModel_products = dfMetsFromModel[dfMetsFromModel['Id'].isin(['M_' + m for m in lProducts])]### Id_converted
            else:
                dfMetsFromModel_products = dfMetsFromModel[dfMetsFromModel['Id'].isin(lProducts)]### Id_converted

            dfMetsFromModel_products = dfMetsFromModel_products.drop_duplicates(subset = ['Id'])
            lProducts_ids_original = list(dfMetsFromModel_products['lIdentifiers'].dropna())
        print('lReactants\t', lReactants)
        print('lProducts\t', lProducts)

        print('lReactants_ids_original\t', lReactants_ids_original)
        print('lProducts_ids_original\t', lProducts_ids_original, '\n')

        lReactants_ids_original_tuple = map(tuple, lReactants_ids_original)
        lProducts_ids_original_tuple = map(tuple, lProducts_ids_original)
        ltransportedMets = list(set(lReactants_ids_original_tuple).intersection(lProducts_ids_original_tuple))
        print('ltransportedMets\t', ltransportedMets)

        if len(ltransportedMets) != 0:
            lReactants_ids = []
            for l in ltransportedMets:
                l = [el for el in l if el.isdigit() == True]
                chebiIdentifiers = gL.intersect(l,list(dfChebiCompounds['Id']))
                print('chebiIdentifiers\t', chebiIdentifiers)

                dfChebiExplodedFilter = dfChebiCompounds_exploded[dfChebiCompounds_exploded['ParentalChebiIds_exploded'].isin(chebiIdentifiers)]
                dfChebiFilter = dfChebiCompounds[dfChebiCompounds['Id'].isin(chebiIdentifiers)]

                lAllChebi2Search = []
                if dfChebiExplodedFilter.empty is False:
                    for riga in list(dfChebiExplodedFilter['ParentalChebiIds']):
                        lAllChebi2Search += riga
                    for riga in list(dfChebiExplodedFilter['AllChebiIds']):
                        lAllChebi2Search += riga

                if dfChebiFilter.empty is False:
                    for riga in list(dfChebiFilter['ParentalChebiIds']):
                        lAllChebi2Search += riga
                    for riga in list(dfChebiFilter['AllChebiIds']):
                        lAllChebi2Search += riga
                lAllChebi2Search = gL.unique(lAllChebi2Search)
                # print('lAllChebi2Search\t', lAllChebi2Search, '\n')
                ## filtro lAllChebi2Search su quelli che hanno formula cosi da non avere un livello troppo generico di descrizione del metabolita
                dfChebiFormula_filtered = dfChebiFormula[dfChebiFormula['COMPOUND_ID'].isin(lAllChebi2Search)]
                if dfChebiFormula_filtered.empty is False:
                    dfChebiFormula_filtered = dfChebiFormula_filtered.reset_index(drop = True)
                    dfChebiFormula_filtered_toKeep = dfChebiFormula_filtered[dfChebiFormula_filtered['TYPE'] == 'FORMULA']
                    if dfChebiFormula_filtered_toKeep.empty is False:
                        # print('Con formula\t', list(dfChebiFormula_filtered_toKeep['COMPOUND_ID']), '\n')
                        l += list(dfChebiFormula_filtered_toKeep['COMPOUND_ID'])
                        for riga in list(dfChebiFormula_filtered_toKeep['COMPOUND_ID']):
                            dfChebiExplodedFilter2 = dfChebiCompounds_exploded[dfChebiCompounds_exploded['ParentalChebiIds_exploded'].isin(list(dfChebiFormula_filtered_toKeep['COMPOUND_ID']))]
                            for riga2 in list(dfChebiExplodedFilter2['ParentalChebiIds']):
                                l += riga2
                                # print('Eventuali alternativi \t', riga2, '\n')

                l = gL.unique(l)
                lReactants_ids.append(l)
            print('lReactants_ids\n', lReactants_ids, '\n')

            if len(lReactants_ids) != 0:
                dfTc2Subs = dfTc2Subs_original.copy()
                dfTc2Subs['Find_allItems'] = dfTc2Subs.apply(lambda row: findTC_transports(lReactants_ids, row['Substrates_converted']), axis=1)
                putativeTransportRxn = dfTc2Subs[dfTc2Subs['Find_allItems'] == True]
                if putativeTransportRxn.empty is False:
                    lPutativeTransports = list(putativeTransportRxn['TCnumber'])
                    # print('lPutativeTransports\t', lPutativeTransports)
                    dfUniprotId = dfTc2Uniprot[dfTc2Uniprot['TCnumber'].isin(lPutativeTransports)]
                    if dfUniprotId.empty is False:
                        dfUniprotId = dfUniprotId[dfUniprotId['UniprotId'].notna()]
                        print('dfUniprotId\t', dfUniprotId.shape)
                        dfUniprotId_grouped = dfUniprotId.groupby('TCnumber')['UniprotId'].apply(list)
                        print(dfUniprotId_grouped, '\n')
                        dUniprotId_grouped = dfUniprotId_grouped.to_dict()
                        # print('dUniprotId_grouped\t', dUniprotId_grouped)
        elif len(ltransportedMets) == 0 and ((len(lReactants_ids_original) != 0 and len(lProducts_ids_original) == 0) or (len(lProducts_ids_original) != 0 and len(lReactants_ids_original) == 0)):
            print('EXCHANGE')

            lReactants_ids = []
            if len(lReactants_ids_original) != 0:
                lAll = lReactants_ids_original.copy()
            elif len(lProducts_ids_original) != 0:
                lAll = lProducts_ids_original.copy()
            for l in lAll:
                l = [el for el in l if el.isdigit() == True]
                chebiIdentifiers = gL.intersect(l,list(dfChebiCompounds['Id']))
                # print('chebiIdentifiers\t', chebiIdentifiers)

                dfChebiExplodedFilter = dfChebiCompounds_exploded[dfChebiCompounds_exploded['ParentalChebiIds_exploded'].isin(chebiIdentifiers)]
                dfChebiFilter = dfChebiCompounds[dfChebiCompounds['Id'].isin(chebiIdentifiers)]

                lAllChebi2Search = []
                if dfChebiExplodedFilter.empty is False:
                    for riga in list(dfChebiExplodedFilter['ParentalChebiIds']):
                        lAllChebi2Search += riga
                    for riga in list(dfChebiExplodedFilter['AllChebiIds']):
                        lAllChebi2Search += riga

                if dfChebiFilter.empty is False:
                    for riga in list(dfChebiFilter['ParentalChebiIds']):
                        lAllChebi2Search += riga
                    for riga in list(dfChebiFilter['AllChebiIds']):
                        lAllChebi2Search += riga
                lAllChebi2Search = gL.unique(lAllChebi2Search)
                # print('lAllChebi2Search\t', lAllChebi2Search, '\n')
                dfChebiFormula_filtered = dfChebiFormula[dfChebiFormula['COMPOUND_ID'].isin(lAllChebi2Search)]
                if dfChebiFormula_filtered.empty is False:
                    dfChebiFormula_filtered = dfChebiFormula_filtered.reset_index(drop = True)
                    dfChebiFormula_filtered_toKeep = dfChebiFormula_filtered[dfChebiFormula_filtered['TYPE'] == 'FORMULA']
                    if dfChebiFormula_filtered_toKeep.empty is False:
                        # print('Con formula\t', list(dfChebiFormula_filtered_toKeep['COMPOUND_ID']), '\n')
                        l += list(dfChebiFormula_filtered_toKeep['COMPOUND_ID'])
                        for riga in list(dfChebiFormula_filtered_toKeep['COMPOUND_ID']):
                            dfChebiExplodedFilter2 = dfChebiCompounds_exploded[dfChebiCompounds_exploded['ParentalChebiIds_exploded'].isin(list(dfChebiFormula_filtered_toKeep['COMPOUND_ID']))]
                            for riga2 in list(dfChebiExplodedFilter2['ParentalChebiIds']):
                                l += riga2
                                # print('Eventuali alternativi \t', riga2, '\n')

                l = gL.unique(l)
                lReactants_ids.append(l)
            print('lReactants_ids\n', lReactants_ids, '\n')

            if all(r != [] for r in lAll) == True:
                dfTc2Subs = dfTc2Subs_original.copy()
                dfTc2Subs['Find_allItems'] = dfTc2Subs.apply(lambda row: findTC_transports(lReactants_ids, row['Substrates_converted']), axis=1)
                putativeTransportRxn = dfTc2Subs[dfTc2Subs['Find_allItems'] == True]
                if putativeTransportRxn.empty is False:
                    print('putativeTransportRxn\n',putativeTransportRxn.shape)
                    lPutativeTransports = list(putativeTransportRxn['TCnumber'])
                    # print('lPutativeTransports\t', len(lPutativeTransports))
                    dfUniprotId = dfTc2Uniprot[dfTc2Uniprot['TCnumber'].isin(lPutativeTransports)]
                    if dfUniprotId.empty is False:
                        dfUniprotId = dfUniprotId[dfUniprotId['UniprotId'].notna()]
                        print('dfUniprotId\t', dfUniprotId.shape)
                        dfUniprotId_grouped = dfUniprotId.groupby('TCnumber')['UniprotId'].apply(list)
                        print(dfUniprotId_grouped, '\n')
                        dUniprotId_grouped = dfUniprotId_grouped.to_dict()
                        # print('dUniprotId_grouped\t', dUniprotId_grouped)

    # print('dUniprotId_grouped\t', dUniprotId_grouped)
    dUniprotId_grouped_all.append(dUniprotId_grouped)

dfRxns['Identifiers_fromTCDB'] = dUniprotId_grouped_all
dfRxns.to_csv(os.path.join(OUTDIR, dfrxnsInfo + '_enriched_tcdb_' + timeStamp + '.csv'), sep = '\t', index = False)

end = time.time()
print('elapsed time\t', end-start, '\tseconds')
