import os
import sys
import cobra as cb
import pandas as pd
import genericLib as gL
import reactionsLib as rxnL
from ast import literal_eval

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]

# setting input data
testModel = sys.argv[1]
if testModel == 'recon3':
    ## Recon 3
    modelXml = 'Recon3D_301_20200923'
    dfmetsInfo = 'recon3D_metabolites'
    dfmetsIds = 'recon3_mappingFuzzyAndClassic'
    dfrxnsInfo = 'recon3D_reactions'
elif testModel == 'y7':
    ## Yeast 7
    modelXml = 'yeast_7.6_cobra'
    dfmetsInfo = 'y7_metabolites'
    dfmetsIds = 'yeast7_mappingFuzzyAndClassic'
    dfrxnsInfo = 'y7_reactions'
elif testModel == 'y8':
    ## Yeast 8
    modelXml = 'yeast8'
    dfmetsInfo = 'y8_metabolites'
    dfmetsIds = 'yeast8_mappingFuzzyAndClassic'
    dfrxnsInfo = 'y8_reactions'
elif testModel == 'hmr':
    ## HMRcore
    modelXml = 'HMRcore_20200328_wReconNames'
    dfmetsInfo = 'hmrCore_metabolites'
    dfmetsIds = 'hmrCore_mappingFuzzyAndClassic'
    dfrxnsInfo = 'hmrCore_reactions'
elif testModel == 'ownData':
    ## specify your input data
    modelXml = ''
    dfmetsInfo = ''
    dfmetsIds = ''
    dfrxnsInfo = ''

model = cb.io.read_sbml_model(os.path.join(RAWDIR, modelXml + '.xml'))
dfRxns = pd.read_csv(os.path.join(OUTDIR, dfrxnsInfo + '.csv'), sep = '\t', dtype = {'Rxn': str, 'KeggId': str, 'GPR': str, 'Name': str, 'IsTransport': bool, 'IsExchange': bool, 'GPRrule': str})
dfRxns['EC number'] = dfRxns['EC number'].apply(literal_eval)
dfRxns['trasportedMets'] = dfRxns['trasportedMets'].apply(literal_eval)

dfTc2Subs_original = pd.read_csv(os.path.join(RAWDIR, 'tcdb_tc2substrates.csv'), sep = '\t', dtype=str, names= ['TCnumber', 'Substrates'])
dfTc2Uniprot = pd.read_csv(os.path.join(RAWDIR, 'tcdb_tc2Uniprot.csv'), sep = '\t', dtype=str, names= ['UniprotId', 'TCnumber'])

dfTc2Subs_original['Substrates_splt'] = dfTc2Subs_original['Substrates'].str.split('|')
dfTc2Subs_original['Substrates_converted'] = dfTc2Subs_original.apply(lambda row: rxnL.extractChebiIds(row['Substrates_splt']), axis=1)

dfMetsFromModel = pd.read_csv(os.path.join(OUTDIR, dfmetsInfo + '_enriched.csv'), sep = '\t', dtype=str)
dfMetsFromModel['lIdentifiers'] = dfMetsFromModel['lIdentifiers'].apply(literal_eval)

if testModel == 'recon3':
    dfMetsFromModel['Id_converted'] = dfMetsFromModel['Id'].str.strip('M_')
    dfMetsFromModel['Id_converted'] = dfMetsFromModel['Id_converted'].str.replace('__91__', '[')
    dfMetsFromModel['Id_converted'] = dfMetsFromModel['Id_converted'].str.replace('__93', ']')
    dfMetsFromModel['Id_converted2'] = dfMetsFromModel['Id'].str.strip('M_')
    dfMetsFromModel['Id_converted2'] = dfMetsFromModel['Id_converted2'] + '__'


## Read from ChEBI all the parental identifiers of each metabolite
dfChebiFormula = pd.read_csv(os.path.join(RAWDIR, 'chebi_chemical_data_20201216.tsv'), sep = '\t', dtype=str)
dfChebiCompounds_exploded = pd.read_csv(os.path.join(RAWDIR, 'chebi_compounds_20201216153117_exploded.csv.bz2'), compression='bz2', sep = '\t', dtype=str)
dfChebiCompounds_exploded['ParentalChebiIds'] = dfChebiCompounds_exploded['ParentalChebiIds'].apply(literal_eval)
dfChebiCompounds_exploded['AllChebiIds'] = dfChebiCompounds_exploded['AllChebiIds'].apply(literal_eval)
dfChebiCompounds_exploded['ParentalKeggC'] = dfChebiCompounds_exploded['ParentalKeggC'].apply(literal_eval)
dfChebiCompounds_exploded['ParentalKeggG'] = dfChebiCompounds_exploded['ParentalKeggG'].apply(literal_eval)
dfChebiCompounds_exploded['ParentalMetacyc'] = dfChebiCompounds_exploded['ParentalMetacyc'].apply(literal_eval)
dfChebiCompounds_exploded['AllKeggC'] = dfChebiCompounds_exploded['AllKeggC'].apply(literal_eval)
dfChebiCompounds_exploded['AllKeggG'] = dfChebiCompounds_exploded['AllKeggG'].apply(literal_eval)
dfChebiCompounds_exploded['AllMetacyc'] = dfChebiCompounds_exploded['AllMetacyc'].apply(literal_eval)

dfChebiCompounds = pd.read_csv(os.path.join(RAWDIR, 'chebi_compounds_20201216153117.csv.bz2'), compression='bz2',  sep = '\t', dtype=str)
dfChebiCompounds['ParentalChebiIds'] = dfChebiCompounds['ParentalChebiIds'].apply(literal_eval)
dfChebiCompounds['AllChebiIds'] = dfChebiCompounds['AllChebiIds'].apply(literal_eval)
dfChebiCompounds['ParentalKeggC'] = dfChebiCompounds['ParentalKeggC'].apply(literal_eval)
dfChebiCompounds['ParentalKeggG'] = dfChebiCompounds['ParentalKeggG'].apply(literal_eval)
dfChebiCompounds['ParentalMetacyc'] = dfChebiCompounds['ParentalMetacyc'].apply(literal_eval)
dfChebiCompounds['AllKeggC'] = dfChebiCompounds['AllKeggC'].apply(literal_eval)
dfChebiCompounds['AllKeggG'] = dfChebiCompounds['AllKeggG'].apply(literal_eval)
dfChebiCompounds['AllMetacyc'] = dfChebiCompounds['AllMetacyc'].apply(literal_eval)


dfMetsFromModel['Name'] = dfMetsFromModel['Name'].str.strip()

# Retrieve from TCDB database the transport reactions identifiers
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
        lIds_all.append(rxn.id)
        lIdentifiersRxn = []
        lReactants = []
        if len(rowRxn.trasportedMets) != 0:
            lReactants = rowRxn.trasportedMets
            dfMetsFromModel_reactants = dfMetsFromModel[dfMetsFromModel['Name'].isin(lReactants)]
            dfMetsFromModel_reactants = dfMetsFromModel_reactants.drop_duplicates(subset = ['Name'])
            lReactants_ids_original = list(dfMetsFromModel_reactants['lIdentifiers'].dropna())
        else:
            for r in rxn.reactants:
                lReactants.append(r.id)
            if testModel == 'recon3' or testModel == 'hmr':
                dfMetsFromModel_reactants = dfMetsFromModel[dfMetsFromModel['Id'].isin(['M_' + m for m in lReactants])]
            else:
                dfMetsFromModel_reactants = dfMetsFromModel[dfMetsFromModel['Id'].isin(lReactants)]
            dfMetsFromModel_reactants = dfMetsFromModel_reactants.drop_duplicates(subset = ['Id'])
            lReactants_ids_original = list(dfMetsFromModel_reactants['lIdentifiers'].dropna())
        lProducts = []
        if len(rowRxn.trasportedMets) != 0:
            lProducts = rowRxn.trasportedMets
            dfMetsFromModel_products = dfMetsFromModel[dfMetsFromModel['Name'].isin(lProducts)]
            dfMetsFromModel_products = dfMetsFromModel_products.drop_duplicates(subset = ['Name'])
            lProducts_ids_original = list(dfMetsFromModel_products['lIdentifiers'].dropna())
        else:
            lProducts = []
            for r in rxn.products:
                lProducts.append(r.id)
            if testModel == 'recon3' or testModel == 'hmr':
                dfMetsFromModel_products = dfMetsFromModel[dfMetsFromModel['Id'].isin(['M_' + m for m in lProducts])]
            else:
                dfMetsFromModel_products = dfMetsFromModel[dfMetsFromModel['Id'].isin(lProducts)]

            dfMetsFromModel_products = dfMetsFromModel_products.drop_duplicates(subset = ['Id'])
            lProducts_ids_original = list(dfMetsFromModel_products['lIdentifiers'].dropna())

        lReactants_ids_original_tuple = map(tuple, lReactants_ids_original)
        lProducts_ids_original_tuple = map(tuple, lProducts_ids_original)
        ltransportedMets = list(set(lReactants_ids_original_tuple).intersection(lProducts_ids_original_tuple))

        if len(ltransportedMets) != 0:
            lReactants_ids = []
            for l in ltransportedMets:
                l = [el for el in l if el.isdigit() == True]
                chebiIdentifiers = gL.intersect(l,list(dfChebiCompounds['Id']))

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
                dfChebiFormula_filtered = dfChebiFormula[dfChebiFormula['COMPOUND_ID'].isin(lAllChebi2Search)]
                if dfChebiFormula_filtered.empty is False:
                    dfChebiFormula_filtered = dfChebiFormula_filtered.reset_index(drop = True)
                    dfChebiFormula_filtered_toKeep = dfChebiFormula_filtered[dfChebiFormula_filtered['TYPE'] == 'FORMULA']
                    if dfChebiFormula_filtered_toKeep.empty is False:
                        l += list(dfChebiFormula_filtered_toKeep['COMPOUND_ID'])
                        for riga in list(dfChebiFormula_filtered_toKeep['COMPOUND_ID']):
                            dfChebiExplodedFilter2 = dfChebiCompounds_exploded[dfChebiCompounds_exploded['ParentalChebiIds_exploded'].isin(list(dfChebiFormula_filtered_toKeep['COMPOUND_ID']))]
                            for riga2 in list(dfChebiExplodedFilter2['ParentalChebiIds']):
                                l += riga2

                l = gL.unique(l)
                lReactants_ids.append(l)

            if len(lReactants_ids) != 0:
                dfTc2Subs = dfTc2Subs_original.copy()
                dfTc2Subs['Find_allItems'] = dfTc2Subs.apply(lambda row: rxnL.findTC_transports(lReactants_ids, row['Substrates_converted']), axis=1)
                putativeTransportRxn = dfTc2Subs[dfTc2Subs['Find_allItems'] == True]
                if putativeTransportRxn.empty is False:
                    lPutativeTransports = list(putativeTransportRxn['TCnumber'])
                    dfUniprotId = dfTc2Uniprot[dfTc2Uniprot['TCnumber'].isin(lPutativeTransports)]
                    if dfUniprotId.empty is False:
                        dfUniprotId = dfUniprotId[dfUniprotId['UniprotId'].notna()]
                        dfUniprotId_grouped = dfUniprotId.groupby('TCnumber')['UniprotId'].apply(list)
                        dUniprotId_grouped = dfUniprotId_grouped.to_dict()
        elif len(ltransportedMets) == 0 and ((len(lReactants_ids_original) != 0 and len(lProducts_ids_original) == 0) or (len(lProducts_ids_original) != 0 and len(lReactants_ids_original) == 0)):
            lReactants_ids = []
            if len(lReactants_ids_original) != 0:
                lAll = lReactants_ids_original.copy()
            elif len(lProducts_ids_original) != 0:
                lAll = lProducts_ids_original.copy()
            for l in lAll:
                l = [el for el in l if el.isdigit() == True]
                chebiIdentifiers = gL.intersect(l,list(dfChebiCompounds['Id']))
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
                dfChebiFormula_filtered = dfChebiFormula[dfChebiFormula['COMPOUND_ID'].isin(lAllChebi2Search)]
                if dfChebiFormula_filtered.empty is False:
                    dfChebiFormula_filtered = dfChebiFormula_filtered.reset_index(drop = True)
                    dfChebiFormula_filtered_toKeep = dfChebiFormula_filtered[dfChebiFormula_filtered['TYPE'] == 'FORMULA']
                    if dfChebiFormula_filtered_toKeep.empty is False:
                        l += list(dfChebiFormula_filtered_toKeep['COMPOUND_ID'])
                        for riga in list(dfChebiFormula_filtered_toKeep['COMPOUND_ID']):
                            dfChebiExplodedFilter2 = dfChebiCompounds_exploded[dfChebiCompounds_exploded['ParentalChebiIds_exploded'].isin(list(dfChebiFormula_filtered_toKeep['COMPOUND_ID']))]
                            for riga2 in list(dfChebiExplodedFilter2['ParentalChebiIds']):
                                l += riga2

                l = gL.unique(l)
                lReactants_ids.append(l)

            if all(r != [] for r in lAll) == True:
                dfTc2Subs = dfTc2Subs_original.copy()
                dfTc2Subs['Find_allItems'] = dfTc2Subs.apply(lambda row: rxnL.findTC_transports(lReactants_ids, row['Substrates_converted']), axis=1)
                putativeTransportRxn = dfTc2Subs[dfTc2Subs['Find_allItems'] == True]
                if putativeTransportRxn.empty is False:
                    lPutativeTransports = list(putativeTransportRxn['TCnumber'])
                    dfUniprotId = dfTc2Uniprot[dfTc2Uniprot['TCnumber'].isin(lPutativeTransports)]
                    if dfUniprotId.empty is False:
                        dfUniprotId = dfUniprotId[dfUniprotId['UniprotId'].notna()]
                        dfUniprotId_grouped = dfUniprotId.groupby('TCnumber')['UniprotId'].apply(list)
                        dUniprotId_grouped = dfUniprotId_grouped.to_dict()

    dUniprotId_grouped_all.append(dUniprotId_grouped)

dfRxns['Identifiers_fromTCDB'] = dUniprotId_grouped_all
dfRxns.to_csv(os.path.join(OUTDIR, dfrxnsInfo + '_enriched_tcdb.csv'), sep = '\t', index = False)
