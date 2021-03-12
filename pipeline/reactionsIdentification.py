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


# Enrich dfMetsFromModel with the information coming from dfMetsIdentifiers
dfmetacyc = pd.read_csv(os.path.join(RAWDIR, 'metacyc_compounds_20201216152513.csv'), sep = '\t', dtype=str)
dfkeggC = pd.read_csv(os.path.join(RAWDIR, 'kegg_compounds_20201216150802.csv'), sep = '\t', dtype=str)
dfkeggC['ChebiId'] = dfkeggC['ChebiId'].apply(literal_eval)
dfkeggG = pd.read_csv(os.path.join(RAWDIR, 'kegg_glycans_20201216150802.csv'), sep = '\t', dtype=str)
dfkeggG['ChebiId'] = dfkeggG['ChebiId'].apply(literal_eval)

dfMetsFromModel = pd.read_csv(os.path.join(OUTDIR, dfmetsInfo + '.csv'), sep = '\t', dtype=str)

dfMetsIdentifiers = pd.read_csv(os.path.join(OUTDIR, dfmetsIds + '.tsv'), sep = '\t')
dfMetsIdentifiers['Name'] = dfMetsIdentifiers['Name'].str.strip()
dfMetsIdentifiers['Identifiers'] = dfMetsIdentifiers['Identifiers'].apply(literal_eval)
dizMetsIdentifiers = dfMetsIdentifiers.set_index('Name')['Identifiers'].to_dict()

for k in dizMetsIdentifiers:
    lCompleteIds = []
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

lIdentifiers = []
for row in dfMetsFromModel.itertuples():
    lInferredIds = dizMetsIdentifiers[row.Name.strip()]
    lIdentifiers.append(lInferredIds)

dfMetsFromModel['lIdentifiers'] = lIdentifiers
dfMetsFromModel.to_csv(os.path.join(OUTDIR, dfmetsInfo + '_enriched.csv'), sep = '\t', index = False)


## Read from ChEBI all the parental identifiers of each metabolite
dfChebiFormula = pd.read_csv(os.path.join(RAWDIR, 'chebi_chemical_data_20201216.tsv'), sep = '\t', dtype=str)

dfChebiCompounds_exploded = pd.read_csv(os.path.join(RAWDIR, 'chebi_compounds_20201216153117_exploded.csv.bz2'),compression='bz2', sep = '\t', dtype=str) ## 1
dfChebiCompounds_exploded['ParentalChebiIds'] = dfChebiCompounds_exploded['ParentalChebiIds'].apply(literal_eval)
dfChebiCompounds_exploded['AllChebiIds'] = dfChebiCompounds_exploded['AllChebiIds'].apply(literal_eval)
dfChebiCompounds_exploded['ParentalKeggC'] = dfChebiCompounds_exploded['ParentalKeggC'].apply(literal_eval)
dfChebiCompounds_exploded['ParentalKeggG'] = dfChebiCompounds_exploded['ParentalKeggG'].apply(literal_eval)
dfChebiCompounds_exploded['ParentalMetacyc'] = dfChebiCompounds_exploded['ParentalMetacyc'].apply(literal_eval)
dfChebiCompounds_exploded['AllKeggC'] = dfChebiCompounds_exploded['AllKeggC'].apply(literal_eval)
dfChebiCompounds_exploded['AllKeggG'] = dfChebiCompounds_exploded['AllKeggG'].apply(literal_eval)
dfChebiCompounds_exploded['AllMetacyc'] = dfChebiCompounds_exploded['AllMetacyc'].apply(literal_eval)

dfChebiCompounds = pd.read_csv(os.path.join(RAWDIR, 'chebi_compounds_20201216153117.csv.bz2'), compression='bz2', sep = '\t', dtype=str) ## 1
dfChebiCompounds['ParentalChebiIds'] = dfChebiCompounds['ParentalChebiIds'].apply(literal_eval)
dfChebiCompounds['AllChebiIds'] = dfChebiCompounds['AllChebiIds'].apply(literal_eval)
dfChebiCompounds['ParentalKeggC'] = dfChebiCompounds['ParentalKeggC'].apply(literal_eval)
dfChebiCompounds['ParentalKeggG'] = dfChebiCompounds['ParentalKeggG'].apply(literal_eval)
dfChebiCompounds['ParentalMetacyc'] = dfChebiCompounds['ParentalMetacyc'].apply(literal_eval)
dfChebiCompounds['AllKeggC'] = dfChebiCompounds['AllKeggC'].apply(literal_eval)
dfChebiCompounds['AllKeggG'] = dfChebiCompounds['AllKeggG'].apply(literal_eval)
dfChebiCompounds['AllMetacyc'] = dfChebiCompounds['AllMetacyc'].apply(literal_eval)

dfMetsFromModel = pd.read_csv(os.path.join(OUTDIR, dfmetsInfo + '_enriched.csv'), sep = '\t', dtype=str)
dfMetsFromModel['lIdentifiers'] = dfMetsFromModel['lIdentifiers'].apply(literal_eval)

if testModel == 'recon3':
    dfMetsFromModel['Id_converted'] = dfMetsFromModel['Id'].str.strip('M_')
    dfMetsFromModel['Id_converted'] = dfMetsFromModel['Id_converted'].str.replace('__91__', '[')
    dfMetsFromModel['Id_converted'] = dfMetsFromModel['Id_converted'].str.replace('__93', ']')
    dfMetsFromModel['Id_converted2'] = dfMetsFromModel['Id'].str.strip('M_')
    dfMetsFromModel['Id_converted2'] = dfMetsFromModel['Id_converted2'] + '__'

# Read the macrodatabase
dfAllDBs = pd.read_csv(os.path.join(RAWDIR, 'dfJoin_metacyc_kegg_rhea_20201218164915.csv'), sep = '\t', dtype=str) ## 3
dfAllDBs['ec_number'] = dfAllDBs['ec_number'].apply(literal_eval)
dfAllDBs['enzymatic_reaction'] = dfAllDBs['enzymatic_reaction'].fillna({i: '[]' for i in dfAllDBs.index})
dfAllDBs['enzymatic_reaction'] = dfAllDBs['enzymatic_reaction'].apply(literal_eval)
dfAllDBs['left'] = dfAllDBs['left'].apply(literal_eval)
dfAllDBs['right'] = dfAllDBs['right'].apply(literal_eval)
dfAllDBs['rxn_locations'] = dfAllDBs['rxn_locations'].apply(literal_eval)
dfAllDBs['subreactions'] = dfAllDBs['subreactions'].apply(literal_eval)
dfAllDBs['synonyms'] = dfAllDBs['synonyms'].apply(literal_eval)
dfAllDBs['enzymes_not_used'] = dfAllDBs['enzymes_not_used'].apply(literal_eval)
dfAllDBs['compartments_of_reaction'] = dfAllDBs['compartments_of_reaction'].apply(literal_eval)
dfAllDBs['enzymes_of_reaction'] = dfAllDBs['enzymes_of_reaction'].apply(literal_eval)
dfAllDBs['genes_of_reaction'] = dfAllDBs['genes_of_reaction'].apply(literal_eval)
dfAllDBs['specific_forms_of_rxn'] = dfAllDBs['specific_forms_of_rxn'].apply(literal_eval)
dfAllDBs['Subs_fromKegg'] = dfAllDBs['Subs_fromKegg'].apply(literal_eval)
dfAllDBs['Prods_fromKegg'] = dfAllDBs['Prods_fromKegg'].apply(literal_eval)

dfAllDBs['left_metacyc'] = dfAllDBs.apply(lambda row: rxnL.removeProtons(row['left'], row['reaction_type']), axis=1)
dfAllDBs['right_metacyc'] = dfAllDBs.apply(lambda row: rxnL.removeProtons(row['right'], row['reaction_type']), axis=1)
dfAllDBs['lSubs_fromKegg'] = dfAllDBs.apply(lambda row: rxnL.removeProtonsKegg(row['Subs_fromKegg'], row['reaction_type']), axis=1)
dfAllDBs['lProds_fromKegg'] = dfAllDBs.apply(lambda row: rxnL.removeProtonsKegg(row['Prods_fromKegg'], row['reaction_type']), axis=1)

dfAllDBs['OtherRheaId_fromKegg'] = dfAllDBs['OtherRheaId_fromKegg'].apply(literal_eval)
lotherrhea = []
for el in dfAllDBs['OtherRheaId_fromKegg']:
    if type(el) == list:
        lotherrhea.append(el)
    elif type(el) == int:
        lotherrhea.append([el])
dfAllDBs['OtherRheaId_fromKegg'] = lotherrhea

dfAllDBs['RheaId_master'] = dfAllDBs['RheaId_master'].fillna({i: '[]' for i in dfAllDBs.index})
dfAllDBs['RheaId_master'] = dfAllDBs['RheaId_master'].apply(literal_eval)
dfAllDBs['RheaId_lr'] = dfAllDBs['RheaId_lr'].fillna({i: '[]' for i in dfAllDBs.index})
dfAllDBs['RheaId_lr'] = dfAllDBs['RheaId_lr'].apply(literal_eval)
dfAllDBs['RheaId_rl'] = dfAllDBs['RheaId_rl'].fillna({i: '[]' for i in dfAllDBs.index})
dfAllDBs['RheaId_rl'] = dfAllDBs['RheaId_rl'].apply(literal_eval)
dfAllDBs['RheaId_bi'] = dfAllDBs['RheaId_bi'].fillna({i: '[]' for i in dfAllDBs.index})
dfAllDBs['RheaId_bi'] = dfAllDBs['RheaId_bi'].apply(literal_eval)
dfAllDBs['KeggId_y'] = dfAllDBs['KeggId_y'].fillna({i: '[]' for i in dfAllDBs.index})
dfAllDBs['KeggId_y'] = dfAllDBs['KeggId_y'].apply(literal_eval)
dfAllDBs['Outgoings'] = dfAllDBs['Outgoings'].fillna({i: '[]' for i in dfAllDBs.index})
dfAllDBs['Outgoings'] = dfAllDBs['Outgoings'].apply(literal_eval)
dfAllDBs['Reactants_name_id'] = dfAllDBs['Reactants_name_id'].fillna({i: '[]' for i in dfAllDBs.index})
dfAllDBs['Reactants_name_id'] = dfAllDBs['Reactants_name_id'].apply(literal_eval)
dfAllDBs['Products_name_id'] = dfAllDBs['Products_name_id'].fillna({i: '[]' for i in dfAllDBs.index})
dfAllDBs['Products_name_id'] = dfAllDBs['Products_name_id'].apply(literal_eval)

dfAllDBs['lReactants_name_id'] = dfAllDBs.apply(lambda row: rxnL.removeProtonsRhea(row['Reactants_name_id'], row['reaction_type']), axis=1)
dfAllDBs['lProducts_name_id'] = dfAllDBs.apply(lambda row: rxnL.removeProtonsRhea(row['Products_name_id'], row['reaction_type']), axis=1)

list_cols = {'lReactants_name_id','lProducts_name_id'}
other_cols = list(set(dfAllDBs.columns) - set(list_cols))
exploded = [dfAllDBs[col].explode() for col in list_cols]
dfAllDBs_explode = pd.DataFrame(dict(zip(list_cols, exploded)))
dfAllDBs_explode = dfAllDBs[other_cols].merge(dfAllDBs_explode, how="right", left_index=True, right_index=True)
dfAllDBs_explode = dfAllDBs_explode.reset_index(drop = True)

dfAllDBs_explode['lReactants_name_id'] = dfAllDBs_explode['lReactants_name_id'].fillna({i: [] for i in dfAllDBs_explode.index})
dfAllDBs_explode['lProducts_name_id'] = dfAllDBs_explode['lProducts_name_id'].fillna({i: [] for i in dfAllDBs_explode.index})

dfAllDBs_explode['lReactants_name_id_flat'] = dfAllDBs_explode.apply(lambda row: rxnL.flatList(row['lReactants_name_id']), axis=1)
dfAllDBs_explode['lProducts_name_id_flat'] = dfAllDBs_explode.apply(lambda row: rxnL.flatList(row['lProducts_name_id']), axis=1)

dfAllDBs_explode['left_metacyc_l'] = dfAllDBs_explode.left_metacyc.str.len()
dfAllDBs_explode['right_metacyc_l'] = dfAllDBs_explode.right_metacyc.str.len()

dfAllDBs_explode['lSubs_fromKegg_l'] = dfAllDBs_explode.lSubs_fromKegg.str.len()
dfAllDBs_explode['lProds_fromKegg_l'] = dfAllDBs_explode.lProds_fromKegg.str.len()

dfAllDBs_explode['lReactants_name_id_flat_l'] = dfAllDBs_explode.lReactants_name_id_flat.str.len()
dfAllDBs_explode['lProducts_name_id_flat_l'] = dfAllDBs_explode.lProducts_name_id_flat.str.len()

dgb_left_metacyc = dfAllDBs_explode.groupby(by=['left_metacyc_l']).indices
dgb_right_metacyc = dfAllDBs_explode.groupby(by=['right_metacyc_l']).indices
dgb_lSubs_fromKegg = dfAllDBs_explode.groupby(by=['lSubs_fromKegg_l']).indices
dgb_lProds_fromKegg = dfAllDBs_explode.groupby(by=['lProds_fromKegg_l']).indices
dgb_lReactants_name_id_flat = dfAllDBs_explode.groupby(by=['lReactants_name_id_flat_l']).indices
dgb_lProducts_name_id_flat = dfAllDBs_explode.groupby(by=['lProducts_name_id_flat_l']).indices

dfAllDBs_equalSP_M = dfAllDBs_explode.copy()
dfAllDBs_equalSP_M = dfAllDBs_equalSP_M[dfAllDBs_equalSP_M['reaction_type'] != '|TRANSPORT|']
dfAllDBs_equalSP_M['equalSP'] = dfAllDBs_equalSP_M.apply(lambda row: rxnL.checkEqualSubsProds(row['left'], row['right']), axis=1)
dfAllDBs_equalSP_filter_M = dfAllDBs_equalSP_M[dfAllDBs_equalSP_M['equalSP'] == True]
dfAllDBs_equalSP_filter_M['equalSP_mets'] = dfAllDBs_equalSP_filter_M.apply(lambda row: rxnL.detectEqualSubsProds(row['left'], row['right']), axis=1)
dfAllDBs_equalSP_filter_M['left_woEqualMet'] = dfAllDBs_equalSP_filter_M.apply(lambda row: rxnL.removeEqualCompound(row['equalSP_mets'], row['left']), axis=1)
dfAllDBs_equalSP_filter_M['right_woEqualMet'] = dfAllDBs_equalSP_filter_M.apply(lambda row: rxnL.removeEqualCompound(row['equalSP_mets'], row['right']), axis=1)

dfAllDBs_equalSP_K = dfAllDBs_explode.copy()
dfAllDBs_equalSP_K = dfAllDBs_equalSP_K[dfAllDBs_equalSP_K['reaction_type'] != '|TRANSPORT|']
dfAllDBs_equalSP_K['equalSP'] = dfAllDBs_equalSP_K.apply(lambda row: rxnL.checkEqualSubsProds(row['lSubs_fromKegg'], row['lProds_fromKegg']), axis=1)
dfAllDBs_equalSP_filter_K = dfAllDBs_equalSP_K[dfAllDBs_equalSP_K['equalSP'] == True]
dfAllDBs_equalSP_filter_K['equalSP_mets'] = dfAllDBs_equalSP_filter_K.apply(lambda row: rxnL.detectEqualSubsProds(row['lSubs_fromKegg'], row['lProds_fromKegg']), axis=1)
dfAllDBs_equalSP_filter_K['left_woEqualMet'] = dfAllDBs_equalSP_filter_K.apply(lambda row: rxnL.removeEqualCompound(row['equalSP_mets'], row['lSubs_fromKegg']), axis=1)
dfAllDBs_equalSP_filter_K['right_woEqualMet'] = dfAllDBs_equalSP_filter_K.apply(lambda row: rxnL.removeEqualCompound(row['equalSP_mets'], row['lProds_fromKegg']), axis=1)


dfAllDBs_equalSP_R = dfAllDBs_explode.copy()
dfAllDBs_equalSP_R = dfAllDBs_equalSP_R[dfAllDBs_equalSP_R['reaction_type'] != '|TRANSPORT|']
dfAllDBs_equalSP_R['equalSP'] = dfAllDBs_equalSP_R.apply(lambda row: rxnL.checkEqualSubsProds(row['lReactants_name_id_flat'], row['lProducts_name_id_flat']), axis=1)
dfAllDBs_equalSP_filter_R = dfAllDBs_equalSP_R[dfAllDBs_equalSP_R['equalSP'] == True]
dfAllDBs_equalSP_filter_R['equalSP_mets'] = dfAllDBs_equalSP_filter_R.apply(lambda row: rxnL.detectEqualSubsProds(row['lReactants_name_id_flat'], row['lProducts_name_id_flat']), axis=1)
dfAllDBs_equalSP_filter_R['left_woEqualMet'] = dfAllDBs_equalSP_filter_R.apply(lambda row: rxnL.removeEqualCompound(row['equalSP_mets'], row['lReactants_name_id_flat']), axis=1)
dfAllDBs_equalSP_filter_R['right_woEqualMet'] = dfAllDBs_equalSP_filter_R.apply(lambda row: rxnL.removeEqualCompound(row['equalSP_mets'], row['lProducts_name_id_flat']), axis=1)
print('dfAllDBs_equalSP_filter_R\t', dfAllDBs_equalSP_filter_R.shape)


# Extract reactions info from the input model
dfRxns = rxnL.getRxnsInfoGEM(os.path.join(RAWDIR, modelXml + ".xml"))
model = cb.io.read_sbml_model(os.path.join(RAWDIR, modelXml + '.xml'))

# Check which reactions are transport reactions
lNames = []
lTransport = []
lMetsTrasportati = []
for row in dfRxns.itertuples():
    if row.Rxn.startswith('R_'):
        rxnId = row.Rxn[2:]
    else:
        rxnId = row.Rxn
    rxn = model.reactions.get_by_id(rxnId)
    lNames.append(rxn.name)
    dReactants = {}
    for reactant in rxn.reactants:
        if testModel == 'y7' or testModel == 'y8':
            dReactants[reactant.id] = reactant.name[:reactant.name.rfind('[')].strip()
        else:
            dReactants[reactant.id] = reactant.name

    lReactants = list(dReactants.values())
    lReactants.sort()
    dProducts = {}
    for product in rxn.products:
        if testModel == 'y7' or testModel == 'y8':
            dProducts[product.id] = product.name[:product.name.rfind('[')].strip()
        else:
            dProducts[product.id] = product.name
    lProducts = list(dProducts.values())
    lProducts.sort()
    if len(gL.intersect(lReactants, lProducts)) != 0:
        lTransport.append(True)
        metTrasportati = gL.intersect(lReactants, lProducts)
        lMetsTrasportati.append(metTrasportati)
    else:
        lTransport.append(False)
        lMetsTrasportati.append([])

dfRxns['Name'] = lNames
dfRxns['IsTransport'] = lTransport
dfRxns['trasportedMets'] = lMetsTrasportati

initialLetter = all(el.startswith('M_') for el in list(dfMetsFromModel['Id']))

# Check which reactions are exchange reactions
lExchange = []
for row in dfRxns.itertuples():
    if row.Rxn.startswith('R_'):
        rxnId = row.Rxn[2:]
    else:
        rxnId = row.Rxn
    rxn = model.reactions.get_by_id(rxnId)
    dReactants = {}
    for reactant in rxn.reactants:
        if initialLetter is True:
            outMet = dfMetsFromModel.loc[dfMetsFromModel['Id'] == 'M_' + reactant.id]
        else:
            outMet = dfMetsFromModel.loc[dfMetsFromModel['Id'] == reactant.id]
        outMet = outMet.reset_index(drop=True)
        dReactants[reactant.id] = outMet.iloc[0]['boundaryCondition']
    dProducts = {}
    for product in rxn.products:
        if initialLetter is True:
            outMet = dfMetsFromModel.loc[dfMetsFromModel['Id'] == 'M_' + product.id]
        else:
            outMet = dfMetsFromModel.loc[dfMetsFromModel['Id'] == product.id]
        outMet = outMet.reset_index(drop=True)
        dProducts[product.id] = outMet.iloc[0]['boundaryCondition']
    if len(dReactants) == 0 or len(dProducts) == 0:
        lExchange.append(True)
    elif any(v == True for v in list(dReactants.values())) is True or any(v == True for v in list(dReactants.values())) is 'True' or any(v == True for v in list(dReactants.values())) is 'true' or any(v == True for v in list(dProducts.values())) is True or any(v == True for v in list(dProducts.values())) is 'True' or any(v == True for v in list(dProducts.values())) is 'true':
        lExchange.append(True)
    else:
        lExchange.append(False)

dfRxns['IsExchange'] = lExchange

# Extract the original GPR rule of the model
lRules = []
for row in dfRxns.itertuples():
    if row.Rxn.startswith('R_'):
        rxnId = row.Rxn[2:]
    else:
        rxnId = row.Rxn
    rxn = model.reactions.get_by_id(rxnId)
    rule = rxn.gene_reaction_rule
    if '__46__' in rule:
        s = '__46__'
        pos = rule.find(s)
        while pos != -1:
            rule = rule.replace(s, '.')
            pos = rule.find(s)
    lRules.append(rule)

dfRxns['GPRrule'] = lRules

dfRxns.to_csv(os.path.join(OUTDIR, dfrxnsInfo + '.csv'), sep = '\t', index = False)

# Infer reaction identifiers
model = cb.io.read_sbml_model(os.path.join(RAWDIR, modelXml + '.xml'))
dfRxns = pd.read_csv(os.path.join(OUTDIR, dfrxnsInfo + '.csv'), sep = '\t', dtype = {'Rxn': str, 'KeggId': str, 'GPR': str, 'Name': str, 'IsTransport': bool, 'IsExchange': bool, 'GPRrule': str})
dfRxns['EC number'] = dfRxns['EC number'].apply(literal_eval)
dfRxns['trasportedMets'] = dfRxns['trasportedMets'].apply(literal_eval)

if testModel == 'recon3':
    dfRxns['Rxn_converted'] = dfRxns['Rxn'].str.replace('__91__', '[')
    dfRxns['Rxn_converted'] = dfRxns['Rxn_converted'].str.replace('__93__', ']')

lcofactors = ['44215', '|NAD|', 'C00003', '15846', '57540', 'C00001', '|WATER|', '15377', '16908', 'C00004', '|NADH|', '57945', '25805', '|OXYGEN-MOLECULE|', '15379', 'C00007', '|NADP|', '58349', '|B-HEP-1:5|', 'C00667', '18009', 'C00006',
            '16474', '77312', '|NADPH|', '|CPD-16005|', 'C20745', '77177', '57783', 'C00005', 'CPD-16005', '15996', 'C00044', '37565', '|GTP|', 'C00010', '57287', '|CO-A|', '15346', 'CARBON-DIOXIDE', 'C00011', '|CARBON-DIOXIDE|', '16526',
            'AMMONIUM', 'C01342', '|AMMONIUM|', '|AMMONIA|', 'C00014', 'AMMONIA', '16134', '28938', 'C00002', '15422', '22258', '|ATP|', '30616', '456216', '73342', '|ADP|', '22251', 'C00008', 'G11113', '22252', '16761']

if testModel == 'y7' or testModel == 'y8':
    dfMetsFromModel['Name'] = dfMetsFromModel['Name'].str.strip()

outFile = open(os.path.join(OUTDIR, dfrxnsInfo + '_wIds.csv'), mode='w')
gL.writeLineByLineToFile(outFile, ['RxnId', 'PutativeIdentifiers'], '\t')

lcofactors = gL.unique(lcofactors)
for rowRxn in dfRxns.itertuples():
    isTransport = rowRxn.IsTransport
    isExchange = rowRxn.IsExchange
    if isTransport == True:
        if rowRxn.Rxn.startswith('R_'):
            rxnId = rowRxn.Rxn[2:]
        else:
            rxnId = rowRxn.Rxn
        rxn = model.reactions.get_by_id(rxnId)
        lIdentifiersRxn = [] # list where all the retrieved identifiers of the current reactions are saved
        lReactants = []
        if isTransport == True:
            lReactants = rowRxn.trasportedMets
        else:
            for r in rxn.reactants:
                lReactants.append(r.id)
        lProducts = []
        if isTransport == True:
            lProducts = rowRxn.trasportedMets
        else:
            for r in rxn.products:
                lProducts.append(r.id)
        if isTransport == True:
            dfMetsFromModel_reactants = dfMetsFromModel[dfMetsFromModel['Name'].isin(lReactants)]
        else:
            if testModel == 'recon3' or testModel == 'hmr':
                dfMetsFromModel_reactants = dfMetsFromModel[dfMetsFromModel['Id'].isin(['M_' + m for m in lReactants])]
            else:
                dfMetsFromModel_reactants = dfMetsFromModel[dfMetsFromModel['Id'].isin(lReactants)]
        lReactants_ids_original = list(dfMetsFromModel_reactants['lIdentifiers'].dropna())

        if isTransport == False and isExchange == False:
            lReactants_ids_original = [l for l in lReactants_ids_original if 'C00080' not in l]
        lReactants_ids = []
        for lIdReactant in lReactants_ids_original:
            lCompsFound = []
            lCompsFound += lIdReactant
            chebiIdentifiers = gL.intersect(lIdReactant,list(dfChebiCompounds['Id']))
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
                    lCompsFound += list(dfChebiFormula_filtered_toKeep['COMPOUND_ID'])

                    for riga in list(dfChebiFormula_filtered_toKeep['COMPOUND_ID']):
                        dfChebiExplodedFilter2 = dfChebiCompounds_exploded[dfChebiCompounds_exploded['ParentalChebiIds_exploded'].isin(list(dfChebiFormula_filtered_toKeep['COMPOUND_ID']))]
                        for riga2 in list(dfChebiExplodedFilter2['ParentalChebiIds']):
                            lCompsFound += riga2

            lCompsFound += [element for el in list(dfChebiFilter['ParentalKeggC']) for element in el] + [element for el in list(dfChebiFilter['ParentalKeggG']) for element in el] + [element for el in list(dfChebiFilter['ParentalMetacyc']) for element in el] + [element for el in list(dfChebiFilter['AllKeggC']) for element in el] + [element for el in list(dfChebiFilter['AllKeggG']) for element in el] + [element for el in list(dfChebiFilter['AllMetacyc']) for element in el]
            lCompsFound = gL.unique(lCompsFound)
            lReactants_ids.append(lCompsFound)
        if isTransport == True:
            dfMetsFromModel_products = dfMetsFromModel[dfMetsFromModel['Name'].isin(lProducts)]
        else:
            if testModel == 'recon3' or testModel == 'hmr':
                dfMetsFromModel_products = dfMetsFromModel[dfMetsFromModel['Id'].isin(['M_' + m for m in lProducts])]
            else:
                dfMetsFromModel_products = dfMetsFromModel[dfMetsFromModel['Id'].isin(lProducts)]
        lProducts_ids_original = list(dfMetsFromModel_products['lIdentifiers'].dropna())
        if isTransport == False and isExchange == False:
            lProducts_ids_original = [l for l in lProducts_ids_original if 'C00080' not in l]
        lProducts_ids = []
        for lIdProduct in lProducts_ids_original:
            lCompsFound = []
            lCompsFound += lIdProduct
            chebiIdentifiers = gL.intersect(lIdProduct,list(dfChebiCompounds['Id']))
            dfChebiFilter = dfChebiCompounds[dfChebiCompounds['Id'].isin(chebiIdentifiers)]
            dfChebiExplodedFilter = dfChebiCompounds_exploded[dfChebiCompounds_exploded['ParentalChebiIds_exploded'].isin(chebiIdentifiers)]

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
                    lCompsFound += list(dfChebiFormula_filtered_toKeep['COMPOUND_ID'])
                    for riga in list(dfChebiFormula_filtered_toKeep['COMPOUND_ID']):
                        dfChebiExplodedFilter2 = dfChebiCompounds_exploded[dfChebiCompounds_exploded['ParentalChebiIds_exploded'].isin(list(dfChebiFormula_filtered_toKeep['COMPOUND_ID']))]
                        for riga2 in list(dfChebiExplodedFilter2['ParentalChebiIds']):
                            lCompsFound += riga2

            lCompsFound += [element for el in list(dfChebiFilter['ParentalKeggC']) for element in el] + [element for el in list(dfChebiFilter['ParentalKeggG']) for element in el] + [element for el in list(dfChebiFilter['ParentalMetacyc']) for element in el] + [element for el in list(dfChebiFilter['AllKeggC']) for element in el] + [element for el in list(dfChebiFilter['AllKeggG']) for element in el] + [element for el in list(dfChebiFilter['AllMetacyc']) for element in el]
            lCompsFound = gL.unique(lCompsFound)
            lProducts_ids.append(lCompsFound)

        dfR = pd.DataFrame(lReactants_ids).transpose()
        dfP = pd.DataFrame(lProducts_ids).transpose()

        if all(r != [] for r in lReactants_ids) == True and all(p != [] for p in lProducts_ids) == True:
            lReactants_ids_cofactors = []
            lReactants_ids_internal = []
            for lReagenti in lReactants_ids:
                if len(gL.intersect(lcofactors, lReagenti)) != 0:
                    lReactants_ids_cofactors += lReagenti
                else:
                    lReactants_ids_internal += lReagenti

            lReactants_ids_cofactors = gL.unique(lReactants_ids_cofactors)
            lReactants_ids_internal = gL.unique(lReactants_ids_internal)

            lProducts_ids_cofactors = []
            lProducts_ids_internal = []
            for lProdotti in lProducts_ids:
                if len(gL.intersect(lcofactors, lProdotti)) != 0:
                    lProducts_ids_cofactors += lProdotti
                else:
                    lProducts_ids_internal += lProdotti
            lProducts_ids_cofactors = gL.unique(lProducts_ids_cofactors)
            lProducts_ids_internal = gL.unique(lProducts_ids_internal)
            dfAllDBs_copy = dfAllDBs_explode.copy()
            if isTransport == True:
                dfAllDBs_copy = dfAllDBs_copy[dfAllDBs_copy['reaction_type'] == '|TRANSPORT|']
            else:
                dfAllDBs_copy = dfAllDBs_copy[dfAllDBs_copy['reaction_type'] != '|TRANSPORT|']

            lIdxs = []
            if len(lReactants_ids) in dgb_left_metacyc and len(lProducts_ids) in dgb_right_metacyc:
                lCommonIdxs = gL.intersect(dgb_left_metacyc[len(lReactants_ids)].tolist(), dgb_right_metacyc[len(lProducts_ids)].tolist())
                lIdxs += lCommonIdxs

            if len(lReactants_ids) in dgb_right_metacyc and len(lProducts_ids) in dgb_left_metacyc:
                lCommonIdxs = gL.intersect(dgb_right_metacyc[len(lReactants_ids)].tolist(), dgb_left_metacyc[len(lProducts_ids)].tolist())
                lIdxs += lCommonIdxs

            if len(lReactants_ids) in dgb_lSubs_fromKegg and len(lProducts_ids) in dgb_lProds_fromKegg:
                lCommonIdxs = gL.intersect(dgb_lSubs_fromKegg[len(lReactants_ids)].tolist(), dgb_lProds_fromKegg[len(lProducts_ids)].tolist())
                lIdxs += lCommonIdxs

            if len(lReactants_ids) in dgb_lProds_fromKegg and len(lProducts_ids) in dgb_lSubs_fromKegg:
                lCommonIdxs = gL.intersect(dgb_lProds_fromKegg[len(lReactants_ids)].tolist(), dgb_lSubs_fromKegg[len(lProducts_ids)].tolist())
                lIdxs += lCommonIdxs

            if len(lReactants_ids) in dgb_lReactants_name_id_flat and len(lProducts_ids) in dgb_lProducts_name_id_flat:
                lCommonIdxs = gL.intersect(dgb_lReactants_name_id_flat[len(lReactants_ids)].tolist(), dgb_lProducts_name_id_flat[len(lProducts_ids)].tolist())
                lIdxs += lCommonIdxs

            if len(lReactants_ids) in dgb_lProducts_name_id_flat and len(lProducts_ids) in dgb_lReactants_name_id_flat:
                lCommonIdxs = gL.intersect(dgb_lProducts_name_id_flat[len(lReactants_ids)].tolist(), dgb_lReactants_name_id_flat[len(lProducts_ids)].tolist())
                lIdxs += lCommonIdxs

            lIdxs = gL.unique(lIdxs)
            lDfIdx = dfAllDBs_copy.index.tolist()
            dfAllDBs_copy_filter = dfAllDBs_copy.loc[gL.intersect(lIdxs, lDfIdx)]

            # Query MetaCyc
            lDfs2Concat = rxnL.filterRows(dfAllDBs_copy_filter, 'left_metacyc', 'right_metacyc', lReactants_ids_cofactors, lReactants_ids_internal, lProducts_ids_cofactors, lProducts_ids_internal)
            if len(lDfs2Concat) != 0:
                lPutativeRxns = rxnL.findRxnsAfterFilter(lDfs2Concat, 'left_metacyc', 'right_metacyc', 'metacyc', dfR, dfP, lReactants_ids, lProducts_ids)
                lIdentifiersRxn += lPutativeRxns

            dfAllDBs_equalSP_filter_M_copy = dfAllDBs_equalSP_filter_M.copy()
            lPutativeRxns = rxnL.findPutativeRxns(dfAllDBs_equalSP_filter_M_copy, 'left_woEqualMet', 'right', 'metacyc')
            lIdentifiersRxn += lPutativeRxns
            dfAllDBs_equalSP_filter_M_copy = dfAllDBs_equalSP_filter_M.copy()
            lPutativeRxns = rxnL.findPutativeRxns(dfAllDBs_equalSP_filter_M_copy, 'left', 'right_woEqualMet', 'metacyc')
            lIdentifiersRxn += lPutativeRxns

            dfAllDBs_equalSP_filter_K_copy = dfAllDBs_equalSP_filter_K.copy()
            lPutativeRxns = rxnL.findPutativeRxns(dfAllDBs_equalSP_filter_K_copy, 'left_woEqualMet', 'lProds_fromKegg', 'metacyc')
            lIdentifiersRxn += lPutativeRxns
            dfAllDBs_equalSP_filter_K_copy = dfAllDBs_equalSP_filter_K.copy()
            lPutativeRxns = rxnL.findPutativeRxns(dfAllDBs_equalSP_filter_K_copy, 'lSubs_fromKegg', 'right_woEqualMet', 'metacyc')
            lIdentifiersRxn += lPutativeRxns

            dfAllDBs_equalSP_filter_R_copy = dfAllDBs_equalSP_filter_R.copy()
            lPutativeRxns = rxnL.findPutativeRxns(dfAllDBs_equalSP_filter_R_copy, 'left_woEqualMet', 'lProducts_name_id_flat', 'metacyc')
            lIdentifiersRxn += lPutativeRxns
            dfAllDBs_equalSP_filter_R_copy = dfAllDBs_equalSP_filter_R.copy()
            lPutativeRxns = rxnL.findPutativeRxns(dfAllDBs_equalSP_filter_R_copy, 'lReactants_name_id_flat', 'right_woEqualMet', 'metacyc')
            lIdentifiersRxn += lPutativeRxns

            # Query KEGG
            lDfs2Concat = rxnL.filterRows(dfAllDBs_copy_filter, 'lSubs_fromKegg', 'lProds_fromKegg', lReactants_ids_cofactors, lReactants_ids_internal, lProducts_ids_cofactors, lProducts_ids_internal)
            if len(lDfs2Concat) != 0:
                lPutativeRxns = rxnL.findRxnsAfterFilter(lDfs2Concat, 'lSubs_fromKegg', 'lProds_fromKegg', 'kegg', dfR, dfP, lReactants_ids, lProducts_ids)
                lIdentifiersRxn += lPutativeRxns

            dfAllDBs_equalSP_filter_M_copy = dfAllDBs_equalSP_filter_M.copy()
            lPutativeRxns = rxnL.findPutativeRxns(dfAllDBs_equalSP_filter_M_copy, 'left_woEqualMet', 'right', 'kegg')
            lIdentifiersRxn += lPutativeRxns
            dfAllDBs_equalSP_filter_M_copy = dfAllDBs_equalSP_filter_M.copy()
            lPutativeRxns = rxnL.findPutativeRxns(dfAllDBs_equalSP_filter_M_copy, 'left', 'right_woEqualMet', 'kegg')
            lIdentifiersRxn += lPutativeRxns

            dfAllDBs_equalSP_filter_K_copy = dfAllDBs_equalSP_filter_K.copy()
            lPutativeRxns = rxnL.findPutativeRxns(dfAllDBs_equalSP_filter_K_copy, 'left_woEqualMet', 'lProds_fromKegg', 'kegg')
            lIdentifiersRxn += lPutativeRxns
            dfAllDBs_equalSP_filter_K_copy = dfAllDBs_equalSP_filter_K.copy()
            lPutativeRxns = rxnL.findPutativeRxns(dfAllDBs_equalSP_filter_K_copy, 'lSubs_fromKegg', 'right_woEqualMet', 'kegg')
            lIdentifiersRxn += lPutativeRxns

            dfAllDBs_equalSP_filter_R_copy = dfAllDBs_equalSP_filter_R.copy()
            lPutativeRxns = rxnL.findPutativeRxns(dfAllDBs_equalSP_filter_R_copy, 'left_woEqualMet', 'lProducts_name_id_flat', 'kegg')
            lIdentifiersRxn += lPutativeRxns
            dfAllDBs_equalSP_filter_R_copy = dfAllDBs_equalSP_filter_R.copy()
            lPutativeRxns = rxnL.findPutativeRxns(dfAllDBs_equalSP_filter_R_copy, 'lReactants_name_id_flat', 'right_woEqualMet', 'kegg')
            lIdentifiersRxn += lPutativeRxns

            # Query Rhea
            lDfs2Concat = rxnL.filterRows(dfAllDBs_copy_filter, 'lReactants_name_id_flat', 'lProducts_name_id_flat', lReactants_ids_cofactors, lReactants_ids_internal, lProducts_ids_cofactors, lProducts_ids_internal)
            if len(lDfs2Concat) != 0:
                lPutativeRxns = rxnL.findRxnsAfterFilter(lDfs2Concat, 'lReactants_name_id_flat', 'lProducts_name_id_flat', 'rhea', dfR, dfP, lReactants_ids, lProducts_ids)
                lIdentifiersRxn += lPutativeRxns

            dfAllDBs_equalSP_filter_M_copy = dfAllDBs_equalSP_filter_M.copy()
            lPutativeRxns = rxnL.findPutativeRxns(dfAllDBs_equalSP_filter_M_copy, 'left_woEqualMet', 'right', 'rhea')
            lIdentifiersRxn += lPutativeRxns
            dfAllDBs_equalSP_filter_M_copy = dfAllDBs_equalSP_filter_M.copy()
            lPutativeRxns = rxnL.findPutativeRxns(dfAllDBs_equalSP_filter_M_copy, 'left', 'right_woEqualMet', 'rhea')
            lIdentifiersRxn += lPutativeRxns

            dfAllDBs_equalSP_filter_K_copy = dfAllDBs_equalSP_filter_K.copy()
            lPutativeRxns = rxnL.findPutativeRxns(dfAllDBs_equalSP_filter_K_copy, 'left_woEqualMet', 'lProds_fromKegg', 'rhea')
            lIdentifiersRxn += lPutativeRxns
            dfAllDBs_equalSP_filter_K_copy = dfAllDBs_equalSP_filter_K.copy()
            lPutativeRxns = rxnL.findPutativeRxns(dfAllDBs_equalSP_filter_K_copy, 'lSubs_fromKegg', 'right_woEqualMet', 'rhea')
            lIdentifiersRxn += lPutativeRxns

            dfAllDBs_equalSP_filter_R_copy = dfAllDBs_equalSP_filter_R.copy()
            lPutativeRxns = rxnL.findPutativeRxns(dfAllDBs_equalSP_filter_R_copy, 'left_woEqualMet', 'lProducts_name_id_flat', 'rhea')
            lIdentifiersRxn += lPutativeRxns
            dfAllDBs_equalSP_filter_R_copy = dfAllDBs_equalSP_filter_R.copy()
            lPutativeRxns = rxnL.findPutativeRxns(dfAllDBs_equalSP_filter_R_copy, 'lReactants_name_id_flat', 'right_woEqualMet', 'rhea')
            lIdentifiersRxn += lPutativeRxns

        lIdentifiersRxn = gL.unique(lIdentifiersRxn)
        gL.writeLineByLineToFile(outFile, [rxn.id, lIdentifiersRxn], '\t')

outFile.close()

dfRxnWId = pd.read_csv(os.path.join(OUTDIR, dfrxnsInfo + '_wIds.csv'), sep = '\t', dtype=str)
dfMerge = pd.merge(dfRxns, dfRxnWId, on = 'RxnId')
dfMerge.to_csv(os.path.join(OUTDIR, dfrxnsInfo + '_enriched.csv'), sep = '\t', index = False)
