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
# from libchebipy._chebi_entity import ChebiEntity
import chebiLib as cL
import time
# import libchebipy
import re
import ast
import RESTmoduleModified as RESTmod
import bioservices.kegg as kegg
import cobra as cb

def checkNadNadpDependencies_and(nadp, nad, id2Search, id2Search_complete, lAnd, dGenesFromKegg):
    if nadp == True:
        if id2Search in dGenesFromKegg:
            dGeneK = dGenesFromKegg[id2Search]
        else:
            dGeneK = getKeggInfo(id2Search_complete)
            dGenesFromKegg[id2Search] = dGeneK
        if dGeneK != 400 and dGeneK != 404 and 'DEFINITION' in dGeneK:
            if '(NAD(+))' not in dGeneK['DEFINITION'] and '(NAD+)' not in dGeneK['DEFINITION'] and '(NADP(+))' not in dGeneK['DEFINITION'] and '(NADP+)' not in dGeneK['DEFINITION']:
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            elif ('(NAD(+))' in dGeneK['DEFINITION'] or '(NAD+)' in dGeneK['DEFINITION']) and ('(NADP(+))' in dGeneK['DEFINITION'] or '(NADP+)' in dGeneK['DEFINITION']):
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            elif ('(NAD(+))' not in dGeneK['DEFINITION'] and '(NAD+)' not in dGeneK['DEFINITION']) and ('(NADP(+))' in dGeneK['DEFINITION'] or '(NADP+)' in dGeneK['DEFINITION']):
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            else:
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
        if dGeneK != 400 and dGeneK != 404 and 'ORTHOLOGY' in dGeneK:
            if '(NAD(+))' not in dGeneK['ORTHOLOGY'] and '(NAD+)' not in dGeneK['ORTHOLOGY'] and '(NADP(+))' not in dGeneK['ORTHOLOGY'] and '(NADP+)' not in dGeneK['ORTHOLOGY']:
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            elif ('(NAD(+))' in dGeneK['ORTHOLOGY'] or '(NAD+)' in dGeneK['ORTHOLOGY']) and ('(NADP(+))' in dGeneK['ORTHOLOGY'] or '(NADP+)' in dGeneK['ORTHOLOGY']):
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            elif ('(NAD(+))' not in dGeneK['ORTHOLOGY'] and '(NAD+)' not in dGeneK['ORTHOLOGY']) and ('(NADP(+))' in dGeneK['ORTHOLOGY'] or '(NADP+)' in dGeneK['ORTHOLOGY']):
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            else:
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
    if nad == True:
        if id2Search in dGenesFromKegg:
            dGeneK = dGenesFromKegg[id2Search]
        else:
            dGeneK = getKeggInfo(id2Search_complete)
            dGenesFromKegg[id2Search] = dGeneK
        if dGeneK != 400 and dGeneK != 404 and 'DEFINITION' in dGeneK:
            if '(NAD(+))' not in dGeneK['DEFINITION'] and '(NAD+)' not in dGeneK['DEFINITION'] and '(NADP(+))' not in dGeneK['DEFINITION'] and '(NADP+)' not in dGeneK['DEFINITION']:
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            elif ('(NAD(+))' in dGeneK['DEFINITION'] or '(NAD+)' in dGeneK['DEFINITION']) and ('(NADP(+))' in dGeneK['DEFINITION'] or '(NADP+)' in dGeneK['DEFINITION']):
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            elif ('(NAD(+))' in dGeneK['DEFINITION'] or '(NAD+)' in dGeneK['DEFINITION']) and ('(NADP(+))' not in dGeneK['DEFINITION'] and '(NADP+)' not in dGeneK['DEFINITION']):
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            else:
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
        if dGeneK != 400 and dGeneK != 404 and 'ORTHOLOGY' in dGeneK:
            if '(NAD(+))' not in dGeneK['ORTHOLOGY'] and '(NAD+)' not in dGeneK['ORTHOLOGY'] and '(NADP(+))' not in dGeneK['ORTHOLOGY'] and '(NADP+)' not in dGeneK['ORTHOLOGY']:
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            elif ('(NAD(+))' in dGeneK['ORTHOLOGY'] or '(NAD+)' in dGeneK['ORTHOLOGY']) and ('(NADP(+))' in dGeneK['ORTHOLOGY'] or '(NADP+)' in dGeneK['ORTHOLOGY']):
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            elif ('(NAD(+))' in dGeneK['ORTHOLOGY'] or '(NAD+)' in dGeneK['ORTHOLOGY']) and ('(NADP(+))' not in dGeneK['ORTHOLOGY'] and '(NADP+)' not in dGeneK['ORTHOLOGY']):
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
            else:
                if id2Search not in lAnd:
                    lAnd.append(id2Search)
    if nadp == False and nad == False:
        if id2Search not in lAnd:
            lAnd.append(id2Search)
    return lAnd, dGenesFromKegg

def checkNadNadpDependencies_or(nadp, nad, id2Search, id2Search_complete, lMetaEnzOR, dGenesFromKegg):
    if nadp == True:
        if id2Search in dGenesFromKegg:
            dGeneK = dGenesFromKegg[id2Search]
        else:
            dGeneK = getKeggInfo(id2Search_complete)
            dGenesFromKegg[id2Search] = dGeneK
        if dGeneK != 400 and dGeneK != 404 and 'DEFINITION' in dGeneK:
            if '(NAD(+))' not in dGeneK['DEFINITION'] and '(NAD+)' not in dGeneK['DEFINITION'] and '(NADP(+))' not in dGeneK['DEFINITION'] and '(NADP+)' not in dGeneK['DEFINITION']:
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            elif ('(NAD(+))' in dGeneK['DEFINITION'] or '(NAD+)' in dGeneK['DEFINITION']) and ('(NADP(+))' in dGeneK['DEFINITION'] or '(NADP+)' in dGeneK['DEFINITION']):
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            elif ('(NAD(+))' not in dGeneK['DEFINITION'] and '(NAD+)' not in dGeneK['DEFINITION']) and ('(NADP(+))' in dGeneK['DEFINITION'] or '(NADP+)' in dGeneK['DEFINITION']):
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            else:
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
        if dGeneK != 400 and dGeneK != 404 and 'ORTHOLOGY' in dGeneK:
            if '(NAD(+))' not in dGeneK['ORTHOLOGY'] and '(NAD+)' not in dGeneK['ORTHOLOGY'] and '(NADP(+))' not in dGeneK['ORTHOLOGY'] and '(NADP+)' not in dGeneK['ORTHOLOGY']:
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            elif ('(NAD(+))' in dGeneK['ORTHOLOGY'] or '(NAD+)' in dGeneK['ORTHOLOGY']) and ('(NADP(+))' in dGeneK['ORTHOLOGY'] or '(NADP+)' in dGeneK['ORTHOLOGY']):
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            elif ('(NAD(+))' not in dGeneK['ORTHOLOGY'] and '(NAD+)' not in dGeneK['ORTHOLOGY']) and ('(NADP(+))' in dGeneK['ORTHOLOGY'] or '(NADP+)' in dGeneK['ORTHOLOGY']):
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            else:
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
    if nad == True:
        if id2Search in dGenesFromKegg:
            dGeneK = dGenesFromKegg[id2Search]
        else:
            dGeneK = getKeggInfo(id2Search_complete)
            dGenesFromKegg[id2Search] = dGeneK
        if dGeneK != 400 and dGeneK != 404 and 'DEFINITION' in dGeneK:
            if '(NAD(+))' not in dGeneK['DEFINITION'] and '(NAD+)' not in dGeneK['DEFINITION'] and '(NADP(+))' not in dGeneK['DEFINITION'] and '(NADP+)' not in dGeneK['DEFINITION']:
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            elif ('(NAD(+))' in dGeneK['DEFINITION'] or '(NAD+)' in dGeneK['DEFINITION']) and ('(NADP(+))' in dGeneK['DEFINITION'] or '(NADP+)' in dGeneK['DEFINITION']):
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            elif ('(NAD(+))' in dGeneK['DEFINITION'] or '(NAD+)' in dGeneK['DEFINITION']) and ('(NADP(+))' not in dGeneK['DEFINITION'] and '(NADP+)' not in dGeneK['DEFINITION']):
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            else:
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
        if dGeneK != 400 and dGeneK != 404 and 'ORTHOLOGY' in dGeneK:
            if '(NAD(+))' not in dGeneK['ORTHOLOGY'] and '(NAD+)' not in dGeneK['ORTHOLOGY'] and '(NADP(+))' not in dGeneK['ORTHOLOGY'] and '(NADP+)' not in dGeneK['ORTHOLOGY']:
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            elif ('(NAD(+))' in dGeneK['ORTHOLOGY'] or '(NAD+)' in dGeneK['ORTHOLOGY']) and ('(NADP(+))' in dGeneK['ORTHOLOGY'] or '(NADP+)' in dGeneK['ORTHOLOGY']):
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            elif ('(NAD(+))' in dGeneK['ORTHOLOGY'] or '(NAD+)' in dGeneK['ORTHOLOGY']) and ('(NADP(+))' not in dGeneK['ORTHOLOGY'] and '(NADP+)' not in dGeneK['ORTHOLOGY']):
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
            else:
                if id2Search not in lMetaEnzOR:
                    lMetaEnzOR.append([id2Search])
    if nadp == False and nad == False:
        if id2Search not in lMetaEnzOR:
            lMetaEnzOR.append([id2Search])
    return lMetaEnzOR, dGenesFromKegg

def getKeggInfo(rxnString):
    try:
        k = kegg.KEGG()
        info = k.get(rxnString)
        dizInfo = k.parse(info)
    except:
        dizInfo = {}
    return dizInfo

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]
FIGUREDIR = workingDirs[4]
DFDIR = workingDirs[5]

## INPUTS:
# - metabolites list

start = time.time()

timeStamp = gL.getTimeStamp()
print('timeStamp\t', timeStamp)

testModel = sys.argv[1]
if testModel == 'recon3':
    ## Recon 3
    modelXml = 'Recon3D_301_20200923'
    rxns = 'recon3D_reactions_enriched_20210208132643'
    transportRxns = 'recon3D_reactions_enriched_tcdb_20210122084205'
    outputName = 'recon3D_reactions_wGenes_'
    metModelFile = 'recon3D_metabolites_20201218172131_enriched'
    orgCode = 'hsa'
    taxId = ['|TAX-9606|']
elif testModel == 'y7':
    ## Yeast 7
    modelXml = 'yeast_7.6_cobra'
    rxns = 'y7_reactions_enriched_20210208130139'
    transportRxns = 'y7_reactions_enriched_tcdb_20210201134957'
    outputName = 'y7_reactions_wGenes_'
    metModelFile = 'y7_metabolites_20201219101924_enriched'
    orgCode = 'sce'
    taxId = ['|TAX-4932|', '|TAX-559292|', '|TAX-580239|', '|TAX-658763|', '|TAX-1294310|']
elif testModel == 'y8':
    ## Yeast 8
    modelXml = 'yeast8'
    rxns = 'y8_reactions_enriched_20210208130145'
    transportRxns = 'y8_reactions_enriched_tcdb_20210201135027'
    outputName = 'y8_reactions_wGenes_'
    metModelFile = 'y8_metabolites_20201219102046_enriched'
    orgCode = 'sce'
    taxId = ['|TAX-4932|', '|TAX-559292|', '|TAX-580239|', '|TAX-658763|', '|TAX-1294310|']
elif testModel == 'hmr':
    ## HMRcore
    modelXml = 'HMRcore_20200328_wReconNames'
    rxns = 'hmrCore_reactions_enriched_20210208093612'
    transportRxns = 'hmrCore_reactions_enriched_tcdb_20210118143656'
    outputName = 'hmrCore_reactions_wGenes_'
    metModelFile = 'hmrCore_metabolites_20210112092737_enriched'
    orgCode = 'hsa'
    taxId = ['|TAX-9606|']



start = time.time()

dfAllDBs = pd.read_csv(os.path.join(OUTDIR, 'dfJoin_metacyc_kegg_rhea_20201218164915.csv'), sep = '\t', dtype=str) ## 3
dfAllDBs['ec_number'] = dfAllDBs['ec_number'].apply(literal_eval)
dfAllDBs['enzymes_of_reaction'] = dfAllDBs['enzymes_of_reaction'].apply(literal_eval)
dfAllDBs['OtherRheaId_fromKegg'] = dfAllDBs['OtherRheaId_fromKegg'].apply(literal_eval)
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
dfAllDBs['ECnumber'] = dfAllDBs['ECnumber'].fillna({i: '[]' for i in dfAllDBs.index})
dfAllDBs['ECnumber'] = dfAllDBs['ECnumber'].apply(literal_eval)
dfAllDBs['UniprotId'] = dfAllDBs['UniprotId'].fillna({i: '[]' for i in dfAllDBs.index})
dfAllDBs['UniprotId'] = dfAllDBs['UniprotId'].apply(literal_eval)



dfMetacyc_proteins = pd.read_csv(os.path.join(OUTDIR, 'metacyc_proteins_20201216152513.csv'), sep = '\t', dtype=str)
dfMetacyc_proteins['lSpecies'] = dfMetacyc_proteins['lSpecies'].apply(literal_eval)
dfMetacyc_proteins['lGenes'] = dfMetacyc_proteins['lGenes'].apply(literal_eval)
dfMetacyc_proteins['Components'] = dfMetacyc_proteins['Components'].apply(literal_eval)
dfMetacyc_proteins['Component_of'] = dfMetacyc_proteins['Component_of'].apply(literal_eval)
dfMetacyc_proteins = dfMetacyc_proteins.explode('Component_of')

dfMetacyc_genes = pd.read_csv(os.path.join(OUTDIR, 'metacyc_genes_20201216152513.csv'), sep = '\t', dtype=str)
dfMetacyc_genes['Accessions'] = dfMetacyc_genes['Accessions'].apply(literal_eval)
dfMetacyc_genes['Names'] = dfMetacyc_genes['Names'].apply(literal_eval)
dfMetacyc_genes['Product'] = dfMetacyc_genes['Product'].apply(literal_eval)
dfMetacyc_genes = dfMetacyc_genes.explode('Product')
dfMetacyc_genes = dfMetacyc_genes.reset_index(drop = True)

dfrxns = pd.read_csv(os.path.join(OUTDIR, rxns + '.csv'), sep = '\t', dtype = str)
dfrxns['PutativeIdentifiers'] = dfrxns['PutativeIdentifiers'].apply(literal_eval)
dfrxns['IsTransport'] = dfrxns['IsTransport'].apply(literal_eval)
dfrxns['IsExchange'] = dfrxns['IsExchange'].apply(literal_eval)
print('dfrxns\t', dfrxns.shape)
dftransportRxns = pd.read_csv(os.path.join(OUTDIR, transportRxns + '.csv'), sep = '\t')
print('dftransportRxns\t', dftransportRxns.shape)

dfmerge = pd.merge(dfrxns, dftransportRxns, on='Rxn')
print('dfmerge\t', dfmerge.shape)
# dfmerge.to_csv(os.path.join(OUTDIR, 'merge_' + orgCode + '.csv'), sep='\t', index = False)


## conversione kegg a ncbi gene id
ncbi2Org = RESTmod.kegg_conv(orgCode, "ncbi-geneid").readlines()
ncbi = []
gene = []
for line in ncbi2Org:
    separo = line.strip().split('\t')
    ncbi.append(separo[0])
    gene.append(separo[1])

dfncbi2Org = pd.DataFrame({'ncbi': ncbi, 'keggGeneId': gene})
# dfncbi2Org.to_csv(os.path.join(OUTDIR, 'ConversionKegg2NcbiGene_' + orgCode + '.csv'), sep='\t', index = False)

## conversione kegg a uniprot
uniprot2Org = RESTmod.kegg_conv(orgCode, "uniprot").readlines()
unip = []
gene = []
for line in uniprot2Org:
    separo = line.strip().split('\t')
    unip.append(separo[0])
    gene.append(separo[1])

dfuniprot2Org = pd.DataFrame({'uniprot': unip, 'keggGeneId': gene})
# dfuniprot2Org.to_csv(os.path.join(OUTDIR, 'ConversionKegg2Uniprot_' + orgCode + '.csv'), sep='\t', index = False)


dfMetEnriched = pd.read_csv(os.path.join(OUTDIR, metModelFile + '.csv'), sep = '\t', dtype=str)
dfMetEnriched['lIdentifiers'] = dfMetEnriched['lIdentifiers'].apply(literal_eval)

model = cb.io.read_sbml_model(os.path.join(RAWDIR, modelXml + '.xml'))

# lModelRxns = []
# for reactions in model.reactions:
#     lModelRxns.append(reactions.id)
#
# print('lModelRxns\t', len(lModelRxns))


lMetaEnzOR_all = []
lEc_all = []


if testModel == 'recon3' or testModel == 'hmr':
    dfmerge['Rxn_conv'] = dfmerge.Rxn.str[2:]
else:
    dfmerge['Rxn_conv'] = dfmerge.Rxn.values


# print(gL.difference(lModelRxns, list(dfmerge['Rxn_conv'])))
# dfmerge_filtered = dfmerge[dfmerge['Rxn_conv'].isin(lModelRxns)]

dGenesFromKegg = {}
dEcFromKegg = {}
dOrthFromKegg = {}

for row in dfmerge.itertuples():
    print('Rxn\t', row.Rxn, '\t', 'convert\t', row.Rxn_conv)
    nadp = False ## flag per indicare se la rxn e' nad o nadp dipendente e quindi devo verificare se i geni che estraggo per la reazione sono nad o nadp dipendenti
    nad = False
    rxnInModel = model.reactions.get_by_id(row.Rxn_conv)
    print(rxnInModel.id, '\n', rxnInModel.reactants)
    # if testModel == 'recon3' or testModel == 'hmr':
    #     rxnInModel = model.reactions.get_by_id(row.Rxn[2:])
    #     print(rxnInModel.id, '\n', rxnInModel.reactants)
    # else:
    #     rxnInModel = model.reactions.get_by_id(row.Rxn)
    #     print(rxnInModel.id, '\n', rxnInModel.reactants)

    ###########################################################
    ## per sapere se rxn è nad o nadp dipendente
    for reactant in rxnInModel.reactants:
        print('reactant\t', reactant)
        if testModel == 'recon3' or testModel == 'hmr':
            metSearch = dfMetEnriched[dfMetEnriched['Id'] == 'M_' + reactant.id]

        else:
            metSearch = dfMetEnriched[dfMetEnriched['Id'] == reactant.id]
        print('metSearch\n', metSearch)
        if metSearch.empty is False:
            for metRow in metSearch.itertuples():
                if 'C00005' in metRow.lIdentifiers or 'C00006' in metRow.lIdentifiers:
                    nadp = True
                elif 'C00004' in metRow.lIdentifiers or 'C00003' in metRow.lIdentifiers:
                    nad = True
    ###########################################################

    lEc = []
    lMetaEnzOR = []
    if (row.IsTransport_x == True or row.IsExchange_x == True) or (row.IsTransport_x == False and row.IsExchange_x == False):
        if len(row.PutativeIdentifiers) != 0:
            print('row.PutativeIdentifiers\t', row.PutativeIdentifiers)

            m_m = dfAllDBs[dfAllDBs['MetaCycId'].isin(row.PutativeIdentifiers)]
            # print('m_m\n', m_m)
            m_r = dfAllDBs[dfAllDBs['MetaCycId_fromRhea'].isin(row.PutativeIdentifiers)]
            # print('m_r\n', m_r)
            k_m = dfAllDBs[dfAllDBs['KeggId_x'].isin(row.PutativeIdentifiers)]
            # print('k_m\n', k_m)
            k_k = dfAllDBs[dfAllDBs['KeggId_fromKegg'].isin(row.PutativeIdentifiers)]
            # print('k_k\n', k_k)

            lDfs = [m_m, m_r, k_m, k_k]

            lRhea = ['RheaId_master', 'RheaId_lr','RheaId_rl','RheaId_bi', 'OtherRheaId_fromKegg']
            for r in lRhea:
                tmp = dfAllDBs.explode(r)
                tmp = tmp.reset_index(drop = True)
                # print('tmp\n', tmp[tmp[r].isin(row.PutativeIdentifiers)])
                lDfs.append(tmp[tmp[r].isin(row.PutativeIdentifiers)])
            r_m = dfAllDBs[dfAllDBs['RheaId'].isin(row.PutativeIdentifiers)]
            # print('r_m\n', r_m)
            lDfs.append(r_m)
            dfSearch = pd.concat(lDfs)
            dfSearch = dfSearch.reset_index(drop = True)
            # dfSearch = dfSearch.drop_duplicates(subset = ['MetaCycId']) ## tolgo questa riga di codice altrimenti cancella righe in cui metacyc id non c'e' ma ci sono altre info che ci servono
            # dfSearch = dfSearch.reset_index(drop = True)

            print('ECCC\n', dfSearch['ec_number'])
            for lec in list(dfSearch['ec_number'].dropna()):
                # print('Leccccccc\t', lec)
                for ec in lec:
                    # print('ec_number\t', ec[4:-1])
                    lEc.append(ec[4:-1])
            print('lista completa lEc\t', lEc)

            allmetacycEnzymes = list(dfSearch['enzymes_of_reaction'].dropna()) ## verranno legati da OR e ciascuno lo sostituisco con la lista dei geni che devono stare tra loro in AND se annotati nell'organismo target
            print('allmetacycEnzymes\n', allmetacycEnzymes)
            for l in allmetacycEnzymes:
                print('l\t', l)
                for el in l: ## el puo' essere un complesso o un monomero
                    print('el\t', el)
                    lAnd = []
                    dfProt1 = dfMetacyc_proteins[dfMetacyc_proteins['MetaCycId'] == el]
                    dfProt2 = dfMetacyc_proteins[dfMetacyc_proteins['Component_of'] == el]
                    dfProt = pd.concat([dfProt1, dfProt2])
                    # dfProt = dfProt.reset_index(drop = True)
                    # print('dfProt\n', dfProt)
                    # print('species\n', dfProt.iloc[0]['lSpecies'])
                    # if dfProt.empty is False and len(gL.intersect(dfProt.iloc[0]['lSpecies'], taxId)) != 0:
                    if dfProt.empty is False:
                        for rowdfProt in dfProt.itertuples():
                            print('rowdfProt.lSpecies\t', rowdfProt.lSpecies)
                            print('rowdfProt.lGenes\t', rowdfProt.lGenes)
                            print('Component_of\t', rowdfProt.Component_of)
                            if len(gL.intersect(rowdfProt.lSpecies, taxId)) != 0:
                                print('UNIPROT\t', rowdfProt.Uniprot, '\t', pd.isna(rowdfProt.Uniprot))
                                if pd.isna(rowdfProt.Uniprot) is False:
                                    dfCorrespondingGenes = dfuniprot2Org[dfuniprot2Org['uniprot'] == 'up:' + rowdfProt.Uniprot]
                                    print('dfCorrespondingGenes\n', dfCorrespondingGenes)
                                    for foundGene in list(dfCorrespondingGenes['keggGeneId']):
                                        print('foundGene\t', foundGene.split(':')[1])
                                        if [foundGene.split(':')[1]] not in lMetaEnzOR:
                                            # lMetaEnzOR.append([foundGene.split(':')[1]])
                                            geneId2search = foundGene.split(':')[1]
                                            print('geneId2search\t', geneId2search)
                                            print('nadp\t', nadp)
                                            print('nad\t', nad)
                                            lMetaEnzOR, dGenesFromKegg = checkNadNadpDependencies_or(nadp, nad, geneId2search, orgCode + ':' + geneId2search, lMetaEnzOR, dGenesFromKegg)


                                lgenes = rowdfProt.lGenes
                                print('lgenes\t', lgenes)
                                lComponents = rowdfProt.Components
                                print('lComponents\t', lComponents)
                                foundGenes = dfMetacyc_genes[dfMetacyc_genes['MetaCycId'].isin(lgenes)]
                                foundComponents = dfMetacyc_genes[dfMetacyc_genes['Product'].isin(lComponents)]
                                found = pd.concat([foundGenes,foundComponents])
                                found = found.drop_duplicates(subset = ['MetaCycId'])
                                # print('found\n', found)
                                print('MetaCycId\t', list(found['MetaCycId']))
                                print('NcbiId\n',list(found['NcbiId']))
                                print('Accessions\n',list(found['Accessions']))
                                print('Uniprot\n',list(found['Uniprot'].dropna()))

                                for g in list(found['NcbiId'].dropna()):
                                    print('g\t', g)
                                    if testModel  == '1':
                                        gId = orgCode + ':' + g
                                        print('gId\t', gId)
                                        if gId in list(dfncbi2Org['keggGeneId']) and g not in lAnd:
                                            print('g found in kegg conv in keggGeneId')
                                            print('nadp\t', nadp)
                                            print('nad\t', nad)
                                            lAnd, dGenesFromKegg = checkNadNadpDependencies_and(nadp, nad, g, gId, lAnd, dGenesFromKegg)
                                    else:
                                        gId = 'ncbi-geneid:' + g
                                        print('Other gId\t', gId)
                                        if gId in list(dfncbi2Org['ncbi']):
                                            searchNcbi = dfncbi2Org[dfncbi2Org['ncbi'] == gId]
                                            print('searchNcbi\n', searchNcbi['keggGeneId'])
                                            for s in list(searchNcbi['keggGeneId']):
                                                print('found\t', s.split(':')[1])
                                                if s.split(':')[1] not in lAnd:

                                                    geneId2search = s.split(':')[1]
                                                    print('geneId2search\t', geneId2search)
                                                    print('nadp\t', nadp)
                                                    print('nad\t', nad)
                                                    lAnd, dGenesFromKegg = checkNadNadpDependencies_and(nadp, nad, geneId2search, orgCode + ':' + geneId2search, lAnd, dGenesFromKegg)

                                for lAcc in list(found['Accessions'].dropna()):
                                    print('lAcc\t', lAcc)
                                    for a in lAcc:
                                        gId = orgCode + ':' + a
                                        print('gId\t', gId)
                                        if gId in list(dfncbi2Org['keggGeneId']) and a not in lAnd:
                                            print('gId found in kegg conv in keggGeneId')
                                            print('nadp\t', nadp)
                                            print('nad\t', nad)
                                            lAnd, dGenesFromKegg = checkNadNadpDependencies_and(nadp, nad, a, gId, lAnd, dGenesFromKegg)

                                for uniprotFound in list(found['Uniprot'].dropna()):
                                    print('uniprotFound\t', uniprotFound)
                                    dfCorrespondingGenes = dfuniprot2Org[dfuniprot2Org['uniprot'] == 'up:' + uniprotFound]
                                    print('dfCorrespondingGenes\n', dfCorrespondingGenes)
                                    for foundGene in list(dfCorrespondingGenes['keggGeneId']):
                                        print('foundGene\t', foundGene.split(':')[1])
                                        if [foundGene.split(':')[1]] not in lMetaEnzOR:
                                            # lMetaEnzOR.append([foundGene.split(':')[1]])
                                            geneId2search = foundGene.split(':')[1]
                                            print('geneId2search\t', geneId2search)
                                            print('nadp\t', nadp)
                                            print('nad\t', nad)
                                            lMetaEnzOR,dGenesFromKegg = checkNadNadpDependencies_or(nadp, nad, geneId2search, orgCode + ':' + geneId2search, lMetaEnzOR,dGenesFromKegg)
                    if len(lAnd) != 0:
                        lAnd.sort()
                        if lAnd not in lMetaEnzOR:
                            lMetaEnzOR.append(lAnd)

            lOrths = list(dfSearch['Orthology_fromKegg'].dropna())
            for orth in lOrths:
                # print('orth\t', orth)
                lOId = list(gL.extractRegexFromItem(orth, r"([K0-9]{6})")[0])
                print('lOId\t', lOId)
                for oId in lOId:
                    print('oId\t', oId)
                    if oId in dOrthFromKegg:
                        dOrth = dOrthFromKegg[oId]
                    else:
                        dOrth = getKeggInfo('ko:' + oId)
                        dOrthFromKegg[oId] = dOrth
                    if dOrth != 400 and dOrth != 404 and 'GENES' in dOrth and orgCode.upper() in dOrth['GENES']:
                        for item in dOrth['GENES'][orgCode.upper()].split():
                            par = item.find('(')
                            if par != -1:
                                print('item[par]\t', item[:par])
                                if [item[:par]] not in lMetaEnzOR:
                                    # lMetaEnzOR.append([item[:par]])
                                    geneId2search = item[:par]
                                    print('geneId2search\t', geneId2search)
                                    print('nadp\t', nadp)
                                    print('nad\t', nad)
                                    lMetaEnzOR,dGenesFromKegg = checkNadNadpDependencies_or(nadp, nad, geneId2search, orgCode + ':' + geneId2search, lMetaEnzOR,dGenesFromKegg)
                            else:
                                if [item] not in lMetaEnzOR:
                                    # lMetaEnzOR.append([item])
                                    print('geneId2search\t', item)
                                    print('nadp\t', nadp)
                                    print('nad\t', nad)
                                    lMetaEnzOR,dGenesFromKegg = checkNadNadpDependencies_or(nadp, nad, item, orgCode + ':' + item, lMetaEnzOR, dGenesFromKegg)

            for el in list(dfSearch['EC_fromKegg'].dropna()):
                print('EC_fromKegg\t',el.split())
                lEc += el.split()

            for llecRhea in list(dfSearch['ECnumber'].dropna()):
                print('llecRhea\t', llecRhea)
                for lecRhea in llecRhea:
                    print('lecRhea\t', lecRhea)
                    lEc += lecRhea

            for llUnipRhea in list(dfSearch['UniprotId'].dropna()):
                # print('llUnipRhea\t', llUnipRhea)
                for lUnipRhea in llUnipRhea:
                    print('lUnipRhea\t', lUnipRhea)
                    # for el in lUnipRhea:
                    #     print('up in df\t', 'up:' + el in list(dfuniprot2Org['uniprot']))
                    if any('up:' + el in list(dfuniprot2Org['uniprot']) for el in lUnipRhea) is True:
                        dfCorrespondingGenes = dfuniprot2Org[dfuniprot2Org['uniprot'].isin(['up:' + el for el in lUnipRhea])]
                        print('dfCorrespondingGenes\n', dfCorrespondingGenes)
                        for foundGene in list(dfCorrespondingGenes['keggGeneId']):
                            print('foundGene\t', foundGene.split(':')[1])
                            if [foundGene.split(':')[1]] not in lMetaEnzOR:
                                # lMetaEnzOR.append([foundGene.split(':')[1]])

                                geneId2search = foundGene.split(':')[1]
                                print('geneId2search\t', geneId2search)
                                print('nadp\t', nadp)
                                print('nad\t', nad)
                                lMetaEnzOR,dGenesFromKegg = checkNadNadpDependencies_or(nadp, nad, geneId2search, orgCode + ':' + geneId2search, lMetaEnzOR, dGenesFromKegg)

    if (row.IsTransport_x == True or row.IsExchange_x == True):
        # print('row.Identifiers_fromTCDB\t', row.Identifiers_fromTCDB)
        dTC2Unip = ast.literal_eval(row.Identifiers_fromTCDB)
        # print('dTC2Unip\t', dTC2Unip)
        if len(dTC2Unip) != 0:
            for tc in dTC2Unip:
                lUnips = dTC2Unip[tc]
                print('tc\t', tc)
                print('lUnips\t', lUnips)
                if any('up:' + el in list(dfuniprot2Org['uniprot']) for el in lUnips) is True:
                    dfCorrespondingGenes = dfuniprot2Org[dfuniprot2Org['uniprot'].isin(['up:' + el for el in lUnips])]
                    print('dfCorrespondingGenes\n', dfCorrespondingGenes)
                    for foundGene in list(dfCorrespondingGenes['keggGeneId']):
                        print('foundGene\t', foundGene.split(':')[1])
                        if [foundGene.split(':')[1]] not in lMetaEnzOR:
                            # lMetaEnzOR.append([foundGene.split(':')[1]])

                            geneId2search = foundGene.split(':')[1]
                            print('geneId2search\t', geneId2search)
                            print('nadp\t', nadp)
                            print('nad\t', nad)
                            lMetaEnzOR, dGenesFromKegg = checkNadNadpDependencies_or(nadp, nad, geneId2search, orgCode + ':' + geneId2search, lMetaEnzOR, dGenesFromKegg)


    lEc = gL.unique(lEc)
    print('lEc FINALE\t', lEc, '\t')
    # if len(lMetaEnzOR) == 0:
    for ec in lEc:
        print('ec\t', ec)
        if ec in dEcFromKegg:
            dEcs = dEcFromKegg[ec]
        else:
            dEcs = getKeggInfo('ec:' + ec)
            dEcFromKegg[ec] = dEcs
        if dEcs != 400 and dEcs != 404 and 'GENES' in dEcs and orgCode.upper() in dEcs['GENES']:
            for item in dEcs['GENES'][orgCode.upper()].split():
                par = item.find('(')
                if par != -1:
                    if [item[:par]] not in lMetaEnzOR:
                        # lMetaEnzOR.append([item[:par]])
                        # print('item\t', item[:par])

                        geneId2search = item[:par]
                        print('geneId2search\t', geneId2search)
                        print('nadp\t', nadp)
                        print('nad\t', nad)
                        lMetaEnzOR, dGenesFromKegg = checkNadNadpDependencies_or(nadp, nad, geneId2search, orgCode + ':' + geneId2search, lMetaEnzOR, dGenesFromKegg)
                else:
                    if [item] not in lMetaEnzOR:
                        # lMetaEnzOR.append([item])
                        # print('item\t', item)

                        print('geneId2search\t', item)
                        print('nadp\t', nadp)
                        print('nad\t', nad)
                        lMetaEnzOR, dGenesFromKegg = checkNadNadpDependencies_or(nadp, nad, item, orgCode + ':' + item, lMetaEnzOR, dGenesFromKegg)

    print('lMetaEnzOR\t', lMetaEnzOR)
    lMetaEnzOR_all.append(lMetaEnzOR)
    lEc_all.append(lEc)
    print('\n\n')


dfmerge['lGenes'] = lMetaEnzOR_all
dfmerge['lEC'] = lEc_all
dfmerge.to_csv(os.path.join(OUTDIR, outputName + timeStamp + '.csv'), sep = '\t', index = False)
end = time.time()
print('elapsed time\t', end-start, '\tseconds')
