# -*- coding: utf-8 -*-
import os
import sys
import pandas as pd
import genericLib as gL
import reactionsLib as rL
import genesLib as genesL
from ast import literal_eval
import RESTmoduleModified as RESTmod

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]

# setting input data
dfRxns2GenesFile = ''
dfRxnId2Equation = ''
orgCode = ''
taxId = []
modelName = ''

dfAllDBs = pd.read_csv(os.path.join(RAWDIR, 'dfJoin_metacyc_kegg_rhea_20201218164915.csv'), sep = '\t', dtype=str)
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


dfMetacyc_proteins = pd.read_csv(os.path.join(RAWDIR, 'metacyc_proteins_20201216152513.csv'), sep = '\t', dtype=str)
dfMetacyc_proteins['lSpecies'] = dfMetacyc_proteins['lSpecies'].apply(literal_eval)
dfMetacyc_proteins['lGenes'] = dfMetacyc_proteins['lGenes'].apply(literal_eval)
dfMetacyc_proteins['Components'] = dfMetacyc_proteins['Components'].apply(literal_eval)
dfMetacyc_proteins['Component_of'] = dfMetacyc_proteins['Component_of'].apply(literal_eval)
dfMetacyc_proteins = dfMetacyc_proteins.explode('Component_of')

dfMetacyc_genes = pd.read_csv(os.path.join(RAWDIR, 'metacyc_genes_20201216152513.csv'), sep = '\t', dtype=str)
dfMetacyc_genes['Accessions'] = dfMetacyc_genes['Accessions'].apply(literal_eval)
dfMetacyc_genes['Names'] = dfMetacyc_genes['Names'].apply(literal_eval)
dfMetacyc_genes['Product'] = dfMetacyc_genes['Product'].apply(literal_eval)
dfMetacyc_genes = dfMetacyc_genes.explode('Product')
dfMetacyc_genes = dfMetacyc_genes.reset_index(drop = True)

# Conversion of KEGG gene identifiers to NCBI gene identifiers
ncbi2Org = RESTmod.kegg_conv(orgCode, "ncbi-geneid").readlines()
ncbi = []
gene = []
for line in ncbi2Org:
    separo = line.strip().split('\t')
    ncbi.append(separo[0])
    gene.append(separo[1])

dfncbi2Org = pd.DataFrame({'ncbi': ncbi, 'keggGeneId': gene})

# Conversion of KEGG gene identifiers to Uniprot identifiers
uniprot2Org = RESTmod.kegg_conv(orgCode, "uniprot").readlines()
unip = []
gene = []
for line in uniprot2Org:
    separo = line.strip().split('\t')
    unip.append(separo[0])
    gene.append(separo[1])

dfuniprot2Org = pd.DataFrame({'uniprot': unip, 'keggGeneId': gene})

lMetaEnzOR_all = []
lEc_all = []
dGenesFromKegg = {}
dEcFromKegg = {}
dOrthFromKegg = {}

## Loop over each reaction to add genes from the macrodatabase
dfRxns2Genes = pd.read_csv(os.path.join(OUTDIR, dfRxns2GenesFile + '.csv'), sep = '\t', dtype=str)
dfRxns2Genes['Genes'] = dfRxns2Genes['Genes'].apply(literal_eval)

dfRxnId2Equation = pd.read_csv(os.path.join(OUTDIR, dfRxnId2EquationFile + '.csv'), sep = '\t', dtype=str)

for rxn in dfRxns2Genes.itertuples():
    nadp = False
    nad = False
    rxnTarget = dfRxnId2Equation[dfRxnId2Equation['RxnId'] == rxn]
    equation = rxnTarget.iloc[0]['Equation']
    # determine if the reaction is NADH or NADPH dependent to filter genes accordingly
    if 'C00005' in equation or 'C00006' in equation:
        nadp = True
    elif 'C00004' in equation or 'C00003' in equation:
        nad = True

    dfSearch = pd.concat(dfAllDBs[dfAllDBs['KeggId_x'] == rxn], dfAllDBs[dfAllDBs['KeggId_fromKegg'] == rxn], dfAllDBs[dfAllDBs['KeggId_y'] == rxn])

    lEc = []
    lMetaEnzOR = []

    dfSearch = dfSearch.reset_index(drop = True)
    for lec in list(dfSearch['ec_number'].dropna()):
        for ec in lec:
            lEc.append(ec[4:-1])
    allmetacycEnzymes = list(dfSearch['enzymes_of_reaction'].dropna())
    for l in allmetacycEnzymes:
        for el in l:
            lAnd = []
            dfProt1 = dfMetacyc_proteins[dfMetacyc_proteins['MetaCycId'] == el]
            dfProt2 = dfMetacyc_proteins[dfMetacyc_proteins['Component_of'] == el]
            dfProt = pd.concat([dfProt1, dfProt2])
            if dfProt.empty is False:
                for rowdfProt in dfProt.itertuples():
                    if len(gL.intersect(rowdfProt.lSpecies, taxId)) != 0:
                        if pd.isna(rowdfProt.Uniprot) is False:
                            dfCorrespondingGenes = dfuniprot2Org[dfuniprot2Org['uniprot'] == 'up:' + rowdfProt.Uniprot]
                            for foundGene in list(dfCorrespondingGenes['keggGeneId']):
                                if [foundGene.split(':')[1]] not in lMetaEnzOR:
                                    geneId2search = foundGene.split(':')[1]
                                    lMetaEnzOR, dGenesFromKegg = rxnL.checkNadNadpDependencies_or(nadp, nad, geneId2search, orgCode + ':' + geneId2search, lMetaEnzOR, dGenesFromKegg)
                        lgenes = rowdfProt.lGenes
                        lComponents = rowdfProt.Components
                        foundGenes = dfMetacyc_genes[dfMetacyc_genes['MetaCycId'].isin(lgenes)]
                        foundComponents = dfMetacyc_genes[dfMetacyc_genes['Product'].isin(lComponents)]
                        found = pd.concat([foundGenes,foundComponents])
                        found = found.drop_duplicates(subset = ['MetaCycId'])
                        for g in list(found['NcbiId'].dropna()):
                            if testModel  == '1':
                                gId = orgCode + ':' + g
                                if gId in list(dfncbi2Org['keggGeneId']) and g not in lAnd:
                                    lAnd, dGenesFromKegg = rxnL.checkNadNadpDependencies_and(nadp, nad, g, gId, lAnd, dGenesFromKegg)
                            else:
                                gId = 'ncbi-geneid:' + g
                                if gId in list(dfncbi2Org['ncbi']):
                                    searchNcbi = dfncbi2Org[dfncbi2Org['ncbi'] == gId]
                                    for s in list(searchNcbi['keggGeneId']):
                                        if s.split(':')[1] not in lAnd:
                                            geneId2search = s.split(':')[1]
                                            lAnd, dGenesFromKegg = rxnL.checkNadNadpDependencies_and(nadp, nad, geneId2search, orgCode + ':' + geneId2search, lAnd, dGenesFromKegg)

                        for lAcc in list(found['Accessions'].dropna()):
                            for a in lAcc:
                                gId = orgCode + ':' + a
                                if gId in list(dfncbi2Org['keggGeneId']) and a not in lAnd:
                                    lAnd, dGenesFromKegg = rxnL.checkNadNadpDependencies_and(nadp, nad, a, gId, lAnd, dGenesFromKegg)

                        for uniprotFound in list(found['Uniprot'].dropna()):
                            dfCorrespondingGenes = dfuniprot2Org[dfuniprot2Org['uniprot'] == 'up:' + uniprotFound]
                            for foundGene in list(dfCorrespondingGenes['keggGeneId']):
                                if [foundGene.split(':')[1]] not in lMetaEnzOR:
                                    geneId2search = foundGene.split(':')[1]
                                    lMetaEnzOR,dGenesFromKegg = rxnL.checkNadNadpDependencies_or(nadp, nad, geneId2search, orgCode + ':' + geneId2search, lMetaEnzOR,dGenesFromKegg)
            if len(lAnd) != 0:
                lAnd.sort()
                if lAnd not in lMetaEnzOR:
                    lMetaEnzOR.append(lAnd)

    lOrths = list(dfSearch['Orthology_fromKegg'].dropna())
    for orth in lOrths:
        lOId = list(gL.extractRegexFromItem(orth, r"([K0-9]{6})")[0])
        for oId in lOId:
            if oId in dOrthFromKegg:
                dOrth = dOrthFromKegg[oId]
            else:
                dOrth = rxnL.getKeggInfo('ko:' + oId)
                dOrthFromKegg[oId] = dOrth
            if dOrth != 400 and dOrth != 404 and 'GENES' in dOrth and orgCode.upper() in dOrth['GENES']:
                for item in dOrth['GENES'][orgCode.upper()].split():
                    par = item.find('(')
                    if par != -1:
                        if [item[:par]] not in lMetaEnzOR:
                            geneId2search = item[:par]
                            lMetaEnzOR,dGenesFromKegg = rxnL.checkNadNadpDependencies_or(nadp, nad, geneId2search, orgCode + ':' + geneId2search, lMetaEnzOR,dGenesFromKegg)
                    else:
                        if [item] not in lMetaEnzOR:
                            lMetaEnzOR,dGenesFromKegg = rxnL.checkNadNadpDependencies_or(nadp, nad, item, orgCode + ':' + item, lMetaEnzOR, dGenesFromKegg)

    for el in list(dfSearch['EC_fromKegg'].dropna()):
        lEc += el.split()

    for llecRhea in list(dfSearch['ECnumber'].dropna()):
        for lecRhea in llecRhea:
            lEc += lecRhea

    for llUnipRhea in list(dfSearch['UniprotId'].dropna()):
        for lUnipRhea in llUnipRhea:
            if any('up:' + el in list(dfuniprot2Org['uniprot']) for el in lUnipRhea) is True:
                dfCorrespondingGenes = dfuniprot2Org[dfuniprot2Org['uniprot'].isin(['up:' + el for el in lUnipRhea])]
                for foundGene in list(dfCorrespondingGenes['keggGeneId']):
                    if [foundGene.split(':')[1]] not in lMetaEnzOR:
                        geneId2search = foundGene.split(':')[1]
                        lMetaEnzOR,dGenesFromKegg = rxnL.checkNadNadpDependencies_or(nadp, nad, geneId2search, orgCode + ':' + geneId2search, lMetaEnzOR, dGenesFromKegg)

    lEc = gL.unique(lEc)
    for ec in lEc:
        if ec in dEcFromKegg:
            dEcs = dEcFromKegg[ec]
        else:
            dEcs = rxnL.getKeggInfo('ec:' + ec)
            dEcFromKegg[ec] = dEcs
        if dEcs != 400 and dEcs != 404 and 'GENES' in dEcs and orgCode.upper() in dEcs['GENES']:
            for item in dEcs['GENES'][orgCode.upper()].split():
                par = item.find('(')
                if par != -1:
                    if [item[:par]] not in lMetaEnzOR:
                        geneId2search = item[:par]
                        lMetaEnzOR, dGenesFromKegg = rxnL.checkNadNadpDependencies_or(nadp, nad, geneId2search, orgCode + ':' + geneId2search, lMetaEnzOR, dGenesFromKegg)
                else:
                    if [item] not in lMetaEnzOR:
                        lMetaEnzOR, dGenesFromKegg = rxnL.checkNadNadpDependencies_or(nadp, nad, item, orgCode + ':' + item, lMetaEnzOR, dGenesFromKegg)

    lMetaEnzOR_all.append(lMetaEnzOR)
    lEc_all.append(lEc)


dfRxns2Genes['Genes_fromMacroDb'] = lMetaEnzOR_all
dfRxns2Genes['lEC'] = lEc_all
dfRxns2Genes['Genes_fromKEGG'] = dfRxns2Genes.Genes.tolist()
dfRxns2Genes['Genes'] = dfRxns2Genes['Genes_fromKEGG'] + dfRxns2Genes['Genes_fromMacroDb']
dfRxns2Genes.to_csv(os.path.join(OUTDIR, modelName + '_Rxns2Genes.csv'), sep = '\t', index = False)

# Convert the list of all KEGG genes identifiers of the target organism to their corresponding Uniprot identifiers
genesL.kegg2UniprotGenesId(orgCode, modelName, OUTDIR)
