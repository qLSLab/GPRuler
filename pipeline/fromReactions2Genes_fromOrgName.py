import sys
import genericLib as gL
import os
import pandas as pd
from ast import literal_eval
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

# setting input data
testModel = sys.argv[1]
if testModel == 'recon3':
    ## Recon 3
    modelXml = 'Recon3D_301_20200923'
    rxns = 'recon3D_reactions_enriched'
    transportRxns = 'recon3D_reactions_enriched_tcdb'
    outputName = 'recon3D_reactions_wGenes'
    metModelFile = 'recon3D_metabolites_enriched'
    orgCode = 'hsa'
    taxId = ['|TAX-9606|']
elif testModel == 'y7':
    ## Yeast 7
    modelXml = 'yeast_7.6_cobra'
    rxns = 'y7_reactions_enriched'
    transportRxns = 'y7_reactions_enriched_tcdb'
    outputName = 'y7_reactions_wGenes'
    metModelFile = 'y7_metabolites_enriched'
    orgCode = 'sce'
    taxId = ['|TAX-4932|', '|TAX-559292|', '|TAX-580239|', '|TAX-658763|', '|TAX-1294310|']
elif testModel == 'y8':
    ## Yeast 8
    modelXml = 'yeast8'
    rxns = 'y8_reactions_enriched'
    transportRxns = 'y8_reactions_enriched_tcdb'
    outputName = 'y8_reactions_wGenes'
    metModelFile = 'y8_metabolites_enriched'
    orgCode = 'sce'
    taxId = ['|TAX-4932|', '|TAX-559292|', '|TAX-580239|', '|TAX-658763|', '|TAX-1294310|']
elif testModel == 'hmr':
    ## HMRcore
    modelXml = 'HMRcore_20200328_wReconNames'
    rxns = 'hmrCore_reactions_enriched'
    transportRxns = 'hmrCore_reactions_enriched_tcdb'
    outputName = 'hmrCore_reactions_wGenes'
    metModelFile = 'hmrCore_metabolites_enriched'
    orgCode = 'hsa'
    taxId = ['|TAX-9606|']
elif testModel == 'ownData':
    ## specify your input data
    modelXml = ''
    rxns = ''
    transportRxns = ''
    outputName = ''
    metModelFile = ''
    orgCode = ''
    taxId = []
    ####
    dfRxns2GenesFile

dfAllDBs = pd.read_csv(os.path.join(OUTDIR, 'dfJoin_metacyc_kegg_rhea_20201218164915.csv'), sep = '\t', dtype=str)
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

# dfrxns = pd.read_csv(os.path.join(OUTDIR, rxns + '.csv'), sep = '\t', dtype = str)
# dfrxns['PutativeIdentifiers'] = dfrxns['PutativeIdentifiers'].apply(literal_eval)
# dfrxns['IsTransport'] = dfrxns['IsTransport'].apply(literal_eval)
# dfrxns['IsExchange'] = dfrxns['IsExchange'].apply(literal_eval)
# print('dfrxns\t', dfrxns.shape)
# dftransportRxns = pd.read_csv(os.path.join(OUTDIR, transportRxns + '.csv'), sep = '\t')
# print('dftransportRxns\t', dftransportRxns.shape)
#
# dfmerge = pd.merge(dfrxns, dftransportRxns, on='Rxn')

## ciclo su ogni reazione per agg geni dal macrodb
dfRxns2Genes = pd.read_csv(os.path.join(OUTDIR, dfRxns2GenesFile + '.csv'), sep = '\t', dtype=str)
dfRxns2Genes['Genes'] = dfRxns2Genes['Genes'].apply(literal_eval)

for rxn in dfRxns2Genes.itertuples():
    dfAllDBs[dfAllDBs['KeggId_x'] == rxn]
    dfAllDBs[dfAllDBs['KeggId_fromKegg'] == rxn]
    dfAllDBs[dfAllDBs['KeggId_y'] == rxn]

##############


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

dfMetEnriched = pd.read_csv(os.path.join(OUTDIR, metModelFile + '.csv'), sep = '\t', dtype=str)
dfMetEnriched['lIdentifiers'] = dfMetEnriched['lIdentifiers'].apply(literal_eval)

model = cb.io.read_sbml_model(os.path.join(RAWDIR, modelXml + '.xml'))

lMetaEnzOR_all = []
lEc_all = []
if testModel == 'recon3' or testModel == 'hmr':
    dfmerge['Rxn_conv'] = dfmerge.Rxn.str[2:]
else:
    dfmerge['Rxn_conv'] = dfmerge.Rxn.values

dGenesFromKegg = {}
dEcFromKegg = {}
dOrthFromKegg = {}

for row in dfmerge.itertuples():
    nadp = False
    nad = False
    rxnInModel = model.reactions.get_by_id(row.Rxn_conv)
    # determine if the reaction is NADH or NADPH dependent to filter genes accordingly
    for reactant in rxnInModel.reactants:
        if testModel == 'recon3' or testModel == 'hmr':
            metSearch = dfMetEnriched[dfMetEnriched['Id'] == 'M_' + reactant.id]
        else:
            metSearch = dfMetEnriched[dfMetEnriched['Id'] == reactant.id]
        if metSearch.empty is False:
            for metRow in metSearch.itertuples():
                if 'C00005' in metRow.lIdentifiers or 'C00006' in metRow.lIdentifiers:
                    nadp = True
                elif 'C00004' in metRow.lIdentifiers or 'C00003' in metRow.lIdentifiers:
                    nad = True
    lEc = []
    lMetaEnzOR = []
    if (row.IsTransport_x == True or row.IsExchange_x == True) or (row.IsTransport_x == False and row.IsExchange_x == False):
        if len(row.PutativeIdentifiers) != 0:
            m_m = dfAllDBs[dfAllDBs['MetaCycId'].isin(row.PutativeIdentifiers)]
            m_r = dfAllDBs[dfAllDBs['MetaCycId_fromRhea'].isin(row.PutativeIdentifiers)]
            k_m = dfAllDBs[dfAllDBs['KeggId_x'].isin(row.PutativeIdentifiers)]
            k_k = dfAllDBs[dfAllDBs['KeggId_fromKegg'].isin(row.PutativeIdentifiers)]
            lDfs = [m_m, m_r, k_m, k_k]
            lRhea = ['RheaId_master', 'RheaId_lr','RheaId_rl','RheaId_bi', 'OtherRheaId_fromKegg']
            for r in lRhea:
                tmp = dfAllDBs.explode(r)
                tmp = tmp.reset_index(drop = True)
                lDfs.append(tmp[tmp[r].isin(row.PutativeIdentifiers)])
            r_m = dfAllDBs[dfAllDBs['RheaId'].isin(row.PutativeIdentifiers)]
            lDfs.append(r_m)
            dfSearch = pd.concat(lDfs)
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
                                            lMetaEnzOR, dGenesFromKegg = checkNadNadpDependencies_or(nadp, nad, geneId2search, orgCode + ':' + geneId2search, lMetaEnzOR, dGenesFromKegg)
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
                                            lAnd, dGenesFromKegg = checkNadNadpDependencies_and(nadp, nad, g, gId, lAnd, dGenesFromKegg)
                                    else:
                                        gId = 'ncbi-geneid:' + g
                                        if gId in list(dfncbi2Org['ncbi']):
                                            searchNcbi = dfncbi2Org[dfncbi2Org['ncbi'] == gId]
                                            for s in list(searchNcbi['keggGeneId']):
                                                if s.split(':')[1] not in lAnd:
                                                    geneId2search = s.split(':')[1]
                                                    lAnd, dGenesFromKegg = checkNadNadpDependencies_and(nadp, nad, geneId2search, orgCode + ':' + geneId2search, lAnd, dGenesFromKegg)

                                for lAcc in list(found['Accessions'].dropna()):
                                    for a in lAcc:
                                        gId = orgCode + ':' + a
                                        if gId in list(dfncbi2Org['keggGeneId']) and a not in lAnd:
                                            lAnd, dGenesFromKegg = checkNadNadpDependencies_and(nadp, nad, a, gId, lAnd, dGenesFromKegg)

                                for uniprotFound in list(found['Uniprot'].dropna()):
                                    dfCorrespondingGenes = dfuniprot2Org[dfuniprot2Org['uniprot'] == 'up:' + uniprotFound]
                                    for foundGene in list(dfCorrespondingGenes['keggGeneId']):
                                        if [foundGene.split(':')[1]] not in lMetaEnzOR:
                                            geneId2search = foundGene.split(':')[1]
                                            lMetaEnzOR,dGenesFromKegg = checkNadNadpDependencies_or(nadp, nad, geneId2search, orgCode + ':' + geneId2search, lMetaEnzOR,dGenesFromKegg)
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
                        dOrth = getKeggInfo('ko:' + oId)
                        dOrthFromKegg[oId] = dOrth
                    if dOrth != 400 and dOrth != 404 and 'GENES' in dOrth and orgCode.upper() in dOrth['GENES']:
                        for item in dOrth['GENES'][orgCode.upper()].split():
                            par = item.find('(')
                            if par != -1:
                                if [item[:par]] not in lMetaEnzOR:
                                    geneId2search = item[:par]
                                    lMetaEnzOR,dGenesFromKegg = checkNadNadpDependencies_or(nadp, nad, geneId2search, orgCode + ':' + geneId2search, lMetaEnzOR,dGenesFromKegg)
                            else:
                                if [item] not in lMetaEnzOR:
                                    lMetaEnzOR,dGenesFromKegg = checkNadNadpDependencies_or(nadp, nad, item, orgCode + ':' + item, lMetaEnzOR, dGenesFromKegg)

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
                                lMetaEnzOR,dGenesFromKegg = checkNadNadpDependencies_or(nadp, nad, geneId2search, orgCode + ':' + geneId2search, lMetaEnzOR, dGenesFromKegg)

    if (row.IsTransport_x == True or row.IsExchange_x == True):
        dTC2Unip = literal_eval(row.Identifiers_fromTCDB)
        if len(dTC2Unip) != 0:
            for tc in dTC2Unip:
                lUnips = dTC2Unip[tc]
                if any('up:' + el in list(dfuniprot2Org['uniprot']) for el in lUnips) is True:
                    dfCorrespondingGenes = dfuniprot2Org[dfuniprot2Org['uniprot'].isin(['up:' + el for el in lUnips])]
                    for foundGene in list(dfCorrespondingGenes['keggGeneId']):
                        if [foundGene.split(':')[1]] not in lMetaEnzOR:
                            geneId2search = foundGene.split(':')[1]
                            lMetaEnzOR, dGenesFromKegg = checkNadNadpDependencies_or(nadp, nad, geneId2search, orgCode + ':' + geneId2search, lMetaEnzOR, dGenesFromKegg)

    lEc = gL.unique(lEc)
    for ec in lEc:
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
                        geneId2search = item[:par]
                        lMetaEnzOR, dGenesFromKegg = checkNadNadpDependencies_or(nadp, nad, geneId2search, orgCode + ':' + geneId2search, lMetaEnzOR, dGenesFromKegg)
                else:
                    if [item] not in lMetaEnzOR:
                        lMetaEnzOR, dGenesFromKegg = checkNadNadpDependencies_or(nadp, nad, item, orgCode + ':' + item, lMetaEnzOR, dGenesFromKegg)

    lMetaEnzOR_all.append(lMetaEnzOR)
    lEc_all.append(lEc)

dfmerge['lGenes'] = lMetaEnzOR_all
dfmerge['lEC'] = lEc_all
dfmerge.to_csv(os.path.join(OUTDIR, outputName + '.csv'), sep = '\t', index = False)
