from Bio.KEGG.REST import *
import sys
import genericLib as gL
import os
import pandas as pd
import time
import RESTmoduleModified as RESTmod

timeStamp = gL.getTimeStamp()
print('timeStamp\t', timeStamp)

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

testModel = sys.argv[1]
if testModel == 'recon3':
    ## Recon 3
    rxnswGenesFileName = 'recon3D_genes2Compartments_wFilter_20210209102550'
    outputFileName = 'Recon3D'
    organismCode = 'hsa'
elif testModel == 'y7':
    ## Yeast 7
    rxnswGenesFileName = 'y7_genes2Compartments_wFilter_20210208172543'
    outputFileName = 'Yeast7'
    organismCode = 'sce'
elif testModel == 'y8':
    ## Yeast 8
    rxnswGenesFileName = 'y8_genes2Compartments_wFilter_20210208172609'
    outputFileName = 'Yeast8'
    organismCode = 'sce'
elif testModel == 'hmr':
    ## HMRcore
    rxnswGenesFileName = 'hmrCore_genes2Compartments_wFilter_20210208133221'
    outputFileName = 'HMRcore'
    organismCode = 'hsa'


## primo file: rxn --> lista geni associata
dfRxnswGenes = pd.read_csv(os.path.join(OUTDIR, rxnswGenesFileName + '.csv'), sep = '\t', dtype=str)

dfRxnswGenes_filter = dfRxnswGenes[['Rxn', 'lGenes_filtered']]
# lRxns = []
# lGenes = []
# for row in dfRxnswGenes.itertuples():
#     lRxns.append(row.Rxn)
#     lGenes.append(row.lGenes_filtered)

# df = pd.DataFrame({'Rxn': lRxns, 'Genes': lGenes})
dfRxnswGenes_filter = dfRxnswGenes_filter.rename(columns={'lGenes_filtered': 'Genes'})
dfRxnswGenes_filter.to_csv(os.path.join(OUTDIR, outputFileName + '_Rxns2Genes.csv'), sep = '\t', index = False)


## secondo file: per ogni kegg gene dell'organismo target --> suo uniprot id
kegg2uniprot = RESTmod.kegg_conv('uniprot', organismCode)
dOrgKegg2Uniprot = {}
for el in kegg2uniprot.readlines():
    elSplt = el.strip().split('\t')
    if elSplt[0].split(':')[1] not in dOrgKegg2Uniprot.keys():
        dOrgKegg2Uniprot[elSplt[0].split(':')[1]] = [elSplt[1].split(':')[1]]
    else:
        dOrgKegg2Uniprot[elSplt[0].split(':')[1]] += [elSplt[1].split(':')[1]]

outFile = open(os.path.join(OUTDIR, outputFileName + '_Kegg2UniprotGenes.csv'), mode='w')
gL.writeLineByLineToFile(outFile, ['keggId', 'uniprotId'], '\t')

for k, v in dOrgKegg2Uniprot.items():
    for vv in v:
        gL.writeLineByLineToFile(outFile, [k, vv], '\t')
outFile.close()

end = time.time()
print('Elapsed time\t', end-start, '\tseconds')
