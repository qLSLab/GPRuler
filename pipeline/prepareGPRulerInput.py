import sys
import genericLib as gL
import os
import pandas as pd
import RESTmoduleModified as RESTmod

# setting working dirs
workingDirs = gL.setWorkingDirs()
OUTDIR = workingDirs[2]

# setting input data
testModel = sys.argv[1]
if testModel == 'recon3':
    ## Recon 3
    rxnswGenesFileName = 'recon3D_genes2Compartments_wFilter'
    outputFileName = 'Recon3D'
    organismCode = 'hsa'
elif testModel == 'y7':
    ## Yeast 7
    rxnswGenesFileName = 'y7_genes2Compartments_wFilter'
    outputFileName = 'Yeast7'
    organismCode = 'sce'
elif testModel == 'y8':
    ## Yeast 8
    rxnswGenesFileName = 'y8_genes2Compartments_wFilter'
    outputFileName = 'Yeast8'
    organismCode = 'sce'
elif testModel == 'hmr':
    ## HMRcore
    rxnswGenesFileName = 'hmrCore_genes2Compartments_wFilter'
    outputFileName = 'HMRcore'
    organismCode = 'hsa'
elif testModel == 'ownData':
    ## specify your input data
    rxnswGenesFileName = 'hmrCore_genes2Compartments_wFilter'
    outputFileName = 'HMRcore'
    organismCode = ''

# Generate the first output: reaction --> list of catalysing genes
dfRxnswGenes = pd.read_csv(os.path.join(OUTDIR, rxnswGenesFileName + '.csv'), sep = '\t', dtype=str)
dfRxnswGenes_filter = dfRxnswGenes[['Rxn', 'lGenes_filtered']]
dfRxnswGenes_filter = dfRxnswGenes_filter.rename(columns={'lGenes_filtered': 'Genes'})
dfRxnswGenes_filter.to_csv(os.path.join(OUTDIR, outputFileName + '_Rxns2Genes.csv'), sep = '\t', index = False)

# Generate the second output: KEGG gene identifier --> corresponding Uniprot identifier
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
