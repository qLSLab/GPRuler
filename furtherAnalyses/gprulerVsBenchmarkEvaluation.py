import sys
import genericLib as gL
import os
import pandas as pd
import re
import numpy as np
import itertools
from sklearn.metrics.pairwise import paired_distances
from scipy.spatial import distance

import matplotlib.pyplot as plt
from matplotlib import font_manager
import matplotlib

def jaccard(list1, list2):
    intersection = len(gL.intersect(list1, list2))
    union = len(gL.union(list1, list2))
    return intersection / union

def generateTruthTable(lGenes1, lGenes2):
    '''
    Generate the truth matrix
    '''
    n1 = len(lGenes1)
    n2 = len(lGenes2)
    lGenes = gL.unique(lGenes1 + lGenes2)
    n = len(lGenes)
    mTruth = np.array([i for i in itertools.product([False, True], repeat=n)])
    return mTruth, lGenes, n


def evaluateBoolExpressions(mTruth, lGenes, rule):
    '''
    Evaluate a given GPR rule for each row of the truth matrix
    '''
    aRes = np.empty([mTruth.shape[0], 1])
    for j in range(0, mTruth.shape[0]):
        row = mTruth[j, :]
        tmp = rule[:]
        dGeneToBoolValue = {}
        i = 0
        while i < len(row):
            dGeneToBoolValue[lGenes[i]] = row[i]
            i += 1
        for g, bool in dGeneToBoolValue.items():
            g_find = tmp.find(g)
            badPositions = []
            while g_find != -1 and g_find not in badPositions:
                if g_find > 0 and (tmp[g_find - 1].isdigit() == False and tmp[g_find - 1].isalpha() == False and tmp[g_find - 1] != '-') \
                        and g_find + len(g) < len(tmp) and (tmp[g_find + len(g)].isdigit() == False and tmp[g_find + len(g)].isalpha() == False and tmp[g_find + len(g)] != '-'):
                    tmp = tmp[:g_find] + str(bool) + tmp[g_find + len(g):]
                elif g_find == 0 and g_find + len(g) < len(tmp) and (tmp[g_find + len(g)].isdigit() == False and tmp[g_find + len(g)].isalpha() == False and tmp[g_find + len(g)] != '-'):
                    tmp = tmp[:g_find] + str(bool) + tmp[g_find + len(g):]
                elif g_find > 0 and g_find + len(g) == len(tmp) and (tmp[g_find - 1].isdigit() == False and tmp[g_find - 1].isalpha() == False and tmp[g_find - 1] != '-'):
                    tmp = tmp[:g_find] + str(bool) + tmp[g_find + len(g):]
                elif g_find == 0 and g_find + len(g) == len(tmp):
                    tmp = str(bool)
                else:
                    badPositions.append(g_find)
                g_find = tmp.find(g)
        res = eval(tmp)
        aRes[j] = [res]
    return aRes


# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]
FIGUREDIR = workingDirs[5]

## load input file
testModel = sys.argv[1]

if testModel == 'hmr':
    ## HMRCore
    modelName = 'HMRcore'
    inputFileName = 'AdditionalFile2_HMRcore_GPRulerExecutionOutput_20210901'
    regexOrgSpecific = r"([0-9]+)"
elif testModel == 'recon3':
    ## Recon3
    inputFileName = 'AdditionalFile3_Recon3D_GPRulerExecutionOutput_20210901'
    regexOrgSpecific = r"([0-9]+)"
    modelName = 'Recon3D'
elif testModel == 'y7':
    ## Yeast7
    inputFileName = 'AdditionalFile4_Yeast7_GPRulerExecutionOutput_20210901'
    regexOrgSpecific = r"([A-Z0-9-]+)"
    modelName = 'Yeast 7'
elif testModel == 'y8':
    ## Yeast8
    inputFileName = 'AdditionalFile5_Yeast8_GPRulerExecutionOutput_20210901'
    regexOrgSpecific = r"([A-Z0-9-]+)"
    modelName = 'Yeast 8'

dfComparison = pd.read_csv(os.path.join(OUTDIR, inputFileName + '.csv'), sep="\t", dtype = str)


## Compute Jaccard score for each reaction and for the entire model
lAllGenes_original = []
lAllGenes_gpruler = []

lJScore_sngRxns = []
for row in dfComparison.itertuples():
    originalRule = row.rule_original
    gprulerRule = row.rule_GPRuler

    if (pd.isna(originalRule) == True or originalRule == ''):
        lGenes_originalRule = []
    else:
        dfGenes_originalRule = gL.extractRegexFromItem(originalRule, regexOrgSpecific)
        lGenes_originalRule = gL.unique(list(dfGenes_originalRule[0]))
        lGenes_originalRule.sort()

    if (pd.isna(gprulerRule) == True or gprulerRule == ''):
        lGenes_gprulerRule = []
    else:
        dfGenes_gprulerRule = gL.extractRegexFromItem(gprulerRule, regexOrgSpecific)
        lGenes_gprulerRule = gL.unique(list(dfGenes_gprulerRule[0]))
        lGenes_gprulerRule.sort()

    lAllGenes_original += lGenes_originalRule
    lAllGenes_gpruler += lGenes_gprulerRule

    if len(lGenes_originalRule) == 0 and len(lGenes_gprulerRule) == 0:
        jScore = 1.0
    elif (len(lGenes_originalRule) == 0 and len(lGenes_gprulerRule) != 0) or (len(lGenes_originalRule) != 0 and len(lGenes_gprulerRule) == 0):
        jScore = 0
    else:
        jScore = jaccard(lGenes_originalRule, lGenes_gprulerRule)
    lJScore_sngRxns.append(jScore)

dfComparison['JaccardScore'] = lJScore_sngRxns
dfComparison['JaccardScore'] = dfComparison['JaccardScore'].round(3)
dfComparison.to_csv(os.path.join(OUTDIR, inputFileName + '.csv'), sep="\t", index = False)

jScore_allModel = jaccard(lAllGenes_original, lAllGenes_gpruler)
print('jScore_allModel\t', jScore_allModel, '\n')


# Normalized Hamming distance on the output vectors of each truth matrix quantifies the extent of similarity in the Boolean logic combination
dfComparison = pd.read_csv(os.path.join(OUTDIR, inputFileName + '.csv'), sep="\t", dtype = {'Rxn':str, 'rule_original': str, 'rule_GPRuler': str, 'Evaluation': str, 'JaccardScore': float})

lDistances = []
lHamming = []

for row in dfComparison.itertuples():
    originalRule = row.rule_original
    gprulerRule = row.rule_GPRuler

    if (pd.isna(originalRule) == True or originalRule == '') and (pd.isna(gprulerRule) == True or gprulerRule == ''):
        distHamming = 0
        lHamming.append(distHamming)
    else:
        if row.JaccardScore == 1:
            dfGenes1 = gL.extractRegexFromItem(originalRule, regexOrgSpecific)
            lGenes1 = gL.unique(list(dfGenes1[0]))
            lGenes1.sort()
            dfGenes2 = gL.extractRegexFromItem(gprulerRule, regexOrgSpecific)
            lGenes2 = gL.unique(list(dfGenes2[0]))
            lGenes2.sort()

            maxGenesNumber = max(len(lGenes1), len(lGenes2))
            if maxGenesNumber >= 20:
                print('DISTANCE FOR MORE THAN 20 GENES')
                distHamming = np.nan
            else:
                mTruth, lGenes, n = generateTruthTable(lGenes1, lGenes2)
                try:
                    aRes1 = evaluateBoolExpressions(mTruth, lGenes, originalRule)
                    aRes2 = evaluateBoolExpressions(mTruth, lGenes, gprulerRule)
                    mTruth1 = np.append(mTruth, aRes1, axis=1)
                    mTruth2 = np.append(mTruth, aRes2, axis=1)
                    distHamming = distance.hamming([el for ll in aRes1.tolist() for el in ll], [el for ll in aRes2.tolist() for el in ll])
                    lHamming.append(distHamming)
                except:
                    print('Error in generation of truth matrix')
        else:
            distHamming = np.nan
    lDistances.append(distHamming)

dfComparison['HammingDistance'] = lDistances
dfComparison.to_csv(os.path.join(OUTDIR, inputFileName + '.csv'), sep="\t", index = False)


## Plot Jaccard index and normalized Hamming distribution
lJaccard = []
lHam = []
for row in dfComparison.itertuples():
    if row.Evaluation == 'Not automatically reconstructed by GPRuler' or row.Evaluation == 'Corrected by GPRuler':
        lJaccard.append(row.JaccardScore)
        if row.JaccardScore == 1:
            lHam.append(row.HammingDistance)

fig = plt.figure(figsize=(13,10))
ax = fig.add_subplot(111)
ax.hist(lJaccard, alpha=0.5, density=False, color = 'steelblue', weights=np.zeros_like(np.array(lJaccard)) + 1. / np.array(lJaccard).size)
ax.set_xlabel('Jaccard index',fontsize= 32)
ax.set_ylabel('Relative frequency',fontsize= 32)
ax.set_title(modelName,fontsize= 37)
plt.xticks(rotation=45, fontsize= 29)
plt.yticks(fontsize= 29)
plt.ylim(0, 1)
plt.savefig(os.path.join(FIGUREDIR, modelName + '_jaccardInAllMismatches.pdf'))
plt.close()

fig = plt.figure(figsize=(13,10))
ax = fig.add_subplot(111)
ax.hist(lHam, alpha=0.5, density=False, color = 'green', weights=np.zeros_like(np.array(lHam)) + 1. / np.array(lHam).size)
ax.set_xlabel('Hamming distance',fontsize= 32)
ax.set_ylabel('Relative frequency',fontsize= 32)
ax.set_title(modelName,fontsize= 37)
plt.xticks(rotation=45, fontsize= 29)
plt.yticks(fontsize= 29)
plt.ylim(0, 1)
plt.xlim(0, 1)
plt.savefig(os.path.join(FIGUREDIR, modelName + '_hammingInMismatchesWJaccard1.pdf'))
plt.close()
