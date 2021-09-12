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
    n1 = len(lGenes1)
    #print('lGenes1\n', lGenes1, '\n', n1, '\n')
    n2 = len(lGenes2)
    #print('lGenes2\n', lGenes2, '\n', n2)

    lGenes = gL.unique(lGenes1 + lGenes2)
    n = len(lGenes)
    #print('n\t', n)
    #print('lGenes\n', lGenes)
    mTruth = np.array([i for i in itertools.product([False, True], repeat=n)])
    # print('mTruth shape\n', mTruth.shape[0], mTruth.shape[1])
    #print('mTruth\n', mTruth)
    return mTruth, lGenes, n


def evaluateBoolExpressions(mTruth, lGenes, rule):
    # print('Evaluate\n')
    # print('lGenes\t', lGenes)
    aRes = np.empty([mTruth.shape[0], 1])
    # print('aRes\n', aRes, '\n')
    for j in range(0, mTruth.shape[0]):  # qui valuto se la riga corrente della mTruth concorda con la regola contenuta in "rule"
        row = mTruth[j, :]
        # print('riga mTruth\t', j)
        # print('ROW\n', row, '\n')
        tmp = rule[:]
        # print('tmp\n', row, '\n')
        dGeneToBoolValue = {}
        i = 0
        while i < len(row):
            dGeneToBoolValue[lGenes[i]] = row[i]
            i += 1
        # print('dGeneToBoolValue\n', dGeneToBoolValue, '\n')

        for g, bool in dGeneToBoolValue.items():  # qua sostituisco i valori contenuti in dGeneToBoolValue nelle corrispondenti chiavi all'interno della regola che stiamo valutando
            g_find = tmp.find(g)
            badPositions = []
            while g_find != -1 and g_find not in badPositions:
                # if g_find > 0 and tmp[g_find-1].isdigit() == False and g_find + len(g) < len(tmp) and tmp[g_find + len(g)].isdigit() == False:
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
        #     print('tmp\n', tmp, '\n')
        # print('\nRULE \n', rule)
        # print('RULE BOOLEANA\n', tmp)

        res = eval(tmp)
        # print('res\n', res, '\n')
        aRes[j] = [res]
    #     print('aRes[j]\n', aRes[j], '\n')
    # print('aRes\n', aRes, '\n')
    return aRes


timeStamp = gL.getTimeStamp()

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]
FIGUREDIR = workingDirs[5]

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

dfComparison = pd.read_csv(os.path.join(OUTDIR, inputFileName + '.csv'), sep="\t", dtype = {'Rxn':str, 'rule_original': str, 'rule_GPRuler': str, 'Evaluation': str, 'JaccardScore': float})


# ## Compute Jaccard score for each reaction and for the entire model
# lAllGenes_original = []
# lAllGenes_gpruler = []
#
# lJScore_sngRxns = []
# for row in dfComparison.itertuples():
#     # print('Rxn:\t', row.Rxn)
#
#     originalRule = row.rule_original
#     gprulerRule = row.rule_GPRuler
#
#     if (pd.isna(originalRule) == True or originalRule == ''):
#         lGenes_originalRule = []
#     else:
#         dfGenes_originalRule = gL.extractRegexFromItem(originalRule, regexOrgSpecific)
#         lGenes_originalRule = gL.unique(list(dfGenes_originalRule[0]))
#         lGenes_originalRule.sort()
#
#     if (pd.isna(gprulerRule) == True or gprulerRule == ''):
#         lGenes_gprulerRule = []
#     else:
#         dfGenes_gprulerRule = gL.extractRegexFromItem(gprulerRule, regexOrgSpecific)
#         lGenes_gprulerRule = gL.unique(list(dfGenes_gprulerRule[0]))
#         lGenes_gprulerRule.sort()
#
#     # print('lGenes_originalRule\t', lGenes_originalRule)
#     # print('lGenes_gprulerRule\t', lGenes_gprulerRule)
#
#     lAllGenes_original += lGenes_originalRule
#     lAllGenes_gpruler += lGenes_gprulerRule
#
#     if len(lGenes_originalRule) == 0 and len(lGenes_gprulerRule) == 0:
#         jScore = 1.0
#     elif (len(lGenes_originalRule) == 0 and len(lGenes_gprulerRule) != 0) or (len(lGenes_originalRule) != 0 and len(lGenes_gprulerRule) == 0):
#         jScore = 0
#     else:
#         jScore = jaccard(lGenes_originalRule, lGenes_gprulerRule)
#     # print('jScore\t', jScore, '\n')
#     lJScore_sngRxns.append(jScore)
#
# dfComparison['JaccardScore'] = lJScore_sngRxns
# dfComparison['JaccardScore'] = dfComparison['JaccardScore'].round(3)
# dfComparison.to_csv(os.path.join(OUTDIR, inputFileName + '.csv'), sep="\t", index = False)
#
# jScore_allModel = jaccard(lAllGenes_original, lAllGenes_gpruler)
# print('jScore_allModel\t', jScore_allModel, '\n')


## Euclidean distance on the output vectors can quantify the extent of similarity in the Boolean logic combination

# print(distance.hamming([1, 0,0,0], [0, 1,1,0]))
# # print(distance.hamming([1, 0, 0], [1, 1, 0]))
# print('dist Norm\t', distance.hamming([1, 0,0,0], [0, 1,1,0]), '\n')
# print('dist2\t', paired_distances(np.array([[1], [0],[0],[0]]), np.array([[0], [1],[1],[0]])), '\n')
#

# sys.exit()

lDistances = []

lHamming = []

hammingNot0 = 0
for row in dfComparison.itertuples():

    originalRule = row.rule_original
    gprulerRule = row.rule_GPRuler

    # print('originalRule\t', originalRule)
    # print('gprulerRule\t', gprulerRule)

    if (pd.isna(originalRule) == True or originalRule == '') and (pd.isna(gprulerRule) == True or gprulerRule == ''):
        distHamming = 0
        lHamming.append(distHamming)
    else:
        if row.JaccardScore == 1:
            # print('JaccardScore\t', row.JaccardScore)
            dfGenes1 = gL.extractRegexFromItem(originalRule, regexOrgSpecific)
            lGenes1 = gL.unique(list(dfGenes1[0]))
            lGenes1.sort()
            # print('lGenes1\t', lGenes1)

            dfGenes2 = gL.extractRegexFromItem(gprulerRule, regexOrgSpecific)
            lGenes2 = gL.unique(list(dfGenes2[0]))
            lGenes2.sort()
            # print('lGenes2\t', lGenes2)

            maxGenesNumber = max(len(lGenes1), len(lGenes2))
            # print('maxGenesNumber\t', maxGenesNumber)
            if maxGenesNumber >= 20:
                print('DISTANCE FOR PIU DI 20 GENI')
                distHamming = np.nan
            else:
                mTruth, lGenes, n = generateTruthTable(lGenes1, lGenes2) # qua genero la truth table con tutte le possibili variabili risultato del'iunione delle due regole che voglio confrontare
                # print('mTruth\n', mTruth, '\n')
                # print('lGenes\t', lGenes, '\n')
                # print('len genes\t', len(lGenes), '\n')
                try:
                    aRes1 = evaluateBoolExpressions(mTruth, lGenes, originalRule)
                    aRes2 = evaluateBoolExpressions(mTruth, lGenes, gprulerRule)
                    # print('aRes1\n', aRes1, '\n')
                    # print('aRes2\n', aRes2, '\n')
                    mTruth1 = np.append(mTruth, aRes1, axis=1)
                    mTruth2 = np.append(mTruth, aRes2, axis=1)
                    # print('mTruth1\n', mTruth1, '\n')
                    # print('mTruth2\n', mTruth2, '\n')
                    # aDistance = paired_distances(mTruth1, mTruth2)
                    # print('aDistance\n', aDistance, '\n')
                    distHamming = distance.hamming([el for ll in aRes1.tolist() for el in ll], [el for ll in aRes2.tolist() for el in ll])
                    # print('HAmming\t', distHamming)
                    lHamming.append(distHamming)
                    if distHamming != 0:
                        print('hamming not 0\t', row.Rxn, '\t', distHamming)
                        hammingNot0 += 1
                except:
                    print('Errore truth matrix')
        else:
            distHamming = np.nan
    # print('distance\n', distance, '\n')
    lDistances.append(distHamming)


# fig = plt.figure(figsize=(13,10))
# ax = fig.add_subplot(111)
#
# font_dirs = ['C:\\Users\\Marzia\\Downloads\\Roboto']
# font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
#
# for font_file in font_files:
#     font_manager.fontManager.addfont(font_file)
#
# from matplotlib.ticker import PercentFormatter
#
#
# matplotlib.rcParams['font.sans-serif'] = "Roboto"
# matplotlib.rcParams['font.family'] = "Roboto"
#
# print('len hamming\t', len(lHamming))
#
# # counts, bins = np.histogram(lHamming)
# # print('counts\n', counts)
# # print('bins\n', bins)
# # freq = counts/float(counts.sum())
# #
# # ax.hist(bins[:-1], freq, weights=counts,
# #         alpha=0.5, density=True, color = 'green')
#
# ax.hist(lHamming, alpha=0.5, density=False, color = 'green', weights=np.zeros_like(np.array(lHamming)) + 1. / np.array(lHamming).size)
#
# ax.set_xlabel('Hamming distance',fontsize= 28)
# ax.set_ylabel('Relative frequency',fontsize= 28)
# ax.set_title(modelName,fontsize= 35)
# plt.xticks(rotation=45, fontsize= 25)
# plt.yticks(fontsize= 25)
# plt.savefig(os.path.join(FIGUREDIR, modelName + '_hammingForJaccard1_2021910.pdf'))
# plt.close()

lJaccard_postCuration = []
for row in dfComparison.itertuples():
    if  row.Evaluation == 'Perfect match' or row.Evaluation == 'Corrected by GPRuler':
        jScore = 1.0
    else:
        jScore = row.JaccardScore
    lJaccard_postCuration.append(jScore)

dfComparison['JaccardScore_postCuration'] = lJaccard_postCuration


dfComparison = dfComparison[['Rxn', 'rule_original', 'rule_GPRuler', 'Evaluation', 'JaccardScore', 'JaccardScore_postCuration']]
dfComparison['HammingDistance'] = lDistances
dfComparison.to_csv(os.path.join(OUTDIR, inputFileName + '.csv'), sep="\t", index = False)

###########
lJaccard = []
lHam = []
for row in dfComparison.itertuples():
    # if row.Evaluation == 'Corrected by GPRuler':
    #     lJaccard.append(1.0)
    if row.Evaluation == 'Not automatically reconstructed by GPRuler' or row.Evaluation == 'Corrected by GPRuler':
        lJaccard.append(row.JaccardScore)
        if row.JaccardScore == 1:
            lHam.append(row.HammingDistance)


print('\nHAMMING\n', lHam, '\n')

## jaccard
fig = plt.figure(figsize=(13,10))
ax = fig.add_subplot(111)
font_dirs = ['C:\\Users\\Marzia\\Downloads\\Roboto']
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
for font_file in font_files:
    font_manager.fontManager.addfont(font_file)

matplotlib.rcParams['font.sans-serif'] = "Roboto"
matplotlib.rcParams['font.family'] = "Roboto"
ax.hist(lJaccard,
        alpha=0.5,
        density=False, color = 'steelblue', weights=np.zeros_like(np.array(lJaccard)) + 1. / np.array(lJaccard).size)

ax.set_xlabel('Jaccard index',fontsize= 32)
ax.set_ylabel('Relative frequency',fontsize= 32)
ax.set_title(modelName,fontsize= 37)
plt.xticks(rotation=45, fontsize= 29)
plt.yticks(fontsize= 29)
plt.ylim(0, 1)
plt.savefig(os.path.join(FIGUREDIR, modelName + '_jaccardAllMismatched_2021911.pdf'))
plt.close()

## Hamming
fig = plt.figure(figsize=(13,10))
ax = fig.add_subplot(111)
font_dirs = ['C:\\Users\\Marzia\\Downloads\\Roboto']
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
for font_file in font_files:
    font_manager.fontManager.addfont(font_file)

from matplotlib.ticker import PercentFormatter
matplotlib.rcParams['font.sans-serif'] = "Roboto"
matplotlib.rcParams['font.family'] = "Roboto"
print('len hamming\t', len(lHamming))
ax.hist(lHam, alpha=0.5, density=False, color = 'green', weights=np.zeros_like(np.array(lHam)) + 1. / np.array(lHam).size)

ax.set_xlabel('Hamming distance',fontsize= 32)
ax.set_ylabel('Relative frequency',fontsize= 32)
ax.set_title(modelName,fontsize= 37)
plt.xticks(rotation=45, fontsize= 29)
plt.yticks(fontsize= 29)
plt.ylim(0, 1)
plt.xlim(0, 1)
plt.savefig(os.path.join(FIGUREDIR, modelName + '_hammingMismatchesWJaccard1_2021911.pdf'))
plt.close()
