import sys
import genericLib as gL
import os
import pandas as pd
import time
import re

# setting working dirs
workingDirs = gL.setWorkingDirs()
RAWDIR = workingDirs[0]
OUTDIR = workingDirs[2]
MODELDIR = workingDirs[3]
FIGUREDIR = workingDirs[4]
DFDIR = workingDirs[5]

# Get Taxa Ids da KEGG
fkeggTaxaIds = open(os.path.join(RAWDIR, 'br08610.keg'), "r") # file get from https://www.biostars.org/p/433244/
match = re.findall(r'\[(.*?)\]', fkeggTaxaIds.read())
lkeggTaxaIds = [m.split('TAX:')[1] for m in match if m.startswith('TAX:')]
lkeggTaxaIds = gL.unique(lkeggTaxaIds)
print('ALL KEGG:\t', len(lkeggTaxaIds))

# Get Taxa Ids da NCBI
dfncbiTaxaIds = pd.read_csv(os.path.join(RAWDIR, 'taxdmp\\names.dmp'), sep = '\t', dtype=str, names = ['tax_id', 'sep1', 'name_txt', 'sep2', 'unique name', 'sep3', 'name class', 'sep4'])
lncbiTaxaIds = gL.unique(dfncbiTaxaIds.tax_id.tolist())
print('ALL NCBI:\t', len(lncbiTaxaIds))

## common and differences
commonIds = gL.intersect(lkeggTaxaIds, lncbiTaxaIds)
print('commonIds\t', len(commonIds))

onlyInKegg = gL.difference(lkeggTaxaIds, lncbiTaxaIds)
print('onlyInKegg\t', len(onlyInKegg), '\t', onlyInKegg)

onlyInNcbi = gL.difference(lncbiTaxaIds, lkeggTaxaIds)
print('onlyInNcbi\t', len(onlyInNcbi))

## coverage
coverage = ((len(commonIds) + len(onlyInKegg))/ len(lncbiTaxaIds))*100
print('coverage\t', coverage)
