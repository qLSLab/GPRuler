# -*- coding: utf-8 -*-
import re
import os
import sys
import pandas as pd
import numpy as np
import re

def setWorkingDirs(dataDir='rawData', outDir=None):
    """Set the working directories dataDir, outDir, mapsDir, reportDir.

    If arguments are omitted current working directory replace them.

    Keyword arguments:
     dataDir: The directory to look for the input files.
     outDir:  The directory to write the output files.

    Return:
    dataDir, outDir, reportDir, modelDir, logDir, figureDir, scriptsDir, mapDir
    """
    cwDir = os.getcwd()
    if dataDir is None:
        dataDir = cwDir + os.sep
    else:
        if not os.path.exists(dataDir):
            # print'the input directory', dataDir, 'does not exists. Fix it.'
            sys.exit()
        else:
            dataDir = os.path.abspath(dataDir)
            if dataDir[-1] != os.sep:
                dataDir += os.sep
    if outDir is None:
        outDir = cwDir
    else:
        outDir = os.path.abspath(outDir) + os.sep
    reportDir = pathJoinCheck('outputs', outDir)
    return dataDir, outDir, reportDir

def unique(a):
	""" return the list with duplicate elements removed """
	return list(set(a))

def difference(a, b):
    """ return the difference of two lists """
    return list(set(a) - set(b))


def extractRegexFromItem(item, regex):
    sItem = pd.Series(data=item)
    if sItem.empty == False:
        dfItem = sItem.str.extractall(regex)
    else:
        dfItem = []
    return dfItem


def writeLineByLineToFile(stream, dataToWrite, spaziatore):
    stream.write(spaziatore.join(str(s) for s in dataToWrite) + '\n')


def intersect(a, b):
	""" return the intersection of two lists """
	return list(set(a) & set(b))

def isCoeff(s):
    """Determine if a string splitted on the spaces the first element is the
    stoichiometric coefficient or not.
    Example: if string is "2 A" return True; if string is "A" it returns False"""
    answer = reModule.match('((\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)$', s)
    if answer is not None:
        # print('\t match ', answer.group(0)
        return True
    else:
        return False

def dizReaProd(s):
    termini = s.strip().split(' + ')
    # print('termini', termini
    diz = {}
    for termine in termini:
        coeffMetabo = termine.split()
        # print('coeffMetabo', coeffMetabo
        coeff = coeffMetabo[0]
        if isCoeff(coeff) is True:
            coeff = float(coeffMetabo[0])
            metabolita = ' '.join(coeffMetabo[1:])
        else:
            metabolita = ' '.join(coeffMetabo)
            coeff = 1.0
        # print('coeff', coeff
        # print('metabolita', metabolita
        diz[metabolita] = coeff
    return diz
