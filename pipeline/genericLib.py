# -*- coding: utf-8 -*-
import re
import os
import sys
import pandas as pd
import numpy as np


def getBaseName(path):
    """The program filename without extension"""
    return os.path.basename(path).split('.')[0]


def beingOrExit(prgPath):
    if not os.path.exists(prgPath):
        print('The file ', prgPath, 'does not exist, check the path')
        sys.exit()


def pathFilename(path, filename):
    """Join path and filename"""
    return os.path.join(path, filename)


def setPath(path):
    """Check the existence of a path and build it if necessary."""
    if not os.path.exists(path):
        os.makedirs(path)


def pathJoinCheck(dir2add, rootPath='.'):
    """Check the existence of a path and build it if necessary."""
    path = os.path.join(rootPath, dir2add)
    if not os.path.exists(path):
        os.makedirs(path)
    return path


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
            print('the input directory', dataDir, 'does not exists. Fix it.')
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


def writeLineByLineToFile(stream, dataToWrite, spacer):
    stream.write(spacer.join(str(s) for s in dataToWrite) + '\n')


def intersect(a, b):
    """ return the intersection of two lists """
    return list(set(a) & set(b))


def isCoeff(s):
    """Determine if a string splitted on the spaces the first element is the
    stoichiometric coefficient or not.
    Example: if string is "2 A" return True; if string is "A" it returns False"""
    answer = re.match('((\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)$', s)
    if answer is not None:
        # print('\t match ', answer.group(0)
        return True
    else:
        return False


def dizReaProd(s):
    termini = s.strip().split(' + ')
    diz = {}
    for termine in termini:
        coeffMetabo = termine.split()
        coeff = coeffMetabo[0]
        if isCoeff(coeff) is True:
            coeff = float(coeffMetabo[0])
            metabolita = ' '.join(coeffMetabo[1:])
        else:
            metabolita = ' '.join(coeffMetabo)
            coeff = 1.0
        diz[metabolita] = coeff
    return diz
