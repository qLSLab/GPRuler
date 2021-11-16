# -*- coding: utf-8 -*-
import re
import os
import sys
import time
import pandas as pd
import numpy as np

dTestModelsParams = {
    'recon': {
        'filename': 'Recon3D_301_20200923.xml',
        'basenameStr': 'Recon3D_301_20200923',
        'prefix': 'recon3D',
        'incComp': False
    },
    'y7': {
        'filename': 'yeast_7.6_cobra.xml',
        'basenameStr': 'yeast_7.6_cobra',
        'prefix': 'y7',
        'incComp': True
    },
    'y8': {
        'filename': 'yeast8.xml',
        'basenameStr': 'yeast8',
        'prefix': 'y8',
        'incComp': True
    },
    'hmrcore': {
        'filename': 'HMRcore_20200328_wReconName.xml',
        'basenameStr': 'HMRcore_20200328_wReconName',
        'prefix': 'hmrCore',
        'incComp': False
    },
    'toy': {
        'filename': 'toymodel.xml',
        'basenameStr': 'toymodel',
        'prefix': 'toy',
        'incComp': False
    },
}

dChebiFiles = {
    'names': 'chebi_names_20201216.tsv.bz2',
    'uniprot': 'chebi_uniprot_20201216.tsv.bz2',
    'compounds': 'chebi_compounds_20201216.tsv.bz2',
    'accession': 'chebi_database_accession_20201216.tsv.bz2',
    'relation': 'chebi_relation_20201216.tsv.bz2',
    'inchi': 'chebiId_inchi_20201216.tsv.bz2'
}

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


def pathJoinOrExit(rootPath='.', filename=None):
    """Check the existence of a path and build it if necessary."""
    path = os.path.join(rootPath, filename)
    beingOrExit(path)
    return path


def setWorkingDirs(wrkDir=None, dizDirNames=None):
    """Set the working directories as recived from input dictionary.

    If arguments are omitted, current working directory is set as
    workingDirectory and no directory is nested.
    If wrkDir does not exists, the script create it.
    wrkDir/dir1
          /dir2
          /dir3

    Keyword arguments:
     wrkDir: The working directory path. if None it is set to cwd
     dizDirNames: The dictionary containing the names of the nested
                  directories, in the form: {'dirName', 'dirPath'}

    Return:
    (wrkDirPath, dizDirPaths) where:
    wrkDirPath: The working dir path
    dizDirPaths: The dictionary containing the paths of the nested directories,
     in the form: {'dirName', 'dirPath'}
    """
    cwDir = os.getcwd()
    dizDirPaths = {}
    if wrkDir is None:
        wrkDirPath = cwDir + os.sep
    else:
        if not os.path.exists(wrkDir):
            print('the working directory', wrkDir, 'does not exist. Fix it.')
            sys.exit()
        else:
            wrkDirPath = os.path.abspath(wrkDir)
            if wrkDirPath[-1] != os.sep:
                wrkDirPath += os.sep
    if dizDirNames is not None:
        for el in dizDirNames:
            dizDirPaths[el] = pathJoinCheck(dizDirNames[el], wrkDirPath)
    return (wrkDirPath, dizDirPaths)


def getTimeStamp():
    return time.strftime('%Y%m%d%H%M%S', time.localtime())


def logFileOpen(logDIR=None, timeStamp=None, aim=None, basename=None):
    if logDIR is None:
        logDIR = os.getcwd() + os.sep
    if timeStamp is None:
        timeStamp = getTimeStamp()
    if aim is None:
        aimStr = ''
    else:
        aimStr = '_' + aim
    if basename is None:
        basenameStr = ''
    else:
        basenameStr = basename + '_'
    logFileName = os.path.join(logDIR,
                               basenameStr + timeStamp + aimStr + '.log')
    logStream = open(logFileName, mode='w')
    return logStream


def toLog(logStream, string):
    logStream.write(string + '\n')


def loadModelParams(sInputLine=['', 'toy'], timeStamp=None, dDirs=None):
    print(sInputLine, timeStamp)
    wrongParams = False
    modelName = sInputLine[1]
    dModelParams = {}
    if modelName in dTestModelsParams.keys():
        dModelParams['filename'] = dTestModelsParams[modelName]['filename']
        dModelParams['basenameStr'] = dTestModelsParams[modelName][
            'basenameStr']
        dModelParams['prefix'] = dTestModelsParams[modelName]['prefix']
        dModelParams['incComp'] = dTestModelsParams[modelName]['incComp']
    else:
        modelXmlFile = modelName  # add checks for existence, xml extension ...
        nParams = len(sInputLine)
        print(nParams)
        dModelParams['filename'] = modelName
        dModelParams['basenameStr'] = getBaseName(modelName)
        if nParams == 2:
            dModelParams[
                'prefix'] = dModelParams['basenameStr'] + '_' + timeStamp
            dModelParams['incComp'] = True
        elif nParams == 3:
            try:
                isinstance(eval(sInputLine[2]), bool)
            except NameError:
                dModelParams['prefix'] = sInputLine[2]
                dModelParams['incComp'] = True
            else:
                dModelParams[
                    'prefix'] = dModelParams['basenameStr'] + '_' + timeStamp
                dModelParams['incComp'] = eval(sInputLine[2])
        elif nParams == 4:
            dModelParams['prefix'] = sInputLine[2]
            dModelParams['incComp'] = sInputLine[3]
        else:
            print('Please check the compliance of your input parameters')
            wrongParams = True
    for param in dModelParams:
        print(param, dModelParams[param])
    if wrongParams:
        sys.exit()
    dModelParams['path'] = pathJoinOrExit(dDirs['raw'], dModelParams['filename'])
    return dModelParams


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
