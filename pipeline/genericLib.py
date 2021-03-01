import re
import os
import sys
import time
import pandas as pd
import numpy as np
import re


def getTimeStamp():
    return time.strftime('%Y%m%d%H%M%S', time.localtime())

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
