import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import time


def createInputFile(fileName, outputDir, inputStrings):
    f = open(fileName,"w")
    for el in inputStrings:
        f.write(el + "\n")
        
    f.close()
    
def plotDOS(fileName):
    data = np.genfromtxt(fileName)
    plt.plot(data[:,0],data[:,1])
    plt.xlabel("E (eV)")
    plt.ylabel("DOS (arb.u.)")
    
def create_tBG_commensurate_cell(integerArray):
    n1 = integerArray[0]
    n2 = integerArray[1]
    n3 = integerArray[2]
    n4 = integerArray[3]
    fileNameIn = "GenMoire.f90_tBG_template.f90"
    fileNameOut = "GenMoire.f90"
    replaceInFile(fileNameIn,fileNameOut,"AAA", n1)
    replaceInFile(fileNameIn,fileNameOut,"BBB", n2)
    replaceInFile(fileNameIn,fileNameOut,"CCC", n3)
    replaceInFile(fileNameIn,fileNameOut,"DDD", n4)
    subprocess.call('gfortran ' + fileNameIn, cwd="./", shell=True)
    time.sleep(3)
    subprocess.call('a.out', cwd="./", shell=True)
    time.sleep(3)
    subprocess.call('python fracToCart.py BLBL.frac', cwd="./", shell=True)
    time.sleep(3)
    subprocess.call('python replaceFrac.py', cwd="./", shell=True)
    time.sleep(3)
    subprocess.call('python createXYZ.py', cwd="./", shell=True)
    
def replaceInFile(fileNameIn, fileNameOut, string, value):
    with open(fileNameIn, "rt") as fin:
        with open(fileNameOut, "wt") as fout:
            for line in fin:
                fout.write(line.replace(string, str(value)))
    os.rename(fileNameOut, fileNameIn)
    
def plotStructure(fileName):
    xVecBottom = []
    yVecBottom = []
    xVecTop = []
    yVecTop = []
    f = open(fileName, "r")
    for i, line in enumerate(f):
        a = line.split()
        if (i>3 and len(a)!=0):
            if (float(a[3])>18):
                xVecTop.append(float(a[1]))
                yVecTop.append(float(a[2]))
            elif (float(a[3])<18):
                xVecBottom.append(float(a[1]))
                yVecBottom.append(float(a[2]))
    plt.scatter(xVecBottom, yVecBottom, c="C0",s=2)
    plt.scatter(xVecTop, yVecTop, c="C1",s=2)
    plt.title("commensurate cell from generate.xyz file")
    
    
    

    
    
    
