#!/usr/bin/env python3

########################################################################
# BME205 Asgn5: problem 25
# Author: James Larbalestier
# Group: Tyler Gaw
########################################################################

import numpy as np

def forwardBackwardAlgorithm(x, emissionStates, pathStates, transMatrix, emissMatrix):
    n = len(x)
    m = len(pathStates)

    # Initialize matrices
    forwardMatrix = np.zeros((n, m))
    backwardMatrix = np.zeros((n, m))
    resultMatrix = np.zeros((n, m))

    # Forward pass
    for i in range(m):
        forwardMatrix[0][i] = emissMatrix[i][emissionStates.index(x[0])] / m

    for j in range(1, n):
        for i in range(m):
            forwardMatrix[j][i] = sum(forwardMatrix[j - 1][k] * transMatrix[k][i] for k in range(m)) * emissMatrix[i][emissionStates.index(x[j])]
    
    # Backward pass
    for i in range(m):
        backwardMatrix[n - 1][i] = 1 / m

    for j in range(n - 2, -1, -1):
        for i in range(m):
            backwardMatrix[j][i] = sum(transMatrix[i][k] * emissMatrix[k][emissionStates.index( x[j + 1])] * backwardMatrix[j + 1][k] for k in range(m))
    
    # Smooth Over
    for i in range(n):
        for j in range(m):
            resultMatrix[i][j] = forwardMatrix[i][j] * backwardMatrix[i][j]
    
    # Normalize the resultMatrix
    rowSums = resultMatrix.sum(axis=1)
    normalizedResult = resultMatrix / rowSums[:, np.newaxis]

    return normalizedResult


def main(inFile = None): 
    with open(inFile) as fh:
        observedEmissions = fh.readline().rstrip()
        junk = fh.readline().rstrip()
        emissionStates = fh.readline().rstrip().split()
        junk = fh.readline().rstrip()
        pathStates = fh.readline().rstrip().split()
        junk = fh.readline().rstrip()
        junk = fh.readline().rstrip()
        transMatrix = []
        for i in range(len(pathStates)):
            line = fh.readline().rstrip()
            vals = [float(value) for value in line.split('\t')[1:]]
            transMatrix.append(vals)
        junk = fh.readline().rstrip()
        junk = fh.readline().rstrip()
        emissMatrix = []
        for i in range(len(pathStates)):
            line = fh.readline().rstrip()
            vals = [float(value) for value in line.split('\t')[1:]]
            emissMatrix.append(vals)

        transMatrix = np.array(transMatrix)
        emissMatrix = np.array(emissMatrix)


    result = forwardBackwardAlgorithm(observedEmissions, emissionStates, pathStates, transMatrix, emissMatrix)

    # Print the result
    print("\t".join(pathStates))
    for j in range(len(observedEmissions)):
        print("\t".join([f"{result[j][i]:.4f}" for i in range(len(pathStates))]))

if __name__ == "__main__":
    main(inFile = 'rosalind_ba10j.txt')
