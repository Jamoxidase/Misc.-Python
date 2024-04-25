#!/usr/bin/env python3

########################################################################
# BME205 Asgn5: problem 21| Compute the Probability of a Hidden Path
# Author: James Larbalestier
# Group: Tyler Gaw
########################################################################

import sys
import numpy as np

class HMMDecoder:
    def __init__(self, alphabet, states, transitionMatrix, emissionMatrix):
        self.alphabet = alphabet
        self.states = states
        self.transitionMatrix = np.array(transitionMatrix, dtype=np.longdouble)
        self.emissionMatrix = np.array(emissionMatrix, dtype=np.longdouble)

    def viterbiAlgorithm(self, x):
        n = len(x)
        m = len(self.states)

        # Initialize the Viterbi matrix
        viterbiMatrix = np.zeros((n, m), dtype=np.longdouble)
        backtrackMatrix = np.empty((n, m), dtype='object')

        # Initialize the first column of the Viterbi matrix
        for i in range(m):
            viterbiMatrix[0][i] = self.emissionMatrix[i][self.alphabet.index(x[0])]

        # Fill in the rest of the Viterbi matrix
        for j in range(1, n):
            for i in range(m):
                maxProb = np.longdouble(0.0)
                maxState = ''
                for k in range(m):
                    prob = viterbiMatrix[j - 1][k] * self.transitionMatrix[k][i] * \
                           self.emissionMatrix[i][self.alphabet.index(x[j])]
                    if prob > maxProb:
                        maxProb = prob
                        maxState = self.states[k]
                viterbiMatrix[j][i] = maxProb
                backtrackMatrix[j][i] = maxState

        # Find the last state with the maximum probability
        lastStateIndex = np.argmax(viterbiMatrix[-1])
        path = [self.states[lastStateIndex]]

        # Backtrack to reconstruct the optimal path
        for j in range(n - 1, 0, -1):
            lastStateIndex = self.states.index(backtrackMatrix[j][lastStateIndex])
            path.insert(0, self.states[lastStateIndex])

        return ''.join(path)


def main(inFile = None): 
    with open(inFile) as fh:
        emissionPath = fh.readline().rstrip()
        junk = fh.readline().rstrip()
        emissionOptions = fh.readline().rstrip().split()
        junk = fh.readline().rstrip()
        stateOptions = fh.readline().rstrip().split()
        junk = fh.readline().rstrip()
        junk = fh.readline().rstrip()
        stateMatrix = []
        for i in range(len(stateOptions)):
            line = fh.readline().rstrip()
            vals = [float(value) for value in line.split('\t')[1:]]
            stateMatrix.append(vals)
        junk = fh.readline().rstrip()
        junk = fh.readline().rstrip()
        emissionMatrix = []
        for i in range(len(stateOptions)):
            line = fh.readline().rstrip()
            vals = [float(value) for value in line.split('\t')[1:]]
            emissionMatrix.append(vals)

    decoder = HMMDecoder(emissionOptions, stateOptions, stateMatrix, emissionMatrix)
    result = decoder.viterbiAlgorithm(emissionPath)
    print(result)
       
if __name__ == "__main__":
    main(inFile = 'rosalind_ba10c.txt')

