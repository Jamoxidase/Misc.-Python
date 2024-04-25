import sys
from collections import Counter
import random
import math
from functools import reduce

'''
fix entropy
'''
class FastAreader:
    ''' Read fasta file from specified fname or STDIN if not given'''
    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''

        self.fname = fname
        self.fileH = None

    def doOpen(self):
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):

        header = ''
        sequence = ''

        with self.doOpen() as self.fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = self.fileH.readline()
            while not line.startswith('>'):
                line = self.fileH.readline()
            header = line[1:].rstrip()

            for line in self.fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence


class Whisper:
    '''
    identify consensus promotor motif allowing for mismatches
    '''
    def __init__(self, k, i, p, g):
        self.k = k
        self.i = i
        self.p = p
        self.g = g
        self.q = Counter('')

        self.nullDict = {}
        self.sequences= []
        self.netLength = 0
        self.workingMotifs = [] #must be reset to handle iterations
        self.countMatrix = []
        self.profileMatrix = []
        self.consensus = ''
        self.profileScore = 0

    def addSequence(self, sequence):
        sequence = sequence.upper()
        self.netLength += len(sequence)
        self.q += self.countBases(sequence) #gross nucleotide count across all added sequences
        self.sequences.append(sequence) #working list of sequences

    def countBases(self, sequence):
        return dict(Counter(sequence))
    
    def randomSampler(self, sequence):
        '''
        return random motif of len k from sequence and append it to working motifs
        '''
        seqLength = len(sequence)
        start = random.randint(0,seqLength - self.k)
        end = start + (self.k) 

        sample = sequence[start:end]
        self.workingMotifs.append(sample)

    def compute(self):
        '''
        organizational, signals that all sequences have been added
        '''
        self.nullDict = {key: value / self.netLength for key, value in  dict(self.q).items()}
        best = ('good luck', -10000000)

        for i in range(self.i):
            self.workingMotifs = []
            for sequence in self.sequences:
                self.randomSampler(sequence) #sample from each saquence added
            while 1 == 1:
                self.makeCountMatrix()
                self.makeProfileMatrix()
                self.makeConsensusMotif()
                self.profileScoring()
                compare = self.fittingMotifs()
                if self.workingMotifs == compare: #if the current working motifs are already the best fitting motifs given their profile, initiate new seed sampling
                    maybeBest = (self.consensus, self.profileScore, self.workingMotifs)
                    if best[1] < maybeBest[1]:
                        best = maybeBest
                    self.countMatrix = []
                    self.profileMatrix = []
                    self.consensus = ''
                    self.profileScore = 0
                    break
                else:                             #if the current working motifs are not the best fitting motifs in the sequence given their profile, replace with the new best fitting motifs
                    self.workingMotifs = compare  #a new profile will be made from these motifs, and the process repeates
                    self.countMatrix = []
                    self.profileMatrix = []
                    self.consensus = ''
                    self.profileScore = 0

        return best

    def fittingMotifs(self):
        '''
        given a profile, this method finds the most probable motifs in each sequence. When the current motifs are 
        already the best given their own profile, new seed sampling is initiated.
        '''
        sequenceMotifs = []
        for sequence in self.sequences:
            bestMotifProb = -1.0
            newMotif = ''
            for startIndex in range(len(sequence)-self.k):
                endIndex = startIndex + self.k #is this right? check
                motif = sequence[startIndex:endIndex]
                zipped_data = zip(motif, self.profileMatrix)
                motifProb = 1.0  # Initialize with 1.0 for multiplication
                # Multiply the values
                for base, dictionary in zipped_data:
                    motifProb *= dictionary[base]
                #motifProb = reduce(lambda x, y: x * y, (self.profileMatrix[base] for base, self.profileMatrix in zip(motif, self.profileMatrix))) error with indexing with base not index

                if motifProb > bestMotifProb: #identify the most probable motif in the sequence given profile
                    newMotif = motif
                    bestMotifProb = motifProb
            sequenceMotifs.append(newMotif) #temp list of these new motifs
        return sequenceMotifs

    def makeCountMatrix(self):
        '''
        call after adding all sequences
        '''
        for index in range(self.k): #all indicies in motifs
            countIndex = []
            
            for seq in self.workingMotifs: 
                #print(index)
                #print(seq)
                base = seq[index]
                countIndex.append(base)
            countDictPartial = dict(self.countBases(countIndex)) #####psu
            countDict = {nt: countDictPartial.get(nt, 0) for nt in 'ATCG'}
            '''for _ in range(self.p):
                for nt in "ACGT":
                    countDict[nt] += 1 ###TH'''
            self.countMatrix.append(countDict)

    def makeProfileMatrix(self):
        for composition in self.countMatrix:
            total = sum(composition.values())
            baseProbabilities = {key: (count + self.p) /( total + (4* self.p)) for key, count in composition.items()}  #####EDITING P
            self.profileMatrix.append(baseProbabilities)

    def makeConsensusMotif(self):
        for indexProbability in self.profileMatrix:
            consensusBase = max(indexProbability, key=indexProbability.get)
            self.consensus += consensusBase

    def profileScoring (self):
        '''
        entropy
        update
        '''
        #Relative Entropy Calculation
        
        for pos in range(self.k):

            for nuc in 'ATCG':
                #print('profile: ', self.profileMatrix[pos][nuc], 'times log2 ', self.profileMatrix[pos][nuc], ' divided by ', (self.nullDict[nuc]))
                try:
                    #print(self.profileMatrix)
                    self.profileScore += self.profileMatrix[pos][nuc] * math.log2((self.profileMatrix[pos][nuc])/(self.nullDict[nuc]))     


                except ValueError:
                    # Handle the math domain error by setting the profile score to 0
                    self.profileScore += 0
                    
class MotifHighlight:
    def __init__(self, motifs, seqData):
        self.motifs = motifs
        self.seqData = seqData

    def highlightMotifsInSequences(self):
        outputList = []
        for i in range(len(self.motifs)):
            header, sequence = self.seqData[i]
            motif = self.motifs[i]

            sequence = sequence.lower()
            motif = motif.lower()
            
            pos = sequence.find(motif)
            
            if pos != -1:
                # Highlight the matching region by making it uppercase
                seq = sequence[:pos] + '  -->' + sequence[pos:pos+len(motif)].upper() + '<--  ' + sequence[pos+len(motif):]
                outputList.append((header, seq))
        return outputList

class CommandLine() :
    """ Constructor implements an arg parser to interpret the command line argv string using argparse. """
    def __init__(self, inOpts=None) :

        import argparse
        self.parser = argparse.ArgumentParser(description = 'Missing Motif',
                                                epilog = '',
                                                add_help = True,
                                                prefix_chars = '-',
                                                usage = 'python3 randomizedMotifSearch.py [options]'
                                                )

        self.parser.add_argument('-g', '--gibbs', default=False, action = 'store_true', help='Use Gibbs sampling to find the optimal consensus motif.')
        self.parser.add_argument('-m', '--print', default=False, action = 'store_true', help='Print the specific motif and the name of the contributing sequence for each of the participating fasta records..')
        self.parser.add_argument('-k', '--length', type=int, default=13, action = 'store', help='Motif length.')
        self.parser.add_argument('-p', '--pseudocounts', type=int, default=1, action = 'store', help='Number of pseudocounts.')
        self.parser.add_argument('-i', '--iterations', type=int, default=100000, action = 'store', help='Number of iterations of descent to complete.')

        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

def main(args=None):
    if args == None:
        commandInput = CommandLine()
    else:
        commandInput = CommandLine(args)

    analyzer = Whisper(commandInput.args.length, commandInput.args.iterations, commandInput.args.pseudocounts, commandInput.args.gibbs)

    seqData = []
    for record in FastAreader().readFasta():
        seqData.append(record)
        seq = record[1]
        analyzer.addSequence(seq)

    #instance is now filled with sequences
    boo = analyzer.compute() #make matrix
    print('Consensus: ', boo[0], ' Profile Score: ', boo[1])
    print()
    if commandInput.args.print:
        highlight = MotifHighlight(boo[2], seqData)
        display = highlight.highlightMotifsInSequences()
        for head, seq in display:
            print(head)
            print(seq)
            print()

if __name__ == "__main__":
    main()

