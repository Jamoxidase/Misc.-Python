'''
python 3.12.0 64bit

Problem 14
Asgn 3 BME 205
James Larbalestier
partners: Tyler Gaw , Fabrice Kurmann, Michael Valdovinos
'''

import sys

#problem 14 
class DeBruijnGraph:
    '''
    a form of an adjacency object that offers comprehensive node data structures.
    '''
    def __init__(self):
        self.Edges = {}  # Edges obj. contains dictionary of the associated nodes for each edge.

    #Each edge is broken into an inNode(prefix) and an outNode(suffix)
    def addEdge(self, kmer):
        inNode = kmer[:-1]
        outNode = kmer[1:]
        if inNode not in self.Edges:
            self.Edges[inNode] = []
        self.Edges[inNode].append(outNode)

    #For each input pattern, we add it as an edge. 
    def constructDeBruijn(self, k, patterns):
        for pattern in patterns:
            for i in range(len(pattern) - k + 1):
                kmer = pattern[i:i + k]
                self.addEdge(kmer)
        return self.Edges


class StringReconstruction:
    '''
    organizes the formation of a pattern polymer
    '''
    def __init__(self, k, patterns):
        self.k = k
        self.patterns = patterns
        self.graph = DeBruijnGraph().constructDeBruijn(k, patterns)
        self.reconstructedString = self.reconstructString()

    def reconstructString(self):
        '''
        Reconstruct the original DNA string from the De Bruijn graph.

        This method reconstructs a DNA string by finding the Eulerian path (via outsourcing to next class) in the De Bruijn graph.
        It starts by finding the starting node (the one with more out-degrees than in-degrees), then follows the Eulerian path.
        The reconstructed string is built by appending to the starting node, the last character of each kmer in the Eulerian path.

        '''
        #identifies the start node
        def findStartNode():  
            inDegrees = {}
            outDegrees = {}
            for inNode, outNode in self.graph.items(): 
                outDegrees[inNode] = len(outNode) #outNode has a list attribute, anyways why
                for edge in outNode: 
                    inDegrees[edge] = inDegrees.get(edge, 0) + 1
            for node in set(outDegrees.keys()) | set(inDegrees.keys()):
                if inDegrees.get(node, 0) < outDegrees.get(node, 0):
                    return node

        startNode = findStartNode()
        eularianPath = self.findEulerianPath(self.graph, startNode)
        reconstructedString = eularianPath[0]
        for kmer in eularianPath[1:]: #add the 
            reconstructedString += kmer[-1]
        return reconstructedString


    def findEulerianPath(self, graph, startNode): #adaped eularianPath method from previous problem to work with string reconstruction
        '''
        Find an Eulerian path in the given graph starting from the specified start node.
        This function performs a depth-first search (DFS) traversal to find the Eulerian path in the graph.
        It starts from the given start node and explores all connected nodes following the Eulerian path algorithm.
        '''
        def dfs(node):
            stack = [node]
            while stack: 
                currentNode = stack[-1]
                if currentNode in graphCopy and graphCopy[currentNode]:
                    neighbor = graphCopy[currentNode].pop(0)
                    stack.append(neighbor)
                else:
                    path.append(stack.pop())

        graphCopy = {inNode: outNode[:] for inNode, outNode in graph.items()}
        path = []
        dfs(startNode)
        eulerianPath = path[::-1]
        return eulerianPath

def main():
    input = [line.strip() for line in sys.stdin] #Read from standard input
    k = int(input[0])
    patterns = input[1:]
    result = StringReconstruction(k, patterns).reconstructedString
    print(result)

if __name__ == "__main__":
    main()
