#!/usr/bin/env python 3.12.0
##############################################################################################################################################
# BME205 Asgn4: problem 15 Rosalind -> DAG highest weighted path
# Author: James Larbalestier
# Group: Tyler Gaw
###############################################################################################################################################
'''
Longest Path in a DAG Problem

Find a longest path between two nodes in an edge-weighted DAG.

Given: An integer representing the source node of a graph, followed by an integer representing the sink node of the graph, followed
by an edge-weighted graph. The graph is represented by a modified adjacency list in which the notation '0->1:7' indicates that an
edge connects node 0 to node 1 with weight 7.

Return: The length of a longest path in the graph, followed by a longest path. (If multiple longest paths exist, you may return any one.)
'''

import argparse
from collections import defaultdict

class LongestPathInDag:
    '''
    Find the longest path in a weighted directed acyclic graph (DAG) from a specified source to sink.

    Args:
        graph (dict): The graph represented as a dictionary, where keys are nodes and values are lists of
                      tuples (neighbor, weight) representing outgoing edges.
        start (int): The starting node.
        end (int): The ending node.

    Returns:
        tuple: A tuple containing two elements. The first element is the length of the longest path. The second
               element is a list representing the longest path.

    Algo overview:
    1. Elute the topological order of the graph.
    2. Initialize scores with negative infinity for all nodes except the start node.
    3. Backtrack to find the path while updating scores.
    4. Reconstruct the longest path from the end node to the start node.

    This implementation preserves the graph structure and allows for backtracking without destroying the graph.
    '''

    def __init__(self, graph, start, end):
        self.graph = graph
        self.start = start
        self.end = end

    def topologicalSort(self):
        # Find in-degrees for each node.
        inDegree = defaultdict(int)
        queue = []
        for node in self.graph:
            for neighbor, _ in self.graph[node]:
                inDegree[neighbor] += 1

        # Initialize the queue with nodes that have no incoming edges
        for node in self.graph:
            if inDegree[node] == 0:
                queue.append(node)

        # Topological sorting
        topOrder = []
        while queue:
            currentNode = queue.pop(0)
            topOrder.append(currentNode)

            for neighbor, _ in self.graph[currentNode]:
                inDegree[neighbor] -= 1
                if inDegree[neighbor] == 0:
                    queue.append(neighbor)

        return topOrder

    def argMax(self, dictionary):
        # Find the key with the maximum value in a given dictionary. Not currently used.
        return max(dictionary, key=dictionary.get)

    def findLongestPath(self):
        # Outsource topological sorting
        topOrder = self.topologicalSort()

        # Initialize scores with negative infinity, except for the start node
        scores = {node: float('-inf') for node in self.graph}
        scores[self.start] = 0

        # Backtrack to find the path
        predecessors = {node: None for node in self.graph}
        for node in topOrder:
            if node == self.end:
                break
            for neighbor, weight in self.graph.get(node, []):
                if scores[node] + weight > scores[neighbor]:
                    scores[neighbor] = scores[node] + weight
                    predecessors[neighbor] = node

        # Reconstruct the longest path
        currentNode = self.end
        path = []
        while currentNode is not None:
            path.append(currentNode)
            currentNode = predecessors[currentNode]

        path.reverse()
        return scores[self.end], path


class ParseGraph:
    '''
    Parses and stores graph data from a file -> provides access to the graph's structure: source, sink, and edges.
    '''
    def __init__(self, start, end, edges):
        '''
        Initializes a ParseGraph instance with start, end, and edges.

        Parameters:
        - start (int): The start node of the graph.
        - end (int): The end node of the graph.
        - edges (defaultdict): A defaultdict representing the graph's edges.

        Returns:
        - ParseGraph: An instance of the ParseGraph class.
        '''
        self.start = start
        self.end = end
        self.edges = edges

    @classmethod
    def fromFile(cls, fileName):
        '''
        Creates a ParseGraph instance by parsing graph data from a file.

        Parameters:
            - cls: The class.
            - fileName (str): The path to the file containing graph data.

        Returns:
            - ParseGraph: An instance of the ParseGraph class.
        '''
        start, end, edges = cls.readGraphFromFile(fileName)
        return cls(start, end, edges)

    @staticmethod
    def parseEdge(line):
        '''
        Parses a line representing an edge in the graph.

        Parameters:
            - line (str): A string representing an edge in the graph.

        Returns:
            - tuple: A tuple containing source node, target node, and weight.
        '''
        parts = line.replace('->',":").split(":")
        source = int(parts[0])
        target = int(parts[1])
        weight = int(parts[2])
        return source, target, weight

    @staticmethod
    def readGraphFromFile(filePath):
        '''
        Reads graph data from a file and returns start, end, and edges.

        Parameters:
            - filePath (str): The path to the file containing graph data.

        Returns:
            - tuple: A tuple containing start node, end node, and edges.
        '''
        with open(filePath, 'r') as file:
            start = int(file.readline().strip())
            end = int(file.readline().strip())
            edges = defaultdict(list)
            for line in file:
                source, target, weight = ParseGraph.parseEdge(line.strip())
                edges[source].append((target, weight))
        return start, end, edges


def main(filePath: str) -> None:
    # Access input file, provides accessible data structure for following process.
    parsedGraph = ParseGraph.fromFile(filePath)
    start = parsedGraph.start
    end = parsedGraph.end
    edges = parsedGraph.edges

    # Instantiate LongestPathInDag
    result = LongestPathInDag(edges, start, end).findLongestPath()
    length = result[0]
    pathList = list(result[1])

    print(length)
    print('->'.join(str(dat) for dat in pathList))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Longest Path in DAG') #note that fPath should be passed without <
    parser.add_argument('filePath', type=str, help='Input file path')
    args = parser.parse_args()
    main(args.filePath)