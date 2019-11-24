#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Trevor Ridgley (tridgley)
# BME-205: DAGs - a generalization of the alignment problem, Problem 15
"""
Created on Mon Oct 21 20:55:38 2019

@author: tridgley

problem15.py - Find the Longest Path in a DAG

Find a longest path between two nodes in an edge-weighted DAG.

Given: An integer representing the source node of a graph, followed by an integer 
       representing the sink node of the graph, followed by an edge-weighted graph. 
       The graph is represented by a modified adjacency list in which the notation "0->1:7" 
       indicates that an edge connects node 0 to node 1 with weight 7.

Return: The length of a longest path in the graph, followed by a longest path. 
        (If multiple longest paths exist, you may return any one.)
"""
import sys

class DAG:
    """Create a directed acyclic graph object."""
    
    def __init__(self, start, end, adjList):
        """Initialize a DAG with start, end, and adjacency info
        
        Params:
            start = name of the src node in the adj list.
            end = name of the sync node in the adj list.
            adjList = a text encoded DAG with form: fromNode -> toNode
        """
        self.adjList = adjList
        self.start = start
        self.end = end
        self.DAG = dict()
            
    def makeGraph2(self):
        """Make a graph using node objects from the text-encoded adjList."""
        for record in self.adjList:
            items = record.split('->')
            fromAdj = items[0]
            connections = items[1].split(':')
            toAdj = connections[0]
            edgeWeight = int(connections[1])
            
            # Never connect start with incoming edges
            if toAdj not in self.DAG.keys():
                toNode = Node(toAdj)
                self.DAG[toAdj] = toNode
            if toAdj != self.start:
                self.DAG[toAdj].addData(fromAdj, edgeWeight)
                
            # Only connect start with its outgoing edges
            if fromAdj not in self.DAG.keys():
                fromNode = Node(fromAdj)
                self.DAG[fromAdj] = fromNode
            if fromAdj == self.start:
                self.DAG[fromAdj].addData(None, 0)
                self.DAG[fromAdj].total = 0
            if toAdj != self.start:
                self.DAG[fromAdj].addNext(self.DAG[toAdj])
                self.DAG[toAdj].addPrev(self.DAG[fromAdj])
            
    def findLongestPath4(self):
        """This version is based on the topological sort problem.
        
        Note that node.next and node.prev fields are node objects
        To get their node number, use node.name
        To get their incoming nodes and weights, use node.weights.keys()
        To get their cumulative score, use node.total
        """
        queue = self.noIncomingEdges()
        done = list()
        while len(queue) > 0:
            node = queue[0]
            done.append(node)
            del queue[0]
            for nextNode in self.DAG[node].next:
                ready = True
                for j in nextNode.weights.keys():
                    score = self.DAG[j].total + nextNode.weights[j]
                    if j not in done:
                        ready = False
                    elif score > nextNode.total:
                        nextNode.total = score
                if ready:
                    queue.append(nextNode.name)
                    
        
    def noIncomingEdges(self):
        """Return the nodes from our graph that have no incoming edges."""
        noIncoming = list()
        for node in self.DAG.keys():
            if not len(self.DAG[node].prev):
                noIncoming.append(node)
        return noIncoming
    
    def backTrack(self):
        """Return the longest path of nodes using a reverse walk."""
        node = self.DAG[self.end]
        pathList = [node.name]
        while node.name != self.start:
            incomingDict = node.weights
            for i in incomingDict.keys():
                inNode = self.DAG[i]
                subtotal = node.total - incomingDict[i]
                if subtotal >= 0 and subtotal == inNode.total:
                    pathList.insert(0, inNode.name)
                    node = self.DAG[inNode.name]
                    break
        return pathList
                
class Node():
    """Store k-1mer data and links to overlapping k-mers"""
    
    def __init__(self, name):
        """Initialize a node for use in a doubly linked list."""
        self.name = name
        self.weights = dict()
        self.total = 1-sys.maxsize
        self.next = list()
        self.prev = list()
        
    def addData(self, backNode, length):
        """Update the weights dict with another node that points to this one.
        
        Params:
            backNode = the name of an incoming node .
            length = the edge weight from incoming node.
        """
        self.weights[backNode] = length
        
    def addPrev(self, prevNode):
        """Hook up the previous (suffix) k-mer to the current node.
        
        Params:
            prevNode = an incoming node object to the current node.
        """
        self.prev.append(prevNode)
        
    def addNext(self, nextNode):
        """Hook up the next node (prefix) k-mer to the current node.
        
        Params:
            nextNode = an outgoing node object from the current node.
        """
        self.next.append(nextNode)

class FastAreader :
    """Class for parsing fasta sequence files."""
	
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        """Handle fasta data via file object or std input stream"""
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
 
    def readFasta (self):
        """Parse a fasta file
        
        Returns: Tuple with record header at index 0, sequence at index 1
        """
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
			
            header = ''
            sequence = ''
 
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
						
        yield header,sequence
        
    def readKmersWithHeader(self):
        """Open and parse a Rosalind file of kmers, except k length is at the top.
        
        Modified from readFasta() above (author David Bernick)
        """

        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            header = ''
            sequence = list()
 
            # skip to first fasta header
            line = fileH.readline()
            while set(line.strip()) <= set('ACGT') :
                line = fileH.readline()
            header = int(line.strip())
            for line in fileH:
                if not set(line.strip()) <= set('ACGT'):
                    yield header,sequence
                    header = int(line.strip())
                    sequence = list()
                else :
                    sequence.append(''.join(line.rstrip().split()).upper())
						
        yield header,sequence  
        
    def readString(self):
        """Open and parse a Rosalind string file with k in header and string thereafter.
        
        Modified from readFasta() above (author David Bernick)
        """
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            header = ''
            sequence = ''
 
            line = fileH.readline()
            while set(line.strip()) <= set('ACGT') :
                line = fileH.readline()
            header = int(line.strip())

            for line in fileH:
                if not set(line.strip()) <= set('ACGT'):
                    yield header,sequence
                    header = int(line.strip())
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
						
        yield header,sequence  
        
    def readKmers(self):
        """Open and parse a Rosalind file of kmers.
        
        Modified from readFasta() above (author David Bernick)
        """
        sequence = ''
        with self.doOpen() as fileH:
            for line in fileH:
                sequence = ''.join(line.strip().split()).upper()
                yield sequence 
                
    def readAdjList(self):
        """Open and parse a Rosalind file containing adjacency list.
        
        Modified from readFasta() above (author David Bernick)
        """
        with self.doOpen() as fileH:
            for line in fileH:
                item = line.strip().split()
                yield item
                
    def readDag(self):
        """Open and parse a Rosalind file with text-encoded weighted DAG with start and end node headers.
        
        Modified from readFasta() above (author David Bernick)
        """
        header1, header2, sequence = (str(), str(), list())
        
        with self.doOpen() as fileH:
            for line in fileH:
                if not header1:
                    header1 = line.strip()
                elif not header2:
                    header2 = line.strip()
                else:
                    sequence.append(line.strip())
        yield header1, header2, sequence
                
def main():
    """Parse a Rosalind file with text-encoded DAG and find its longest path."""
    myReader = FastAreader()
    myDAG = None
    for start,end,adjList in myReader.readDag():
        myDAG = DAG(start, end, adjList)
    myDAG.makeGraph2()
    for key in sorted(myDAG.DAG.keys()):
        prevData, nextData = (list(), list())
        for i in range(len(myDAG.DAG[key].prev)):
            prevData.append(myDAG.DAG[key].prev[i].name)
        for i in range(len(myDAG.DAG[key].next)):
            nextData.append(myDAG.DAG[key].next[i].name)
    myDAG.findLongestPath4()
    print(myDAG.DAG[myDAG.end].total)
    print('->'.join(myDAG.backTrack()))

main()

