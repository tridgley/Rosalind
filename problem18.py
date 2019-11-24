#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Trevor Ridgley
# BME-205: HMM Problems 16-19
"""
Created on Tue Nov 12 20:31:18 2019

@author: tridgley
Usage: python3 problem18.py < rosalind_ba10c.txt > result18.txt

Decoding Problem

Given: A string x, followed by the alphabet Σ from which x was constructed, 
       followed by the states States, transition matrix Transition, and emission 
       matrix Emission of an HMM (Σ, States, Transition, Emission).

Return: A path that maximizes the (unconditional) probability Pr(x, π) over all possible paths π.
"""
import sys, numpy as np

class Viterbi:
    """Viterbi representation of an HMM that finds the most probable emission states."""
    
    def __init__(self, string, alpha, states, transition, emission):
        """Create a viterbi object and initialize its dependencies.
                
        Params:
            string = the emitted data
            alpha = the alphabet that the emitted string is composed of
            states = the possible viterbi states that produced each emission
            transition = the matrix of probabilities between states
            emission = the matrix of probabilities for each emission
        """
        self.transMatrix = dict()
        self.emitMatrix = dict()
        self.data = string
        self.alpha = alpha
        self.states = states
        self.probs = [[np.float64() for j in range(len(self.states))] for i in range(len(self.data))]
        self.maxStates = [[int() for j in range(len(self.states))] for i in range(len(self.data))]
        self.pi = list()
        self.pi2 = [[str() for j in range(len(self.states))] for i in range(len(self.data))]
        for i in range(len(states)):
            convertTrans = list()
            convertEmit = list()
            for j in range(len(alpha)):
                convertEmit.append(float(emission[i][j]))
            for j in range(len(states)):
                convertTrans.append(float(transition[i][j]))
            self.transMatrix[states[i]] = convertTrans
            self.emitMatrix[states[i]] = convertEmit
    
    def calcPathProbs(self):
        """Score all possible paths using Viterbi algorithm."""
        start = self.alpha.index(self.data[0])
        self.probs[0] = [self.emitMatrix[i][start]/len(self.states) for i in self.states]
        for i in range(1, len(self.data)):
            # First iterate every character in the observed sequence
            currentData = self.alpha.index(self.data[i])
            for j in range(len(self.states)):
                # Then iterate every state j that possibly emitted the sequence character i
                currentEmit = self.emitMatrix[self.states[j]][currentData]
                candidateProbList = list()
                candidateProb = float('-inf')
                for k in range(len(self.states)):
                    # Finally, iterate every possible previous state to calc its transition
                    currentTrans = self.transMatrix[self.states[k]][j]
                    cumulative = self.probs[i-1][k]
                    score = np.float64(cumulative*currentEmit*currentTrans)
                    candidateProbList.append(score)
                    if score > candidateProb:
                        candidateProb = score
                        self.pi2[i-1][j] = self.states[k]
                self.probs[i][j] = max(candidateProbList)
                self.maxStates[i][j] = candidateProbList.index(max(candidateProbList))
    
    def getPath(self):
        """Get the hidden state path with max probability over the sequence data."""
        maxIndex = self.probs[-1].index(max(self.probs[-1]))
        maxState = self.states[maxIndex]
        for i in range(len(self.data)-1, -1, -1):
            # Iterate back through the probabilities matrix to accumulate PI path
            maxIndex = self.maxStates[i][maxIndex]
            self.pi.append(maxState)  
            maxState = self.states[maxIndex]
        self.pi.reverse()
        return self.pi

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
        
    def readRosalindFile(self):
        """Parse a custom-formatted Rosalind file.
        
        The format for this exercise is:
        AABBBAABABAAAABBBBAABBABABBBAABBAAAABABAABBABABBAB
        --------
        A   B
        --------
            A   B
        A   0.194   0.806
        B   0.273   0.727
        """
        with self.doOpen() as fileH:
            content = fileH.readlines()
            string = content[0].strip()
            alpha = content[2].strip().split()
            states = content[4].strip().split()
            transition = list()
            emission = list()
            for i in range(len(states)):
                transition.append(content[7+i].strip().split()[1:])
                emission.append(content[9+len(states)+i].strip().split()[1:])
        return string, alpha, states, transition, emission
        
def main():
    """Calculate the most probable state path PI that emitted the sequence."""
    parser = FastAreader()
    string, alpha, states, transition, emission = parser.readRosalindFile()
    v = Viterbi(string, alpha, states, transition, emission)
    v.calcPathProbs()
    print(''.join(v.getPath()))
    
main()

