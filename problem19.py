#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Trevor Ridgley
# BME-205: HMM Problems 16-19
"""
Created on Fri Nov 15 21:03:58 2019

@author: tridgley
Usage: python3 problem19.py < rosalind_ba10d.txt > result19.txt

Outcome Likelihood Problem

Given: A string x, followed by the alphabet Σ from which x was constructed,
       followed by the states States, transition matrix Transition, and emission matrix 
       Emission of an HMM (Σ, States, Transition, Emission).

Return: The probability Pr(x) that the HMM emits x.
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
        self.probs = [[float() for j in range(len(self.states))] for i in range(len(self.data))]
        self.maxStates = [[str() for j in range(len(self.states))] for i in range(len(self.data))]
        self.pi = list()
        self.pi2 = [[str() for j in range(len(self.states))] for i in range(len(self.data))]
        for i in range(len(states)):
            convertTrans = list()
            convertEmit = list()
            # print(emission)
            for j in range(len(alpha)):
                convertEmit.append(float(emission[i][j]))
            for j in range(len(states)):
                convertTrans.append(float(transition[i][j]))
            self.transMatrix[states[i]] = convertTrans
            self.emitMatrix[states[i]] = convertEmit
    
    def calcPathProbs(self):
        """Score all possible paths in the Viterbi graph that emit the data."""
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
                    score = np.longdouble(cumulative*currentTrans*currentEmit)
                    candidateProbList.append(score)
                    if score > candidateProb:
                        candidateProb = score
                        self.pi2[i-1][j] = self.states[k]
                self.probs[i][j] = sum(candidateProbList)
                self.maxStates[i][j] = candidateProbList.index(max(candidateProbList))

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
    """Calculate the probability of an emitted string using Viterbi algorithm."""
    parser = FastAreader()
    string, alpha, states, transition, emission = parser.readRosalindFile()
    v = Viterbi(string, alpha, states, transition, emission)
    v.calcPathProbs()
    print(sum(v.probs[-1]))
    
main()

