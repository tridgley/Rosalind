#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Trevor Ridgley
# BME-205: HMM Problems 16-19
"""
Created on Mon Nov 11 00:34:58 2019

@author: tridgley
Usage: python3 problem16.py < rosalind_ba10a.txt > result16.txt

Probability of a Hidden Path Problem

Given: A hidden path π followed by the states States and transition matrix 
       Transition of an HMM (Σ, States, Transition, Emission).
       
Return: The probability of this path, Pr(π). You may assume that initial 
        probabilities are equal.
"""
import sys

class Markov:
    """A Markovian graph"""
    
    def __init__(self, seq, states, probs):
        """Create a Markovian transition matrix object"""
        self.transMatrix = dict()
        self.pi = seq
        self.states = states
        for i in range(len(states)):
            convertedProbs = list()
            for j in range(len(states)):
                convertedProbs.append(float(probs[i][j]))
            self.transMatrix[states[i]] = convertedProbs
            
    def getProb(self):
        """Return the probability of state sequence using transition matrix."""
        prob = 1/len(self.states)
        for i in range(1, len(self.pi)):
            prev = self.pi[i-1]
            current = self.states.index(self.pi[i])
            prob *= self.transMatrix[prev][current]
        return prob

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
            seq = content[0].strip()
            states = content[2].strip().split()
            probs = list()
            for i in range(len(states)):
                probs.append(content[5+i].strip().split()[1:])
        return seq, states, probs
        
def main():
    """Calculate the probability of a sequence of hidden states using transision matrix."""
    parser = FastAreader()
    seq, states, probs = parser.readRosalindFile()
    mark = Markov(seq, states, probs)
    print(mark.getProb())
    
main()