#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Trevor Ridgley
# BME-205: HMM Problems 16-19
"""
Created on Mon Nov 11 01:18:03 2019

@author: tridgley
Usage: python3 problem17.py < rosalind_ba10b.txt > result17.txt

Probability of an Outcome Given a Hidden Path Problem

Given: A string x, followed by the alphabet Σ from which x was constructed, 
       followed by a hidden path π, followed by the states States and emission 
       matrix Emission of an HMM (Σ, States, Transition, Emission).

Return: The conditional probability Pr(x|π) that string x will be emitted by 
        the HMM given the hidden path π.
"""
import sys

class Markov:
    """A Markovian graph"""
    
    def __init__(self, string, alpha, seq, states, probs):
        """Create a Markovian transition matrix object"""
        self.transMatrix = dict()
        self.emit = string
        self.alpha = alpha
        self.pi = seq
        self.states = states
        for i in range(len(states)):
            convertedProbs = list()
            for j in range(len(alpha)):
                convertedProbs.append(float(probs[i][j]))
            self.transMatrix[states[i]] = convertedProbs
    
    def getCondProb(self):
        """Return the probability of a string using emitted + transition probabilities."""
        prob = 1
        for i in range(0, len(self.emit)):
            transition = self.pi[i]
            emission = self.alpha.index(self.emit[i])
            prob *= self.transMatrix[transition][emission]
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
            string = content[0].strip()
            alpha = content[2].strip().split()
            seq = content[4].strip()
            states = content[6].strip().split()
            probs = list()
            for i in range(len(states)):
                probs.append(content[9+i].strip().split()[1:])
        return string, alpha, seq, states, probs
        
def main():
    """Calculate the probability of a sequence of hidden states using transision matrix."""
    parser = FastAreader()
    string, alpha, seq, states, probs = parser.readRosalindFile()
    mark = Markov(string, alpha, seq, states, probs)
    print(mark.getCondProb())
    
main()

