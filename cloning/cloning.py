# -*- coding: utf-8 -*-

"""
Created on Sun Mar 20 19:13:56 2022
Created by @Kung-Fuzi

A script to find restriction sites in a list of sequences and generate corresponding 
digest fragments. Sequences should be provided as a tab-delimited text file formatted 
as follows:
    ID\tRestrictionEnzymes\tSequence\n

Sequences will be output as a tab-delimited text file that is formatted as follows:
    ID\tRestrictionEnzymes\tSequence\tFragmentSequences...\n

Usage:
    python restrictiondigest.py <sequence file>
    
Example:
    python restrictiondigest.py sequences.txt
"""


