# -*- coding: utf-8 -*-

"""
Created on Thu Mar 17 21:03:43 2022
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


import os
import sys
from Bio.Seq import Seq
from Bio.Restriction import *


def check_input(args):
    """Validates user input/options"""
    
    # Defaults
    fh = sys.stdin # File handle
    
    if not len(args):
        # Reading from pipe
        if sys.stdin.isatty():
            sys.stderr.write(__doc__)
            sys.exit(1)
            
    elif len(args) == 1:
        # Reading from file
        if not os.path.isfile(args[0]):
            errmsg = 'File not found or not readable: \'{}\'\n'
            sys.stderr.write(errmsg.format(args[0]))
            sys.stderr.write(__doc__)
            sys.exit(1)
            
        fh = open(args[0],'r')
        
    else:
        # For any other input
        sys.stderr.write(__doc__)
        sys.exit(1)
        
    return (fh)


def run(fhandle):
    """Generator that searches for cutsites and yields digest fragments"""
    
    header = True
    for line in fhandle:
        if header == True:
            header = False
            line = line.rstrip()
            line = line + '\tFragmentSequences\n'
            yield line
            
        elif header == False:
            line = line.rstrip()
            seqid,seqre,seqdna = line.split('\t')
            
            # Search input sequences for restriction sites
            restriction_seqdna = Seq(seqdna.upper())
            restriction_seqre = seqre.replace('"','').replace(' ','').split(',')
            restriction_batch = RestrictionBatch(restriction_seqre)
            restriction_cuts = restriction_batch.search(restriction_seqdna)
            
            # Gather all restriction sites
            cutsites = []
            for cuts in restriction_cuts.values():
                cutsites.extend(cuts)
            cutsites.sort()
            
            # Create digest fragments using restriction sites
            fiveprime = 0
            for i in cutsites:
                threeprime = int(i)
                fragmentseq = seqdna[fiveprime:threeprime]
                line = line + f'\t{fragmentseq}'
                fiveprime = int(i)
            fragmentseq = seqdna[fiveprime:]
            line = line + f'\t{fragmentseq}\n'
            
            yield line


def main():
    # Check input
    seqfh = check_input(sys.argv[1:])
    
    # Run
    seqout = run(seqfh)
    
    # Write output
    try:
        with open('output.txt','w') as f:
            for line in seqout:
                f.write(line)
    except IOError:
        pass
    
    # Close
    seqfh.close()
    sys.exit(0)


if __name__ == '__main__':
    main()