# primerDesign.py
'''
Allow user to pick a cDNA sequence, vector, find the open frame, decide which 
restriction enzyme to use, and finally design the primer

Author: John Cooper Hopkin
Git Repo: coopman66/BioPython
6-25-2019
Version: 0.1.0
'''

import openFrame
import sqlite3
from tkinter import filedialog
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

if __name__ == '__main__':
    #TODO: allow for command line input (or data file)

    #cDNA file
    # cdnaFile = str(input("Enter the cDNA filename: "))
    cdnaFile = 'insulin.gb'

    #vector file
    # vectorFile = str(input("Enter the cDNA filename: "))
    vectorFile = 'pet32a.gb'

    '''i = 1
    for frame in openFrame.findFrames(100, dataFile=cdnaFile):
        print(f'{i}:\tcDNA:\t\t{frame[0][:30]}...{frame[0][-15:]}')
        print(f'\tProtSeq:\t{frame[1][:30]}...{frame[1][-15:]}')
        print(f'\tStart: {frame[2]+1}\t\tEnd: {frame[3]+1}\t\tLength: {frame[4]}')
        i += 1'''

    vector = SeqIO.read(vectorFile, 'genbank')
    # print(vector.seq)

    for seq_record in SeqIO.parse(vectorFile, 'genbank'):
        print(seq_record.id)
        print(repr(seq_record.seq))
        print(len(seq_record))
