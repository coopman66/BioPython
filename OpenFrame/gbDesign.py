# gbDesign.py
'''
Finds primers for genbank files

1) create vector and cDNA sequences
2) Find promoters/terminators
3) Determine promoter/terminator that user wants to use
4) Determine strand, reverse_comp if neccessary, and slice sequence for primer design

TODO:
    *If user can't find Promoter/Terminator, allow for manual creation
    *Slice sequence
    *Find Cloning Sites
    *Have user pick cloning sites
    *Find Open Frame of cDNA
    *Create Primer

    *Make a pretty printout of Design process

Author: John Cooper Hopkin
Git Repo: coopman66/BioPython/OpenFrame
6-25-2019
Version: 0.0.2
'''

import openFrame
import sqlite3
from tkinter import filedialog
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def prom_term_select(promoters, terminators):
    # checks that vector has valid promoters or terminators
    if len(promoters) < 1:
        print('No promoters found.\nPlease choose a new vector.')
        return False, False
    if len(terminators) < 1:
        print('No terminators found.\nPlease choose a new vector.')
        return False, False

    #default, stays this way if only one option for each list
    promoter = promoters[0]
    terminator = terminators[0]

    # prompts user for promoter choice
    if len(promoters) > 1:
        print('There are several different promoters.')
        for x in range(len(promoters)):
            print(f'{x+1}: {promoters[x].qualifiers["label"]}')
        try:
            num = int(input('Please choose a number:'))
            if num > len(promoters):
                raise ValueError
            promoter = promoters[num-1]
        except:
            print('Invalid Input.')
            while True:
                try:
                    num = int(input('Please choose a number:'))
                    if num > len(promoters):
                        raise ValueError
                    promoter = promoters[num-1]
                    break
                except:
                    print('Invalid Input.')

    #prompts user for choice of terminator
    if len(terminators) > 1:
        print('There are several different terminators.')
        for x in range(len(terminators)):
            print(f'{x+1}: {terminators[x].qualifiers["label"]}')
        try:
            num = int(input('Please choose a number:'))
            if num > len(terminators):
                raise ValueError
            terminator = terminators[num-1]
        except:
            print('Invalid Input.')
            while True:
                try:
                    num = int(input('Please choose a number:'))
                    if num > len(terminator):
                        raise ValueError
                    terminator = terminators[num-1]
                    break
                except:
                    print('Invalid Input.')

    # double checks user choice
    print(f'Your promotor is: {promoter.qualifiers["label"]}')
    print(f'Your terminator is: {terminator.qualifiers["label"]}')

    answer = str(input('Does this look correct? (y or n): ')).lower().strip()
    if answer[0] == 'y':
        return (promoter, terminator)
    elif answer[0] == 'n':
        print('No more promoters/terminators found.\nPlease choose a new vector.')
        return False, False
    else:
        print('Invalid Input')
        while True:
            answer = str(input('Does this look correct? (0 to exit): ')).lower().strip()
            if answer[0] == 'y':
                return (promoter, terminator)
            elif answer[0] == 'n':
                print('No more promoters/terminators found.\nPlease choose a new vector.')
                return False, False
            elif answer[0] == '0':
                return False, False
            else:
                print('Invalid Input')

def genBankprimerDesign(cdnaFile, vectorFile):
    cdna = SeqIO.read(cdnaFile, 'genbank')
    vector = SeqIO.read(vectorFile, 'genbank')
    # print(vector.seq)

    promoters = []
    terminators = []
    for feature in vector.features:
        if feature.type == 'promoter':
            #add the promoter to a list
            promoters.append(feature)

        if feature.type == 'terminator':
            # add terminator to a list
            terminators.append(feature)

    # print('length of promoters:', len(promoters))
    # print('length of terminators:', len(terminators))

    promoter, terminator = prom_term_select(promoters, terminators)

    #basically same as error code from prom_term_select
    if not promoter:
        return

    # Check that promoters are on the same strand
    if promoter.location.strand != terminator.location.strand:
        print('Something is wrong...')
        return

    # determine strand that term/prom is on
    if promoter.location.strand == 1 or promoter.location.strand == 0:
        print(promoter.location.start)
        print(terminator.location.end)
    elif promoter.location.strand == -1:
        print(terminator.location.start)
        print(promoter.location.end)
    
    # Slice the vector from the beginning of the terminator to the end of the promoter

    
    

if __name__ == '__main__':
    #cDNA file
    # cdnaFile = str(input("Enter the cDNA filename: "))
    cdnaFile = 'insulin.gb'

    #vector file
    # vectorFile = str(input("Enter the cDNA filename: "))
    vectorFile = 'pet32a.gb'

    genBankprimerDesign(cdnaFile, vectorFile)