# gbDesign.py
'''
Finds primers for genbank files

1) create vector and cDNA sequences
2) Find promoters/terminators
3) Determine promoter/terminator that user wants to use
4) Determine strand, reverse_comp if neccessary, and slice sequence for primer design

TODO:
    *If user can't find Promoter/Terminator, allow for manual creation
    *Find Open Frame of cDNA
    *Create Primer
    *Allow user to select RE vendor
    *Warn user if RE is in an important location (his tag, lacI operon, etc.)

    *Make a pretty printout of Design process

Author: John Cooper Hopkin
Git Repo: coopman66/BioPython/OpenFrame
6-25-2019
Version: 0.0.2
'''

import openFrame
import sqlite3
from tkinter import filedialog
from Bio import SeqIO, Restriction
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_dna

class ExitError(Exception):
    '''
    Simple Exception class for when the whole program need to quit
    '''
    def __init__(self, value='Exiting...'):
        self.value = value

    def __str__(self):
        return repr(self.value)


class GenBankPrimer:

    def __init__(self, vectorgb, cdnagb):
        self.cdna = SeqIO.read(cdnagb, 'genbank')
        self.vector = SeqIO.read(vectorgb, 'genbank')

        self.promoter, self.terminator = self.find_prom_term()
        self.VectorSeq = Seq(str(self.vector.reverse_complement().seq), IUPAC.ambiguous_dna)

        # Slice the vector from the beginning of the terminator to the end of the promoter
        if self.promoter.location.strand ==-1:
            self.codingvector = Seq(str(self.vector[self.terminator.location.start:self.promoter.location.end].reverse_complement().seq), IUPAC.ambiguous_dna)
        else:
            self.codingvector = Seq(str(self.vector[self.promoter.location.start:self.terminator.location.end].seq), IUPAC.ambiguous_dna)

        # find restriction sites and allow user to select restriction Enzymes
        self.startEnzyme, self.endEnzyme = self.restriction_select()
        


    def manual_prom(self):
        #allow user to define a promoter and terminator
        return 'Dummy'

    def manual_term(self):
        #manually define a terminator
        return 'Dummy'

    def restriction_select(self):
        self.rb = Restriction.RestrictionBatch([], ['B'])
        codingStrandAna = Restriction.Analysis(self.rb, self.codingvector)
        codingStrandAna.print_as('number')
        codingStrandAna.print_that(codingStrandAna.with_N_sites(1))
        print()

        first = False
        while not first:
            print("Enzyme names are case sensitive.")
            firstEnzyme = str(input("Enter the name of the first restriction enzyme you want to use (q to quit): "))
            if firstEnzyme == 'q' or firstEnzyme == 'Q':
                raise ExitError
            for enzyme in self.rb:
                    if str(enzyme) == firstEnzyme:
                        tempEnzyme = enzyme
                        if len(codingStrandAna.full()[tempEnzyme]) == 1:
                            firstEnzyme = enzyme
                            first = True
                            break
            else: 
                print('That is not a valid restriction enzyme for this vector. Did you misspell the name?')

        print()
        second = False
        while not second:
            print("Enzyme names are case sensitive.")
            secondEnzyme = str(input("Enter the name of the second restriction enzyme you want to use (q to quit): "))
            if secondEnzyme == 'q' or secondEnzyme == 'Q':
                raise ExitError
            for enzyme in self.rb:
                    if str(enzyme) == secondEnzyme:
                        tempEnzyme = enzyme
                        if len(codingStrandAna.full()[tempEnzyme]) == 1:
                            secondEnzyme = enzyme
                            second = True
                            break
            else: 
                print('That is not a valid restriction enzyme for this vector. Did you misspell the name?')

        while True:
            print(f'RestrictionEnzymes are: \n\t{str(self.startEnzyme)} at locus: {self.rb.search(self.codingvector)[self.startEnzyme][0]}')
            print(f'\t{str(self.endEnzyme)} at locus: {self.rb.search(self.codingvector)[self.endEnzyme][0]}')
            answer = str(input('Does this look correct? (y or n, q to quit): ')).lower()
            if answer[0] == 'y':
                return firstEnzyme, secondEnzyme
            elif answer[0] == 'n':
                return self.restriction_select()
            elif answer == 'q':
                raise ExitError
            else:
                print('Invalid Input.')

        return firstEnzyme, secondEnzyme

    def prom_select(self, promoters):
        # Default Selection
        promoter = promoters[0]

        # prompts user for promoter choice
        if len(promoters) > 1:
            print('There are several different promoters.')
            for x in range(len(promoters)):
                print(f'{x+1}: {promoters[x].qualifiers["label"]}')
            try:
                num = int(input('Please choose a promoter (0 if none):'))
                if num > len(promoters) or num < 0:
                    raise ValueError
                #check if user wants to define own Promoter
                elif num == 0:
                    answer = input('No promoters found.\nWould you like to define your own? (y or n): ').lower().strip()
                    if answer[0] == 'y':
                        promoter = self.manual_prom()
                    #user Exits
                    elif answer == 'n':
                        raise ExitError
                    else:
                        raise ValueError
                else:
                    promoter = promoters[num-1]

            except ExitError:
                raise ExitError
            except:
                print('Invalid Input.')
                while True:
                    try:
                        num = int(input('Please choose a number (0 if none):'))
                        #check for invalid option
                        if num > len(promoters) or num < 0:
                            raise ValueError
                        #check if user wants to define own Promoter
                        elif num == 0:
                            answer = input('No promoters found.\nWould you like to define your own? (y or n): ').lower().strip()
                            if answer[0] == 'y':
                                promoter = self.manual_prom()
                            # User Exits
                            elif answer == 'n':
                                raise ExitError
                            else:
                                raise ValueError
                        else:
                            promoter = promoters[num-1]
                        break
                    except ExitError:
                        raise ExitError
                    except:
                        print("Invalid Input.")

        return promoter

    def term_select(self, terminators):
        # Default Selection
        terminator = terminators[0]

        #prompts user for choice of terminator
        if len(terminators) > 1:
            print('There are several different terminators.')
            for x in range(len(terminators)):
                print(f'{x+1}: {terminators[x].qualifiers["label"]}')
            try:
                num = int(input('Please choose a terminator (0 if none):'))
                if num > len(terminators) or num < 0:
                    raise ValueError
                #check if user wants to define own terminator
                elif num == 0:
                    answer = input('No terminators found.\nWould you like to define your own? (y or n): ').lower().strip()
                    if answer[0] == 'y':
                        terminator = self.manual_term()
                    #user Exits
                    elif answer == 'n':
                        raise ExitError
                    else:
                        raise ValueError
                else:
                    terminator = terminators[num-1]

            except ExitError:
                raise ExitError
            except:
                print('Invalid Input.')
                while True:
                    try:
                        num = int(input('Please choose a number (0 if none):'))
                        #check for invalid option
                        if num > len(terminators) or num < 0:
                            raise ValueError
                        #check if user wants to define own terminator
                        elif num == 0:
                            answer = input('No terminator found.\nWould you like to define your own? (y or n): ').lower().strip()
                            if answer[0] == 'y':
                                promoter = self.manual_term()
                            # User Exits
                            elif answer == 'n':
                                raise ExitError
                            else:
                                raise ValueError
                        else:
                            terminator = terminatorss[num-1]
                        break
                    except ExitError:
                        raise ExitError
                    except:
                        print("Invalid Input.")

        return terminator

    def prom_term_select(self, promoters, terminators):
        # checks that vector has valid promoters or terminators
        if len(promoters) < 1:
            answer = input('No promoters found.\nWould you like to define your own? (y or n): ').lower().strip()
            if answer[0] == 'y':
                promoter = self.manual_prom()
            else:
                raise ExitError
        else:
            #modularize promoter selection
            promoter = self.prom_select(promoters)

        if len(terminators) < 1:
            print('No terminators found.\nWould you like to define your own? (y or n): ').lower().strip()
            if answer[0] == 'y':
                promoter = self.manual_prom()
            else:
                raise ExitError
        else:
            #modularize terminator selection
            terminator = self.term_select(terminators)

        # double checks user choice
        print(f'Your promotor is: {promoter.qualifiers["label"]}')
        print(f'Your terminator is: {terminator.qualifiers["label"]}')

        answer = str(input('Does this look correct? (y or n): ')).lower().strip()
        if answer[0] == 'y':
            return (promoter, terminator)
        elif answer[0] == 'n':
            print('No more promoters/terminators found.')
            answer = str(input('Would you like to start Promoter/Terminator selection again? (y or n): ')).lower().strip()
            if answer[0] == 'y':
                self.prom_term_select(promoters, terminators)
            elif answer[0] == 'n':
                raise ExitError
        else:
            print('Invalid Input')
            print(f'Your promotor is: {promoter.qualifiers["label"]}')
            print(f'Your terminator is: {terminator.qualifiers["label"]}')
            while True:
                answer = str(input('Does this look correct? (y or n): ')).lower().strip()
                if answer[0] == 'y':
                    return (promoter, terminator)
                elif answer[0] == 'n':
                    print('No more promoters/terminators found.\nPlease choose a new vector.')
                    answer = str(input('Would you like to start Promoter/Terminator selection again? (y or n): ')).lower().strip()
                    if answer[0] == 'y':
                        self.prom_term_select(promoters, terminators)
                    elif answer[0] == 'n':
                        raise ExitError
                else:
                    print('Invalid Input')
                    print(f'Your promotor is: {promoter.qualifiers["label"]}')
                    print(f'Your terminator is: {terminator.qualifiers["label"]}')
        return promoter, terminator

    def find_prom_term(self):
        promoters = []
        terminators = []
        for feature in self.vector.features:
            # print(feature.type)
            if feature.type == 'promoter':
                #add the promoter to a list
                promoters.append(feature)
            if feature.type == 'terminator':
                # add terminator to a list
                terminators.append(feature)

        promoter, terminator = self.prom_term_select(promoters, terminators)

        # Check that promoters are on the same strand
        if promoter.location.strand != terminator.location.strand:
            print('Something is wrong...')
            raise ExitError
        else:
            return promoter, terminator
        


if __name__ == '__main__':
    #cDNA file
    # cdnaFile = str(input("Enter the cDNA filename: "))
    cdnaFile = 'insulin.gb'

    #vector file
    # vectorFile = str(input("Enter the cDNA filename: "))
    vectorFile = 'pet32a.gb'

    try:
        GenBankPrimer(vectorFile, cdnaFile)
    except ExitError:
        print('Exiting....')
