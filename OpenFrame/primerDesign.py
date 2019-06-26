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

if __name__ == '__main__':
    #TODO: allow for command line input (or data file)

    

    #cDNA file
    cdnaFile = filedialog.askopenfilename(title="cDNA File", filetypes=())

    #vector file
    vectorFile = filedialog.askopenfile(title='Vector File', filetypes=())
