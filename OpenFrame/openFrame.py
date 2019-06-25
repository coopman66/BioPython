# openFrame.py
'''
Determine the open frame of a cDNA sequence, starting with a MET AA, and 
ending with a stop codon.

Author: John Cooper Hopkin
Git Repo: coopman66/BioPython
6-25-2019
Version: 0.1.0
'''

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re


def findFrames(NTLength, dataFile='insulin.gb', fileType='genbank'):
    cdna = SeqIO.read(dataFile, fileType)
    # print(cdna)
    # print(len(cdna)//3)

    openFrames = []
    for match in re.finditer(r'ATG(\w{3})*?(TAA|TAG|TGA)', str(cdna.seq), re.IGNORECASE):
        frame = Seq(match.group(0), IUPAC.unambiguous_dna)
        if len(match.group(0)) > NTLength:
            openFrames.append((frame, frame.translate(1, to_stop=True), match.start(), match.end(), len(frame)))

    return list(sorted(openFrames, key=lambda n: n[4], reverse=True))


if __name__ == '__main__':

    i = 1
    for frame in findFrames(100):
        print(f'{i}:\tcDNA:\t\t{frame[0][:30]}...{frame[0][-15:]}')
        print(f'\tProtSeq:\t{frame[1][:30]}...{frame[1][-15:]}')
        print(f'\tStart: {frame[2]+1}\t\tEnd: {frame[3]+1}\t\tLength: {frame[4]}')
        i += 1





