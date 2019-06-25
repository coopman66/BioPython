from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# parsing a fasta file into sequence objects
'''for seq_record in SeqIO.parse('ls_orchid.fasta','fasta'):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))'''

# parsing a genbank record file into sequence objects
'''for seq_record in SeqIO.parse('insulin.gb', 'genbank'):
    print(seq_record.id)
    print(seq_record.seq)
    print(len(seq_record))'''

'''
cdna = SeqIO.read('insulin.gb', 'genbank')
# print(cdna)
# print(len(cdna)//3)

i = 1
for strand, nuc in [(+1, cdna.seq)]:
    for frame in range(3):
        length = 3 * ((len(cdna)-frame) // 3)  # Multiple of three
        for pro in nuc[frame:frame+length].translate(1).split("*"):
            if len(pro) >= 100 and re.match(r'^\w*?M', str(pro)):
                # find the frame that starts with methianine
                pro = Seq(re.sub(r'^\w*?M', 'M', str(pro), re.IGNORECASE), IUPAC.protein)

                print("%i: %s...%s - length %i, strand %i, frame %i"
                        % (i, pro[:30], pro[-3:], len(pro), strand, frame))
                i += 1
'''

# Seq record objects from a fasta file
'''record = SeqIO.read("NC_005816.fna", 'fasta')
print(record.seq)'''
