import os
import re
#Reading user inut file
fasta = input("Enter fasta file: ")

def get_filehandle(infile):
    seqs = []
    headers = []
    with open(infile) as f:
        sequence = ""
        header = None
        for line in f:
            if line.startswith('>'):
                headers.append(line[1:-1])
                if header:
                    seqs.append([sequence])
                sequence = ""
                header = line[1:]
            else:
                sequence += line.rstrip()
        seqs.append(sequence)
    return headers, seqs
def get_fasta_lists():
    get_fasta_lists.headers, get_fasta_lists.seqs = get_filehandle(fasta)
    

get_fasta_lists()
headers = get_fasta_lists.headers
seqs = get_fasta_lists.seqs
seqss = [seqs]
flat_seqs = [item for sublist in seqss for item in sublist]
sequence_flat = ''.join(flat_seqs)
import re
start_ind =[m.start() for m in re.finditer('ATG', sequence_flat)]
#print(start_ind)
openRF = []
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
seq = sequence_flat
seqs_rev = "".join(complement.get(base, base) for base in reversed(seq))
start_ind_rev =[m.start() for m in re.finditer('ATG', seqs_rev)]
#print(reverse_complement)
number = input("Enter minimum ORF length: ")
orf_len = int(number)
f = 1
while(f <= 6):
    frame = f
    if frame <= 3:
        for i in start_ind[:len(start_ind) - 6]:
            j = i+frame
            #print(seqs[0][i])
            #print(seqs[0][i:i+50])
            ex = seqs[0][j+orf_len:]
            #print(seqs[0][0:seqs[0].index('TAA')+3])
            stop_codon = []
            stop_codon.append(ex.index('TAA'))
            stop_codon.append(ex.index('TGA'))
            stop_codon.append(ex.index('TAG'))
            stop_codon.sort()
            #print(stop_codon)
            #print("STOP codon list",ex[0:stop_codon[0]+3])
            ORF = seqs[0][j:j+orf_len] + ex[0:stop_codon[0]+3]
            orf_trip = list(map(''.join, zip(*[iter(ORF)]*3)))
            if 'TAA' == orf_trip[-1] or 'TGA' == orf_trip[-1] or 'TAG' == orf_trip[-1]:
                #print(orf_trip)
                openRF.append(f">{headers[0]}| FRAME = {frame} POS = {j} LEN = {len(ORF)}")
                openRF.append(' '.join(orf_trip))
    elif frame >3 and frame <= 6:
        for i in start_ind_rev[:len(start_ind_rev) - 6]:
            j = i+frame
            #print(seqs[0][i])
            #print(seqs[0][i:i+50])
            ex = seqs_rev[j+orf_len:]
            #print(seqs[0][0:seqs[0].index('TAA')+3])
            stop_codon = []
            stop_codon.append(ex.index('TAA'))
            stop_codon.append(ex.index('TGA'))
            stop_codon.append(ex.index('TAG'))
            stop_codon.sort()
            #print(stop_codon)
            #print("STOP codon list",ex[0:stop_codon[0]+3])
            ORF = seqs_rev[j:j+orf_len] + ex[0:stop_codon[0]+3]
            orf_trip = list(map(''.join, zip(*[iter(ORF)]*3)))
            if 'TAA' == orf_trip[-1] or 'TGA' == orf_trip[-1] or 'TAG' == orf_trip[-1]:
                #print(orf_trip)
                openRF.append(f">{headers[0]}| FRAME = {frame} POS = {-j} LEN = {len(ORF)}")
                openRF.append(' '.join(orf_trip))
    f += 1
#print(openRF)
print('\n'.join(openRF))
with open('ORFs.fasta', 'w') as wri:
    wri.write('\n'.join(openRF))

