#!/usr/bin/env python3
"""Translate 3 frames of a genomic sequence and search for a protein sequence.

Usage: python3 find_geneX.py sequence.fasta geneX.fasta

Outputs: frame (1/2/3), nucleotide start position (1-based), start codon check (ATG?)
"""
import sys

GENETIC_CODE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}


def read_fasta(path):
    header = None
    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_lines)
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            yield header, ''.join(seq_lines)


def translate(seq, offset=0):
    seq = seq.upper()
    aa = []
    for i in range(offset, len(seq)-2, 3):
        codon = seq[i:i+3]
        aa.append(GENETIC_CODE.get(codon, 'X'))
    return ''.join(aa)


def main():
    if len(sys.argv) != 3:
        print('Usage: python3 find_geneX.py sequence.fasta geneX.fasta')
        sys.exit(1)

    seq_path = sys.argv[1]
    prot_path = sys.argv[2]

    # read genomic sequence (first record)
    seqs = list(read_fasta(seq_path))
    if not seqs:
        print('No sequences found in', seq_path); sys.exit(1)
    seq_header, seq = seqs[0]

    prots = list(read_fasta(prot_path))
    if not prots:
        print('No protein sequences found in', prot_path); sys.exit(1)
    prot_header, prot = prots[0]
    prot = prot.strip()

    found = []
    for frame in range(3):
        aa = translate(seq, offset=frame)
        idx = aa.find(prot)
        if idx != -1:
            # compute nucleotide start (1-based)
            nt_start = frame + idx*3 + 1
            start_codon = seq[nt_start-1:nt_start-1+3].upper()
            found.append((frame+1, nt_start, start_codon, len(prot)))

    if not found:
        print('Protein not found in any frame.')
        # Optionally, show best partial matches
        return

    # report all matches (usually one)
    for framenum, nt_start, start_codon, plen in found:
        print(f'Found protein "{prot_header}" (len {plen}) in frame {framenum}')
        print(f'Nucleotide start (1-based): {nt_start}')
        print(f'Start codon at that position: {start_codon} (ATG?) -> {start_codon == "ATG"}')


if __name__ == '__main__':
    main()
