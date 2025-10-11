#!/usr/bin/env python3
"""Simple FASTA reader and pretty-printer.

Usage: python3 read_fq.py path/to/file.fasta

This will print each sequence header and its sequence (wrapped at 60 chars).
"""

import sys


def read_fasta(path):
	"""Yield (header, sequence) tuples from a FASTA file."""
	header = None
	seq_lines = []
	with open(path, 'r') as f:
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


def pretty_print_fasta(path, wrap=60):
	for hdr, seq in read_fasta(path):
		print(f'>{hdr}')
		for i in range(0, len(seq), wrap):
			print(seq[i:i+wrap])


def main():
	if len(sys.argv) != 2:
		print('Usage: python3 read_fq.py path/to/file.fasta')
		sys.exit(1)
	pretty_print_fasta(sys.argv[1])


if __name__ == '__main__':
	main()

