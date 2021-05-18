#!/usr/bin/env python
from Bio import SeqIO
import sys

#This would be nicer with argparse
if(len(sys.argv) != 2) :
    sys.exit("Usage: %s file.fa" % sys.argv[0])

f=SeqIO.parse(sys.argv[1], "fasta") #We should ensure this exists
for rec in f :
    of=open("%s.fa" % (rec.id), "w")
    SeqIO.write(rec, of, "fasta")
    of.close()
f.close()
