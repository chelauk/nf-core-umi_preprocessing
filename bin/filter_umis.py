#!/usr/bin/env python
# coding: utf-8
import gzip
import argparse
import os
import sys 

def stream_fastq(fqfile):
    """Read a fastq file and provide an iterable of the sequence ID, the
    full header, the sequence, and the quaity scores.
    Note that the sequence ID is the header up until the first space,
    while the header is the whole header.
    """
    qin = gzip.open(fqfile, 'rb')

    while True:
        header = qin.readline()
        if not header:
            break
        header = header.decode().strip()
        seqidparts = header.split(' ')
        seqid = seqidparts[0]
        seq = qin.readline()
        qualheader = qin.readline()
        qualscores = qin.readline()
        qualscores = qualscores
        yield seqid, header, seq, qualscores

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="validate umi fastq file, remove those reads that are not in R1.fastq.gz")
    parser.add_argument('-u', help='UMI reads file')
    parser.add_argument('-v', help='The shorter validated R1 file')
    args = parser.parse_args()
    qout = gzip.open("{}".format(args.u.replace("fq.gz","val_umi.fq.gz")),'wb')
    R1_generator = stream_fastq(args.v) # the shorter list 
    R2_generator = stream_fastq(args.u) # the (longer) umi list

    def list_cycle(seqid, longer_list):
        for (seqid2, header2, seq2, qualscores2) in longer_list:
            if seqid == seqid2:
                qout.write(header2.encode('utf-8') + seq2 + b"+\n" + qualscores2)
                return True
        return False

    for (seqid, header, seq, qualscores) in R1_generator:
        list_cycle(seqid, R2_generator)

    qout.close()
