#!/usr/bin/env python
"""
This script will transform the riboseq mapping bam to reads P-site bed file. Default reads length is [26,34], and the default p-site is 12 for [26,29], and 13 for [30,34].
Usage:
python P_site_position_of_riboseq_reads.py -i uniq.sorted.bam -o Psite.bed

Author: XKLuan
Date: 20171004
"""
import sys
import pybedtools
import argparse
import os
def main(i,minn,maxn,o):
    ifile = pybedtools.example_bedtool(os.path.abspath('.')+'/'+i)
    ibed = ifile.bam_to_bed(stream=True,bed12=True,split=True)
    for li in ibed:
        ppos = {'26':12,'27':12,'28':12,'29':12,'30':13,'31':13,'32':13,'33':13,'34':13}
        line = str(li)
        chrom = line.split()[0]
        strand = line.split()[5]
        reads_start = int(line.split()[1])
        reads_end = int(line.split()[2])
        reads_exon_count = int(line.split()[9])
        reads_len = line.split()[10]
        reads_name = line.split()[3]
        reads_exon_start = line.split()[11]
        if reads_exon_count == 1:
            reads_length = int(reads_len)
        else:
            reads_length = int(reads_len.split(',')[0]) + int(reads_len.split(',')[1])
#You can select the reads length from this words:        
        if reads_length >=int(minn) and reads_length <= int(maxn):
            if strand == '+' and reads_exon_count == 1:
                p_site = reads_start + ppos[reads_len]
            elif strand == '+' and reads_exon_count ==2:
                reads_length_key = str(int(reads_len.split(',')[0])+int(reads_len.split(',')[1]))
                if int(reads_len.split(',')[0]) <=ppos[reads_length_key]:
                    p_site = reads_start + ppos[reads_length_key] + int(reads_exon_start.split(',')[1]) - int(reads_len.split(',')[0])
                elif int(reads_len.split(',')[0]) >ppos[reads_length_key]:
                    p_site = reads_start + ppos[reads_length_key]
            elif strand == '-' and reads_exon_count == 1:
                p_site = reads_end - (ppos[reads_len]+1)
            elif strand == '-' and reads_exon_count == 2:
                reads_length_key = str(int(reads_len.split(',')[0])+int(reads_len.split(',')[1]))
                if int(reads_len.split(',')[1]) <=ppos[reads_length_key]:
                    p_site = reads_start + int(reads_len.split(',')[0]) + int(reads_len.split(',')[1]) - ppos[reads_length_key]
                elif int(reads_len.split(',')[1]) >ppos[reads_length_key]:
                    p_site = reads_end - ppos[reads_length_key]
            o.write(chrom+'\t'+str(p_site)+'\t'+str(p_site+1)+'\t'+reads_name+'\t'+'50'+'\t'+strand+'\n') 

if __name__=="__main__":
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-i",default="uniq.sorted.bam",help="The riboseq mapping bed file.")
    parser.add_argument('-minn',type=int,default=26,help="The minimum reads length you select.")
    parser.add_argument('-maxn',type=int,default=34,help="The maximum reads length you select.")
    parser.add_argument('-o',nargs="?",type=argparse.FileType('w'),default=sys.stdout,help="The output file.")
    args = parser.parse_args()
    main(i=args.i,minn=args.minn,maxn=args.maxn,o=args.o)
