#!/usr/bin/env python
"""
We only selected the 27,28,29,30,31,32,33nt reads and defined the 12nt as the P-site.
"""
import argparse
import os
import sys

def main(gpe,bed,minn,maxn,o):
#def main(gpe,bed,nu,op):
    orf0=[]
    orf1=[]
    orf2=[]
    bed_pos = {}
    for line in bed:
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
        if reads_length >int(minn) and reads_length < int(maxn):
            if strand == '+' and reads_exon_count == 1:
                p_site = reads_start + 13
            elif strand == '+' and reads_exon_count ==2:
                if int(reads_len.split(',')[0]) <=13:
                    p_site = reads_start + 13 + int(reads_exon_start.split(',')[1]) - int(reads_len.split(',')[0])
                elif int(reads_len.split(',')[0]) >13:
                    p_site = reads_start + 13
            elif strand == '-' and reads_exon_count == 1:
                p_site = reads_end - 14
            elif strand == '-' and reads_exon_count == 2:
                if int(reads_len.split(',')[1]) <=13:
                    p_site = reads_start + int(reads_len.split(',')[0]) + int(reads_len.split(',')[1]) - 13
                elif int(reads_len.split(',')[1]) >13:
                    p_site = reads_end - 13
            key_bed = chrom + ':' + strand
            bed_pos.setdefault(key_bed,[]).append(p_site)
            
    for line2 in gpe:
        chrom = line2.split()[1]
        strand = line2.split()[2]
        cds_start = int(line2.split()[5])
        cds_end = int(line2.split()[6])
        exon_nu = int(line2.split()[7])
        exon_start = line2.split()[8]
        exon_end = line2.split()[9]
        inframe = line2.split()[14]
        search_key = chrom + ':' + strand
        if search_key in bed_pos:
            for p_site in bed_pos[search_key]:
                if p_site >= cds_start and p_site <= cds_end:
                    if exon_nu ==1 and strand =='+':
                        orf_type = (p_site - cds_start)%3
                        if orf_type == 0:
                            orf0.append(p_site)
                        elif orf_type == 1:
                            orf1.append(p_site)
                        elif orf_type == 2:
                            orf2.append(p_site)
                    elif exon_nu >1 and strand == "+":
                        for i in range(exon_nu):
                            if int(exon_start.split(',')[i]) < p_site and int(exon_end.split(',')[i]) > p_site:
                                orf_type = (p_site - (int(exon_start.split(',')[i]))%3 + int(inframe.split(',')[i]))%3
                                if orf_type == 0:
                                    orf0.append(p_site)
                                elif orf_type == 1:
                                    orf1.append(p_site)
                                elif orf_type == 2:
                                    orf2.append(p_site)
                    elif exon_nu ==1 and strand == '-':
                        orf_type = (cds_end - 1 - p_site)%3
                        if orf_type == 0:
                            orf0.append(p_site)
                        elif orf_type == 1:
                            orf1.append(p_site)
                        elif orf_type == 2:
                            orf2.append(p_site)
                    elif exon_nu >1 and strand == "-":
                        for i in range(exon_nu):
                            if int(exon_start.split(',')[i]) < p_site and int(exon_end.split(',')[i]) > p_site:
                                orf_type = ((int(exon_end.split(',')[i]) - 1 - p_site)%3 + int(inframe.split(',')[i]))%3
                                if orf_type == 0:
                                    orf0.append(p_site)
                                elif orf_type == 1:
                                    orf1.append(p_site)
                                elif orf_type == 2:
                                    orf2.append(p_site)
#11th as P-site:
    o.write("11th as P-site:\n")
    o.write(str(len(orf2))+'\t'+str(len(orf0))+'\t'+str(len(orf1))+'\n')
#12th as P-site:
    o.write("12th as P-site:\n")
    o.write(str(len(orf1))+'\t'+str(len(orf2))+'\t'+str(len(orf0))+'\n')
#13th as P-site:
    o.write("13th as P-site:\n")
    o.write(str(len(orf0))+'\t'+str(len(orf1))+'\t'+str(len(orf2))+'\n')




if __name__=="__main__":
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('-g',type=argparse.FileType('r'),default="/mnt/share/luanxk/database/GENCODE/hg19.gcCodingGene.gpe",help="The gpe file.")
    parser.add_argument("-b",type=argparse.FileType('r'),default="uniq.sorted.bed12",help="The riboseq mapping bed file.")
#    parser.add_argument("-l",help="The comma separated reads length.")
    parser.add_argument('-minn',type=int,help="The minimum reads length you select.(Not include)")
    parser.add_argument('-maxn',type=int,help="The maximum reads length you select.(Not include)")
    parser.add_argument('-o',nargs="?",type=argparse.FileType('w'),default=sys.stdout,help="The output file.")
    args = parser.parse_args()
    main(gpe=args.g,bed=args.b,minn=args.minn,maxn=args.maxn,o=args.o)
