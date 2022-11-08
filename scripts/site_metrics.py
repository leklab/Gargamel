import pysam
import argparse


def reset_base_counts():
    base_counts = {}
    base_counts['A'] = 0
    base_counts['C'] = 0
    base_counts['T'] = 0
    base_counts['G'] = 0
    base_counts['N'] = 0

    return base_counts

def site_metrics(args):
    infile = pysam.AlignmentFile(args.input,"rb")
    
    print("chr\tpos\tA\tC\tG\tT")
    for pileupcolumn in infile.pileup(args.contig, int(args.start), int(args.end),min_base_quality=int(args.base_quality),
        truncate=True,max_depth=100000000,stepper='nofilter'):
        
        #pileupcolumn.set_min_base_quality(0)
        #print("\ncoverage at base %d = %d %d" % (pileupcolumn.pos, pileupcolumn.nsegments,pileupcolumn.get_num_aligned()))
        base_counts = reset_base_counts()

        for pileupread in pileupcolumn.pileups:

            if not pileupread.is_del and not pileupread.is_refskip:
            #if True:

                #Can take either query_sequence OR query_alignment_sequence (ignores soft clipped bases)
                #query_alignment_sequence is most likely going to be out of sync with pileupread.query_position?

                base = pileupread.alignment.query_sequence[pileupread.query_position]
                base_quality = pileupread.alignment.query_qualities[pileupread.query_position]
                
                #print(pileupread.query_position)

                '''print('\tbase in read %s = %s %d' % (pileupread.alignment.query_name, 
                    pileupread.alignment.query_sequence[pileupread.query_position],
                    pileupread.alignment.query_qualities[pileupread.query_position]))
                '''

                #if base_quality > args.base_quality:
                    #base_counts[base] += 1
                base_counts[base] += 1
        
        #print(base_counts)
        print("%s\t%d\t%d\t%d\t%d\t%d" %(args.contig,pileupcolumn.pos+1, base_counts['A'], base_counts['C'], base_counts['G'], base_counts['T']))

        #break


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--input', '-i', help='input bam file', required=True)
    parser.add_argument('--contig', '-c', help='start position', required=True)
    parser.add_argument('--start', '-s', help='start position', required=True)
    parser.add_argument('--end', '-e', help='end position', required=True)
    parser.add_argument('--base_quality', '-q', help='base quality threshold', default=0)

    #parser.add_argument('--out', '-o', help='output sam file', required=True)
    #parser.add_argument('--overlap', '-l', help='read overlap size', default=4)
    #"LARGE1_CDS"

    args = parser.parse_args()
    site_metrics(args)
