import pysam
import argparse

#samfile = pysam.AlignmentFile("downlow1.bam","rb")
#samfile = pysam.AlignmentFile("downhigh1_byname_pa.bam","rb")
#samfile = pysam.AlignmentFile("downlow1_byname_pa.bam","rb")
#clean_reads.py -i downlow2_byname_pa.bam -o test.sam


def check_overlap_region(md_str1, md_str2,overlap_size):

	mismatches = 0
	len_str = ''

	for c in list(md_str2):

		if c.isdigit():
			len_str += c
		else:
			break
			#mismatches += 1
		#print(c)

	#print("R2 position: %d" % int(len_str))

	R2_var_pos = int(len_str)

	return R2_var_pos < overlap_size



def count_mismatches(md_str):

	mismatches = 0

	for c in list(md_str):
		if c.isalpha():
			mismatches += 1
		#print(c)


	return mismatches


def filter_read(read):
	mapq = read.mapping_quality < 50
	#read_chr = read.reference_id != 18
	#template_len = abs(read.template_length) < 250
	template_len = abs(read.template_length) < 240

	clipped_bases = 0


	if read.cigartuples is not None:
		#print(read.cigartuples)
		for cigar_op,cigar_len in read.cigartuples:
			if (cigar_op == 4) or (cigar_op == 5):
				clipped_bases += cigar_len



	clipped_len = clipped_bases > 40

	'''
	if high_clipped_bases:
		print("Filtering read based on clipping length")
		print(read)

	'''

	#print("clipped_len: %d template_len: %d mapping quality: %d" % (clipped_len,read.template_length,read.mapping_quality))

	#return (mapq or read_chr or template_len or clipped_len)
	return (mapq or template_len or clipped_len)


def filter_bam_file(args):

	infile = pysam.AlignmentFile(args.input,"rb")
	outfile = pysam.AlignmentFile(args.out,"w",template=infile)
	overlap_size = int(args.overlap)

	i = 1

	for read in infile:
		if i % 2 == 0:
			try:
				if not (filter_read(read) or filter_read(prev_read)):
					md_tag_R1 = prev_read.get_tag('MD')
					md_tag_R2 = read.get_tag('MD')
					#mismatches = count_mismatches(md_tag_R1) + count_mismatches(md_tag_R2)
					R1_mismatches = count_mismatches(md_tag_R1)
					R2_mismatches = count_mismatches(md_tag_R2)

					mismatches = R1_mismatches + R2_mismatches


					#print("MD tags: %s %s\tRead has mismatches: %d" % (md_tag_R1,md_tag_R2,mismatches))
					#print(mismatches)
					if mismatches < 2:
						#pass
						outfile.write(prev_read)
						outfile.write(read)

					elif R1_mismatches == 1 and R2_mismatches == 1 and check_overlap_region(md_tag_R1, md_tag_R2,overlap_size):
						'''
						print("R1")
						print(prev_read)
						print(prev_read.cigartuples)
						print("R2")
						print(read)
						print(read.cigartuples)
						'''
						#pass

						outfile.write(prev_read)
						outfile.write(read)

			except KeyError:
				pass
				#print("Missing MD tag")

		prev_read = read
		i += 1

	infile.close()
	outfile.close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--input', '-i', help='input bam file', required=True)
    parser.add_argument('--out', '-o', help='output sam file', required=True)
    parser.add_argument('--overlap', '-l', help='read overlap size', default=4)

    args = parser.parse_args()

    #print("Using overlap size: %d" % int(args.overlap))
    filter_bam_file(args)
