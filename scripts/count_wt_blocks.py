import pysam
import argparse

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
	template_len = abs(read.template_length) < 240

	clipped_bases = 0


	if read.cigartuples is not None:
		#print(read.cigartuples)
		for cigar_op,cigar_len in read.cigartuples:
			if (cigar_op == 4) or (cigar_op == 5):
				clipped_bases += cigar_len

	

	clipped_len = clipped_bases > 40
	
	#return (mapq or read_chr or template_len or clipped_len)
	return (mapq or template_len or clipped_len)


def parse_bam_file(args):

	infile = pysam.AlignmentFile(args.input,"rb")

	i = 1
	wt_blocks = 0
	variant_blocks = 0

	for read in infile:
		if i % 2 == 0:
			try:
				if not (filter_read(read) or filter_read(prev_read)):
					md_tag_R1 = prev_read.get_tag('MD')
					md_tag_R2 = read.get_tag('MD')
					mismatches = count_mismatches(md_tag_R1) + count_mismatches(md_tag_R2)

					#print(mismatches)

					if mismatches == 0:
						wt_blocks += 1
					else:
						variant_blocks += 1
						
			except KeyError:
				pass
				#print("Missing MD tag")

		prev_read = read
		i += 1

	infile.close()

	#print("WT_blocks:\t%d" % wt_blocks)
	#print("Variant_blocks:\t%d" % variant_blocks)
	print("%s\t%d\t%d" % (args.sample,wt_blocks,variant_blocks))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--input', '-i', help='input bam file', required=True)
    parser.add_argument('--sample', '-s', help='sample name', required=True)

    args = parser.parse_args()

    parse_bam_file(args)

