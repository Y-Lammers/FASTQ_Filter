#!/usr/bin/env python

# Filter a set of fastq read based on minimum read length and mean quality score

# import the argparse module to handle the input commands
import argparse, itertools, os, sys

# get the commandline arguments for the files and settings
parser = argparse.ArgumentParser(description = 'Filter the input file')

parser.add_argument('--fastq_file', metavar='', type=str,
                        help='Enter the fastq sequence file.')
parser.add_argument('--min_length', metavar='', type=int,
                        help='Enther the minimum sequence length.', default = 0)
parser.add_argument('--median_qual', metavar='', type=int,
                        help='The minimum median quality score for a read.', default = 0)
parser.add_argument('--trim_5', metavar='', type=int,
                        help='Trim x bases from the 5\' end.', default = 0)
parser.add_argument('--trim_3', metavar='', type=int,
                        help='Trim x bases from the 3\' end.', default = 0)
parser.add_argument('--min_quality', metavar='', type=int,
			help='Discard reads with basecalls lower than this threshold.', default = 0)
parser.add_argument('--trim_qual_5', metavar='', type=int,
			help='Trim 5\' bases lower than this threshold.', default = 0)
parser.add_argument('--trim_qual_3', metavar='', type=int,
			help='Trim 3\' bases lower than this threshold.', default = 0)

args = parser.parse_args()

def extract_sequences():

	# open the sequence file submitted by the user, get 
	# the file format and rewind the file
	sequence_file = open(args.fastq_file)

	# create a iterative index of all the headers
        lines = (x[1] for x in itertools.groupby(sequence_file, key=lambda line: line[0] == '@'))

	# walk through the header and obtain the sequences (and quality score if applicable)
        for headers in lines:
		header = headers.next().strip()
		temporary_list, sequence, quality = [line.strip() for line in lines.next()], [], []
	
		# get the multi line sequences and break at the sequence - quality
		# seperator symbol (+)
		while len(temporary_list) > 0:
			line = temporary_list.pop(0)
			if line[0] == '+':
				break
			sequence.append(line)
		quality = temporary_list

		# if the length of the sequences differs from the length of the
		# quality scores (because the quality line starts with a '@') get
		# the next quality line and append it to the previous one
		while len(quality) < len(sequence):
			if len(quality) == 0 and len(sequence) == 1:
				quality.append(headers.next().strip())
			else:
				quality += [line.strip() for line in lines.next()]
		
		# yield the header + sequence
		yield [header, [''.join(sequence), ''.join(quality), [(ord(let)-33) for let in ''.join(quality)]]]


def trim_left_right(sequence, five, three):

	count = 0
	# trim left right
	for item in sequence[1]:	
		sequence[1][count] = item[five:three]
		count += 1
	
	# retun the trimmed sequence
	return sequence


def find_low_qual(qual_string, quality_thresh):


	position_count = 0
	for base in qual_string:
		if base < quality_thresh: position_count += 1
		else: break
	
	return position_count


def main():

	# obtain the sequences and filter them	

	output_file = open(os.path.splitext(args.fastq_file)[0] + '_filtered.fastq','w')

	for sequence in extract_sequences():

		try:

			if args.trim_5 > 0 or args.trim_3 > 0:
				sequence = trim_left_right(sequence, args.trim_5, (len(sequence[1][1])-args.trim_3))

			if args.trim_qual_5 > 0 or args.trim_qual_3 > 0:
				qual_string = sequence[1][2]
				forward = find_low_qual(qual_string, args.trim_qual_5)
				reverse = find_low_qual(qual_string[::-1], args.trim_qual_3)
				sequence = trim_left_right(sequence, forward, (len(sequence[1][1])-reverse))

			if len(sequence[1][0]) < args.min_length or min(sequence[1][2]) < args.min_quality or sorted(sequence[1][2])[(len(sequence[1][2])/2)] < args.median_qual:
				continue

			output_file.write('{0}\n{1}\n+\n{2}\n'.format(sequence[0], sequence[1][0], sequence[1][1]))

		except:
			pass

if __name__ == '__main__':
	main()

