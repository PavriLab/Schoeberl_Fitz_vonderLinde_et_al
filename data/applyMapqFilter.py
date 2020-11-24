import argparse as ap
import pysam as ps
import logging

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('--input', '-i', required=True,
                    help='input samfile')
parser.add_argument('--mapq', '-q', default=30, type=int,
                    help='mapping quality cutoff')
parser.add_argument('--output', '-o', required=True,
                    help='output bamfile')
args = parser.parse_args()

with ps.AlignmentFile(args.input, 'r') as in_sam:
    out_bam = ps.AlignmentFile(args.output,
                               'wb',
                               threads=2,
                               header=in_sam.header)

    while True:
        try:
            alnseg = in_sam.__next__()

        except StopIteration:
            break

        if alnseg.mapq < args.mapq:
            tmp = ps.AlignedSegment(header=out_bam.header)
            tmp.query_name = alnseg.query_name
            tmp.query_sequence = alnseg.query_sequence
            tmp.flag = 4
            tmp.query_qualities = alnseg.query_qualities
            tmp.set_tags([('YT', 'UU')])
            out_bam.write(tmp)

        else:
            out_bam.write(alnseg)

out_bam.close()
