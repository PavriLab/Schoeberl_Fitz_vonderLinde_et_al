import argparse as ap
import pysam as ps
import logging

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('--input', '-i', required=True,
                    help='input sam/bam file')
parser.add_argument('--output', '-o', required=True,
                    help='output bam file')
args = parser.parse_args()

# def checkDirections(readL, merged_reads, header):
#     for i in range(0, len(readL)-1):
#         if not readL[i].is_unmapped and not readL[i+1].is_unmapped:
#             if (not readL[i].is_reverse and not readL[i+1].is_reverse \
#             and readL[i].reference_end - 4 == readL[i+1].reference_start) \
#             or (readL[i].is_reverse and readL[i+1].is_reverse \
#             and readL[i+1].reference_end - 4 == readL[i].reference_start):
#                 # print('Fragments with equal directionality found, second read is marked as unaligned.')
#                 merged_reads +=1
#                 # continue with the readblanking and proper flagging
#                 replacement = ps.AlignedSegment(header=header)
#                 replacement.query_name = readL[i+1].query_name
#                 replacement.query_sequence = readL[i+1].query_sequence
#                 replacement.flag = 4
#                 replacement.query_qualities = readL[i+1].query_qualities
#                 replacement.set_tags([('YD', 'DD'), ('YT', 'UU')]) # the YD DD tag marks reads as part of one continuous fragment
#                 readL[i+1] = replacement
#     return readL, merged_reads

def checkDirections(readL, merged_reads, header):
    for i in range(0, len(readL)-1):
        if not readL[i].is_unmapped and not readL[i+1].is_unmapped:
            if (not readL[i].is_reverse and not readL[i+1].is_reverse) \
            or (readL[i].is_reverse and readL[i+1].is_reverse):
                # print('Fragments with equal directionality found, second read is marked as unaligned.')
                merged_reads +=1
                # continue with the readblanking and proper flagging
                replacement = ps.AlignedSegment(header=header)
                replacement.query_name = readL[i].query_name
                replacement.query_sequence = readL[i].query_sequence
                replacement.flag = 4
                replacement.query_qualities = readL[i].query_qualities
                replacement.set_tags([('YD', 'DD'), ('YT', 'UU')]) # the YD DD tag marks reads as part of one continuous fragment
                readL[i] = replacement
    return readL, merged_reads


with ps.AlignmentFile(args.input, 'r') as in_sam:
    # print('bam read')
    out_bam = ps.AlignmentFile(args.output,
                               'wb',
                               threads=2,
                               header=in_sam.header)

    old_seqn = ['', '']
    readList = []
    counter = 0
    while True:
        try:
            alnseg = in_sam.__next__()

        except StopIteration:
            break
        
        seqn = alnseg.query_name
        seqn = seqn.split(':')[-5:-3]
        if seqn == old_seqn:
            readList.append(alnseg)
        else:
            if len(readList) >= 3:
                # print('3Way found')
                readList, counter = checkDirections(readList, counter, out_bam.header)
            for read in readList: out_bam.write(read)
            readList = [alnseg]
        
        old_seqn = seqn

    logging.info(f'a total of {counter} fragments were discarded as secondary fragments')





out_bam.close()
