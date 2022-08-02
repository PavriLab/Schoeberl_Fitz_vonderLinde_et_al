import argparse as ap
import pysam as ps
import logging

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('--input', '-i', required=True,
                    help='input sam/bam file')
parser.add_argument('--output', '-o', required=True,
                    help='output bam file')
parser.add_argument('--filter', '-f', required=True, choices=['conv', 'div', 'fwd', 'rev', 'fwdrev'],
                    help='filter type')
parser.add_argument('--mapq', '-q', required=False, default=0, type=int,
                    help='mapping quality filter value. 0 means no mapq filter.')
args = parser.parse_args()

def checkDirections(readL, merged_reads, low_mapq, header, mapq, filter):
    for i in range(0, len(readL)-1):
        if not readL[i].is_unmapped and not readL[i+1].is_unmapped:
            if readL[i].get_reference_positions()[0] > readL[i+1].get_reference_positions()[0]:
                # switch places
                A = readL[i+1]
                B = readL[i]
            else:
                A = readL[i]
                B = readL[i+1]      
            if not A.is_reverse and B.is_reverse:
                directionality = 'conv'
            elif A.is_reverse and not B.is_reverse:
                directionality = 'div'
            elif not A.is_reverse and not B.is_reverse:
                directionality = 'fwd'
            elif A.is_reverse and B.is_reverse:
                directionality = 'rev'
            if directionality in filter:
                # print('Fragments with specific directionality that should be filtered out found, first read is marked as unaligned.')
                merged_reads +=1
                # continue with the readblanking and proper flagging
                replacement = ps.AlignedSegment(header=header)
                replacement.query_name = readL[i].query_name
                replacement.query_sequence = readL[i].query_sequence
                replacement.flag = 4
                replacement.query_qualities = readL[i].query_qualities
                replacement.set_tags([('YD', 'DD'), ('YT', 'UU')]) # the YD DD tag marks reads as part of one continuous fragment
                readL[i] = replacement
            elif readL[i].mapq < mapq:
                # print('Fragment with poor mapping quality found, is now marked as unaligned.')
                low_mapq +=1
                # continue with the readblanking and proper flagging
                replacement = ps.AlignedSegment(header=header)
                replacement.query_name = readL[i].query_name
                replacement.query_sequence = readL[i].query_sequence
                replacement.flag = 4
                replacement.query_qualities = readL[i].query_qualities
                replacement.set_tags([('YM', 'MM'), ('YT', 'UU')]) # the YM MM tag marks reads as below mapq threshold
                readL[i] = replacement
    return readL, merged_reads, low_mapq


with ps.AlignmentFile(args.input, 'r') as in_sam:
    # print('bam read')
    out_bam = ps.AlignmentFile(args.output,
                               'wb',
                               threads=2,
                               header=in_sam.header)

    old_seqn = ['', '']
    readList = []
    counter = 0
    mapqcounter = 0
    while True:
        try:
            alnseg = in_sam.__next__()

        except (StopIteration, OSError):
            break
        
        seqn = alnseg.query_name
        seqn = seqn.split(':')[-5:-3]
        if seqn == old_seqn:
            readList.append(alnseg)
        else:
            if len(readList) >= 3:
                # print('3Way found')
                readList, counter, mapqcounter = checkDirections(readList, counter, mapqcounter, out_bam.header, args.mapq, ['fwd', 'rev'] if args.filter == 'fwdrev' else [args.filter])
            for read in readList: out_bam.write(read)
            readList = [alnseg]
        
        old_seqn = seqn

    logging.info(f'a total of {counter} fragments were discarded as secondary fragments \na total of {mapqcounter} fragments were further discarded because of low mapq values')





out_bam.close()
