import argparse as ap
import pysam as ps
import logging

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('-b', '--bam', required = True,
                    help = 'COMBINED_reported_capture_reads_CS5.bam from CCseq pipeline')
parser.add_argument('-n', '--numberOfFragments', default = 2, type = int,
                    help = 'specifies how many fragments the reported reads are allowed to have')
parser.add_argument('--larger', action = 'store_true', default = False,
                    help = 'if set also includes reads that have more than --numberOfFragments fragments')
parser.add_argument('-o', '--outputFile', required = True,
                    help = 'name of the outputfile')
args = parser.parse_args()

inbam = ps.AlignmentFile(args.bam, 'rb')
outbam = ps.AlignmentFile(args.outputFile, 'wb', template = inbam)


readNames = set()
fragmentList = []
fragmentCounts = 0
while True:
    try:
        inputFragment = inbam.__next__()

    except StopIteration:
        if not args.larger:
            if fragmentCounts == args.numberOfFragments:
                for outputFragment in fragmentList:
                    outbam.write(outputFragment)

        else:
            if fragmentCounts >= args.numberOfFragments:
                for outputFragment in fragmentList:
                    outbam.write(outputFragment)

        break

    readname = ':'.join(inputFragment.query_name.split(':')[:-3])

    if readname in readNames:
        fragmentCounts += 1
        fragmentList.append(inputFragment)

    else:
        if not args.larger:
            if fragmentCounts == args.numberOfFragments:
                for outputFragment in fragmentList:
                    outbam.write(outputFragment)

        else:
            if fragmentCounts >= args.numberOfFragments:
                for outputFragment in fragmentList:
                    outbam.write(outputFragment)

        readNames.add(readname)
        fragmentCounts = 1
        fragmentList = [inputFragment]


inbam.close()
outbam.close()