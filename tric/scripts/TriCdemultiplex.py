import argparse as ap
import regex
import pysam
import logging
import multiprocessing as mp
from io import StringIO


def bam2fq(read, rn, bc):
    '''
    takes a bamread as dictionary (i.e. read.to_dict) and returns a complementary fastq entry

    :param read:    bam read as dictionary as returned by read.to_dict()
    :param rn:      read number
    :param bc:      barcode sequence

    :return:        bam read in fastq format
    '''

    name = '@' + read['name'] + ' ' + ':'.join([str(rn), 'N', '0', bc])
    seq = read['seq']
    qual = read['qual']

    return '\n'.join([name, seq, '+', qual]) + '\n'


def getKeyAssociatedToPattern(barcodeOfRead1, barcodeOfRead2, regularExpressions, patterns):
    match = False
    key = None
    for re1, re2 in regularExpressions:
        key = '+'.join([patterns[re1], patterns[re2]])
        if re1.match(barcodeOfRead1) and re2.match(barcodeOfRead2):
            match = True
            break

    return key if match else None


def fastqWriter(outhandles, writerQueue, regularExpressions, patterns, lock, processnumber):
    counter = {k: 0 for k in outhandles.keys()}
    filechunks = {k: (StringIO(), StringIO()) for k in outhandles.keys()}
    logging.info('writeProcess-{0} waiting for input'.format(processnumber))
    while True:
        readlist = writerQueue.get()
        if readlist:
            logging.info('writeProcess-{0} processing input'.format(processnumber))
            while len(readlist) > 0:
                BC, B2, r1, r2 = readlist.pop(-1)
                key = getKeyAssociatedToPattern(BC, B2, regularExpressions, patterns)

                if key:
                    for stream, read, readnumber in zip(filechunks[key], (r1, r2), [1, 2]):
                        stream.write(bam2fq(read, readnumber, key))
                    counter[key] += 1

                else:
                    for stream, read, readnumber in zip(filechunks['unknown'], (r1, r2), [1, 2]):
                        stream.write(bam2fq(read, readnumber, 'unknown'))
                    counter['unknown'] += 1

                del BC, B2, r1, r2

            del readlist

            lock.acquire()
            logging.info('writeProcess-{0} writing output'.format(processnumber))
            for key, streams in filechunks.items():
                for filehandle, stream in zip(outhandles[key], streams):
                    try:
                        filehandle.write(stream.getvalue())
                        # needed to empty the file buffer which would
                        # otherwise result in not all reads being written
                        # to files
                        filehandle.flush()

                    except OSError as e:
                        logging.info(repr(e))

                filechunks[key] = (StringIO(), StringIO())

            logging.info('writeProcess-{0} wrote output'.format(processnumber))
            lock.release()

        else:
            logging.info('writeProcess-{0} has nothing left to process'.format(processnumber))
            writerQueue.put([])
            break

    for s1, s2 in filechunks.values():
        s1.close()
        s2.close()

    resultQueue.put(counter)
    logging.info('writeProcess-{0} terminating'.format(processnumber))


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('-i', '--input', required = True,
                    help = 'bamfile containing paired reads as adjacent entries')
parser.add_argument('-b', '--barcodes', required = True,
                    help = '''tab separated file containing the name of the sample, the barcode of read 1 and 
                              the barcode of read 2 in one line e.g. barcode_r1\tbarcode_r2\tid\toutname''')
parser.add_argument('-m', '--mismatches', default = 1, type = int,
                    help = 'number of allowed mismatches in the barcode')
parser.add_argument('-p', '--prefix', default = '.',
                    help = 'prefix for the files to write, can also include directory')
parser.add_argument('-s', '--stats', required = True,
                    help = 'file to write count stats to')
parser.add_argument('-t', '--writerThreads', default = 1, type = int,
                    help = 'number of processes to spawn for writing fastq files')
args = parser.parse_args()

# setting this to large will result in a deadlock of the program due to the thread waiting for an item in
# the queue which will never happen because the pipe used to feed data in the queue is full see also
# https://stackoverflow.com/questions/10028809/maximum-size-for-multiprocessing-queue-item/10029074
BUFFERSIZE = 200000

# maximum memory and average read package size in byte usage to compute maxsize of queue
# the size of one element in the readlist was computed with pympler.asizeof.asizeof((BC, B2, r1, r2))
# AVGSIZEOFQUEUEPACKAGE = 3368 * BUFFERSIZE
# MAXQUEUESIZE = (args.maxMemory*1000000000)//AVGSIZEOFQUEUEPACKAGE - ((args.writerThreads + 1) * 2)
MAXQUEUESIZE = args.writerThreads * 6

logging.info('maximum size of Queue is: {0}'.format(MAXQUEUESIZE))
logging.info('initializing...')
statcounter, names = {}, {}
res, outhandles, patterns = [], {}, {}
with open(args.barcodes, 'r') as bfile:
    for line in bfile:
        b1, b2, i, outname = line.rstrip().split('\t')

        key = '+'.join([b1, b2])
        outhandles[key] = (open('_'.join([args.prefix, outname, '1.fastq']), 'w'),
                           open('_'.join([args.prefix, outname, '2.fastq']), 'w'))

        names[key] = outname
        statcounter[key] = 0

        re1 = regex.compile('(' + b1 + '){s<=' + str(args.mismatches) + '}')
        re2 = regex.compile('(' + b2 + '){s<=' + str(args.mismatches) + '}')
        res.append((re1, re2))

        for re, b in zip(res[-1], [b1, b2]):
            patterns[re] = b

outhandles['unknown'] = (open('_'.join([args.prefix, 'unknown_1.fastq']), 'w'),
                         open('_'.join([args.prefix, 'unknown_2.fastq']), 'w'))

names['unknown'] = 'unknown'
statcounter['unknown'] = 0

writerQueue = mp.Queue(maxsize = MAXQUEUESIZE)
resultQueue = mp.Queue()
writerLock = mp.Lock()
writerProcesses = []
for processnumber in range(args.writerThreads):
    writerProcesses.append(mp.Process(target = fastqWriter,
                           args = [outhandles, writerQueue, res.copy(), patterns.copy(), writerLock, processnumber]))

    writerProcesses[-1].start()


logging.info('demultiplexing...')
with pysam.AlignmentFile(args.input, 'rb', check_sq = False) as infile:
    totalNumberOfReads = 0
    numberOfReadsInBuffer = 0
    readlist = []
    while True:
        try:
            r1 = infile.__next__()
            r2 = infile.__next__().to_dict()

        except StopIteration:
            writerQueue.put(readlist)
            del readlist
            writerQueue.put([])
            break

        BC = r1.get_tag('BC')
        B2 = r1.get_tag('B2')
        r1 = r1.to_dict()

        r1.pop('tags')
        r2.pop('tags')

        readlist.append((BC, B2, r1, r2))
        totalNumberOfReads += 1
        numberOfReadsInBuffer += 1

        if numberOfReadsInBuffer == BUFFERSIZE:
            writerQueue.put(readlist)
            del readlist
            numberOfReadsInBuffer = 0
            readlist = []

        if totalNumberOfReads%1000000 == 0:
            logging.info('processed %i read pairs' % totalNumberOfReads)


logging.info('waiting for writerprocess to finish')
for process in writerProcesses:
    process.join()

del writerQueue, writerProcesses


logging.info('retrieving demux statistics')
resultQueue.put({})
while True:
    writerProcessCounter = resultQueue.get()
    if writerProcessCounter:
        for key, readCount in writerProcessCounter.items():
            statcounter[key] += readCount

    else:
        break

logging.info('writing statsfile...')
with open(args.stats, 'w') as sfile:
    sfile.write('file\t#readpairs\n')

    for key, count in statcounter.items():
        sfile.write('\t'.join([names[key], str(count)]) + '\n')

for handles in outhandles.values():
    for handle in handles:
        handle.close()

logging.info('done')
