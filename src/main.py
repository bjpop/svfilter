import sys
import re
from bx.intervals.intersection import IntervalTree
from argparse import ArgumentParser
import vcf
from version import version
import logging
import csv
import re

DEFAULT_LOG_FILE = "svfilter.log"

def parse_args():
    'Parse the command line arguments for the program.'
    parser = ArgumentParser(
        description="Filter structural variants to a set of genomic coordinates")
    parser.add_argument(
        '--version', action='version', version='%(prog)s ' + version)
    parser.add_argument(
        '--coords', required=False, type=str,
        help='TSV coordinates file for region of interest')
    parser.add_argument(
        '--type', required=True, choices=['vcf', 'socrates'],
        help='Type of structural variant file to filter')
    parser.add_argument(
        '--log', metavar='FILE', type=str, default=DEFAULT_LOG_FILE,
        help='Log progress in FILENAME.')
    parser.add_argument(
        '--custom', metavar='FILE', type=str,
        help='Custom variant row filter as a Python script')
    parser.add_argument(
        '--sample', metavar='STR', type=str,
        help='Sample identifier')
    parser.add_argument(
        'variants', metavar='FILE', type=str,
        help='Input file containing variant calls')
    return parser.parse_args()


def start_log(log):
    '''Initiate program logging. If no log file is specified then
    log output goes to DEFAULT_LOGFILE.'''
    logging.basicConfig(
        filename=log,
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S')
    logging.info('program started')
    # Log the command line that was used to run the program
    logging.info('command line: {0}'.format(' '.join(sys.argv)))


def get_target_coords(coords_filename):
    '''Read target genomic coordinates from input file.
    Format is: chrom start_pos end_pos [optional_annotation]
    The optional annotation will typically be the gene name or
    exon id.
    Return a dictionary of interval trees. Dictionary is
    indexed by chromosome. Interval trees store (start, end, annotation)
    for each coordinate span in the input file.''' 
    result = {}
    with open(coords_filename) as coords_file:
        for row in coords_file:
            fields = row.split()
            if len(fields) >= 3:
                chrom, start, end = fields[0:3]
                start = int(start)
                end = int(end)
                if len(fields) >= 4:
                    annotation = fields[3]
                else:
                    annotation = None
                if chrom in result:
                    result[chrom].add(start, end, annotation)
                else:
                    result[chrom] = IntervalTree()
                    result[chrom].add(start, end, annotation)
    return result

def parse_bnd_alt(alt):
    '''Parse a break end (BND) entry'''
    chrom = alt.chr
    connecting_sequence = alt.connectingSequence
    orientation = alt.orientation
    pos = alt.pos
    remote_orientation = alt.remoteOrientation
    type = alt.type
    within_main_assembly = alt.withinMainAssembly
    return chrom, pos 

chrom_pos_regex = re.compile(r'(chr.+):(\d+)')

def parse_chrom_colon_pos(text):
    matches = chrom_pos_regex.match(text)
    if matches is not None: 
        groups = matches.groups()
        if len(groups) == 2:
            return groups[0], int(groups[1])
    # Exit the program if we can't parse the coordinate
    exit("Badly formated genomic coordinate: {}".format(text))


def find_intersections(target_coords, record, chrom, pos1, pos2):
    '''Filter a record based on its coordinates.'''
    if chrom in target_coords:
        chrom_targets = target_coords[chrom]
        intersections = chrom_targets.find(pos1, pos2)
        if len(intersections) > 0:
            return set(intersections) 
    return set() 

class Socrates(object):
    def __init__(self, target_coords):
        self.target_coords = target_coords


    def sample_annotate(self, sample, record):
        '''Append the sample id to a record'''
        record.append(sample)
        yield record


    def coord_filter(self, record):
        c1_realign = record[0]
        c1_anchor = record[3]
        c2_realign = record[12]
        c2_anchor = record[15]
        chrom1, pos1 = parse_chrom_colon_pos(c1_realign)
        chrom2, pos2 = parse_chrom_colon_pos(c1_anchor)
        chrom3, pos3 = parse_chrom_colon_pos(c2_realign)
        chrom4, pos4 = parse_chrom_colon_pos(c2_anchor)
        intersections = set()
  
        if chrom1 == chrom2:
            pos_low = min(pos1, pos2)
            pos_high = max(pos1, pos2)
            intersections.update(find_intersections(self.target_coords, record, chrom1, pos_low, pos_high))
        else:
            intersections.update(find_intersections(self.target_coords, record, chrom1, pos1, pos1))
            intersections.update(find_intersections(self.target_coords, record, chrom2, pos2, pos2))

        if len(intersections) > 0:
            record.append(','.join(intersections))
            yield record


class VCF(object):
    def __init__(self, target_coords):
        self.target_coords = target_coords


    def sample_annotate(self, sample, record):
        record.INFO['sample'] = sample
        yield record


    def coord_filter(self, record):
        intersections = set()
        info = record.INFO
        svtype = info['SVTYPE']
        if svtype == 'BND':
            # Break end. These can happen on different chromosomes.
            chrom1 = record.CHROM 
            pos1 = record.POS
            alt = record.ALT
            for item in alt:
                chrom2, pos2 = parse_bnd_alt(item)
                if chrom1 == chrom2:
                    # Both break ends are on the same chromosome, we assume one interval
                    # for the pair. 
                    pos_low = min(pos1, pos2)
                    pos_high = max(pos1, pos2)
                    intersections.update(find_intersections(self.target_coords, record, chrom1, pos_low, pos_high))
                else:
                    # Break ends are on different chromosomes.
                    intersections.update(find_intersections(self.target_coords, record, chrom1, pos1, pos1))
                    intersections.update(find_intersections(self.target_coords, record, chrom2, pos2, pos2))
        else:
            # INV, DEL, INS, DUP
            # These all happen on the same chromosome
            chrom = record.CHROM
            start_pos = record.POS
            if 'END' in info:
                end_pos = info['END']
            else:
                end_pos = start_pos
            intersections.update(find_intersections(self.target_coords, record, chrom, start_pos, end_pos))

        if len(intersections) > 0:
            record.INFO['hits'] = ','.join(intersections) 
            yield record


def identity_filter(record):
    '''This is a filter which just yields its argument; it is 
    a filter which accepts any input'''
    yield record 

def get_custom_filter(args):
    # default filter is the identity function,
    # will include all rows in the output
    # WARNING: this code execs arbitrary Python code. Do not use this in
    # an untrusted environment, such as a web application!
    custom_filter = identity_filter 
    if args.custom:
        try:
            with open(args.custom) as custom_filter_file:
                contents = custom_filter_file.read()
                exec(contents)
        except Exception as e:
            print(e)
            exit(1)
    return custom_filter


def run_filter(reader, writer, custom_filter, coord_filter, sample_annotate):
    # apply the coord filter first then the custom filter, then sample annotator
    for record_1 in reader:
        for record_2 in coord_filter(record_1):
            for record_3 in custom_filter(record_2): 
                for record_4 in sample_annotate(record_3):
                    writer(record_4)

def main():
    args = parse_args()
    start_log(args.log)
    if args.coords:
        target_coords = get_target_coords(args.coords)
    else:
        target_coords = None
    custom_filter = get_custom_filter(args)
    coord_filter = identity_filter
    sample_annotate = identity_filter
    with open(args.variants) as variants_file:
        if args.type == 'vcf':
            reader = vcf.Reader(variants_file)
            writer = vcf.Writer(sys.stdout, reader).write_record
            vcf_filter = VCF(target_coords)
            if target_coords is not None:
                coord_filter = vcf_filter.coord_filter
            if args.sample:
                sample_annotate = lambda record: vcf_filter.sample_annotate(args.sample, record)
        elif args.type == 'socrates':
            reader = csv.reader(variants_file, delimiter='\t')
            # skip the header row
            reader.next()
            writer = csv.writer(sys.stdout, delimiter='\t').writerow
            socrates_filter = Socrates(target_coords)
            if target_coords is not None:
                coord_filter = socrates_filter.coord_filter
            if args.sample:
                sample_annotate = lambda record: socrates_filter.sample_annotate(args.sample, record)
        run_filter(reader, writer, custom_filter, coord_filter, sample_annotate)

if __name__ == "__main__":
    main()
