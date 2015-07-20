import sys
import re
from bx.intervals.intersection import IntervalTree
from argparse import ArgumentParser
import vcf
from version import version
import logging

DEFAULT_LOG_FILE = "svfilter.log"

def parse_args():
    'Parse the command line arguments for the program.'
    parser = ArgumentParser(
        description="Filter structural variants to a set of genomic coordinates")
    parser.add_argument(
        '--version', action='version', version='%(prog)s ' + version)
    parser.add_argument(
        '--coords', required=True, type=str,
        help='TSV coordinates file for region of interest')
    parser.add_argument(
        '--type', required=True, choices=['lumpy', 'manta', 'socrates'],
        help='Type of structural variant file to filter')
    parser.add_argument(
        '--log', metavar='FILE', type=str, default=DEFAULT_LOG_FILE,
        help='Log progress in FILENAME.')
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
    '''Read target genomic coordinates from input file.''' 
    result = {}
    with open(coords_filename) as coords_file:
        for row in coords_file:
            chrom, start, end = row.split()
            start = int(start)
            end = int(end)
            if chrom in result:
                result[chrom].add(start, end, None)
            else:
                result[chrom] = IntervalTree()
                result[chrom].add(start, end, None)
    return result


'''
            info = record.INFO
            support = info['SU'][0]
            paired_end_support = info['PE'][0]
            split_read_support = info['SR'][0]
            #print((support, paired_end_support, split_read_support))
            if support > 10 and paired_end_support > 5 and split_read_support > 5:
                vcf_writer.write_record(record)
'''

def parse_bnd_alt(alt):
    chrom = alt.chr
    connecting_sequence = alt.connectingSequence
    orientation = alt.orientation
    pos = alt.pos
    remote_orientation = alt.remoteOrientation
    type = alt.type
    within_main_assembly = alt.withinMainAssembly
    return chrom, pos 

def find_intersections(target_coords, vcf_writer, record, chrom, pos1, pos2):
    if chrom in target_coords:
        chrom_targets = target_coords[chrom]
        intersections = chrom_targets.find(pos1, pos2)
        if len(intersections) > 0:
            vcf_writer.write_record(record)


def filter_variants_lumpy(target_coords, variants_filename):
    output_file = variants_filename + '.filtered'
    with open(variants_filename) as variants_file:
        vcf_reader = vcf.Reader(variants_file)
        vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
        for record in vcf_reader:
            info = record.INFO
            svtype = info['SVTYPE']
            if svtype == 'BND':
                chrom1 = record.CHROM 
                pos1 = record.POS
                find_intersections(target_coords, vcf_writer, record, chrom1, pos1, pos1)
                alt = record.ALT
                for item in alt:
                    chrom2, pos2 = parse_bnd_alt(item)
                    find_intersections(target_coords, vcf_writer, record, chrom2, pos2, pos2)
            else:
                chrom = record.CHROM
                start_pos = record.POS
                if 'END' in info:
                    end_pos = info['END']
                else:
                    end_pos = start_pos
                find_intersections(target_coords, vcf_writer, record, chrom, start_pos, end_pos)


def main():
    args = parse_args()
    start_log(args.log)
    target_coords = get_target_coords(args.coords)
    #print(target_coords)
    filter_variants_lumpy(target_coords, args.variants)

if __name__ == "__main__":
    main()
