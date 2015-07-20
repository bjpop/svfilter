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
    '''Read target genomic coordinates from input file.''' 
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

def find_intersections(target_coords, sample_id, vcf_writer, record, chrom, pos1, pos2):
    if target_coords is not None:
        if chrom in target_coords:
            chrom_targets = target_coords[chrom]
            intersections = chrom_targets.find(pos1, pos2)
            if len(intersections) > 0:
                annotations = ','.join(intersections)
                record.INFO['hits'] = annotations
                if sample_id is not None:
                    record.INFO['sample'] = sample_id
                vcf_writer.write_record(record)
    else:
        vcf_writer.write_record(record)


def filter_variants_vcf(target_coords, sample_id, variants_filename, custom_filter):
    output_file = variants_filename + '.filtered'
    with open(variants_filename) as variants_file:
        vcf_reader = vcf.Reader(variants_file)
        vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
        for record in vcf_reader:
            filtered_record = custom_filter(record)
            if filtered_record is not None:
                info = filtered_record.INFO
                svtype = info['SVTYPE']
                if svtype == 'BND':
                    # Break end. These can happen on different chromosomes.
                    chrom1 = filtered_record.CHROM 
                    pos1 = filtered_record.POS
                    alt = filtered_record.ALT
                    for item in alt:
                        chrom2, pos2 = parse_bnd_alt(item)
                        if chrom1 == chrom2:
                            # Both break ends are on the same chromosome, we assume one interval
                            # for the pair. 
                            pos_low = min(pos1, pos2)
                            pos_high = max(pos1, pos2)
                            find_intersections(target_coords, sample_id, vcf_writer, filtered_record, chrom1, pos_low, pos_high)
                        else:
                            # Break ends are on different chromosomes.
                            find_intersections(target_coords, sample_id, vcf_writer, filtered_record, chrom1, pos1, pos1)
                            find_intersections(target_coords, sample_id, vcf_writer, filtered_record, chrom2, pos2, pos2)
                else:
                    # INV, DEL, INS, DUP
                    # These all happen on the same chromosome
                    chrom = filtered_record.CHROM
                    start_pos = filtered_record.POS
                    if 'END' in info:
                        end_pos = info['END']
                    else:
                        end_pos = start_pos
                    find_intersections(target_coords, sample_id, vcf_writer, filtered_record, chrom, start_pos, end_pos)

def get_custom_filter(args):
    custom_filter = lambda row: row 
    if args.custom:
        try:
            with open(args.custom) as custom_filter_file:
                contents = custom_filter_file.read()
                exec(contents)
        except Exception as e:
            print(e)
            exit(1)
    return custom_filter


def main():
    args = parse_args()
    start_log(args.log)
    if args.coords:
        target_coords = get_target_coords(args.coords)
    else:
        target_coords = None
    custom_filter = get_custom_filter(args)
    if args.type == 'vcf':
        filter_variants_vcf(target_coords, args.sample, args.variants, custom_filter)

if __name__ == "__main__":
    main()
