import sys, getopt
import os, argparse

DEFAULT_DAT = os.path.join('data', 'fra.dev')
DEFAULT_OUT = os.path.join('data', 'out', 'fra.dev.out')
DEFAULT_DIF = os.path.join('data', 'diff', 'fra.dev.diff')

def get_diffs(datfile, outfile, difffile, incl_all=False):

    with open(datfile) as datfile, open(outfile) as outfile, open(difffile, 'w', encoding='utf8') as difffile:
        outlines = [line.strip().split('\t') for line in outfile if line.strip() != '']
        includelines = [vals[:2] for vals in outlines]
        datlines = [line.strip().split('\t') for line in datfile if line.strip().split('\t')[:2] in includelines]

        for dat, out in zip(datlines, outlines):
            out_guess = out[-1]

            if dat != out:
                mark_bad = 'BAD'
            else:
                mark_bad = ''
            
            if incl_all:
                difffile.write('\t'.join(dat + [out_guess, mark_bad]) + '\n')
            elif dat != out:
                difffile.write('\t'.join(dat + [out_guess]) + '\n')

def main(parsed_args):
    DAT,OUT,DIFF,ALL = parsed_args.datafile,parsed_args.outfile,parsed_args.diff,parsed_args.all

    if DIFF == DEFAULT_DIF and OUT != DEFAULT_OUT:
        data_path = os.path.commonpath([DAT, OUT])
        file_name = os.path.split(OUT)[1].replace('.out', '.diff')
        DIFF = os.path.join(data_path, 'diff', file_name)
    
    get_diffs(DAT, OUT, DIFF, incl_all=ALL)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate difference file between French data file and output from nonneural.py. '
                                     + 'Default behavior is to generate difference file for data/fra.dev and data/out/fra.dev.out')
    parser.add_argument('datafile',
                        nargs='?',
                        default=DEFAULT_DAT,
                        type=str,
                        help='Path to UTF-8 encoded file containing correct French data.')
    parser.add_argument('outfile',
                        nargs='?',
                        default=DEFAULT_OUT,
                        type=str,
                        help='Path to UTF-8 encoded file containing French produced by nonneural.py')
    parser.add_argument('-D', '--diff',
                        dest='diff',
                        default=DEFAULT_DIF,
                        type=str,
                        help='Filepath for difference file to generate. '
                        + 'Defaults to ' + os.path.join('data', 'diff', '<filename>.diff')
                        + ', where <filename> is the outfile name before .out')
    parser.add_argument('-a', '--all',
                        dest='all',
                        action='store_true',
                        help='Difference file includes correctly produced lines and marks incorrect ones with BAD.')

    main(parser.parse_args())