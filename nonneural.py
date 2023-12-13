#!/usr/bin/env python3
"""
Non-neural baseline system for the SIGMORPHON 2020 Shared Task 0.
Author: Mans Hulden
Modified by: Tiago Pimentel
Modified by: Jordan Kodner
Modified by: Omer Goldman
Last Update: 22/03/2021
"""

import os, argparse
from functools import wraps
from warnings import warn

# global constants
PRE = 'pre'
SUF = 'suf'

def hamming(s, t):
    """
    Arguments:
        s,t -- two strings to calculate hamming distance for

    Returns:
        Hamming distance between input strings (int)
    """

    return sum(1 for x,y in zip(s,t) if x != y)

def halign(s, t):
    """
    Arguments:
        s,t -- strings to align by hamming distance

    Returns:
        Newly aligned strings with leading and trailing _ characters added as needed.

    Aligns two input strings by 'sliding' one along the other to find the alignment with
    the lowest hamming distance.
    """

    mindist = len(s) + len(t) + 1 # init mindist high enough to ensure best_s and best_t are set

    for i in range(len(s)+len(t)+1): # loop through all alignments using one variable
        left_pad = min(i, len(t)) # increase left pad of s until i >= len(t)
        right_pad = max(0, i-len(t)) # increase right pad of t once i >= len(t)

        padded_s = '_'*left_pad + s + '_'*(len(t)-left_pad) # pad strings
        padded_t = '_'*(len(s)-right_pad) + t + '_'*right_pad

        curdist = hamming(padded_s, padded_t) # check hamming dist of current alignment
        if curdist < mindist: # save it if better than previous ones
            best_s = padded_s
            best_t = padded_t
            mindist = curdist
    
    zipped = [z for z in zip(best_s, best_t) if z != ('_', '_')] # remove extra leading/trailing pads
    new_s = ''.join([z[0] for z in zipped]) # join tuple elements back into strings
    new_t = ''.join([z[1] for z in zipped])
    return new_s, new_t

def levenshtein(s, t, inscost = 1.0, delcost = 1.0, substcost = 1.0):
    """
    Arguments:
        s,t -- input strings to align with levenshtein algorithm

        inscost, delcost, substcost
            -- (default values 1.0) costs to use respectively for insertion, deletion and
               substitution (in input string s)

    Returns:
        Newly aligned input s (str), newly aligned output t (str), cost required to make get
        from s to t (int)

    Recursive implementation of the levenshtein algorithm to find the least expensive way to
    get from s to t using given costs for insertion, deletion, and substitution.
    """

    @memolrec
    def lrec(spast, tpast, srem, trem, cost):
        # Base case:
        #   either s or t has been fully processed, insert/delete remaining characters and return
        if len(srem) == 0:
            return spast + '_'*len(trem), tpast + trem, '', '', cost + len(trem)*inscost
        if len(trem) == 0:
            return spast + srem, tpast + '_'*len(srem), '', '', cost + len(srem)*delcost
        
        # Recursive case:

        curcharcost = 0.0 # 'substitution' cost 0 if chars at current pos are the same
        if srem[0] != trem[0]:
            curcharcost = substcost
        
        # Return minimum cost of substitution, inserting, and deleting at current position
        return min(
            (
                lrec(spast+srem[0], tpast+trem[0], srem[1:], trem[1:], cost+curcharcost),
                lrec(spast+'_', tpast+trem[0], srem, trem[1:], cost+inscost),
                lrec(spast+srem[0], tpast+'_', srem[1:], trem, cost+delcost)
            ),
            key=lambda x: x[4]
        )
    
    res = lrec('', '', s, t, 0) # make starting call to recursive function
    return res[0], res[1], res[4] # only need aligned strings + cost

def memolrec(func):
    """Wrapper function/memoizer for recursive levenshtein implementation. Returns 'decorated' version
    of levenshtein."""

    cache = {} # dictionary to keep track of previously calculated values

    @wraps(func)
    def wrap(spast, tpast, srem, trem, cost):
        # check whether values for aligning srem and trem previously cached
        if (srem, trem) not in cache: # if not, calculate and add to dictionary
            res = func('', '', srem, trem, 0)
            cache[(srem, trem)] = (res[0], res[1], res[4])
        
        aln_srem, aln_trem, rem_cost = cache[(srem, trem)] # retrieve values for srem and trem from dictionary
        return spast+aln_srem, tpast+aln_trem, '', '', cost + rem_cost # add them to previous values + return
    
    return wrap # return decorated function

def alignprs(lemma, form):
    """
    Arguments:
        lemma -- string representing a 'root'/'basic' form of a word
        form  -- string representing derived form of lemma with added morphological features

    Returns:
        6 tuple of strings representing approximate breakdowns of lemma and form into prefix, suffix, and root

              p   r   s
              _________
        lemma 0 | 1 | 2
              ---------
        form  3 | 4 | 5

    Uses trailing and leading 'blank' characters in aligned lemma and form to approximately separate strings into
    prefix, root, and suffix substrings
    """

    # align with preference for insertion/deletion
    aln_lemma, aln_form, _ = levenshtein(lemma, form, substcost=1.1) # unpack only aligned forms

    pre_len = max(numleadingsyms(aln_lemma), numleadingsyms(aln_form))
    suf_len = max(numtrailingsyms(aln_lemma), numtrailingsyms(aln_form))

    lem_pr, lem_rt, lem_sf = aln_lemma[:pre_len], aln_lemma[pre_len:len(aln_lemma)-suf_len], aln_lemma[len(aln_lemma)-suf_len:]
    frm_pr, frm_rt, frm_sf = aln_form[:pre_len], aln_form[pre_len:len(aln_form)-suf_len], aln_form[len(aln_form)-suf_len:]

    return lem_pr, lem_rt, lem_sf, frm_pr, frm_rt, frm_sf

def prefix_suffix_rules_get(lemma, form):
    """
    Arguments:
        lemma -- string representing a 'root'/'basic' form of a word
        form  -- string representing derived form of lemma with added morphological features

    Returns:
        rules -- length-2 list of sets, where index 0 is the prefix rule set and 1 is the suffix set.
                 the rule sets are sets of 2-tupes containing inputs and outputs for possible affix
                 rules used to turn lemma into form

    Aligns lemma and form using alignprs, and identifies all possible prefixing/suffixing rules that
    convert lemma to form. Assumes suffixing to be a slightly more complex/elaborate process than
    prefixing (can be applied to reverse strings for languages which prefer prefixing).
    """

    lp, lr, ls, fp, fr, fs = alignprs(lemma, form) # get approximate prefix, root, and suffix of lemma+form

    rules = {PRE: set(), SUF: set()} # store rules in dictionary so redundant code can be in for-loops later on
    
    inpre = '<' + lp.replace('_', '') # removed if-statement checking whether prefix strings were at least
    outpre = '<' + fp.replace('_', '') # length 0, seemed unnecessary
    for i in range(len(lr)):
        rules[PRE].add((inpre + lr[:i].replace('_',''), outpre + fr[:i].replace('_','')))
    
    insuf = lr + ls + '>' # suffixing rules loop through entire string, not just the root
    outsuf = fr + fs + '>'
    for i in range(len(insuf)):
        in_clean = insuf[i:].replace('_','')
        out_clean = outsuf[i:].replace('_','')

        if len(in_clean) == 5 and is_ecer(lemma):
            rules[SUF].add((in_clean[0] + 'C' + in_clean[2:], replace_first_inst(out_clean, in_clean[1], 'C')))

        rules[SUF].add((in_clean, out_clean))
    
    return rules

def apply_best_rule(lemma, msd, allrules, keys=None, debug=False):
    """
    Arguments:
        lemma -- 'root' or 'base' form of word to transform

        msd   -- unimorph string representing desired features in derived form of lemma

        allrules -- 2-length list of dictionaries mapping all msds to possible affixing rules found
                    in the language. prefix rules at 0 index and suffix rules at 1

        keys -- (default: None) the alternative functions to use when comparing rules against one
                another. must be a dictionary, where the 'pre' key maps to the prefix rule func and
                'suf' to the suffix func. when this argument is None, compares prefixing rules
                using frequency and suffixing rules using input length, with ties broken by frequency,
                and then output length.

    Applies the longest-matching suffix-changing rule given an input form and the MSD.
    """

    if keys == None:
        keys = {PRE: lambda x: x[2], SUF: lambda x: (len(x[0]), x[2], len(x[1]))}

    affixes = ('pre', 'suf')
    
    if debug: print(f'Lemma: {lemma}\nFeatures: {msd}')
    
    base = '<' + lemma + '>' # surrounding characters anchor prefix and suffix rules to the ends of the word

    for fix in (PRE, SUF):
        if msd not in allrules[fix]:
            continue # do nothing if no X-fixing rules exist for this form

        applicablerules = [(rule[0],rule[1],freq) for rule,freq in allrules[fix][msd].items() if rule[0] in base]
        if fix == SUF and is_ecer(lemma):
            for rule,freq in allrules[fix][msd].items():
                if rule[0].replace('C', lemma[-3]) in base:
                    applicablerules.append((rule[0].replace('C', lemma[-3]), rule[1].replace('C', lemma[-3]), freq))

        if len(applicablerules) > 0:
            bestrule = max(applicablerules, key=keys[fix])
            base = base.replace(bestrule[0], bestrule[1])
            if debug: print(f'Applicable {fix}fixing rules:\n{applicablerules}\nUsing: {bestrule}')
    
    return base[1:-1] # trim boundary chars

def is_ecer(lemma):
    return lemma[-4] in 'ée' and lemma[-3] not in 'aâeéèêiïîoôuû' and lemma[-2:] == 'er'

def numleadingsyms(s, sym='_'):
    return len(s) - len(s.lstrip(sym))

def numtrailingsyms(s, sym='_'):
    return len(s) - len(s.rstrip(sym))

def is_in_subset(lemma, subset, wordmap):
    subset_opts_map = {             # Use dictionary to define which categories are included
        'mod': ['-', 'dtd', 'obs'], # in which subsets
        'no_old': ['-', 'dtd', 'obs', 'M'],
        'no_dtd': ['-', 'obs'],
        'no_obs': ['-', 'dtd'],
        'mod_strict': ['-']
    }

    if subset is None or subset == 'all': # If there was no subset argument, do not filter the data
        return True
    
    if '.' in subset: # Include option to run on certain prefixes or suffixes
        fix, substr = subset.split('.', maxsplit=1)
        if fix == PRE:
            return lemma[:len(substr)] == substr
        if fix == SUF:
            return lemma[-len(substr):] == substr
        warn(f'{fix} is not a valid affix type. Running script on unfiltered data set.')
        return True
    
    if subset in subset_opts_map.keys():
        return wordmap[lemma] in subset_opts_map[subset]
    
    warn('Invalid argument passed to -T -V or -S. Running script on unfiltered data set.')
    return True # Run as though no subset was passed if the subset arg invalid

def load_word_map():
    # programmatically generate list of CSVs
    files = [ 'marked-fra.{}.csv'.format(x) for x in ('dev', 'trn', 'tst') ]

    wordmap = {} # Initialize empty dictionary
    for f in files: # Loop through all the files
        try:
            with open(os.path.join('data', f)) as mapfile:
                for line in mapfile:
                    if line.strip() != '': # Lines are all <word>,<category>
                        word, fam = line.strip().split(',')
                        wordmap[word] = fam # Dictionary maps words to their categories
        except FileNotFoundError: # Handle error if one of the files does not exist
            print('File {} not yet generated.'.format(f))
    
    return wordmap

def replace_first_inst(s, to_replace, replace_with):
    t = ''
    for c in s:
        if c == to_replace:
            t += replace_with
            rem = len(s)-len(t)
            t += s[-rem:]
            return t
        else:
            t += c
    return t

def main(parsed):
    # Look how fun! the argparse library lets you access your options/arguments as object attributes!
    # Isn't it so much nicer than the ugly, ugly getopt.gnu_getopt() line? and now the function isn't all
    # cluttered with argument parsing and setting defaults

    out, evl_ext, debug = parsed.out, parsed.eval, parsed.debug
    trn_set, eval_set, subset = parsed.trn_set, parsed.eval_set, parsed.subset
    path = parsed.path

    if debug:
        evl_ext = '.dbg'
    if not subset is None:
        trn_set = subset
        eval_set = subset

    totalavg, numlang = 0.0, 0
    wordmap = load_word_map()

    for lang in [os.path.splitext(d)[0] for d in os.listdir(path) if d[-4:] == '.trn']:

        lang_fname = os.path.join(path, lang) # makes some lines shorter/cleaner

        if not os.path.isfile(lang_fname+'.trn') or not os.path.isfile(lang_fname+evl_ext):
            continue # don't do anything if trn or dev/tst files don't exist in the path

        with open(lang_fname+'.trn', 'r', encoding='utf8') as trnfile:
            # directly extract lines as 3-tuples of (lemma, msd, form)
            trnlines = [line.strip().split('\t') for line in trnfile if line != '\n']

        prefbias,suffbias = 0,0
        for lemma,_,form in trnlines:
            hal_lemma, hal_form = halign(lemma, form)

            prefbias += numleadingsyms(hal_lemma) + numleadingsyms(hal_form)
            suffbias += numtrailingsyms(hal_lemma) + numtrailingsyms(hal_form)
        
        slicedir = 1 # instead of several if-statements, slice every string before and after applying algorithm
        if prefbias > suffbias:
            slicedir = -1
        
        allrules = {PRE: {}, SUF: {}} # init empty dictionaries to store prefix and suffix rules
        for lemma,msd,form in trnlines:

            if is_in_subset(lemma, trn_set, wordmap):
                lemma = lemma[::slicedir] # slice string forward/backward depending on affixing bias
                form = form[::slicedir]

                rules = prefix_suffix_rules_get(lemma, form)

                for fix in (PRE, SUF):

                    if msd not in allrules[fix]: # init empty freq dictionary for new msd
                        allrules[fix][msd] = {}
                
                    for rule in rules[fix]:
                        if rule not in allrules[fix][msd]: # add freq entry for new rule
                            allrules[fix][msd][rule] = 0
                    
                        allrules[fix][msd][rule] += 1 # increase freq by 1 for every repeat rule
        
        # run model on dev file, calculate accuracy
        with open(lang_fname+evl_ext, 'r', encoding='utf8') as evalfile:
            evallines = [line.strip().split('\t') for line in evalfile if line != '\n']

        numcorrect = 0
        numguesses = 0

        if out: outfile = open(lang_fname+'.out', 'w', encoding='utf8')

        for lemma,msd,correct in evallines:
            if is_in_subset(lemma, eval_set, wordmap):
                lemma = lemma[::slicedir]
                outform = apply_best_rule(lemma, msd, allrules, debug=debug)

                lemma = lemma[::slicedir]
                outform = outform[::slicedir]

                if outform == correct:
                    numcorrect += 1
                numguesses += 1
            
                if out:
                    outfile.write('\t'.join((lemma,msd,correct,outform))+'\n')
        
        if out: outfile.close()
        
        pctcorrect = numcorrect / float(numguesses)
        print(f'{lang} accuracy: {pctcorrect:.2%}')

        totalavg += pctcorrect
        numlang += 1
    
    totalavg /= float(numlang)
    print(f'average accuracy: {totalavg:.2%}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='NonNeuralMyTeam',
                                     description='Our cleaned and very slightly edited version'
                                        + ' of the original nonneural.py code.',
                                     epilog='Generating this help message was done with the argparse library'
                                        + ' instead of manually with a bunch of print-statements.')
    
    parser.add_argument('-o', '--output',
                        dest='out',
                        action='store_true',
                        help='generate output files with guesses. files are written to the same place as'
                            + ' the path argument (or default value if no path was specified) under <lang>.out')
    parser.add_argument('-t', '--test',
                        dest='eval',
                        action='store_const',
                        const='.tst',
                        default='.dev',
                        help='evaluate models on test instead of dev data.')
    parser.add_argument('-d', '--debug',
                        dest='debug',
                        action='store_true',
                        help='evaluate on debug (or subset when specified) instead of dev and print debug statements')
    parser.add_argument('-T', '--trn-set',
                        dest='trn_set',
                        default='all',
                        help='subset of data to train on. must specify one of the following: all, mod, no_old,'
                        + 'no_obs, no_dtd, pre.<prefix string>, suf.<suffix string>. defaults to \'all\'.')
    parser.add_argument('-V', '--eval-set',
                        dest='eval_set',
                        help='subset of data to evaluate on. categories are the same as for the -T arg.'
                        + ' defaults to \'all\'')
    parser.add_argument('-S', '--subset',
                        dest='subset',
                        help='subset of data to work with for both training and evaluation. overwrites arguments passed'
                        + ' to -T and -V. when nothing is passed, uses the values of -T and -V, so default behavior is'
                        + ' to run on the unfiltered data, but \'mod\' is recommended for the best results.')
    parser.add_argument('-p', '--path',
                        dest='path',
                        default='data',
                        help='path to the directory containing data files. defaults to \'data\'.')
    
    main(parser.parse_args())