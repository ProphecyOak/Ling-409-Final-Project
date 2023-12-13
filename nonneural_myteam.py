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

# global constants
PRE = 0
SUF = 1

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

    rules = [set(), set()] # store rules in length-2 list so redundant code can be in for-loops later on
    
    inpre = '<' + lp.replace('_', '') # removed if-statement checking whether prefix strings were at least
    outpre = '<' + fp.replace('_', '') # length 0, seemed unnecessary
    for i in range(len(lr)):
        rules[PRE].add((inpre + lr[:i].replace('_',''), outpre + fr[:i].replace('_','')))
    
    insuf = lr + ls + '>' # suffixing rules loop through entire string, not just the root
    outsuf = fr + fs + '>'
    for i in range(len(insuf)):
        rules[SUF].add((insuf[i:].replace('_',''), outsuf[i:].replace('_','')))
    
    return rules

def apply_best_rule(lemma, msd, allrules, keys=None):
    """
    Arguments:
        lemma -- 'root' or 'base' form of word to transform

        msd   -- unimorph string representing desired features in derived form of lemma

        allrules -- 2-length list of dictionaries mapping all msds to possible affixing rules found
                    in the language. prefix rules at 0 index and suffix rules at 1

        keys -- (default: None) the alternative functions to use when comparing rules against one
                another. must be at least length-2 list, where 0 index is the prefix rule key and
                index 1 is the suffix rule key. when this argument is None, compares prefixing rules
                using frequency and suffixing rules using input length, with ties broken by frequency,
                and then output length.

    Applies the longest-matching suffix-changing rule given an input form and the MSD.
    """

    if keys == None:
        keys = [lambda x: x[2], lambda x: (len(x[0]), x[2], len(x[1]))]
    
    base = '<' + lemma + '>' # surrounding characters anchor prefix and suffix rules to the ends of the word

    for fix in (PRE, SUF):
        if msd not in allrules[fix]:
            continue # do nothing if no X-fixing rules exist for this form

        applicablerules = [(rule[0],rule[1],freq) for rule,freq in allrules[fix][msd].items()]

        if len(applicablerules) > 0:
            bestrule = max(applicablerules, key=keys[fix])
            base = base.replace(bestrule[0], bestrule[1])
    
    return base[1:-1] # trim boundary chars

def numleadingsyms(s, sym='_'):
    return len(s) - len(s.lstrip(sym))

def numtrailingsyms(s, sym='_'):
    return len(s) - len(s.rstrip(sym))

# =========================== MAIN PROGRAM STARTS HERE =============================

def main(parsed):
    out = parsed.out
    evl_ext = parsed.eval
    path = parsed.path

    totalavg, numlang = 0.0, 0
    for lang in [os.path.splitext(d)[0] for d in os.listdir(path) if d[-4:] == '.trn']:

        lang_fname = os.path.join(path, lang) # makes some lines shorter/cleaner

        if not os.path.isfile(lang_fname+'.trn') or not os.path.isfile(lang_fname+'.dev'):
            continue # don't do anything if trn or dev files don't exist in the path

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
        
        allrules = [{}, {}] # init list of empty dictionaries to store prefix and suffix rules
        for lemma,msd,form in trnlines:

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

        if out: outfile = open(lang_fname+'.out', 'w', encoding='utf8')

        for lemma,msd,correct in evallines:
            lemma = lemma[::slicedir]
            outform = apply_best_rule(lemma, msd, allrules)

            lemma = lemma[::slicedir]
            outform = outform[::slicedir]

            if outform == correct:
                numcorrect += 1
            
            if out:
                outfile.write('\t'.join((lemma,msd,outform))+'\n')
        
        if out: outfile.close()
        
        pctcorrect = numcorrect / float(len(evallines))
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
    parser.add_argument('-p', '--path',
                        dest='path',
                        default='data',
                        help='path to the directory containing data files. defaults to \'data\'.')
    
    main(parser.parse_args())