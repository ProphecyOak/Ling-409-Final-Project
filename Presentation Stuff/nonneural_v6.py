#!/usr/bin/env python3
"""
Non-neural baseline system for the SIGMORPHON 2020 Shared Task 0.
Author: Mans Hulden
Modified by: Tiago Pimentel
Modified by: Jordan Kodner
Modified by: Omer Goldman
Last Update: 22/03/2021
"""

import sys, os, getopt, re
from functools import wraps
from glob import glob
from warnings import warn


def hamming(s,t):
    """
    Arguments:
        s,t -- two strings to calculate hamming distance for

    Returns:
        Hamming distance between input strings (int)
    """

    return sum(1 for x,y in zip(s,t) if x != y)


def halign(s,t):
    """
    Arguments:
        s,t -- strings to align by hamming distance

    Returns:
        Newly aligned strings with leading and trailing _ characters added as needed.

    Aligns two input strings by 'sliding' one along the other to find the alignment with
    the lowest hamming distance.
    """

    minscore = len(s) + len(t) + 1 # Starting "minimum" score 1 above maximum possible score

    # 'Slide' forward string s along string t
    for upad in range(0, len(t)+1):
        upper = '_' * upad + s + (len(t) - upad) * '_'
        lower = len(s) * '_' + t
        score = hamming(upper, lower) # Calculate hamming distance between current alignment
        if score < minscore: # Save best alignment so far if lower score is found
            best_upper = upper
            best_lower = lower
            minscore = score

    # Repeat above algorithm 'sliding' t backwards along s
    for lpad in range(0, len(s)+1):
        upper = len(t) * '_' + s
        lower = (len(s) - lpad) * '_' + t + '_' * lpad
        score = hamming(upper, lower)
        if score < minscore:
            best_upper = upper
            best_lower = lower
            minscore = score

    # Remove unnecessary trailing and leading _ characters from both strings
    zipped = list(zip(best_upper,best_lower))
    newin  = ''.join(i for i,o in zipped if i != '_' or o != '_')
    newout = ''.join(o for i,o in zipped if i != '_' or o != '_')
    return newin, newout


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
        # Base cases:
        #     If either s or t has been fully 'processed', the remaining distance is simply
        #     the cost of inserting/deleting remaining characters in the input string.
        if len(srem) == 0:
            return spast + len(trem) * '_', tpast + trem, '', '', cost + len(trem) * inscost
        if len(trem) == 0:
            return spast + srem, tpast + len(srem) * '_', '', '', cost + len(srem)

        addcost = 0
        if srem[0] != trem[0]: # Cost of 'substituting' a character with itself is 0
            addcost = substcost

        # Return the minimum of substituting, inserting, or deleting at the current character position
        return min((lrec(spast + srem[0], tpast + trem[0], srem[1:], trem[1:], cost + addcost),
                   lrec(spast + '_', tpast + trem[0], srem, trem[1:], cost + inscost),
                   lrec(spast + srem[0], tpast + '_', srem[1:], trem, cost + delcost)),
                   key = lambda x: x[4])

    answer = lrec('', '', s, t, 0) # Starting call to recursive function
    return answer[0],answer[1],answer[4] # Indices 2 and 3 of the returned tuple are empty strings


def memolrec(func):
    """Wrapper function/memoizer for recursive levenshtein implementation. Returns 'decorated' version
    of levenshtein."""
    cache = {} # Initialize empty memoization dictionary
    @wraps(func)
    def wrap(sp, tp, sr, tr, cost):
        if (sr,tr) not in cache: # Add new entry to cache dictionary if cost for sr and tr not previously calclated
            result = func(sp, tp, sr, tr, cost)
            cache[(sr,tr)] = (result[0][len(sp):], result[1][len(tp):], result[4] - cost)
        # Use previously calculated 'sub-values' to get total cost
        return sp + cache[(sr,tr)][0], tp + cache[(sr,tr)][1], '', '', cost + cache[(sr,tr)][2]
    return wrap


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

    al = levenshtein(lemma, form, substcost = 1.1) # Slight preference for insertion/deletion over substitution
    alemma, aform = al[0], al[1] # Only unpack aligned strings, not cost
    # leading spaces
    lspace = max(len(alemma) - len(alemma.lstrip('_')), len(aform) - len(aform.lstrip('_')))
    # trailing spaces
    tspace = max(len(alemma[::-1]) - len(alemma[::-1].lstrip('_')), len(aform[::-1]) - len(aform[::-1].lstrip('_')))
    return alemma[0:lspace], alemma[lspace:len(alemma)-tspace], alemma[len(alemma)-tspace:], aform[0:lspace], aform[lspace:len(alemma)-tspace], aform[len(alemma)-tspace:]


def prefix_suffix_rules_get(lemma, form):
    """
    Arguments:
        lemma -- string representing a 'root'/'basic' form of a word
        form  -- string representing derived form of lemma with added morphological features

    Returns:
        prules -- set of 2-tuples containing inputs and outputs for possible prefix rules used to
            turn lemma into form
        srules -- set of 2-tuples containing inputs and outputs for possible suffix rules used to
            turn lemma into form

    Aligns lemma and form using alignprs, and identifies all possible prefixing/suffixing rules that
    convert lemma to form. Assumes suffixing to be a slightly more complex/elaborate process than
    prefixing (can be applied to reverse strings for languages which prefer prefixing).
    """
    lp,lr,ls,fp,fr,fs = alignprs(lemma, form) # Get six parts, three for in three for out

    # Suffix rules
    ins  = lr + ls + ">"
    outs = fr + fs + ">"

    srules = set()
    for i in range(min(len(ins), len(outs))):
        clean = ins[i:].replace('_', '')
        if len(clean) == 5 and clean[-3:] == 'er>' and clean[-5] in 'eé' and clean[-4] not in 'aâeéèêiïîoôuû':
            srules.add((clean[:-4]+'C'+clean[-3:], replace_first_inst(outs[i:],clean[-4],'C')))
        srules.add((ins[i:], outs[i:])) # Each possible substring pair is added as a suffixing rule
    srules = {(x[0].replace('_',''), x[1].replace('_','')) for x in srules} # Remove _ alignment characters

    # Prefix rules
    prules = set()
    if len(lp) >= 0 or len(fp) >= 0: # Assumes simple prefixing, does not attempt to loop over prefix substrings
        inp = "<" + lp
        outp = "<" + fp
        for i in range(0,len(fr)):
            prules.add((inp + fr[:i],outp + fr[:i]))
            prules = {(x[0].replace('_',''), x[1].replace('_','')) for x in prules}

    return prules, srules

def apply_best_rule(lemma, msd, allprules, allsrules, debug=False, no_pref=False, subset=None, wordmap=None):
    """
    Arguments:
        lemma -- 'root' or 'base' form of word to transform
        msd   -- unimorph string representing desired features in derived form of lemma
        allprules -- dictionary mapping all msds to possible prefixing rules found in the language
        allsrules -- dictionary mapping all msds to possible suffixing rules found in the language
        debug -- boolean whether print statements should be run

    Applies the longest-matching suffix-changing rule given an input
    form and the MSD. Length ties in suffix rules are broken by frequency.
    For prefix-changing rules, only the most frequent rule is chosen.
    """

    if not subset is None and not is_in_subset(lemma, subset, wordmap):
        return

    if debug: print("Lemma: %s\nFeatures: %s" % (lemma, msd))

    base = "<" + lemma + ">"

    if msd not in allprules and msd not in allsrules:
        return lemma # Haven't seen this inflection, so bail out

    if msd in allsrules:
        # One applicable rule is a 3-tuple containing the input, output, and frequency
        applicablerules = [(x[0],x[1],y) for x,y in allsrules[msd].items() if x[0] in base]
        if base[-4] not in 'aâeéèêiïîoôuû':
            for x,y in allsrules[msd].items():
                if len(x[0]) == 5 and x[0] == base[-5] + 'C' + base[-3:]:
                    applicablerules.append((x[0].replace('C',base[-4]), x[1].replace('C',base[-4]), y))

        if applicablerules: # If there are applicable rules, find the best one
            bestrule = max(applicablerules, key = lambda x: (len(x[0]), x[2], len(x[1])))
            if debug: print("Applicable suffix rules:\n%s\nUsing: %s" % (applicablerules, bestrule))
            base = base.replace(bestrule[0], bestrule[1]) # Apply best rule to base form

    # Use above method to apply relevant prefixing rule to base form
    if msd in allprules:
        applicablerules = [(x[0],x[1],y) for x,y in allprules[msd].items() if x[0] in base]
        if applicablerules:
            bestrule = max(applicablerules, key = lambda x: (x[2])) 
            if debug and not no_pref: print("Applicable prefix rules:\n%s\nUsing: %s" % (applicablerules, bestrule))
            base = base.replace(bestrule[0], bestrule[1])
    
    if debug: print()

    # Remove extra characters
    base = base.replace('<', '')
    base = base.replace('>', '')
    return base


def numleadingsyms(s, symbol):
    """Get number of leading symbol characters in s"""
    return len(s) - len(s.lstrip(symbol))


def numtrailingsyms(s, symbol):
    """Get number of trailing symbol characters in s"""
    return len(s) - len(s.rstrip(symbol))

def is_in_subset(lemma, subset, wordmap):
    subset_opts_map = {
        'mod': ['-', 'dtd', 'obs'],
        'no_old': ['-', 'dtd', 'obs', 'M'],
        'no_dtd': ['-', 'obs'],
        'no_obs': ['-', 'dtd'],
        'mod_strict': ['-']
    }

    if subset is None:
        return True
    
    if '.' in subset:
        fix, substr = subset.split('.', maxsplit=1)
        if fix == 'pre':
            return lemma[:len(substr)] == substr
        if fix == 'suf':
            return lemma[-len(substr):] == substr
        warn('{} is not a valid affix type. Running script on unfiltered data set.'.format(fix))
        return True
    
    if subset in subset_opts_map.keys():
        return wordmap[lemma] in subset_opts_map[subset]
    
    warn('Invalid argument passed to --subset option. Running script on unfiltered data set.')
    return True

def load_word_map():
    files = [ 'marked-fra.{}.csv'.format(x) for x in ('dev', 'trn', 'tst') ]

    wordmap = {}
    for f in files:
        try:
            with open(os.path.join('data', f)) as mapfile:
                for line in mapfile:
                    if line.strip() != '':
                        word, fam = line.strip().split(',')
                        wordmap[word] = fam
        except FileNotFoundError:
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

###############################################################################


def main(argv):
    options, remainder = getopt.gnu_getopt(argv[1:], 'odSX:Y:thp:', ['output','debug','suffix-only','x-set=','y-set=','test','help','path='])
    DEBUG, NO_PREF, XSET, YSET, TEST, OUTPUT, HELP, path = False,False,None,None,False,False, False, './Data/'
    for opt, arg in options:
        if opt in ('-o', '--output'):
            OUTPUT = True
        if opt in ('-d', '--debug'):
            DEBUG = True
        if opt in ('-S', '--suffix-only'):
            NO_PREF = True
        if opt in ('-X', '--x-set'):
            XSET = arg
        if opt in ('-Y', '--y-set'):
            YSET = arg
        if opt in ('-t', '--test'):
            TEST = True
        if opt in ('-h', '--help'):
            HELP = True
        if opt in ('-p', '--path'):
            path = arg

    if HELP:
            print("\n*** Baseline for the SIGMORPHON 2020 shared task ***\n")
            print("By default, the program runs all languages only evaluating accuracy.")
            print("To create output files, use -o")
            print("The training and dev-data are assumed to live in ./part1/development_languages/")
            print("Options:")
            print(" -o         create output files with guesses (and don't just evaluate)")
            print(" -t         evaluate on test instead of dev")
            print(" -d         evaluate on debug (or subset when specified) and print debug statements instead of dev")
            print(" -S         when -d flag is set, only print information for suffixing rules")
            print(" -X [set]   trains on subset of data. Must specify one of the following: mod, no_old,")
            print("            no_obs, no_dtd, pre.<prefix string>, suf.<suffix string>")
            print(" -Y [set]   tests on subset of data. Must specify one of the following: mod, no_old,")
            print("            no_obs, no_dtd, pre.<prefix string>, suf.<suffix string>")
            print(" -p [path]  data files path. Default is ../data/")
            quit()

    totalavg, numlang = 0.0, 0
    wordmap = load_word_map()

    for lang in [os.path.splitext(d)[0] for d in os.listdir(path) if '.trn' in d]:
        allprules, allsrules = {}, {}
        if not os.path.isfile(path + lang +  ".trn"):
            continue
        lines = [line.strip() for line in open(path + lang + ".trn", "r", encoding='utf8') if line != '\n']

        # First, test if language is predominantly suffixing or prefixing
        # If prefixing, work with reversed strings
        prefbias, suffbias = 0,0
        for l in lines:
            lemma, _, form = l.split(u'\t')
            aligned = halign(lemma, form)
            if ' ' not in aligned[0] and ' ' not in aligned[1] and '-' not in aligned[0] and '-' not in aligned[1]:
                prefbias += numleadingsyms(aligned[0],'_') + numleadingsyms(aligned[1],'_')
                suffbias += numtrailingsyms(aligned[0],'_') + numtrailingsyms(aligned[1],'_')

        for l in lines: # Read in lines and extract transformation rules from pairs
            lemma, msd, form = l.split(u'\t')

            if XSET is None or is_in_subset(lemma, XSET, wordmap):
                if prefbias > suffbias:
                    lemma = lemma[::-1]
                    form = form[::-1]
                prules, srules = prefix_suffix_rules_get(lemma, form)

                if msd not in allprules and len(prules) > 0:
                    allprules[msd] = {}
                if msd not in allsrules and len(srules) > 0:
                    allsrules[msd] = {}

                for r in prules:
                    if (r[0],r[1]) in allprules[msd]:
                        allprules[msd][(r[0],r[1])] = allprules[msd][(r[0],r[1])] + 1
                    else:
                        allprules[msd][(r[0],r[1])] = 1

                for r in srules:
                    if (r[0],r[1]) in allsrules[msd]:
                        allsrules[msd][(r[0],r[1])] = allsrules[msd][(r[0],r[1])] + 1
                    else:
                        allsrules[msd][(r[0],r[1])] = 1

        if not YSET is None:
            outname = lang + '-' + YSET.replace('.', '-').replace('_', '-')
        else:
            outname = lang

        # Run eval on dev
        
        devlines = [line.strip() for line in open(path + lang + ".dev", "r", encoding='utf8') if line != '\n']
        if TEST:
            devlines = [line.strip() for line in open(path + lang + ".tst", "r", encoding='utf8') if line != '\n']
        if DEBUG and ('pre' not in YSET or 'suf' not in YSET):
            devlines = [line.strip() for line in open(path + lang + ".dbg", "r", encoding='utf8') if line != '\n']
        numcorrect = 0
        numguesses = 0
        if OUTPUT:
            if TEST:
                outfile = open(path + outname + ".tst.out", 'w', encoding='utf8')
            else:
                outfile = open(path + outname + ".dev.out", "w", encoding='utf8')
        
        for l in devlines:
            lemma, msd, correct = l.split(u'\t')
#                    lemma, msd, = l.split(u'\t')
            if prefbias > suffbias:
                lemma = lemma[::-1]
            if not YSET is None:
                outform = apply_best_rule(lemma, msd, allprules, allsrules,
                                          debug=DEBUG,
                                          no_pref=NO_PREF,
                                          subset=YSET,
                                          wordmap=wordmap)
            else:
                outform = apply_best_rule(lemma, msd, allprules, allsrules, debug=DEBUG)

            if not outform is None:
                if prefbias > suffbias:
                    outform = outform[::-1]
                    lemma = lemma[::-1]
                if outform == correct:
                    numcorrect += 1
                numguesses += 1
                if OUTPUT:
                    outfile.write(lemma + "\t" + msd + "\t" + outform + "\n")

        if OUTPUT:
            outfile.close()

        numlang += 1
        totalavg += numcorrect/float(numguesses)

        print(lang + ": " + str(str(numcorrect/float(numguesses)))[0:7])

    print("Average accuracy", totalavg/float(numlang))


if __name__ == "__main__":
    main(sys.argv)