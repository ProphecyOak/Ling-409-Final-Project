import os

path = "./Data/"

CONSONANTS = ['m','n','ɲ','ŋ','p','t','k','b','d','ɡ','f','s','ʃ','v','z','ʒ','ʁ','l','j','ɥ','w']
VOWELS = ['i','y','u','e','ø','ə','o','ɛ','ɛː','œ','ɔ','a','ɑ','ɛ̃','œ̃','ɔ̃','ɑ̃']

def getLines(lang):
    return [line.strip() for line in open(path + lang + ".ipa", "r", encoding='utf8') if line != '\n']

def syllabifyWord(line, debug=False):
    lineSplit = line.strip().split("\t")
    orthography, ipa = lineSplit[0], lineSplit[1].split(" ")
    spot = len(ipa) - 1
    out = []
    while spot >= 0:
        #CODA
        while ipa[spot] in CONSONANTS and spot > 0:
            out.insert(0,ipa[spot])
            spot -= 1
        #NUCLEUS
        out.insert(0,ipa[spot])
        spot -= 1
        #GLIDE
        if ipa[spot] in 'jɥw':
            out.insert(0,ipa[spot])
            spot -= 1
        #ONSET
        if spot >= 0:
            out.insert(0,ipa[spot])
            spot -= 1
            if ipa[spot+1] in 'lʁ' and ipa[spot] not in 'lʁ' and ipa[spot] in CONSONANTS:
                out.insert(0,ipa[spot])
                spot -= 1
        out.insert(0,'.')
    out.remove('.')
    outString = ''.join(out)
    if debug: print("%-15s  %-15s %-15s" % (orthography,''.join(ipa), outString))
    return outString
