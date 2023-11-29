import os

diffLineCounter = 0
problemWords = {}
problemWordErrors = {}

with open('./Data/fra.dev','r') as inFile, open('./Data/fra.out','r') as outFile:
    for inLine, outLine in zip(inFile, outFile):
        if inLine != outLine:
            inLineSplit, outLineSplit = inLine.strip().split("\t"), outLine.strip().split("\t")
            curWord = inLineSplit[0]
            diffLineCounter += 1
            if curWord not in problemWords.keys():
                problemWords[curWord] = 1
            else:
                problemWords[curWord] += 1

problemWordsList = sorted(problemWords.items(), key=lambda item: item[1], reverse=True)
with open("./WordsWithIssues.txt",'w') as WordsWithIssues:
    WordsWithIssues.write("Count\tWord\n")
    for x in problemWordsList:
        WordsWithIssues.write("%s\t%s\n" % (x[1],x[0]))
