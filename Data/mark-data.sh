#!/bin/bash

base_url='https://en.wiktionary.org/wiki/'
marked_name="marked-$1.csv"

if [[ $(ls "$marked_name") ]] ; then
    rm "$marked_name"
fi

cut -f 1 $1 | sort | uniq | while read -r line ; do
    curl -s "$base_url$line" > tmp
    grep -e '<span class=\"mw-headline[^>]*>[^<]*<' tmp | sed -r 's;.*<span class=\"mw-headline[^>]*>([^<]*)<.*;\1;g' > tmp2
    if [[ $(grep '^French' tmp2) ]] ; then
        if [[ $(ggrep -P '<span class=\"usage-label-sense((?<!<a).)*<a[^>]*>dated<' tmp) ]] ; then
            echo "$line,dtd" >> "tmp-marked.csv"
        elif [[ $(grep 'Category:French obsolete forms' tmp) ]] ; then
            echo "$line,obs" >> "tmp-marked.csv"
        else
            echo "$line,-" >> "tmp-marked.csv"
        fi
    else
        if [[ $(grep 'Middle French' tmp2) ]] && [[ $(grep 'Old French' tmp2) ]] ; then
            echo "$line,M-O" >> "tmp-marked.csv"
        elif [[ $(grep 'Middle French' tmp2) ]] ; then
            echo "$line,M" >> "tmp-marked.csv"
        else
            echo "$line,O" >> "tmp-marked.csv"
        fi
    fi
done

cat tmp-marked.csv >> "$marked_name"

rm tmp*