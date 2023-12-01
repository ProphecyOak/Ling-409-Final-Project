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
            echo "$line,dtd" >> "$marked_name"
        elif [[ $(grep 'Category:French obsolete forms' tmp) ]] ; then
            echo "$line,obs" >> "$marked_name"
        else
            echo "$line,-" >> "$marked_name"
        fi
    else
        if [[ $(grep 'Middle French' tmp2) ]] && [[ $(grep 'Old French' tmp2) ]] ; then
            echo "$line,M-O" >> "$marked_name"
        elif [[ $(grep 'Middle French' tmp2) ]] ; then
            echo "$line,M" >> "$marked_name"
        else
            echo "$line,O" >> "$marked_name"
        fi
    fi
done

rm tmp*