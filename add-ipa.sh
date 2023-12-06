#!/bin/bash

baseurl="https://en.wiktionary.org/wiki/"
lastword=""
infile="data/fra.$1"
ipafile="data/fra.ipa.$1"

rm "$ipafile"
touch "$ipafile"

while read -r line ; do
    lemma=$(echo "$line" | cut -f 1)
    if [ "$lemma" != "$lastword" ] ; then
        curl "$baseurl$lemma" > wikt.html
    fi
    lemma_ipa=$(grep ">IPA<.*French pronunciation.*<span class=\"IPA\">[^>]*" wikt.html | sed -r "s/.*<span class=\"IPA\">([^<]*).*/\1/g")

    msd=$(echo "$line" | cut -f 2)
    form=$(echo "$line" | cut -f 3)
    form_ipa=$(grep ">$form<.*<span class=\"IPA\">[^<]*" wikt.html | sed -r "s/.*>$form<.*<span class=\"IPA\">([^<]*)<.*/\1/g")
    if [ "$lemma" == "$form" ] ; then
        form_ipa="$lemma_ipa"
    fi
    if [[ -z $lemma_ipa ]] ; then
        lemma_ipa="NONE"
    fi
    if [[ -z $form_ipa ]] ; then
        form_ipa="NONE"
    fi
    form_ipa=$(echo $form_ipa | cut -w -f 1)
    echo -e "$lemma\t$lemma_ipa\t$msd\t$form\t$form_ipa" >> "$ipafile"
    lastword=$lemma
done < "$infile"

rm wikt.html