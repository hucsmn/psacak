#!/bin/bash

cd "$(dirname "$(realpath "$(which "$0")")")"

mkdir data &> /dev/null

echo "Pizza&Chili Corpus downloader"
echo "============================="
for url in {texts/code/sources.gz,texts/music/pitches.gz,texts/protein/proteins.gz,texts/dna/dna.gz,texts/nlang/english.gz,texts/xml/dblp.xml.gz}; do
  corpus="$(basename "$url" | sed 's/.gz$//')"
  if [[ -f "data/$corpus" ]]; then
    echo "already downloaded '$corpus'."
  else
    echo "downloading '$corpus'..."
    tmp="$(mktemp)"
    curl -o - "http://pizzachili.dcc.uchile.cl/$url" | gunzip > "$tmp"
    mv "$tmp" "data/$corpus"
  fi
done
