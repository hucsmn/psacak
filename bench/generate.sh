#!/bin/bash

cd "$(dirname "$(realpath "$(which "$0")")")"

random_data() {
  local length="$1"
  local alphabet="$2"

  local set0="$(for i in {0..255}; do printf "\\%o" "$i"; done)"
  local set1="$(for i in {0..255}; do printf "\\%o" "$((i % alphabet))"; done)"
  
  dd status=none if=/dev/urandom bs=1M count="$length" | tr "$set0" "$set1"
}

mkdir data &> /dev/null
for length in {16,64,128}; do
  for alphabet in {4,16,128,256}; do
    output="data/random-k${alphabet}-${length}m"
    echo "generate $output (alphabet: ${alphabet}, size: ${length}m)"
    random_data "$length" "$alphabet" > "$output"
  done
done
