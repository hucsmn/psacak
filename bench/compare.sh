#!/bin/bash

cd "$(dirname "$(realpath "$(which "$0")")")"

time_log="result/compare_$(date +"%s")"
rss_log="${time_log}_rss.tsv"
time_log="${time_log}_time.tsv"

echo bin/* | sed 's:\s*bin/:\t:g' | tee "$time_log" "$rss_log" > /dev/null
for data in data/*; do
  echo -ne "${data#data/}\t" >> "$time_log"
  echo -ne "${data#data/}\t" >> "$rss_log"
  for bin in bin/*; do
    echo -n "benchmark ${bin#bin/} ${data#data/} ... "
    "$bin" "$data" 2>&1 | sed -nE '/time:/ { s/^\s*time:\s*([0-9.]+).*$/\1/; p }; /rss:/ { s/^\s*rss:\s*([0-9.]+).*$/\1/; p };' | sed -z 's/\n/ /g' | {
      read time rss;
      echo "${time}s ${rss}MiB"
      echo -ne "${time}\t" >> "$time_log"
      echo -ne "${rss}\t" >> "$rss_log"
    }
  done
  echo >> "$time_log"
  echo >> "$rss_log"
done
