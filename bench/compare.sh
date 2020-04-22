#!/bin/bash

cd "$(dirname "$(realpath "$(which "$0")")")"

time_log="result/compare_$(date +"%s")"
rss_log="${time_log}_rss.tsv"
time_log="${time_log}_time.tsv"

mkdir bin data result &> /dev/null
bin_list=( $(find bin/ -maxdepth 1 -type f -executable | sort | sed 's:^bin/::g') )
data_list=( $(find data/ -maxdepth 1 -type f -exec du -k {} + | sort -n | cut -f2 | sed 's:^data/::g') )
echo ${bin_list[@]} | sed -E 's/^|\s+/\t/g' | tee "$time_log" "$rss_log" > /dev/null
for data in ${data_list[@]}; do
  echo -ne "${data}\t" >> "$time_log"
  echo -ne "${data}\t" >> "$rss_log"
  for bin in ${bin_list[@]}; do
    echo -n "benchmark ${bin} ${data} ... "
    "bin/$bin" "data/$data" 2>&1 | sed -nE '/time:/ { s/^\s*time:\s*([0-9.]+).*$/\1/; p }; /rss:/ { s/^\s*rss:\s*([0-9.]+).*$/\1/; p };' | sed -z 's/\n/ /g' | {
      read time rss;
      echo "${time}s ${rss}MiB"
      echo -ne "${time}\t" >> "$time_log"
      echo -ne "${rss}\t" >> "$rss_log"
    }
  done
  echo >> "$time_log"
  echo >> "$rss_log"
done
