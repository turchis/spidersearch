#!/bin/bash
[[ ! -z "$(ls ?filter/*.fits 2>/dev/null)" ]] && mv ?filter/*.fits ./ && rm -rf ?filter/
mkdir -p {i,g,r}filter
for prefix in g i r ; do while read i ; do [[ -f "$i" ]] &&  mv "$i" "${prefix}filter/" ; done < "${prefix}filter.txt" ; done
for prefix in g i r ; do  mv "median_${prefix}.fits" "${prefix}filter/" ; done
