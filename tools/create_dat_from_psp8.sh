#!/bin/bash

for fname in "$@"
do
    echo "processing $fname"
    name=${fname%.psp8}
    awk '/<INPUT>/{flag=1;next}/<\/INPUT>/{flag=0}flag' $fname > ${name}.dat
    cp ${name}.dat ${name}_upf.dat
    sed -i "s|psp8|upf|g" ${name}_upf.dat
done







