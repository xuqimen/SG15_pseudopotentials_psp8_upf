#!/bin/bash
# syntax: ./convert_psp8_to_upf.sh FILE1 FILE2, ...

EXE="../../run.sh"

if [ $# -lt 1 ]
then
   printf "\n"
   printf "Syntax: ./convert_psp8_to_upf.sh FILE1 FILE2 ...\n\n"
   printf "Requirement: ONCVPSP (http://www.mat-simresearch.com/)\n"
   printf "Converts PSP8 files FILE1, FILE2, ... to UPF format.\n"
   printf "Leaves PSP8 formatted files unmodified.\n\n"
   printf "convert_psp8_to_upf.sh -- Qimen Xu, 07/2021.\n\n"
   exit 1
fi


for fname in "$@"
do
    echo "processing $fname"
    name=${fname%.psp8}

    # create a .dat file by extracting lines between <INPUT> and </INPUT>
    awk '/<INPUT>/{flag=1;next}/<\/INPUT>/{flag=0}flag' $fname > ${name}.dat

    # replace psp8 format to upf in input
    sed -i "s|psp8|upf|g" ${name}.dat
    
    # run oncvpsp with the new .dat for upf
    # since ONCVPSP requires the input file to be in current folder, otherwise there would
    # an "sed -e unknown option to s" error, we copy it to the current folder first
    # ${EXE} ${name} -np
    cp ${name}.dat temp_input_dat.dat
    ${EXE} temp_input_dat -np
    rm temp_input_dat.dat
    mv temp_input_dat.out ${name}.out

    # extract the pseudopotential part between "Begin PSP_UPF" and "END_PSP"
    awk '/Begin PSP_UPF/{flag=1;next}/END_PSP/{flag=0}flag' ${name}.out > ${name}.upf

    rm ${name}.{dat,out}
done







