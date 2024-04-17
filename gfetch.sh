#!/usr/bin/env bash
#
#       gfetch.sh
#
#       Generates correlator input file for bgv* Backus-Gilbert
#       inversion script.
#

gen2l_location="../../gen-2l"      #Location of the gen-2l dataset

outfile="output.dat"    #File to output correlator information
tmpfile="tmp_$RANDOM"   #Temporary arithmetic file (will be destroyed on exit)

> "$outfile"
> "$tmpfile"
> "seed"

channel=$1 #Meson channel (see Don's nomenclature - s01 = Upsilon, p10 = Chi_b1)

# Check for cmd line args
if [ "$#" -eq 4 ]
then
    Nt=$2       #Number of Euclidean time points (N_tau)
    Ns=$3       #Number of sampling slices (=Ns+1 points in omega space)
    Nb=$4       #Size of bootstrap sample
else
    echo "Defaulting to presets..."
    Nt=128      
    Ns=100      
    Nb=-1       #-1 designates full set.
fi

gentype="covgen"              #covgen version ("covgen" or "covgend")

echo "Using $gentype."
echo "Channel: $channel"

if [ "${channel:0:1}" = "s" ]
then
    suffix="s"
    channels=("spp_0" "spp_i" "sxx_0" "sxx_i")
else
    suffix="p"
    channels=("pp_A" "pp_0" "pp_i" "pp_ij" "xx_A" "xx_0" "xx_i" "xx_ij")
fi

datapath="$gen2l_location/32x$Nt/${suffix}onia*"

echo "Collecting data..."

# Calculate channel number ch_no

for ch in "${!channels[@]}"; do
    if [[ "${channels[$ch]}" = "$channel" ]]; then
        ch_no=$ch
    fi
done

echo "Channel number:" $ch_no

confcount=0
for confpath in $datapath
do
    echo -e "$Nt\n$ch_no\n$confpath\n" | ./chread >> "$tmpfile"
    confcount=$((confcount+1))
done

if [ "$Nb" -eq -1 ]
then
    echo "Defaulting to maximum sampling..."
    Nb=$confcount
fi

echo -e "$Nt\n$Ns" >> "$outfile"
echo "Done."

if [ "$gentype" = "covgen" ]
then
    echo -e "$Nt\n$confcount\n$(cat "$tmpfile")" > "$tmpfile"
    cat "$tmpfile" | ./covgen >> "$outfile"
elif [ "$gentype" = "covgend" ]
then
    echo -e "$Nt\n$confcount\n$Nb\n$(cat "$tmpfile")" > "$tmpfile"
    cat "$tmpfile" | ./covgend >> "$outfile"
fi

rm $tmpfile     #Delete arithmetic file
