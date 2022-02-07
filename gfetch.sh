#!/usr/bin/env bash
#
#       gfetch.sh
#
#       Generates correlator input file for bgv* Backus-Gilbert
#       inversion script.
#

gen2l_location="../gen-2l"      #Location of the gen-2l dataset

outfile="output.dat"    #File to output correlator information
tmpfile="tmp_$RANDOM"   #Temporary arithmetic file (will be destroyed on exit)

> "$outfile"
> "$tmpfile"
> "seed"


# Check for cmd line args
if [ "$#" -eq 3 ]
then
    Nt=$1       #Number of Euclidean time points (N_tau)
    Ns=$2       #Number of sampling slices (=Ns+1 points in omega space)
    Nb=$3       #Size of bootstrap sample
else
    echo "Defaulting to presets..."
    Nt=128      
    Ns=500      
    Nb=-1       #-1 designates full set.
fi

channel="s11"           #Meson channel (see Don's nomenclature - s01 = Upsilon, p10 = Chi_b1)
gentype="covgen"              #covgen version ("covgen" or "covgend")

echo "Using $gentype."
echo "Channel: $channel"

datapath="$gen2l_location/32x$Nt/${channel:0:1}onia*"

echo "Collecting data..."

confcount=0
for confpath in $datapath
do
    echo -e "$Nt\n$channel\n$confpath\n" | ./chread >> "$tmpfile"
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
