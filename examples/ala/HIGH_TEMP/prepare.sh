#!/usr/bin/env bash

source params.sh

gps_templ="\
length: 1
type: 4

TEMP.000000
"

params="\
names=(unfolded)
time=20
outfreq=1
radius=10
iterations=50
"

for t in ${temps[@]}; do
    outdir=t_${t}
    test -d $outdir || mkdir -p $outdir
    tpr=$outdir/$(basename $start)
    gps=$outdir/temp.gps
    test -f $outdir/SUCCESS && continue
    cp $start $tpr
    printf "%s" "$gps_templ" >$gps
    sed -i "s/TEMP/$t/" $gps
    echo "Setting up $t K"
    guamps_set -f $tpr -s ref_t -i $gps -O &>/dev/null
    test $? -eq 0 || exit 1
    rm $gps
    printf "%s" "$params" >$outdir/params.sh
done
