#!/usr/bin/env bash

source params.sh
outputfreq_ps=2

get_timestep() {
    local tpr="$1"
    guamps_get -f $tpr -s deltat
}

set_nsteps() {
    local tpr="$1"
    local time_ps="$2"
    local dt=$(get_timestep $tpr)
    local nsteps=$(python -c "print $time_ps / $dt")
    echo $nsteps | guamps_set -f $tpr -s nsteps -O
}

set_outputfreq() {
    local tpr="$1"
    local freq_ps="$2"
    local dt=$(get_timestep $tpr)
    local steps=$(python -c "print $freq_ps / $dt")

    for s in nstlog nstxtcout nstxout nstvout nstfout; do
	echo $steps | guamps_set -f $tpr -s $s -O
    done

}

get_time_ps() {
    local asdir="$1"
    source $asdir/params.sh
    awk "{SUM += \$1} END {print SUM * $time}" $asdir/AS/iteration/*/nwalkers.txt
}

for t in ${temps[@]}; do
    outdir=t_${t}
    mddir=$outdir/MD
    tpr=$mddir/topol.tpr

    test -f $outdir/SUCCESS || continue
    test -f $tpr            && continue

    test -d $mddir || mkdir -p $mddir
    cp $outdir/unfolded.tpr $tpr

    time_ps=$(get_time_ps $outdir)
    set_nsteps $tpr $time_ps
    set_outputfreq $tpr $outputfreq_ps

done
