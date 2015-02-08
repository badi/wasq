#!/usr/bin/env bash






temperatures=(200 273 300 330 360 400 500 600 1000)
md_outfreq=2

for temp in ${temperatures[@]}; do
    echo $temp
    dir=t_${temp}
    source $dir/params.sh
    set -x
    ../CompareToPlainMD.py \
	-r $radius -c $dir/AS/cells -x $dir/MD/traj.xtc -t ../topol.pdb -i $dir/AS/iteration \
	-w $time \
	-o $dir/AS_vs_MD.png \
	-X $dir/AS_vs_MD.trax
    set +x

done
