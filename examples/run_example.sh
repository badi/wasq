

params=params.sh

if ! test -f $params; then
    echo "'$params' needs to exist with the following definitions:"
    echo "  - names: and array of names of the pdb files without the extension"
    echo "           e.x.: names=(unfolded folded)"
    echo "  - time: the time (in ps) of each task. e.x.: time=10"
    echo "  - outfreq: the output frequency (in ps) of the simulation. e.x.: outfreq=1"
    echo "  - radius: the cell radius. e.x.: radius=10"
    echo "  - iterations: the number of AS iterations to run. e.x.: iterations=10"
    exit 1
fi

set -o errexit
set -x

source $params

for name in ${names[@]}; do
    t=$name.tpr
    test -f $t && continue
    python $WASQ_ROOT/wasq/PrepareTpr.py -T $time -O $outfreq -p $name.pdb -o $t -d info2
done

tprs=""
for name in ${names[@]}; do
    tprs="$tprs $name.tpr"
done

if test -e AS; then
    set +x
    echo "Adaptive Sampling has already run."
    echo "To continue, delete the 'AS' directory and rerun this script."
    exit 1
else
    $(which time) -p python $WASQ_ROOT/wasq/AdaptiveSampling.py -r $radius -i $iterations $tprs
fi
