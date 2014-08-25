

radius=10
iterations=10


export PYTHONPATH=$WASQ_ROOT:$PYTHONPATH

names=(unfolded folded)

for name in ${names[@]}; do
    t=$name.tpr
    test -f $t && continue
    python $WASQ_ROOT/wasq/PrepareTpr.py -T 10 -O 1 -p $name.pdb -o $t
done

tprs=""
for name in ${names[@]}; do
    tprs="$tprs $name.tpr"
done

if test -e AS; then
    echo "Adaptive Sampling has already run."
    echo "To continue, delete the 'AS' directory and rerun this script."
    exit 1
else
    python $WASQ_ROOT/wasq/AdaptiveSampling.py -r $radius -i $iterations $tprs
fi
