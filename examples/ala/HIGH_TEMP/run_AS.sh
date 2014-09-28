#!/usr/bin/env bash

runner=$(readlink -f ../../run_example.sh)

for tdir in t_*; do
    test -f $tdir/SUCCESS && continue
    pushd $tdir &>/dev/null
    bash $runner
    test $? -eq 0 && touch SUCCESS || exit 1
    popd &>/dev/null
done
