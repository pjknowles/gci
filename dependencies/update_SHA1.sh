#!/usr/bin/env bash

cwd=$(pwd)
for f in ./*/
do
    cd "${cwd}/${f}"
    fname="${cwd}/${f%%/}_SHA1"
    $(git rev-parse HEAD > $fname)
done
