#!/bin/bash

data_dir=$1

docker run \
    -it \
    --mount type=bind,source="$(pwd)/$data_dir",target=/sql \
    psql13 \
    "./sql/psql_container.sh"
