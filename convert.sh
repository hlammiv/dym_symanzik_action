#!/bin/bash

./kentucky2nersc << EOF >> out_convert.log
dims 24 24 24 24
input $1
output $1.nersc
EOF
