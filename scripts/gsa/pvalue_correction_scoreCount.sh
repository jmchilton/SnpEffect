#!/bin/sh

in="$1"
out="$2"

cut -f 1,2 "$in" > "$out"
