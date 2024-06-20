#!/bin/bash


IN_PEAKS=$1
OUT_BED=$2

cut -f 1,2,3,4,5,6,7,8,9 $IN_PEAKS > $OUT_BED

