#!/bin/csh
sprep96 -M $1 -R -NMOD 10 -PARR PARR.txt 
sdisp96 -v
sregn96 -HS 10 -HR 0 -DE -NOQ
sdpder96 -R -TXT -K 2 -XLEN 3.0 -X0 1.0 -YLEN 4.0
