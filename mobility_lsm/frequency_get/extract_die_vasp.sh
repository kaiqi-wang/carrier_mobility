#!/bin/bash

file=vasprun.xml

echo "extracting dielectric function from " $file
awk '
/<dielectricfunction/ { on=1 }
on==1 && /<r>/ { print $2,($3+$4+$5)/3 } 
#on==1 && /<r>/ { print $2,$3,$4,$5 } 
on==1 && /<real>/ { print " " }
on==1 && /<imag>/ { print " " }
/<\/dielectricfunction/ { on=0 }
' < $file > diel.dat
