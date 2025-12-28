#!/bin/bash


# Declare array with 4 elements
ARRAY=( 'Debian Linux' 'Redhat Linux' 'Ubuntu Linux' )
ELEMENT=${#ARRAY[@]}
echo $ELEMENT
for ((i=0; i<$ELEMENT; i++)); do
    echo ${ARRAY[${i}]}
done


declare -a pmatrix
# parameter set
# d Nk nd tau tp tpp J mup u0 neta dom 

pmatrix[0]="12 16 0.75 0.057 -0.2 0.0 0.17 -0.78 0.69 7 20"
pmatrix[1]="12 16 0.75 0.057 -0.2 0.0 0.17 -0.78 0.69 7 20"
pmatrix[2]="12 16 0.75 0.057 -0.2 0.0 0.17 -0.78 0.69 7 20"





N=${#pmatrix[*]}

echo $N

for ((i=0; i<$N; i++))
do
    echo $i
    params.out ${pmatrix[$i]} # >> log.txt &
    wait $!

done
#
#
