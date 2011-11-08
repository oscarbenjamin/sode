#!/bin/bash
#
# Test time taken by the different implementations of examples

# Integration variables
t1=0
t2=10
dtmax=0.0001
dtout=0.01

echo
echo Performance comparison between pure python and c implementation of the
echo examples here.
echo

for sys in\
    lin2d\
    linmult\
    sinmult\
    tanmult\
    weiner\
    weiner2
    do
        echo ----------- Example $sys -----------
        echo ----- In Pure python
        time ./pyexamples.py solve --system-type=$sys --t1=$t1 --t2=$t2\
                                   --dtmax=$dtmax --dtout=$dtout > tmp
        echo
        echo ----- In Pure c
        time ./cexamples.exe $sys $t1 $t2 $dtmax $dtout > tmp
        echo
        echo
done

rm tmp
