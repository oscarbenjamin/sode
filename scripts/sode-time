#!/bin/bash
#
# Test time taken by the different implementations of examples

set -o errexit

# Integration variables
t1=0
t2=100
t2l=$(($t2*100))
dtmax=0.001
dtout=0.1

echo
echo Performance comparison between pure python and c implementation of the
echo examples here.
echo

for sys in\
    weiner\
    weiner2\
    lin2d\
    linmult\
    sinmult\
    tanmult
    do
        echo
        echo ----------- Example $sys -----------
        echo
        echo ----- In Pure python $t2 secs
        time sode-pyexamples solve --system-type=$sys --t1=$t1 --t2=$t2\
                                   --dtmax=$dtmax --dtout=$dtout > .tmppy
        echo
        echo ----- In Cython      $t2 secs
        time sode-cyexamples solve --system-type=$sys --t1=$t1 --t2=$t2\
                                   --dtmax=$dtmax --dtout=$dtout > .tmpcy
        echo
        echo ----- In Pure c      $t2 secs
        time sode-cexamples $sys $t1 $t2 $dtmax $dtout > .tmpc
        echo

        echo
        echo ----- In Cython      $t2l secs
        time sode-cyexamples solve --system-type=$sys --t1=$t1 --t2=$t2l\
                                   --dtmax=$dtmax --dtout=$dtout > .tmpcy
        echo
        echo ----- In Pure c      $t2l secs
        time sode-cexamples $sys $t1 $t2l $dtmax $dtout > .tmpc
        echo

done

rm .tmpc .tmppy .tmpcy
