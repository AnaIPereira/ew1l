#!/bin/bash

mkdir -p log results

NDIAS1L=$(grep -c "\[ d1l" ../ew_tapir/dia/munuenu-1l.dia)
#NDIAS2L=$(grep -c "\[ d2l" ../tapir/dia/ggH-2l.dia)


#parallel -j 8 "form -M \
#	-d LOOPS=1 -d DIA={} compute-dia.frm > log/d1l{}-dia.log; echo Finished d1l{}" ::: $(seq 1 #$NDIAS1L)

#parallel -j 8 "form -M \
#	-d LOOPS=2 -d DIA={} compute-dia.frm > log/d2l{}-dia.log; echo Finished d2l{}" ::: $(seq 1 $NDIAS2L)

#parallel -j 8 "form -M \
#	-d LOOPS=1 -d DIA={} compute-col.frm > log/d1l{}-col.log; echo Finished d1l{} colour" ::: $(seq 1 $NDIAS1L)

#parallel -j 8 "form -M \
#	-d LOOPS=2 -d DIA={} compute-col.frm > log/d2l{}-col.log; echo Finished d2l{} colour" ::: $(seq 1 $NDIAS2L)


# Without parallel:
for dia in $(seq 1 $NDIAS1L); do
	form -d LOOPS=1 -d DIA=$dia compute-dia.frm > log/d1l$dia-dia.log
	#form -d LOOPS=1 -d DIA=$dia compute-col.frm > log/d1l$dia-col.log
	echo "Finished d1l$dia"
done
#
#for dia in $(seq 1 $NDIAS2L); do
#	form -d LOOPS=2 -d DIA=$dia compute-dia.frm > log/d2l$dia-dia.log
#	form -d LOOPS=2 -d DIA=$dia compute-col.frm > log/d2l$dia-col.log
#	echo "Finished d2l$dia"
#done
