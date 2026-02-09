#!/bin/bash

export PATH=$PATH:/home/ana/Documents/Software/src/tapir

mkdir -p drawings dia topo

for loops in 1; do

	tapir \
		--conf tapir.conf \
		--qlist qlist.${loops}l \
		--topologyart drawings/topologies-${loops}l.tex \
		--diagramart drawings/dias-${loops}l.tex \
		--diaout dia/munuenu-${loops}l.dia \
		--topselout topo/topsel.${loops}l

	mv topo/mapping.h topo/mapping-${loops}l.h
	mv topo/topologyList.yaml topo/topologyList-${loops}l.yaml
	mv topo/topologyList.m topo/topologyList-${loops}l.m

	# Draw diagrams and topologies
	cd drawings
	echo Generating ${loops}l drawings...
	lualatex dias-${loops}l.tex > /dev/null
	lualatex topologies-${loops}l.tex > /dev/null
	cd ..

done
