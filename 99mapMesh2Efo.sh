#!/usr/bin/env bash
/GWD/bioinfo/apps/bin/perl /GWD/bioinfo/projects/cb-software/ontology/bin/NCBO_get_ontology_cross_mappings.pl -a=MESH -b=EFO > ../data/mesh2efo.mappings.txt
