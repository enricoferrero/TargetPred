#!/usr/bin/env bash

## gene_target
curl -H 'Accept:text/tab-separated-values' 'http://us1salx00648.corpnet2.com:9999/api/ds/v1/search/co/sentence/sentencedetail/flat/*/*/*/*?terms=type%3AGENE%20id%3ABIOPROC%24BP71232&inorder=false&slop=1000&limit=100000&from=0&sortby=document_date%3Adesc' > ../data/docstore.gene_target.results.txt
