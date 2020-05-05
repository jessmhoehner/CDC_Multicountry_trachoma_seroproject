# -*- coding: utf-8 -*- 
#
#Author: JR
#Maintainer: JR
#License: 2020, EICC, GPL v2 or later
# 
# multicountrytrachseroproject/clean/Makefile
#----------------------------------------------------
#
HERE := $(shell Rscript --vanilla -e 'write(here::here(), stdout())')

clean_data = \
				$(HERE)/model/input/DRC1Ct694_clean.csv \ 
				$(HERE)/plot/input/DRC1Ct694_agespbin_df.csv
				
.PHONY:	all clean

all: $(clean_data)

clean:
	-rm -r $(HERE)/clean/output/*
	-rm -r $(HERE)/plot/input/*
	-rm -r $(HERE)/model/input/*
		
$(clean_data): \
				src/clean.R \
				$(HERE)/clean/input/DRC1Ct694.csv 
			Rscript --vanilla $<
			
# done.