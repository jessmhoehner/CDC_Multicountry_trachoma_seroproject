# -*- coding: utf-8 -*- 
#
#Author: JR
#Maintainer: JR
#License:  2019, EICC, GPL v2 or later
# 
#NCV19/clean/Makefile
#----------------------------------------------------
#

HERE := $(shell Rscript --vanilla -e 'write(here::here(), stdout())')

clean_data = \
				$(HERE)/clean/output/WHO_sit_rep_clean.csv \
				$(HERE)/model/input/nga_cv19_clean.csv \
				$(HERE)/write/input/nga_cv19_clean.csv \
				$(HERE)/model/input/ncdc_clean.csv \
				$(HERE)/write/input/ncdc_clean.csv \
				$(HERE)/clean/output/state_actions_clean.csv \
				$(HERE)/clean/output/state_size_clean.csv \
				$(HERE)/model/input/nigeria_2020_agedemo_clean.csv \
				$(HERE)/model/input/ncdc_discharges_clean.csv

.PHONY:	all clean

all: $(clean_data)

clean:
		rm $(HERE)/clean/output/*
		rm $(HERE)/write/input/*
		rm $(HERE)/model/input/*
		
		
$(clean_data): \
				src/clean.R \
				$(HERE)/clean/input/WHO_sit_rep.csv \
				$(HERE)/clean/input/NGA_COVID_19_JR.csv \
				$(HERE)/clean/input/ncdc_rep.csv \
				$(HERE)/clean/input/state_actions.csv \
				$(HERE)/clean/input/state_size.csv \
				$(HERE)/clean/input/Nigeria_2020_agedemo.csv \
				$(HERE)/clean/input/ncdc_deaths_discharges.csv
			Rscript --vanilla $<
	
# done.