Code for the manuscript:
<br>
<br>
"Comparison of platforms for testing antibodies to Chlamydia trachomatis 
antigens: data from the Democratic Republic of the Congo and Togo"

Gwyn, S., Awoussi, M.S., Bakhtiari, A. et al. Comparison of platforms for testing antibodies to Chlamydia trachomatis antigens in the Democratic Republic of the Congo and Togo. Sci Rep 11, 7225 (2021). https://doi.org/10.1038/s41598-021-86639-8

Project structure
<br>
  1. clean
      - input (raw spreadsheets of age, titre, and sero_pos data)
      - src (script to clean the data in the input folder)
  2. observed
      - input (files for each sheet from the clean task)
      - src (script to calculate age seroprevalence)
  3. adjust
      - input (files for each sheet from the clean task)
      - output (age adjusted data)
      - src (script to adjust data by age weights)
  4. model
      - input (files for each sheet from the clean task)
      - src (script to run model)
  5. plot
      - input (cleaned and prepared files from the observed and model tasks)
      - output ( all original figures generated from analysis)
      - src (script to make figures)
 6. write
	- contains files from clean task, rest come from plot task
	- output (contains report and image files)
	- src (contains R markdown file)
