Code for the manuscript:
<br>
<br>
"Comparison of platforms for testing antibodies to Chlamydia trachomatis 
antigens: data from the Democratic Republic of the Congo and Togo"

Project structure
<br>
  1. clean
      - input (raw spreadsheets of age, titre, and sero_pos data)
      - src (script to clean the data in the input folder)
  2. observed
      - input (files for each sheet from the clean task)
      - src (script to calculate age seroprevalence)
  2. model
      - input (files for each sheet from the clean task)
      - src (script to run model)
  3. plot
      - input (cleaned and prepared files from the observed and model tasks)
      - output ( all original figures generated from analysis)
      - src (script to make figures)
4. adjust
      - input (files for each sheet from the clean task)
      - output (age adjusted data)
      - src (script to adjust data by age weights)
5. write
	- contains files from clean task, rest come from plot task
	- output (contains report and image files)
	- src (contains R markdown file)
