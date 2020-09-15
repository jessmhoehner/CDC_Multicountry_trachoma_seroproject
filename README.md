Code for the manuscript:
<br>
<br>
"Comparison of multiple tests for determination 
of seroconversion rates to the Chlamydia trachomatis
antigen Pgp3: a multi-country analysis"

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
