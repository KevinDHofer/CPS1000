# CPS1000

This repository contains data and scripts to reproduce the figures and analyses presented in the paper:

**Functional profiling of leukemia and lymphoma by ex vivo drug screening**

Kevin D. Hofer, Sebastian Scheinost, Lena Ben Taarit, Jennifer Hüllein, Tatjana Walther, Kerstin Putzker, Leopold Sellner, Małgorzata Oleś, Sascha Dietrich, Jérôme Huber, Michael W.M. Kühn, Thomas Kindler, Olivier Bernard, Florence Nguyen-Khac, Marta Crespo Maull, Francesc Bosch, Alexandre Theocharides, Markus G. Manz, Beat Bornhauser, Jean-Pierre Bourquin, Joe Lewis, Wolfgang Huber, Junyan Lu, Thorsten Zenz

---

All data of the CPS1000 cohort and linked omcis matrices are stored in /data. The analysis scripts can be found under CPS1000_analysis.Rmd with the associated .html file. 

The respective validation data are stored in /validation. To rerun the analysis, add the respective data from the CLL-map Portal (cllmap.org): download the TSV file and the 41588_2022_1140_MOESM3_ESM.xlsx file to the validation folder.

Please refer to the methods of the manuscript for details on the experimental setup.

If you use data from this work in published research, please cite the paper. 

*Installation guide*

To run the entire analysis, clone the repository, download the files as mentioned above and run the script CPS1000_analysis.Rmd.
