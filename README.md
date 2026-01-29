# CPS1000

---

*Note:*
In the R Markdown document, the following figures are still missing: FigS1, FigS6, FigS7, FigS9, FigS10 + tables (reorder)

These figures not yet working properly: Fig4 (correlation plots), FigS5 (network plot)

---

This repository contains data and scripts to reproduce the figures and analyses presented in the paper:

---

**Drug response of leukemia and lymphoma**

Kevin D. Hofer, Sebastian Scheinost, Lena Ben Taarit, Jennifer Hüllein, Tatjana Walther, Kerstin Putzker, Leopold Sellner, Małgorzata Oleś, Sascha Dietrich, Jérôme Huber, Michael W.M. Kühn, Thomas Kindler, Olivier Bernard, Florence Nguyen-Khac, Marta Crespo Maull, Francesc Bosch, Alexandre Theocharides, Markus G. Manz, Beat Bornhauser, Jean-Pierre Bourquin, Joe Lewis, Wolfgang Huber, Junyan Lu, Thorsten Zenz

---

All data of the CPS1000 cohort is stored in /data/CPS1000_data.RData. The analysis scripts can be found under CPS1000_Analysis.Rmd with the associated .html file.

For the validation cohort (EMBL2016) data are stored in /validation/EMBL2016_data.RData.

Please refer to the methods of the manuscript for details on the experimental setup and quality-control measures.

If you use data from this work in published research, please cite the paper. 

*Installation guide*

To run the entire analysis, clone the repository and run the script CPS1000_Analysis.Rmd.

*System requirements*

To run the analysis, R (at least version 4.1.0) and all dependency libraries are required.

*Output*

The CPS1000 analysis script takes roughly ... minutes to run. The expected output is shown in /figures and /tables.

