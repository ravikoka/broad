# Acceptance Studies for the ALICE Experiment
## Ravi Koka
## UT Austin ALICE Group

This repository contains my analyses for the C. Markert group. In short, we aim to determine how detector acceptance affects hadron-$\Lambda$ angular correlations. 

I've generated $10^7$ proton-proton events with `PYTHIA6`, serving as the base of my analysis. I write the relevant data from these events into *analysisMulti.root*. You can download this file, but the raw data is too large to upload. 

The main workhorse here is *analysisMulti.ipynb*. Here you can find my implementation of the `SparseAnalyzer` class. This class couples my `THnSparse` distributions with kinematic cuts and various useful analysis methods.  

