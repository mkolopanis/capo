This README is about the PAPER pspec pipeline, 2017 version.

The pipeline is meant to be a method for interacting with redundant baseline
type data which has been lstbinned over many nights.

Modules:
capo.oqe
capo.pspec

Scripts:
pspec_2d_to_1d.py
plot_pspec_final.py
pspec_average_seps.py
pspec_batch.sh
plot_trms.py
plot_sigloss.py
frfilter_numbers.py
frf_filter.py
simple_pspec.py
pspec_oqe_2d.py
pspec_final_sigloss_dist.py
pspec_final_sigloss.py
pspec_final_confidence.py

Pipeline components:
pspec_oqe_2d.py   :   computes Pk vs lst.
    many bootstraps over baseline combinations, injects simulated eor signal at given level
    generates many power spectrum channels:
      injected eor, eor+data, simulated noise, simulated noise+eor
      weighted (C) and unweighted (I) are generated for each
pspec_2d_to_1d.py    :  collapses bootstraps and lsts into a single power spectrum with error bars
pspec_average_seps.py : averages multiple kinds of baselines within a single uv "cell". run after 2d_to_1d
pspec_batch.sh    :  runs pspec_oqe_2d and pspec_2d_to_1d across a range of injected eor levels
pspec_final_XXX.py :  uses power spectra with a range of injected power levels to estimate a loss-corrected power spectrum. There are multiple renditions of this step.
frf_filter.py : fringe rate filter data
plot_pspec_final.py : plot the outputs of pspec_final_XXX.py
plot_trms.py : difference even and odd datasets to get an noise spectrum and compare with model
plot_sigloss.py : plot the power spectrum vs injected eor power to visualize signal loss-corrected
