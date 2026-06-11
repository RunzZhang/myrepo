1. Geant4 sims:
    a. macro Cf252_source_iso.mac
    b. Geant 4 file: sbcGeant4.11/RL_branch
        i. in material file, I have included different PE density, CF4 densities on different temperatures
        ii. in construction file, I have added active source volume -> cf_active
        iii. I usually put what information need to saved in the StepAction.cc file, feel free to change it as your need
            Now it only save steps at particular volumes and for particular particles(neutron/Ar isotopes)
2.read root file from simulation: becasue of the ucsb RAM settings, when the output file is too huge, python cannot directly read them.
    To adjust for this, some modification was developed. i.e first chunk files into different smaller files(chunk_root_format) then sort out
    required information(readroot_list_Cf.py) and then analysis them (SN_plotCf_both.py) for a single config simulation result. Because at last, we need to read results from different
    config of simulations, (Wrapup_analysis.py) is to achieve this.

    Some codes are not well structured :( but I am trying to clean it up to make it more sense
    a. chunk_root_format_Cf.py. I made marks A.B.C to indicate where I need to modify
    b. readroot.py. Only class ReadRoot is useful. I also marked A.B.C for easier understanding.
    c. SN_plotCf_both.py Analysis module for single configuration simulation. A.B.C.D comments describe main functions
        to proceed to plot multiple configurations of NR spectrum
    d. Wrapup_Analaysis. This is a code that wrapup both gamma analysis and Cf analysis. A.B state where the Cf simulatin
    plots are.

