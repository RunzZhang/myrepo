1. Geant4 sims:
    a. macro Cf252_source_iso.mac
    b. Geant 4 file: sbcGeant4.11/RL_branch
        i. in material file, I have included different PE density, CF4 densities on different temperatures
        ii. in construction file, I have added active source volume -> cf_active
        iii. I usually put what information need to saved in the StepAction.cc file, feel free to change it as your need
            Now it only save steps at particular volumes and for particular particles(neutron/Ar isotopes)
2.
