# cosmic_rays

Calculates the cosmic ray proton and lepton (electron and positron) spectrum for starbursts galaxies. For protons, includes ionization, pion and wind advection losses. For leptons, includes ionization, bremsstrahlung, wind advection, inverse Compton and synchrotron losses. Also includes the knock-on secondary lepton spectrum. Once complete, the secondary lepton spectrum from pion decay will also be included.

General files: functions.py

Proton spectrum files: protons.py

Lepton spectrum files: 
(1) master - electrons.py
(2) IC losses - ic.py, qint.py, mUph.py, ir.csv
(3) secondaries - secondary_e.py, pion_param.py, pion_sec_e.py, ko_sec_e.py

Synchrotron files (included here for completeness, but they are not necessary since naima calculates synchrotron emission with enough precision): synchrotron.py, bessel.py
