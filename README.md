# tert-Butyl_Hydroperoxide_Dimer
1D local mode model code used to calculate frequencies and intensities used in the published paper "FTIR Detection of the tert-Butyl Hydroperoxide Dimer"

All electronic structure calculations are performed in Gaussian 16

Once an equilibrium structure of the dimer is obtained, make a single point electronic energy .gjf file. The script "displaceOH" takes the .gjf file and creates a set of new .gjf files with the H displaced in increments chosen by the user.

After all the single point calculations have run, open the script "get_E_DM". This makes a textfile collecting the HF energy and dipole moment of all single points. The patterns that is searched for to find the numbers in the files might need to be replaced depending on the format of the output file from gaussian.

The textfile is saved and can be opened by the 3rd script "1D_LM". Now the vibrational schrödinger equation is solved and the vibrational energies and intensities are calculated.

