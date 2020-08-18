# Vargeom-Poisson

* Program calculating the Poisson potential due to charges and dipoles (if present) in model supercapacitors with electrodes of various geometries.

* The equation for the Poisson potential can be found for example in J. Phys. Chem. C, 115, 16613 (2011) - equation 6.

* To compile the program, you can use any c compiler, e.g.:

      gcc -o Vargeom_Poisson.x Vargeom_Poisson.c -lm

* To execute the program, you just have to do:

      ./program < inputfile
      
Some input and output examples are given. These examples can be used to test the program but the systems they correspond to are not described since the number of configurations is not sufficient to obtain a suitable convergence of the Poisson potential.
      
The input file should provide the following information

1. Name of the xyz file with the positions

2. Value to indicate the existence of dipoles: 

      0 if no dipoles, 1 if dipoles

      Supplementary line if there are dipoles: name of the dipoles' file and number of columns (3 if just dipoles, 4 if indices in the first column)

3. Number of configurations and number of independent time zones (put 1 time zone to analyse the full trajectory)

4. Direction along which the ionic densities will be calculated (x, y or z)

5. Length of the box in the three cartesian directions (one value per line)

6. Number of bins for the density calculation

7. Number of different species (all ion/atom types for mobile and fix species)

8. Number of atoms/ions for each species

9. Number of species with non fluctuating charge

10. Type and charge for each species with non fluctuating charge

11. Initial condition for the charge potential psiq(z0) in volts

12. Initial condition for the dipole potential psimu(z0) in volts (0.0 if there are no dipoles)

13. Number of configurations to skip at the beginning of the files (put 0 to analyse the full trajectory)

14. If there are atoms with gaussian charge densities:

      Name of the file containing the fixed atoms charges for each configuration

      Number of lines with comments to skip for each config

      Width of the gaussian for the metallic atoms in Anstroms

      Shift to integrate before and after the physical ends of the cell (related to the Gaussian width of fluctuating charges)
