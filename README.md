# sclsVofFoam-OpenFoam9
Coupled level-set VOF method implementation by Sankar Menon (http://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2015/SankarMenon/Report_SankarMenon.pdf) has been adapted for OpenFoam 9.

Currently the solver is set to compile from user libraries present in (https://github.com/Venky-94/interFoamMod), as I'm trying to implement a change in the property estimation equation to the level set method, based on https://doi.org/10.1115/1.4053548, to improve results in systems with higher viscosity ratios. 

But the rest of the solver has been adapted as is from Sankar Menon's report and anyone should be able to compile the solver standalone by swapping the user libraries to OpenFoam's default libraries in Make/options file.
