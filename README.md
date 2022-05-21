# sclsVofFoam-OpenFoam9
Coupled level-set VOF method implementation by Sankar Menon (http://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2015/SankarMenon/Report_SankarMenon.pdf) has been adapted for OpenFoam 9

Currently the solver is set to compile from user libraries present in (https://github.com/Venky-94/interFoamMod), as I'm testing out same changes to the level set method. But the rest of the solver has been adapted as is from Sankar Menon's report and anyone should be able to compile the solver standalone by swapping the user libraries to OpenFoam's default libraries in Make/options file.
