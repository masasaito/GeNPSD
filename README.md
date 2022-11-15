--------------- READ-ME file for the GeNPSD.package -----------------

 This is a Generalized Nonspherical Particle Size Distribution (GeNPSD) package
 for converting particle size distributions of nonspherical particles in
 atmospheric science applications.


 Lead developer: Dr. Masanori Saito (masa.saito@tamu.edu)


 History
 (MM/DD/YYYY)	(Author)	(Contents)
 11/15/2022	Masa Saito	Version 1.0.0 released.




1. General information

   Reference:

     1). Saito, M., and P. Yang, Generalization of Atmospheric 
         Nonspherical Particle Size: Interconversions of Size Distributions
         and Optical Equivalence, J. Atmos. Sci., in press, 
         https://doi.org/10.1175/JAS-D-22-0086.1.


   GeNPSD.package consists of the following directories:
     
     src			:: Source codes
     src/hparx-20160408  	:: HPARX library
     bin			:: Executable command path
     examples			:: Example files for test runs

   The latest version of the source codes and examples may be 
   available from https://github.com/masasaito/GeNPSD.




 2. Installation

   First, to install the database creation codes, edit the following 
   files:

     ./src/Common.mk
     ./src/hparx/Common.mk  
     
   In these files, users can specify their compiler (FC; e.g., gfortran)
   and associate options (FCFLAGS).

   Then, command as follows:

     $ cd ./src/hparx-20160408
     $ make
     $ cd ../
     $ make install
 
   In the end, the following two executable files will be created in the 
   "bin" directory:

     genpsd_generate	:: Create nonspherical PSD for a given size descriptor 
                           and physical variables.
     genpsd_interconv   :: convert nonspherical PSD for a set of given size 
                           descriptors and physical variables.
      
   Please read Usage.txt carefully to learn how to run these codes.




  3. Requirements

    The lead developer, Dr. Masanori Saito encourages the research
    community to utilize GeNPSD.package for atmospheric science applications. 
    The only requirement in regards to utilizing this procedure is to 
    acknowledge our contribution in a paper to be published by citing 
    Saito and Yang (2022) in a relevant section of the main text. 
      



  4. Contact Information

     Dr. Masanori Saito (masa.saito@tamu.edu) 
     Department of Atmospheric Sciences, Texas A&M University. 

 
