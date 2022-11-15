#!/bin/bash

## Initialization
mpsd=1
nsiz=5001
sizmin=0.01
sizmax=1000.0
sizlog=1
mgeoX=0
msizP=0
mlogP=1
mlogQ=1
mnor=0
exe='../src/genpsd_interconv'


# Jobs
for msizQ in 0 1 2 3
do

   # Configuration
   rmed=0.8
   s=2.0
   ntot=300.0   #N/cm^3 for Ntot; µm^2/cm^3 for Atot; µm^3/cm^3 for Vtot
   psdtyp=0
   outfile="genpsd_interconversionDust_rmed${rmed}s2.0_msizQ${msizQ}.txt"

   for mgeoY in 0 1 2
   do
      for psi in 0.695 0.785
      do

         if   [ $psi == '0.695' ]; then
            dmax=2.0
            are=1.0483758E+00
            vol=4.6786062E-01 
         elif [ $psi == '0.785' ]; then
            dmax=2.0
            are=1.5984287E+00
            vol=1.0573242E+00 
         fi

         ${exe} $psdtyp $mgeoX $msizP $mlogP $mgeoY $msizQ $mlogQ $mnor $dmax $are $vol $nsiz $sizmin $sizmax $sizlog $mpsd $ntot $rmed $s > tmp${psi}${mgeoY}
      done
   done
   paste tmp0.6950 tmp0.6951 tmp0.6952 tmp0.7850 tmp0.7851 tmp0.7852 > ${outfile}
done 
rm tmp*

