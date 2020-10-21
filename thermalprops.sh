#!/bin/bash

echo Please identify this structure.
read struc

rm t*
phonopy -p -s mesh.conf > thermalprops
cat thermalprops | tail -n 208 | head -n 201 > thermalprops${struc}.txt
cp thermalprops${struc}.txt /Users/mujanseif/Documents/Research/cathodes/W-surfaces/thermal-properties/surface-energy-plots/
rm thermalprops
