# Unicycle (unified cycles of earthquakes)
Unicycle is a community code for simulating the evolution of fault slip and distributed deformation using the integral method under the radiation damping approximation. The method is computationally efficient and may be used to simulate to simulate earthquake cycles (with radiation damping), afterslip, slow-slip events, viscoelastic relaxation and other aseismic processes.

If you find any bugs or wish to get involved in the development of Unicycle, please contact James D P Moore (earth@jamesdpmoore.com).

## Authors
James D P Moore, 
Sylvain Barbot,
Lujia Feng,
Yu Hang,
Valere Lambert,
Eric Lindsey,
Sagar Masuti,
Takanori Matsuzawa,
Jun Muto,
Priyamvada Nanjundiah,
Rino Salman,
Sharadha Sathiakumar,
and Harpreet Sethi

## Details
Unicycle implements a quasistatic approximation for earthquake cycles on rectangular or triangular fault patches, with both full rate and state frictional laws or rate-strengthening afterslip. We have also implemented off-fault viscoelastic deformation using Greene's functions for volumetric deformation and a range of rheological laws, including: Maxwell viscoelastic, Burgers viscoelastic, power-law viscoelastic, and power-law Burgers viscoelastic. Components from this code have also been used to calculate stress or displacement kernels to be used in inversion strategies.

## Languages
Currently, codes exist in Matlab and Fortran, but we hope to expand that to other languages in the future. A Python version is under active development, for more details or to get involved in developing out community code please email (earth@jamesdpmoore.com).

## Acknowledgements
This software has benefitted from discussions with a wide range of people over the years, leading to new and interesting projects, including:
Adam Switzer, Emma Hill, Judith Hubbard, Jun Muto, Rafael Almeida, Raquel Felix, Rishav Malick, and many others.

This research was supported by the National Research Foundation (NRF) of Singapore under the NRF Fellowship scheme (National Research Fellow Award NRF-NRFF2013-04), the Earth Observatory of Singapore under the Research Centres of Excellence initiative, and the Singapore Ministry of Education (Tier 2 project ARC5/14) the Natural Environment Research Council, UK (NE/R00515X/1) and the Royal Society of New Zealand (14-VUW-085).

## Licence
The intellectual property for this software belongs to Nanyang Technological University, Singapore, and it is distributed under the Creative Commons CC BY-NC-SA 4.0 licence subject to the condition that this readme shall be included in all copies or substantive portions of the Software.
https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode

## References
If you use this code, please cite as James D P Moore, Sylvain Barbot, Lujia Feng, Yu Hang, Valere Lambert, Eric Lindsey, Sagar Masuti, Takanori Matsuzawa, Jun Muto, Priyamvada Nanjundiah, Rino Salman, Sharadha Sathiakumar, and Harpreet Sethi. (2019, September 25). jdpmoore/unicycle: Unicycle (Version 1.0). Zenodo. http://doi.org/10.5281/zenodo.4471162       

Depending on which Green's functions you use, also cite:

Rectangular dislocations: Okada, Y. (1992), Internal deformation due to shear and tensile faults in a half-space, Bull. Seism. Soc. Am., 82, 1018–1040.

Triangular dislocations: Nikkhoo, M., and T. R. Walter (2015), Triangular dislocation: an analytical, artefact-free solution, Geophys. J. Int., 201(2), 1119–1141, https://doi.org/10.1093/gji/ggv035.

Distributed deformation: Barbot S, Moore J D P, Lambert V. (2017), Displacements and stress associated with distributed anelastic deformation in a half-space. Bulletin of the Seismological Society of America, 107(2):821, https://doi.org/10.1785/0120160237

## Related papers
A number of papers which have been published using this code, or utilise some of its constituent components. Some of these papers include:

Lindsey E O, Mallick R, Hubbard J, Bradley K E, Almeida R V, Moore J D P, Bürgmann R, Hill E M. Slip rate deficit and earthquake potential on shallow megathrust. (2021), Nature Geoscience. https://doi.org/10.1038/s41561-021-00736-x

Mallick, R., Meltzner, A.J., Tsang, L.L.H. et al. Long-lived shallow slow-slip events on the Sunda megathrust. Nat. Geosci. 14, 327–333 (2021). https://doi.org/10.1038/s41561-021-00727-y

Moore J D P, Yu H, Tang C, Wang T, Barbot S, Peng D, Masuti S, Dauwels J, Hsu Y, Lambert V, Nanjundiah P, Wei S, Lindsey E, Feng L, Shibazaki B. (2017) Imaging the distribution of transient viscosity after the 2016 Mw7.1 Kumamoto earthquake. Science, Vol. 356, Issue 6334, pp. 163-167
https://doi.org/10.1126/science.aal3422

Qiu, Q., Moore, J.D.P., Barbot, S., L Feng L., Hill E. M. Transient rheology of the Sumatran mantle wedge revealed by a decade of great earthquakes. Nat Commun 9, 995 (2018). https://doi.org/10.1038/s41467-018-03298-6

Lamb S, Arnold R, Moore J D P. (2018), Locking on a megathrust as a cause of distributed faulting and ‘fault jumping’ earthquakes. Nature Geoscience 11, 871–875 (2018). https://doi.org/10.1038/s41561-018-0230-5

Tang C, Hsu Y, Barbot S, Moore J D P, Chang W. (2019) Lower-crustal rheology and thermal gradient in the Taiwan orogenic belt illuminated by the 1999 Chi-Chi earthquake. Science Advances, Vol. 5, no. 2, eaav3287. https://doi.org/10.1126/sciadv.aav3287

Muto J, Moore J D P, Barbot S, Iinuma T, Ohta Y, Horiuchi S, Hikaru I. (2019), Coupled afterslip and transient mantle flow after the 2011 Tohoku earthquake.  Science Advances, Vol. 5, no. 9, eaaw1164. https://doi.org/10.1126/sciadv.aaw1164
