# PSH
## Planetary Spherical Harmonics community code

These are Community Codes for Planetary Spherical Harmonic Internal Field codes, in **MATLAB**, **IDL** and **Python 3**.  These are platform independent, and should work on PC, Mac or Linux. These are essentially the same code translated in to the three languages, with our testing so far the 3 languages give the same results to less than 10<sup>-11</sup> nT (rounding errors).

Just download the particular file you want (language, particular model, and if Cartesian (xyz) or Spherical (rtp)) to your local directory, and run.  No 'install' required and each code is independent.  They can be run with scalar inputs, or with 1D vector inputs.  For MATLAB, the 1D vector must be a column vector not a row vector (while the code could check and transpose if neccessary, that would slow it down.) For Python, you must have NumPy installed and the inputs can be a list or NumPy array.

The comments at the top of each code provide citations to papers where the *g* and *h* coefficients were sourced, and other info.
All inputs are in a **right-handed System III (1965)** frame, with distances in units of R<sub>J</sub>, where **1 R<sub>J</sub> = 71,492 km always**, and angles (for spherical) in units of radians.  All these codes are set up to use degree = order.

These codes will take inputs in Cartesian ([*x*,*y*,*z*] in units of R<sub>J</sub>) or Spherical ([*r*, *Co-Lat.*, *East Long.*] in units of R<sub>J</sub> and radians) and return the internal planetary field model values ([*B<sub>x</sub>*, *B<sub>y</sub>*, *B<sub>z</sub>*] or [*B<sub>r</sub>*, *B<sub>theta</sub>*, *B<sub>phi</sub>*] respectively, in units of nT). For older models that used different values of R<sub>J</sub>, our codes below still expect inputs where R<sub>J</sub> = 71,492 km, and will adjust the R<sub>J</sub> for you within the code.

Files can be found under the Planet directory, then language subdirectories.  File names are similar to "*jovian_jrm33_order13_internal_xyz*", where underscores separate out the information, e.g. Planet = Jovian (= Jupiter), model = JRM33, order (= degree) is 13, internal to highlight this is just the internal field (excluding things like magnetic fields due to a current disk), and xyz = Cartesian inputs/outputs (or rtp = Spherical inputs/outputs).

*The directory Mother_Source is not expected to be used by you, but contains the script that writes out all the MATLAB, IDL, Python codes, and is only here for completeness. Ignore it.*

These run quickly, but the higher the order, the longer it takes simply due to more computations required.  The actual speed will depend on your machine, your operating system, and which versions of IDL, Matlab or Python 3 you have installed.  In a test in 2022, running jrm33 order 13 in MATLAB on a Mac for a vector of ~2 million points took around 5s. 

## Spherical Harmonic Models Included

We recommend using JRM33 order 13, but here is the list of existing models.  All were developed in a right-handed System III (1965) system (some earlier models used System III (1957)), but assumed different values for what 1 R<sub>J</sub> was.

| Model | Order | Planet  | Original Model Reference | R<sub>J</sub> used in original work | Reference used for *g* and *h* values | Other references for *g* and *h* values |
| ----- | ----- | ------- |----------------------- | ---- | ---------------------- | ---------------------- |
| JRM33 |   18  | Jupiter | *Connerney et al.*, 2022  | 71,492 km | *Connerney et al.*, 2022 | *Connerney*, 2022 |
| **JRM33** |  **13**   | **Jupiter** | ***Connerney et al., 2022*** | 71,492 km | ***Connerney et al.***, 2022 | *Connerney*, 2022 |
| JRM09 |   10  | Jupiter | *Connerney et al.*, 2018 | 71,492 km | *Connerney et al.*, 2018 | *Connerney*, 2022 |
| ISaAC |   10  | Jupiter | *Hess et al.*, 2017       | 71,492 km | *Hess et al.*, 2017 | |
| VIPAL |    5  | Jupiter | *Hess et al.*, 2011       | 71,492 km | *Hess et al.*, 2011 | |
| VIT4  |    4  | Jupiter | *Connerney et al.*, 1998  | 71,323 km | *Connerney* 2007 | *Hess et al.*, 2011  |
| VIP4  |    4  | Jupiter | *Connerney et al.*, 1998  | 71,323 km | *Connerney* 2007 | *Hess et al.*, 2011 |
| O6    |    3  | Jupiter | *Connerney* 1992          | 71,372 km |*Connerney* 1992 | *Connerney et al.*, 1998 |

The reference papers may provide *g* and *h* values to higher orders, but the authors do not always trust those higher order values (see their papers). Hence the order used here may be lower than given what you can find in publications.  In the case of JRM33, the authors used both order 13 and order 18 for plots in their paper, hene we provide code for both, but we recommend using JRM order 13 for your studies (*personal communication with authors*).

## References

- Connerney, J. E. P., Acuna, M. H., Ness, N. F. (1982). Voyager 1 assessment of Jupiter's planetary magnetic field. *J. Geophys. Res.*, 87 (A5), 3623-3623. doi: [10.1029/JA087iA05p03623](https://doi.org/10.1029/JA087iA05p03623)
- Connerney, J. E. P. (1992). Doing more with Jupiter's magnetic field. In *Planetary radio emissions iii* (p. 13-33).
- Connerney, J. E. P., Acuna, M. H., Ness, N. F.,  Satoh, T. (1998). New models of Jupiter's magnetic field constrained by the Io flux tube footprint. *J. Geophys. Res.*, 103 (A6), 11929-11940. doi: [10.1029/97JA03726](https://doi.org/10.1029/97JA03726)
- Connerney, J. E. P. (2007). Planetary Magnetism. In G. Schubert (Ed.), *Planets and moons* (Vol. 10, p. 243-280). doi: [10.1016/B978-044452748-6.00159-0](https://doi.org/10.1016/B978-044452748-6.00159-0)
- Connerney, J. E. P., et al. (2018). A New Model of Jupiter's Magnetic Field From Juno's First Nine Orbits. *Geophys. Res. Lett.*, 45 (6), 2590-2596. doi: [10.1002/2018GL077312](https://doi.org/10.1002/2018GL077312)
- Connerney, J. E. P., et al. (2022). A New Model of Jupiter's Magnetic Field at the Completion of Juno's Prime Mission. *Journal of Geophysical Research (Planets)*, 127 (2), e07055. doi: [10.1029/2021JE007055](https://doi.org/10.1029/2021JE007055)
- Connerney, J. (2022). Juno magnetometer jupiter archive. NASA Planetary Data System. Retrieved from https://pds.nasa.gov/ds-view/pds/viewDataset.jsp?dsid=JNO-J-3-FGM-CAL-V1.0 doi: [10.17189/1519711](https://doi.org/10.17189/1519711)
- Hess, S. L. G., Bonfond, B., Zarka, P., Grodent, D. (2011). Model of the Jovian magnetic field topology constrained by the Io auroral emissions. *Journal of Geophysical Research (Space Physics)*, 116 (A5), A05217. doi: [10.1029/2010JA016262](https://doi.org/10.1029/2010JA016262)
- Hess, S. L. G., Bonfond, B., Bagenal, F., Lamy, L. (2017). A model of the Jovian internal field derived from in-situ and auroral constraints. In G. Fischer, G. Mann, M. Panchenko, & P. Zarka (Eds.), *Planetary radio emissions viii* (p. 157-167). doi: [10.1553/PRE8s157](https://doi.org/10.1553/PRE8s157)
