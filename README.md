# PSH: Planetary Spherical Harmonics community code

**Authors**:  M. Brennan, M.K. James, G. Provan, M.F. Vogt and R.J. Wilson

**Citation**: *TBD*

**Thanks to**: F. Bagenal for the original idea, Masafumi Imai for sharing their original code we could test against, Stan Cowley for helpful discussions/advice and Jack Connerney for verifcations!

**Table of Contents**:
- [Initial Problem](README.md#initial-problem)
- [Solution](README.md#solution)
- [Spherical Harmonic Models Included](README.md#spherical-harmonic-models-included)
- [Examples](README.md#examples)
  - [MATLAB](README.md#matlab)
  - [IDL](README.md#idl)
  - [Python 3](README.md#python-3)
- [Solution #2: JupiterMag](README.md#solution-2-jupitermag)
- [Speed Tests](README.md#speed-tests)
- [References](README.md#references)


## Initial Problem
There are many past models, with the *g* and *h* coefficients often quoted in multiple places: with different units (G or nT, or stated as G but really nT), or different levels of precision (and occasionally some typos creep in between copies), or were originally in a coordinate system that wasn't right-handed System III (1965), or assumed a different planetary radius (Jovian Radius = R<sub>J</sub>).  Different users making they own codes may get different results from each other depending on which paper they used for their *g* and *h* values, to what precision their $g$ and $h$ values were, and we suspect users may have forgotten to adjust older models for their different R<sub>J</sub>.  The aim here is to have a standard set of efficient codes, that people can just use and cite.  

## Solution
These are Community Codes for Planetary Spherical Harmonic Internal Field codes, in [**MATLAB**](https://www.mathworks.com/products/matlab.html), [**IDL**](https://www.l3harrisgeospatial.com/Software-Technology/IDL) and [**Python 3**](https://www.python.org/).  These are platform independent, and should work on PC, Mac or Linux. These are essentially the same code translated in to the three languages, with our testing so far the 3 languages give the same results to less than 10<sup>-11</sup> nT (rounding errors).

Just download the particular file you want (language, particular model, and if Cartesian (xyz) or Spherical (rtp)) to your local directory, and run.  No 'install' required and each code is independent.  They can be run with scalar inputs, or with 1D vector inputs.  For MATLAB, the 1D vector must be a column vector not a row vector (while the code could check and transpose if neccessary, that would slow it down.) For Python, you must have NumPy installed and the inputs can be a list or NumPy array.

The comments at the top of each code provide citations to papers where the *g* and *h* coefficients were sourced, and other info.
All inputs are in a **right-handed System III (1965)** frame, with distances in units of R<sub>J</sub>, where **1 R<sub>J</sub> = 71,492 km always**, and angles (for spherical) in units of radians.  All these codes are set up to use degree = order, and for convience, we simply refer to both degree and order as order.

These codes will take inputs in Cartesian ([*x*,*y*,*z*] in units of R<sub>J</sub>) or Spherical ([*r*, *Co-Lat.*, *East Long.*] in units of R<sub>J</sub> and radians) and return the internal planetary field model values ([*B<sub>x</sub>*, *B<sub>y</sub>*, *B<sub>z</sub>*] or [*B<sub>r</sub>*, *B<sub>theta</sub>*, *B<sub>phi</sub>*] respectively, in units of nT). For older models that used different values of R<sub>J</sub>, our codes below still expect inputs where R<sub>J</sub> = 71,492 km, and will adjust the R<sub>J</sub> for you within the code.

Files can be found under the Planet directory, then language subdirectories.  File names are similar to "*jovian_jrm33_order13_internal_xyz*", where underscores separate out the information, e.g. Planet = Jovian (= Jupiter), model = JRM33, order (= degree) is 13, internal to highlight this is just the internal field (excluding things like magnetic fields due to a current disk), and xyz = Cartesian inputs/outputs (or rtp = Spherical inputs/outputs).

*The directory Mother_Source is not expected to be used by you, but contains the script that writes out all the MATLAB, IDL, Python codes, and is only here for completeness. Ignore it.*

These run quickly, but the higher the order, the longer it takes simply due to more computations required.  The actual speed will depend on your machine, your operating system, and which versions of IDL, Matlab or Python 3 you have installed.  In a test in 2022, running jrm33 order 13 in MATLAB on a Mac for a vector of ~2 million points took around 5s. 

## Spherical Harmonic Models Included

We recommend using JRM33 order 13, but here is the list of existing models.  All these were developed in a right-handed System III (1965) system (some earlier models had used System III (1957)[^1]), but assumed different values for what 1 R<sub>J</sub> was.
[^1]: The SHA model was originally calculated in the System III (1957) frame, but *Connerney* 2007 converted the *g* and *h* coefficients in to System III (1965), however G<sub>3</sub><sup>0</sup> is listed as 11,300 nT, but we believe this has a sign typo and should be -11,300 nT.  The R<sub>J</sub> value is not specified, but we believe they used 71,398 km.

| Model | Order | Planet  | Original Model Reference | R<sub>J</sub> used in original work | Reference used for *g* and *h* values | Other references for *g* and *h* values |
| ----- | ----- | ------- |----------------------- | ---- | ---------------------- | ---------------------- |
| JRM33 |   18  | Jupiter | *Connerney et al.*, 2022  | 71,492 km | *Connerney et al.*, 2022 | *Connerney*, 20120^4] |
| **JRM33** |  **13**   | **Jupiter** | ***Connerney et al., 2022*** | 71,492 km | ***Connerney et al.***, 2022 | *Connerney*, 2020[^4] |
| JRM09 |   10  | Jupiter | *Connerney et al.*, 2018 | 71,492 km | *Connerney et al.*, 2018 | *Connerney*, 2020[^4] |
| ISaAC |   10  | Jupiter | *Hess et al.*, 2017       | 71,492 km | *Hess et al.*, 2017 | |
| VIPAL |    5  | Jupiter | *Hess et al.*, 2011       | 71,492 km | *Hess et al.*, 2011 | |
| VIT4  |    4  | Jupiter | *Connerney et al.*, 1998  | 71,323 km[^2] | *Connerney* 2007 | *Hess et al.*, 2011[^3]  |
| VIP4  |    4  | Jupiter | *Connerney et al.*, 1998  | 71,323 km[^2] | *Connerney* 2007 | *Hess et al.*, 2011 |
| O6    |    3  | Jupiter | *Connerney* 1992          | 71,372 km |*Connerney* 1992 | *Connerney et al.*, 1998 |
[^2]: R<sub>J</sub> = 71,323 km based on table 1 of the original paper, and thus used here.  However the original paper also states a value of 71,398 km earlier in the text, while the 2007 book suggests 71,372 km earlier in the text then doesn't explicitly list an R<sub>J</sub> with the table for VIP4 and VIT4 coefficients (yet does list specific Rs for some other models).  However, Connerney (private communication) says to use the earliest publication, hence table 1 of original paper.
[^3]: h<sub>4</sub><sup>4</sup> has a typo, probably should be 0.1264 G.
[^4]: The PDS archive is cited as year 2017.  The models were added to this dataset in later years, but same doi for the whole dataset.

The reference papers may provide *g* and *h* values to higher orders, but the authors do not always trust those higher order values (see their papers). Hence the order used here may be lower than given what you can find in publications.  In the case of JRM33, the authors used both order 13 and order 18 for plots in their paper, hene we provide code for both, but we recommend using JRM order 13 for your studies (*personal communication with authors*).

## Examples
For all 3 languages, the output's 1<sup>st</sup> dimension is always number of records (=1 if scalar) and the 2<sup>nd</sup> dimension is always size 3 for the B-vector components, but some langues are row-major, other column-major. This is best seen in the example outputs below that all give the same inputs in each test, but the outputs may be transposed from each other.

The following examples (same sitution for each language) all use *jovian_jrm33_order13_internal_rtp* and *jovian_jrm33_order13_internal_xyz*, but these can be simply swapped with any of the other models files in this collection, they all have the same input and output formats.

### MATLAB
```MATLAB
function  Matlab_test
 
% Spherical coordinate example for scalar at 10 Rj, Colatitude on equator
% at East longitude of 38 degrees (converted to radians)
Brtp_scalar = jovian_jrm33_order13_internal_rtp(10, pi/2, 38*pi/180);
Brtp_scalar % print to screen
 
% Spherical coordinate example.
% 4 Quadrants all on the equator, with increasing r
r = [8;9;10;11];
t = pi/2 *[1;1;1;1];
p = [0; pi/2; pi; 1.5*pi];
Brtp = jovian_jrm33_order13_internal_rtp(r, t, p);
Brtp % print to screen
 
% Same 4 positions but now in Cartesian
x = [  8;  0;-10;  0];
y = [  0;  9;  0;-11];
z = [  0;  0;  0;  0];
Bxyz = jovian_jrm33_order13_internal_xyz(x, y, z);
Bxyz % print to screen
 
whos Brtp_scalar Brtp Bxyz
```
Gives the output:
```
Brtp_scalar =

  -79.8398  399.4828  -53.4823

Brtp =

 -250.0396  779.3628  -48.0068
   38.3079  551.8268  -92.0848
  152.3795  421.1121   17.0471
  -42.3894  313.6507   55.7882

Bxyz =

 -250.0396  -48.0068 -779.3628
   92.0848   38.3079 -551.8268
 -152.3795  -17.0471 -421.1121
   55.7882   42.3894 -313.6507

  Name             Size            Bytes  Class     Attributes

  Brtp             4x3                96  double              
  Brtp_scalar      1x3                24  double              
  Bxyz             4x3                96  double              
```

### IDL
```IDL
PRO IDL_Test
  ; Spherical coordinate example for scalar at 10 Rj, Colatitude on equator
  ; at East longitude of 38 degrees (converted to radians)
  Brtp_scalar = jovian_jrm33_order13_internal_rtp(10d, !DPI/2d, 38d*!DPI/180d)
  PRINT,'Brtp_scalar = ‘
  PRINT, Brtp_scalar

  ; Spherical coordinate example.
  ; 4 Quadrants all on the equator, with increasing r
  r = [     8d,      9d,     10d,       11d]
  t = [!DPI/2d, !DPI/2d, !DPI/2d, !DPI/2d  ]
  p = [     0d, !DPI/2d, !DPI   , !DPI*1.5d]
  Brtp = jovian_jrm33_order13_internal_rtp(r, t, p)
  PRINT,'Brtp = ‘
  PRINT, Brtp

  ; Same 4 positions but now in Cartesian
  x = [  8d,  0d,-10d,  0d]
  y = [  0d,  9d,  0d,-11d]
  z = [  0d,  0d,  0d,  0d]
  Bxyz = jovian_jrm33_order13_internal_xyz(x, y, z)
  PRINT,'Bxyz = ‘
  PRINT, Bxyz

  HELP, Brtp_scalar, Brtp, Bxyz
END
```
Gives the output:
```
Brtp_scalar = 
      -79.839826
       399.48279
      -53.482325
Brtp = 
      -250.03964       38.307890       152.37954      -42.389403
       779.36280       551.82685       421.11209       313.65069
      -48.006775      -92.084823       17.047116       55.788200
Bxyz = 
      -250.03964       92.084823      -152.37954       55.788200
      -48.006775       38.307890      -17.047116       42.389403
      -779.36280      -551.82685      -421.11209      -313.65069
BRTP_SCALAR     DOUBLE    = Array[1, 3]
BRTP            DOUBLE    = Array[4, 3]
BXYZ            DOUBLE    = Array[4, 3]
```

### Python 3
```Python
import numpy as np
import jovian_jrm33_order13_internal_rtp as jrm33o13_rtp
import jovian_jrm33_order13_internal_xyz as jrm33o13_xyz

# Spherical coordinate example for scalar at 10 Rj, Colatitude on equator
# at East longitude of 38 degrees (converted to radians)
Brtp_scalar = jrm33o13_rtp.jovian_jrm33_order13_internal_rtp(10, np.pi/2, 38*np.pi/180);
print('Brtp_scalar = ')
print(Brtp_scalar)

# Spherical coordinate example.
# 4 Quadrants all on the equator, with increasing r
# Can be numpy arrays or not, if you use pi, should be
r =          [      8,       9,      10,        11]
t = np.array([np.pi/2, np.pi/2, np.pi/2, np.pi/2  ])
p = np.array([      0, np.pi/2, np.pi  , np.pi*1.5])
Brtp = jrm33o13_rtp.jovian_jrm33_order13_internal_rtp(r, t, p);
print('Brtp = ')
print(Brtp)

# Same 4 positions but now in Cartesian
# Can be numpy arrays, but since just numbers here, does not have to be
x = [  8,  0,-10,  0]
y = [  0,  9,  0,-11]
z = [  0,  0,  0,  0]
Bxyz = jrm33o13_xyz.jovian_jrm33_order13_internal_xyz(x, y, z);
print('Bxyz = ')
print(Bxyz)

print('Shape of Brtp_scalar: ', Brtp_scalar.shape)
print('Shape of Brtp       : ', Brtp.shape)
print('Shape of Bxyz       : ', Bxyz.shape)
```

Gives the output:
```
Brtp_scalar = 
[-79.8398262  399.48279169 -53.48232536]
Brtp = 
[[-250.03964154  779.36280353  -48.0067748 ]
 [  38.30789003  551.82684557  -92.08482301]
 [ 152.37953988  421.11209401   17.04711562]
 [ -42.38940341  313.65069342   55.78819978]]
Bxyz = 
[[-250.03964154  -48.0067748  -779.36280353]
 [  92.08482301   38.30789003 -551.82684557]
 [-152.37953988  -17.04711562 -421.11209401]
 [  55.78819978   42.38940341 -313.65069342]]
Shape of Brtp_scalar:  (1, 3)
Shape of Brtp       :  (4, 3)
Shape of Bxyz       :  (4, 3)
```

## Solution #2: JupiterMag
There is sister community code that will do the same models here, and give the same results, over at [https://github.com/mattkjames7/JupiterMag](https://github.com/mattkjames7/JupiterMag).  This is a Python 3 package that requires a simple install, and has more flexibility than this code, e.g. you could have Cartesian inputs, but outputs in Spherical.  It also includes code for a current sheet, and field line tracing.

We have tested the JupiterMag codes against the codes here, and for same inputs we still get the same outputs to within the same rounding errors.

## Speed Tests
The following speed tests were done on a Mac in 2022, but speed depends on your physical computer, your operating system, what else you're running and even which version of IDL or Matlab you have.  e.g. IDL 8.4 took 17s to run our code in a test, but the same test on IDL 8.8 took 14s.
For below we test both the spherical (RTP) codes and Cartesian (xyz), when running 75641 test positions once as a vector, or as 75641 scalars in a FOR loop.  We show comparisions of the 3 langauge codes in this respository, and also the sister [*JupiterMag* code](README.md#solution-2-jupitermag).

![speedtest](https://user-images.githubusercontent.com/91491246/178030181-12f68efe-b109-4a75-9d55-bcc26be2ce84.png)


## References

- Connerney, J. E. P., Acuna, M. H., Ness, N. F. (1982). Voyager 1 assessment of Jupiter's planetary magnetic field. *J. Geophys. Res.*, 87 (A5), 3623-3623. doi: [10.1029/JA087iA05p03623](https://doi.org/10.1029/JA087iA05p03623)
- Connerney, J. E. P. (1992). Doing more with Jupiter's magnetic field. In *Planetary radio emissions iii* (p. 13-33).  [*No doi but try here for chapter.*](https://ui.adsabs.harvard.edu/abs/1992pre3.conf...13C/abstract) Book doi: [10.1553/0x0015ce1d](https://doi.org/10.1553/0x0015ce1d)
- Connerney, J. E. P., Acuna, M. H., Ness, N. F.,  Satoh, T. (1998). New models of Jupiter's magnetic field constrained by the Io flux tube footprint. *J. Geophys. Res.*, 103 (A6), 11929-11940. doi: [10.1029/97JA03726](https://doi.org/10.1029/97JA03726)
- Connerney, J. E. P. (2007). Planetary Magnetism. In G. Schubert (Ed.), *Planets and moons* (Vol. 10, p. 243-280). doi: [10.1016/B978-044452748-6.00159-0](https://doi.org/10.1016/B978-044452748-6.00159-0)
- Connerney, J. E. P., et al. (2018). A New Model of Jupiter's Magnetic Field From Juno's First Nine Orbits. *Geophys. Res. Lett.*, 45 (6), 2590-2596. doi: [10.1002/2018GL077312](https://doi.org/10.1002/2018GL077312)
- Connerney, J.E.P. (2020), Juno MAG CALIBRATED DATA J V1.0, JNO-J-3-FGM-CAL-V1.0, NASA Planetary Data System, doi: [10.17189/1519711](https://doi.org/10.17189/1519711)
- Connerney, J. E. P., et al. (2022). A New Model of Jupiter's Magnetic Field at the Completion of Juno's Prime Mission. *Journal of Geophysical Research (Planets)*, 127 (2), e07055. doi: [10.1029/2021JE007055](https://doi.org/10.1029/2021JE007055)
- Hess, S. L. G., Bonfond, B., Zarka, P., Grodent, D. (2011). Model of the Jovian magnetic field topology constrained by the Io auroral emissions. *Journal of Geophysical Research (Space Physics)*, 116 (A5), A05217. doi: [10.1029/2010JA016262](https://doi.org/10.1029/2010JA016262)
- Hess, S. L. G., Bonfond, B., Bagenal, F., Lamy, L. (2017). A model of the Jovian internal field derived from in-situ and auroral constraints. In G. Fischer, G. Mann, M. Panchenko, & P. Zarka (Eds.), *Planetary radio emissions viii* (p. 157-167). doi: [10.1553/PRE8s157](https://doi.org/10.1553/PRE8s157)
