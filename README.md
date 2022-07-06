# PSH
## Planetary Spherical Harmonics community code

These are Community Codes for Planetary Spherical Harmonic Internal Field codes, in **MATLAB**, **IDL** and **Python 3**.  These are platform independent, and should work on PC, Mac or Linux. These are essentially the same code translated in to the three languages, with our testing so far the 3 languages give the same results to less than 10<sup>-11</sup> nT (rounding errors).

Just download the particular file you want (language, particular model, and if Cartesian (xyz) or Spherical (rtp)) to your local directory, and run.  No 'install' required and each code is independent.  They can be run with scalar inputs, or with 1D vector inputs.  For MATLAB, the 1D vector must be a column vector not a row vector (while the code could check and transpose if neccessary, that would slow it down.) For Python, you must have NumPy installed and the inputs can be a list or NumPy array.

The comments at the top of each code provide citations to papers where the g and h coefficients were sourced, and other info.
All inputs are in a **right-handed System III (1965)** frame, with distances in units of R<sub>J</sub>, where **1 R<sub>J</sub> = 71,492 km always**, and angles (for spherical) in units of radians.  All these codes are set up to use degree = order.

These codes will take inputs in Cartesian ([*x*,*y*,*z*] in units of R<sub>J</sub>) or Spherical ([*r*, *Co-Lat.*, *East Long.*] in units of R<sub>J</sub> and radians) and return the internal planetary field model values ([*B<sub>x</sub>*, *B<sub>y</sub>*, *B<sub>z</sub>*] or [*B<sub>r</sub>*, *B<sub>theta</sub>*, *B<sub>phi</sub>*] respectively, in units of nT). For older models that used different values of R<sub>J</sub>, our codes below still expect inputs where R<sub>J</sub> = 71,492 km, and will adjust the R<sub>J</sub> for you within the code.

Files can be found under the Planet directory, then language subdirectories.  File names are similar to "*jovian_jrm33_order13_internal_xyz*", where underscores separate out the information, e.g. Planet = Jovian (= Jupiter), model = JRM33, order (= degree) is 13, internal to highlight this is just the internal field (excluding things like magnetic fields due to a current disk), and xyz = Cartesian inputs/outputs (or rtp = Spherical inputs/outputs).

*The directory Mother_Source is not expected to be used by you, but contains the script that writes out all the MATLAB, IDL, Python codes, and is only here for completeness. Ignore it.*

These run quickly, but the higher the order, the longer it takes simply due to more computations required.  The actual speed will depend on your machine, your operating system, and which versions of IDL, Matlab or Python 3 you have installed.  In a test in 2022, running jrm33 order 13 in MATLAB on a Mac for a vector of ~2 million points took around 5s. 
