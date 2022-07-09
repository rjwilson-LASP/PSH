FUNCTION jovian_vit4_order04_internal_rtp, r_rj, colat_rads, elong_rads
  ; Code to calculate the VIT4_ORDER04 model of Jupiter's internal magnetic field model
  ; with Degree 4 and Order 4.
  ; Reference: Connerney (2007), https://doi.org/10.1016/B978-044452748-6.00159-0
  ;
  ; Required inputs (System III (1965) Spherical, right handed, and assuming 1 Rj = 71492 km):
  ;  r_rj       - radial distance, in Rj.
  ;  colat_rads - colatitude, in radians.                    Value(s) should be 0 <= colat_rads <=  pi.
  ;  elong_rads - East longitude, right handed, in radians.  Value(s) should be 0 <= elong_rads <= 2pi.
  ;
  ; Outputs:
  ;  B - Spherical Magnetic field vector the VIT4_ORDER04 internal magnetic field model, [Br, Btheta, Bphi], units of nT.
  ;
  ; Usage:
  ; For internal field only: B = jovian_vit4_order04_internal_rtp(r_rj, colat_rads, elong_rads)
  ;
  ; This code was written by Marissa Vogt (mvogt@bu.edu) and Rob Wilson (rob.wilson@lasp.colorado.edu).
  ; It is based on a routine originally written by K. Khurana, translated into IDL by Marissa Vogt in 2009.
  ;
  ; Citation Info:
  ;  DOI: 10.5281/zenodo.6814109     This doi links to all versions of code at the Github.
  ;  Individual versions released on https://github.com/rjwilson-LASP/PSH) can have a different DOI,
  ;  See the DOI above for a list of DOIs for each specific Github released version.
  ; Version Info:
  ;  Last update of this file: 2022-07-09 11:37:47.511772 by user wilsonr. 
  ;  This code was re-written/re-formatted by Rob's python code:
  ;   /Users/wilsonr/Documents/JADE/Level2_Processing_Code/IDL/Field_Model/2022/Git_initial/Mother_Source/MOP_spherical.py
  ;   which itself was last updated at UTC 2022-07-09T17:37:36.
  ;
  ;  The Spherical Harmonic g and h values used for this order 4 code are below: 
  ;  
  ;  g[i,j] values (nT) used are:
  ; g[ 1, 0] =     428077, g[ 1, 1] =     -75306, 
  ; g[ 2, 0] =      -4283, g[ 2, 1] =     -59426, g[ 2, 2] =      44386, 
  ; g[ 3, 0] =       8906, g[ 3, 1] =     -21447, g[ 3, 2] =      21130, g[ 3, 3] =      -1190, 
  ; g[ 4, 0] =     -22925, g[ 4, 1] =      18940, g[ 4, 2] =      -3851, g[ 4, 3] =       9926, g[ 4, 4] =       1271, 
  ;
  ;  h[i,j] values (nT) used are:
  ;                        h[ 1, 1] =      24616, 
  ;                        h[ 2, 1] =     -50154, h[ 2, 2] =      38452, 
  ;                        h[ 3, 1] =     -17187, h[ 3, 2] =      40667, h[ 3, 3] =     -35263, 
  ;                        h[ 4, 1] =      16088, h[ 4, 2] =      11807, h[ 4, 3] =       6195, h[ 4, 4] =      12641, 

  ON_ERROR, 2 ; Exit code if an error in main, don't stop in code - no MATLAB equivalent, just delete line in MATLAB

  ; Check inputs are same size.
  N_input = N_ELEMENTS(r_rj)
  scalar_input = (N_input EQ 1)  ; scalar or not

  ; Check inputs r_rj, colat_rads and elong_rads are all numbers,  and same size (scalar or 1D only)
  IF (ISA(r_rj      , NUMBER=1) EQ 0) OR (SIZE(r_rj      , N_DIMENSIONS=1) GT 1) THEN MESSAGE,'ERROR: First  argument    r_rj    must be a scalar number or 1D array of numbers'
  IF (ISA(colat_rads, NUMBER=1) EQ 0) OR (SIZE(colat_rads, N_DIMENSIONS=1) GT 1) THEN MESSAGE,'ERROR: Second argument colat_rads must be a scalar number or 1D array of numbers'
  IF (ISA(elong_rads, NUMBER=1) EQ 0) OR (SIZE(elong_rads, N_DIMENSIONS=1) GT 1) THEN MESSAGE,'ERROR: Third  argument elong_rads must be a scalar number or 1D array of numbers'
  IF (N_input NE N_ELEMENTS(colat_rads)) THEN MESSAGE,'ERROR: First argument r_rj must be the same size as 2nd argument colat_rads'
  IF (N_input NE N_ELEMENTS(elong_rads)) THEN MESSAGE,'ERROR: First argument r_rj must be the same size as 3rd argument elong_rads'

  ; Do this check to be sure that user hasn't got position in km, must be in planetary radii.
  IF scalar_input THEN BEGIN
    IF (    r_rj   LE 0d) OR (    r_rj   GE 200d   ) THEN MESSAGE,'ERROR: First  argument, Position    r_rj   , must be in units of Rj and >0 and <200 only, and not outside that range (did you use km instead?)'
    IF (colat_rads LT 0d) OR (colat_rads GT !DPI   ) THEN MESSAGE,'ERROR: Second argument, Position colat_rads, must be in units of radians and >= 0 and <=   Pi only, and not outside that range (did you use degrees instead?)'
    IF (elong_rads LT 0d) OR (elong_rads GT !DPI*2d) THEN MESSAGE,'ERROR: Third  argument, Position elong_rads, must be in units of radians and >= 0 and <= 2*Pi only, and not outside that range (did you use degress instead?)'
  ENDIF ELSE BEGIN
    min_x = MIN(r_rj, MAX=max_x)
    IF (min_x LE 0d) OR (max_x GE 200d   ) THEN MESSAGE,'ERROR: First  argument, Position    r_rj   , must be in units of Rj and >0 and <200 only, and not outside that range (did you use km instead?)'
    min_y = MIN(colat_rads, MAX=max_y)
    IF (min_y LT 0d) OR (max_y GT !DPI   ) THEN MESSAGE,'ERROR: Second argument, Position colat_rads, must be in units of radians and >= 0 and <=   Pi only, and not outside that range (did you use degrees instead?)'
    min_z = MIN(elong_rads, MAX=max_z)
    IF (min_z LT 0d) OR (max_z GT !DPI*2d) THEN MESSAGE,'ERROR: Third  argument, Position elong_rads, must be in units of radians and >= 0 and <= 2*Pi only, and not outside that range (did you use degress instead?)'
  ENDELSE

  ; Changing inputs to Doubles, and not using input names (so as not to alter inputs, an IDL issue)
  r_rj_dbl       = DOUBLE(      r_rj  )
  colat_rads_dbl = DOUBLE(  colat_rads)
  elong_rads_dbl = DOUBLE(  elong_rads)

  ; Scaling distances since vit4_order04 expects 1Rj to be 71323 km (not the 71492 km that the inputs expect)
  r_rj_dbl = r_rj_dbl * double(71492.000000)  /  double(71323.000000)

  ; ============
  ; Begin hard-coding for VIT4_ORDER04
  ; Values from Connerney (2007), https://doi.org/10.1016/B978-044452748-6.00159-0
  ; Original paper is Connerney et al (1998) [https://doi.org/10.1029/97JA03726], however table 3 of Connerney (2007) [https://doi.org/10.1016/B978-044452748-6.00159-0] provides the g and h values to more significant figures, which are used here. i.e. 4.205 G (1998) -> 420543 nT (2007)
  ; ============

  ; order = 4  ; degree = order for this code 
  ; k     = order + 1
  k       = 5  ; order + 1 

  ; Arrays rec, g and h are processed (depending on degree) but otherwise do not
  ; change. So we calculate them once and use in the code. The initial g and h 
  ; values are given in the comments at the top of this code, and are reformatted
  ; here in to 1D arrays.
  ; g = [            0                         ,        428077.00000000000000000000000 ,        -75306.00000000000000000000000 ,         -4283.00000000000000000000000 , 
  ;             -59426.00000000000000000000000 ,         44386.00000000000000000000000 ,          8906.00000000000000000000000 ,        -21447.00000000000000000000000 , 
  ;              21130.00000000000000000000000 ,         -1190.00000000000000000000000 ,        -22925.00000000000000000000000 ,         18940.00000000000000000000000 , 
  ;              -3851.00000000000000000000000 ,          9926.00000000000000000000000 ,          1271.00000000000000000000000 ]
  ; h = [            0                         ,             0                         ,         24616.00000000000000000000000 ,             0                         , 
  ;             -50154.00000000000000000000000 ,         38452.00000000000000000000000 ,             0                         ,        -17187.00000000000000000000000 , 
  ;              40667.00000000000000000000000 ,        -35263.00000000000000000000000 ,             0                         ,         16088.00000000000000000000000 , 
  ;              11807.00000000000000000000000 ,          6195.00000000000000000000000 ,         12641.00000000000000000000000 ]
  ; These arrays are then extended and manipulated to make larger g and h arrays, and a rec array.
; ######################################################################
; The following is the Python code that was used to expand and process the
; g and h arrays, and create the rec array for pasting the numbers in to
; this source code:
;
; degree = 4      # = order
; g, h, rec = expand_out_g_and_h(degree,order,g,h)
; ----------------------------------------------------------------------
; import numpy as np
; def expand_out_g_and_h(degree,sh_order,g,h):
;
;     # Expand out g and h for later use. i.e. want length = 232 if degree is 20
;     max_gh_len = int( (degree +1)*(degree)/2+1 + degree + 1 )
;     # if g and h arrays aren't long enough, pad them to correct size with zeros
;     if (max_gh_len > len(g)):
;         g = np.append(g,np.zeros(max_gh_len - len(g),dtype='float64'))
;     if (max_gh_len > len(h)):
;         h = np.append(h,np.zeros(max_gh_len - len(h),dtype='float64'))
;
;     one_float = np.float64(1)  # = 1.0
;     two_float = np.float64(2)  # = 2.0
;     rec = np.zeros(max_gh_len,dtype='float64')
;
;     for n in range(1, degree +1 +1):
;         n2 = np.float64( 2*n-1 )
;         n2 = n2 * (n2 - two_float)
;         for m in range(1, n +1):
;             mn = int( n*(n-1)/2 + m )
;             rec[mn] = np.float64( (n-m)*(n+m-2) )/n2
;
;     s = one_float.copy() # = 1.0
;     for n in range(2, degree+1 +1):
;         mn = int( n*(n-1)/2 + 1 )
;         s = s * np.float64( 2*n - 3 )/np.float64( n - 1 )
;         p = s.copy() # = a copy of s, not a pointer to s
;         g[mn] = g[mn] * s
;         h[mn] = h[mn] * s
;         for m in range (2, n +1):
;             if (m == 2):
;                 aa = two_float.copy() # = 2.0
;             else:
;                 aa = one_float.copy() # = 1.0
;             p = p * np.sqrt( aa*np.float64( n-m+1 )/np.float64( n+m-2 ) )
;             mnn = int( mn+m-1 )
;             g[mnn] = g[mnn] * p;
;             h[mnn] = h[mnn] * p;
;
;     # In use, max index called is k*(k-1)/2 + k , where k = order + 1.
;     # so for k = 11, that's index 66, so size 67 (as indexes start at 0 in Python)
;     k = sh_order + 1
;     max_index = int( k*(k-1)/2 + k )
;     if (len(g) > max_index +1 ):  # +1 for index 0
;         g   =   g[0:(max_index +1)]
;     if (len(h) > max_index +1 ):  # +1 for index 0
;         h   =   h[0:(max_index +1)]
;     if (len(rec) > max_index +1 ):  # +1 for index 0
;         rec = rec[0:(max_index +1)]
;
;     # Done, return arrays back to main code
;     return g, h, rec
; ----------------------------------------------------------------------
; ######################################################################

  rec = [     0d                        , $
              0d                        ,             0.33333333333333331482962d,             0d                        ,             0.26666666666666666296592d, $
              0.20000000000000001110223d,             0d                        ,             0.25714285714285711748062d,             0.22857142857142856429142d, $
              0.14285714285714284921269d,             0d                        ,             0.25396825396825395415590d,             0.23809523809523808202115d, $
              0.19047619047619046561692d,             0.11111111111111110494321d,             0d                        ]

  ; This is the modified g array, not the original g coefficients, and will be further modified.
  g = [       0d                        , $
              0d                        ,        428077.00000000000000000000000d,        -75306.00000000000000000000000d,         -6424.50000000000000000000000d, $
        -102928.85129058809252455830574d,         38439.40357237608986906707287d,         22265.00000000000000000000000d,        -65667.75814183852344285696745d, $
          40918.06905268136324593797326d,          -940.77760390009291313617723d,       -100296.87500000000000000000000d,        104813.69304628092504572123289d, $
         -15069.42111736545848543755710d,         20761.71855844308447558432817d,           939.91717553995147227396956d]

  ; This is the modified h array, not the original h coefficients, and will be further modified.
  h = [       0d                        , $
              0d                        ,             0d                        ,         24616.00000000000000000000000d,             0d                        , $
         -86869.27620280947303399443626d,         33300.40882631923159351572394d,             0d                        ,        -52624.22526151809870498254895d, $
          78751.30687010851397644728422d,        -27877.84928262939138221554458d,             0d                        ,         89030.76524438054184429347515d, $
          46202.19556809503410477191210d,         12957.77216094649520528037101d,          9348.14556727028138993773609d]

  ; ============
  ; End parts that are hard-coded for VIT4_ORDER04
  ; ============

  IF scalar_input THEN BEGIN
    a         = [ 0d, 0d, 0d, 0d, 0d, 0d]  ; = DBLARR(k+1)
    DINDGEN_k = [ 0d, 1d, 2d, 3d, 4d, 5d]  ; = DINDGEN(k+1), done manually for speed
  ENDIF ELSE BEGIN
    a         = DBLARR(N_input,k+1)
    DINDGEN_k = a
    FOR i = 0,k DO DINDGEN_k[*,i] = i
  ENDELSE

  da = 1d/r_rj_dbl
  IF scalar_input THEN BEGIN
    a[0] = da
    FOR i=1,k DO a[i] = a[i-1]*da
  ENDIF ELSE BEGIN ; the following vectorized 2 lines works for scalars, but is slower.
    a[*,0] = da
    FOR i=1,k DO a[*,i] = a[*,i-1]*da
  ENDELSE

  b = a  * DINDGEN_k

  cos_phi   = cos(elong_rads_dbl)
  sin_phi   = sin(elong_rads_dbl)
  cos_theta = cos(colat_rads_dbl)
  sin_theta = sin(colat_rads_dbl)
  not_bk = (sin_theta GE 0.00001d)  ; = 1d-5 - also see bk both times below
  IF scalar_input THEN BEGIN
    ; bk = (sin_theta LT 0.00001d)  ; bk not needed for scalar
    zero_array = 0d
    p   = 1d
    d   = 0d
    bbr = 0d
    bbt = 0d
    bbf = 0d
    x = 0d
    y = 1d
  ENDIF ELSE BEGIN
    bk = (sin_theta LT 0.00001d)
    zero_array = DBLARR(N_input)
    p   = zero_array + 1d
    d   = zero_array
    bbr = zero_array
    bbt = zero_array
    bbf = zero_array
    x = zero_array
    y = p  ; 1s
  ENDELSE

  FOR m = 1,k DO BEGIN
    bm  = (m NE 1)
    IF bm THEN BEGIN
      m_minus_1 = DOUBLE(m - 1)
      w = x
      x = w *cos_phi + y *sin_phi
      y = y *cos_phi - w *sin_phi
    ENDIF
    q = p
    z = d
    bi = zero_array
    p2 = zero_array
    d2 = zero_array
    FOR n = m,k DO BEGIN
      mn = n*(n-1)/2 + m
      w  = g[mn]*y + h[mn]*x
      IF scalar_input THEN BEGIN
        bbr += b[  n]*w*q
        bbt -= a[  n]*w*z
        IF bm THEN BEGIN
          IF not_bk THEN bi += a[n] * (g[mn]*x-h[mn]*y) * q  $
          ELSE           bi += a[n] * (g[mn]*x-h[mn]*y) * z
        ENDIF
      ENDIF ELSE BEGIN
        bbr += b[*,n]*w*q
        bbt -= a[*,n]*w*z
        IF bm THEN BEGIN
          qq = q
          ind = WHERE(bk, n_ind, NULL = 1)
          IF (n_ind NE 0) THEN qq[ind] = z[ind]
          bi += a[*,n] * (g[mn]*x-h[mn]*y) * qq
        ENDIF
      END
      xk = rec[mn] ; in IDL 8.4 it's faster to write this to xk, to use below twice.  In Matlab quicker to just use rec(nm)
      dp = cos_theta *z - sin_theta *q - d2*xk
      pm = cos_theta *q                - p2*xk
      d2 = z
      p2 = q
      z = dp
      q = pm
    ENDFOR
    d = sin_theta *d + cos_theta *p
    p = sin_theta *p
    IF bm THEN BEGIN
      bi  *= m_minus_1
      bbf += bi
    ENDIF
  ENDFOR

  ; br = bbr  ; This doesn't change again
  ; bt = bbt  ; This doesn't change again
  IF scalar_input THEN BEGIN
    IF not_bk THEN bf = bbf/sin_theta $
    ELSE BEGIN
      IF (cos_theta GE 0d) THEN bf =  bbf ELSE bf = -bbf
    ENDELSE
  ENDIF ELSE BEGIN
    bf = bbf ; set size of array and do the 3rd case
    ind = WHERE(bk AND (cos_theta LT 0d), NULL = 1)
    IF (N_ELEMENTS(ind) NE 0) THEN bf[ind] = -bbf[ind]
    ind = WHERE(bk EQ 0, NULL = 1)
    IF (N_ELEMENTS(ind) NE 0) THEN bf[ind] =  bbf[ind]/sin_theta[ind]
  ENDELSE

  RETURN,[[bbr],[bbt],[bf]]
END
