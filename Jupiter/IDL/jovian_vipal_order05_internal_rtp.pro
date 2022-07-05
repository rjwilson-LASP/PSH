FUNCTION jovian_vipal_order05_internal_rtp, r_rj, colat_rads, elong_rads
  ; Code to calculate the VIPAL_ORDER05 model of Jupiter's internal magnetic field model
  ; with Degree 5 and Order 5.
  ; Reference: Hess et al. (2011), https://doi.org/10.1029/2010JA016262
  ;
  ; Required inputs (System III (1965) Spherical, right handed, and assuming 1 Rj = 71492 km):
  ;  r_rj       - radial distance, in Rj.
  ;  colat_rads - colatitude, in radians.                    Value(s) should be 0 <= colat_rads <=  pi.
  ;  elong_rads - East longitude, right handed, in radians.  Value(s) should be 0 <= elong_rads <= 2pi.
  ;
  ; Outputs:
  ;  B - Spherical Magnetic field vector the VIPAL_ORDER05 internal magnetic field model, [Br, Btheta, Bphi], units of nT.
  ;
  ; Usage:
  ; For internal field only: B = jovian_vipal_order05_internal_rtp(r_rj, colat_rads, elong_rads)
  ;
  ; This code was written by Marissa Vogt (mvogt@bu.edu) and Rob Wilson (rob.wilson@lasp.colorado.edu).
  ; It is based on a routine originally written by K. Khurana, translated into IDL by Marissa Vogt in 2009.
  ;
  ; Version Info:
  ;  Last update of this file: 2022-07-05 14:33:42.022131 by user wilsonr. 
  ;  This code was re-written/re-formatted by Rob's python code:
  ;   /Volumes/wilsonr/Documents/JADE/Level2_Processing_Code/IDL/Field_Model/2022/Git_initial/Mother_Source/MOP_spherical.py
  ;   which itself was last updated at UTC 2022-07-05T20:33:19.
  ;
  ;  The Spherical Harmonic g and h values used for this order 5 code are below: 
  ;  
  ;  g[i,j] values (nT) used are:
  ; g[ 1, 0] =     420000, g[ 1, 1] =     -69750, 
  ; g[ 2, 0] =      64410, g[ 2, 1] =     -86720, g[ 2, 2] =      95980, 
  ; g[ 3, 0] =     -10580, g[ 3, 1] =     -59000, g[ 3, 2] =      63220, g[ 3, 3] =      46710, 
  ; g[ 4, 0] =     -74660, g[ 4, 1] =      32820, g[ 4, 2] =     -33800, g[ 4, 3] =      18260, g[ 4, 4] =     -14290, 
  ; g[ 5, 0] =      -6600, g[ 5, 1] =       7370, g[ 5, 2] =     -17110, g[ 5, 3] =     -17930, g[ 5, 4] =       -770, g[ 5, 5] =      -7400, 
  ;
  ;  h[i,j] values (nT) used are:
  ;                        h[ 1, 1] =      19730, 
  ;                        h[ 2, 1] =     -40410, h[ 2, 2] =      60300, 
  ;                        h[ 3, 1] =     -23100, h[ 3, 2] =      51600, h[ 3, 3] =     -11310, 
  ;                        h[ 4, 1] =      32830, h[ 4, 2] =     -21310, h[ 4, 3] =      -6060, h[ 4, 4] =      -4860, 
  ;                        h[ 5, 1] =      20650, h[ 5, 2] =     -11670, h[ 5, 3] =      -2880, h[ 5, 4] =       -500, h[ 5, 5] =     -22790, 

  ON_ERROR, 2 ; Exit code if an error in main, don't stop in code - no Matlab equivalent, just delete line in Matlab

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

  ; ============
  ; Begin hard-coding for VIPAL_ORDER05
  ; Values from Hess et al. (2011), https://doi.org/10.1029/2010JA016262
  ; Values are give to whole nT here, but in the original paper were given to 1 to 4 decimal places in units of G.
  ; ============

  ; order = 5  ; degree = order for this code 
  ; k     = order + 1
  k       = 6  ; order + 1 

  ; Arrays rec, g and h are processed (depending on degree) but otherwise do not
  ; change. So we calculate them once and use in the code. The initial g and h 
  ; values are given in the comments at the top of this code, and are reformatted
  ; here in to 1D arrays.
  ; g = [            0                         ,        420000.00000000000000000000000 ,        -69750.00000000000000000000000 ,         64410.00000000000000000000000 , 
  ;             -86720.00000000000000000000000 ,         95980.00000000000000000000000 ,        -10580.00000000000000000000000 ,        -59000.00000000000000000000000 , 
  ;              63220.00000000000000000000000 ,         46710.00000000000000000000000 ,        -74660.00000000000000000000000 ,         32820.00000000000000000000000 , 
  ;             -33800.00000000000000000000000 ,         18260.00000000000000000000000 ,        -14290.00000000000000000000000 ,         -6600.00000000000000000000000 , 
  ;               7370.00000000000000000000000 ,        -17110.00000000000000000000000 ,        -17930.00000000000000000000000 ,          -770.00000000000000000000000 , 
  ;              -7400.00000000000000000000000 ]
  ; h = [            0                         ,             0                         ,         19730.00000000000000000000000 ,             0                         , 
  ;             -40410.00000000000000000000000 ,         60300.00000000000000000000000 ,             0                         ,        -23100.00000000000000000000000 , 
  ;              51600.00000000000000000000000 ,        -11310.00000000000000000000000 ,             0                         ,         32830.00000000000000000000000 , 
  ;             -21310.00000000000000000000000 ,         -6060.00000000000000000000000 ,         -4860.00000000000000000000000 ,             0                         , 
  ;              20650.00000000000000000000000 ,        -11670.00000000000000000000000 ,         -2880.00000000000000000000000 ,          -500.00000000000000000000000 , 
  ;             -22790.00000000000000000000000 ]
  ; These arrays are then extended and manipulated to make larger g and h arrays, and a rec array.
; ######################################################################
; The following is the Python code that was used to expand and process the
; g and h arrays, and create the rec array for pasting the numbers in to
; this source code:
;
; degree = 5      # = order
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
              0.19047619047619046561692d,             0.11111111111111110494321d,             0d                        ,             0.25252525252525254151337d, $
              0.24242424242424243097105d,             0.21212121212121212709967d,             0.16161616161616162989922d,             0.09090909090909091161414d, $
              0d                        ]

  ; This is the modified g array, not the original g coefficients, and will be further modified.
  g = [       0d                        , $
              0d                        ,        420000.00000000000000000000000d,        -69750.00000000000000000000000d,         96615.00000000000000000000000d, $
        -150203.44603237303090281784534d,         83121.11825523042352870106697d,        -26450.00000000000000000000000d,       -180649.86853025938034988939762d, $
         122425.00357361645728815346956d,         36927.49737661625113105401397d,       -326637.50000000000000000000000d,        181625.41741177084622904658318d, $
        -132263.42086911256774328649044d,         38193.53021128055115696042776d,        -10567.59751256168965483084321d,        -51975.00000000000000000000000d, $
          74927.70406156523677054792643d,       -131493.99570417654467746615410d,        -84382.39280112733831629157066d,         -1708.26803737001364424941130d, $
          -5191.54962414884175814222544d]

  ; This is the modified h array, not the original h coefficients, and will be further modified.
  h = [       0d                        , $
              0d                        ,             0d                        ,         19730.00000000000000000000000d,             0d                        , $
         -69992.17313385833404026925564d,         52221.33184820164751727133989d,             0d                        ,        -70729.01632286426320206373930d, $
          99922.97033215136616490781307d,         -8941.34008412609364313539118d,             0d                        ,        181680.75727082381490617990494d, $
         -83388.56505091091094072908163d,        -12675.39940199124430364463478d,         -3594.01846823301684707985260d,             0d                        , $
         209939.90351035579806193709373d,        -89686.43657906138105317950249d,        -13553.89242985202145064249635d,         -1109.26495933117757886066101d, $
         -15988.56972085839333885814995d]

  ; ============
  ; End parts that are hard-coded for VIPAL_ORDER05
  ; ============

  IF scalar_input THEN BEGIN
    a         = [ 0d, 0d, 0d, 0d, 0d, 0d, 0d]  ; = DBLARR(k+1)
    DINDGEN_k = [ 0d, 1d, 2d, 3d, 4d, 5d, 6d]  ; = DINDGEN(k+1), done manually for speed
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
