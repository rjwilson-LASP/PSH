FUNCTION jovian_jrm09_order10_internal_rtp, r_rj, colat_rads, elong_rads
  ; Code to calculate the JRM09_ORDER10 model of Jupiter's internal magnetic field model
  ; with Degree 10 and Order 10.
  ; Reference: Connerney et al. (2018), https://doi.org/10.1002/2018GL077312
  ;
  ; Required inputs (System III (1965) Spherical, right handed, and assuming 1 Rj = 71492 km):
  ;  r_rj       - radial distance, in Rj.
  ;  colat_rads - colatitude, in radians.                    Value(s) should be 0 <= colat_rads <=  pi.
  ;  elong_rads - East longitude, right handed, in radians.  Value(s) should be 0 <= elong_rads <= 2pi.
  ;
  ; Outputs:
  ;  B - Spherical Magnetic field vector the JRM09_ORDER10 internal magnetic field model, [Br, Btheta, Bphi], units of nT.
  ;
  ; Usage:
  ; For internal field only: B = jovian_jrm09_order10_internal_rtp(r_rj, colat_rads, elong_rads)
  ;
  ; This code was written by Marissa Vogt (mvogt@bu.edu) and Rob Wilson (rob.wilson@lasp.colorado.edu).
  ; It is based on a routine originally written by K. Khurana, translated into IDL by Marissa Vogt in 2009.
  ; Thanks to Masafumi Imai for providing code for his version of the JRM09 model, which was used to test and validate this code.
  ;
  ; Citation Info:
  ;  DOI: 10.5281/zenodo.6814109     This DOI links to all versions of code at the Github.
  ;  Github: https://github.com/rjwilson-LASP/PSH
  ;  Individual versions released on the Github repository can have a different DOI,
  ;  See the DOI above for a list of DOIs for each specific Github released version.
  ;
  ; Version Info:
  ;  Last update of this file: 2022-07-09 11:51:05.644133 by user wilsonr. 
  ;  This code was re-written/re-formatted by the Mother_Source python code:
  ;   /Users/wilsonr/Documents/JADE/Level2_Processing_Code/IDL/Field_Model/2022/Git_initial/Mother_Source/MOP_spherical.py
  ;   which itself was last updated at UTC 2022-07-09T17:50:57.
  ;
  ;  The Spherical Harmonic g and h values used for this order 10 code are below: 
  ;  
  ;  g[i,j] values (nT) used are:
  ; g[ 1, 0] =   410244.7, g[ 1, 1] =   -71498.3, 
  ; g[ 2, 0] =    11670.4, g[ 2, 1] =   -56835.8, g[ 2, 2] =    48689.5, 
  ; g[ 3, 0] =     4018.6, g[ 3, 1] =   -37791.1, g[ 3, 2] =    15926.3, g[ 3, 3] =    -2710.5, 
  ; g[ 4, 0] =   -34645.4, g[ 4, 1] =    -8247.6, g[ 4, 2] =    -2406.1, g[ 4, 3] =   -11083.8, g[ 4, 4] =   -17837.2, 
  ; g[ 5, 0] =   -18023.6, g[ 5, 1] =     4683.9, g[ 5, 2] =    16160.0, g[ 5, 3] =   -16402.0, g[ 5, 4] =    -2600.7, g[ 5, 5] =    -3660.7, 
  ; g[ 6, 0] =   -20819.6, g[ 6, 1] =     9992.9, g[ 6, 2] =    11791.8, g[ 6, 3] =   -12574.7, g[ 6, 4] =     2669.7, g[ 6, 5] =     1113.2, g[ 6, 6] =     7584.9, 
  ; g[ 7, 0] =      598.4, g[ 7, 1] =     4665.9, g[ 7, 2] =    -6495.7, g[ 7, 3] =    -2516.5, g[ 7, 4] =    -6448.5, g[ 7, 5] =     1855.3, g[ 7, 6] =    -2892.9, g[ 7, 7] =     2968.0, 
  ; g[ 8, 0] =    10059.2, g[ 8, 1] =     1934.4, g[ 8, 2] =    -6702.9, g[ 8, 3] =      153.7, g[ 8, 4] =    -4124.2, g[ 8, 5] =     -867.2, g[ 8, 6] =    -3740.6, g[ 8, 7] =     -732.4, g[ 8, 8] =    -2433.2, 
  ; g[ 9, 0] =     9671.8, g[ 9, 1] =    -3046.2, g[ 9, 2] =      260.9, g[ 9, 3] =     2071.3, g[ 9, 4] =     3329.6, g[ 9, 5] =    -2523.1, g[ 9, 6] =     1787.1, g[ 9, 7] =    -1148.2, g[ 9, 8] =     1276.5, g[ 9, 9] =    -1976.8, 
  ; g[10, 0] =    -2299.5, g[10, 1] =     2009.7, g[10, 2] =     2127.8, g[10, 3] =     3498.3, g[10, 4] =     2967.6, g[10, 5] =       16.3, g[10, 6] =     1806.5, g[10, 7] =      -46.5, g[10, 8] =     2897.8, g[10, 9] =      574.5, g[10,10] =     1298.9, 
  ;
  ;  h[i,j] values (nT) used are:
  ;                        h[ 1, 1] =    21330.5, 
  ;                        h[ 2, 1] =   -42027.3, h[ 2, 2] =    19353.2, 
  ;                        h[ 3, 1] =   -32957.3, h[ 3, 2] =    42084.5, h[ 3, 3] =   -27544.2, 
  ;                        h[ 4, 1] =    31994.5, h[ 4, 2] =    27811.2, h[ 4, 3] =     -926.1, h[ 4, 4] =      367.1, 
  ;                        h[ 5, 1] =    45347.9, h[ 5, 2] =     -749.0, h[ 5, 3] =     6268.5, h[ 5, 4] =    10859.6, h[ 5, 5] =     9608.4, 
  ;                        h[ 6, 1] =    14533.1, h[ 6, 2] =   -10592.9, h[ 6, 3] =      568.6, h[ 6, 4] =    12871.7, h[ 6, 5] =    -4147.8, h[ 6, 6] =     3604.4, 
  ;                        h[ 7, 1] =    -7626.3, h[ 7, 2] =   -10948.4, h[ 7, 3] =     2633.3, h[ 7, 4] =     5394.2, h[ 7, 5] =    -6050.8, h[ 7, 6] =    -1526.0, h[ 7, 7] =    -5684.2, 
  ;                        h[ 8, 1] =    -2409.7, h[ 8, 2] =   -11614.6, h[ 8, 3] =     9287.0, h[ 8, 4] =     -911.9, h[ 8, 5] =     2754.5, h[ 8, 6] =    -2446.1, h[ 8, 7] =     1207.3, h[ 8, 8] =    -2887.3, 
  ;                        h[ 9, 1] =    -8467.4, h[ 9, 2] =    -1383.8, h[ 9, 3] =     5697.7, h[ 9, 4] =    -2056.3, h[ 9, 5] =     3081.5, h[ 9, 6] =     -721.2, h[ 9, 7] =     1352.5, h[ 9, 8] =     -210.1, h[ 9, 9] =     1567.6, 
  ;                        h[10, 1] =    -4692.6, h[10, 2] =     4445.8, h[10, 3] =    -2378.6, h[10, 4] =    -2204.3, h[10, 5] =      164.1, h[10, 6] =    -1361.6, h[10, 7] =    -2031.5, h[10, 8] =     1411.8, h[10, 9] =     -714.3, h[10,10] =     1676.5, 

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

  ; ============
  ; Begin hard-coding for JRM09_ORDER10
  ; Values from Connerney et al. (2018), https://doi.org/10.1002/2018GL077312
  ; See supplemental online information Table S1, https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2F2018GL077312&file=grl57087-sup-0005-2018GL077312-ds01.txt
  ; ============

  ; order = 10  ; degree = order for this code 
  ; k     = order + 1
  k       = 11  ; order + 1 

  ; Arrays rec, g and h are processed (depending on degree) but otherwise do not
  ; change. So we calculate them once and use in the code. The initial g and h 
  ; values are given in the comments at the top of this code, and are reformatted
  ; here in to 1D arrays.
  ; g = [       0   ,   410244.7 ,   -71498.3 ,    11670.4 , 
  ;        -56835.8 ,    48689.5 ,     4018.6 ,   -37791.1 , 
  ;         15926.3 ,    -2710.5 ,   -34645.4 ,    -8247.6 , 
  ;         -2406.1 ,   -11083.8 ,   -17837.2 ,   -18023.6 , 
  ;          4683.9 ,    16160.0 ,   -16402.0 ,    -2600.7 , 
  ;         -3660.7 ,   -20819.6 ,     9992.9 ,    11791.8 , 
  ;        -12574.7 ,     2669.7 ,     1113.2 ,     7584.9 , 
  ;           598.4 ,     4665.9 ,    -6495.7 ,    -2516.5 , 
  ;         -6448.5 ,     1855.3 ,    -2892.9 ,     2968.0 , 
  ;         10059.2 ,     1934.4 ,    -6702.9 ,      153.7 , 
  ;         -4124.2 ,     -867.2 ,    -3740.6 ,     -732.4 , 
  ;         -2433.2 ,     9671.8 ,    -3046.2 ,      260.9 , 
  ;          2071.3 ,     3329.6 ,    -2523.1 ,     1787.1 , 
  ;         -1148.2 ,     1276.5 ,    -1976.8 ,    -2299.5 , 
  ;          2009.7 ,     2127.8 ,     3498.3 ,     2967.6 , 
  ;            16.3 ,     1806.5 ,      -46.5 ,     2897.8 , 
  ;           574.5 ,     1298.9 ]
  ; h = [       0   ,        0   ,    21330.5 ,        0   , 
  ;        -42027.3 ,    19353.2 ,        0   ,   -32957.3 , 
  ;         42084.5 ,   -27544.2 ,        0   ,    31994.5 , 
  ;         27811.2 ,     -926.1 ,      367.1 ,        0   , 
  ;         45347.9 ,     -749.0 ,     6268.5 ,    10859.6 , 
  ;          9608.4 ,        0   ,    14533.1 ,   -10592.9 , 
  ;           568.6 ,    12871.7 ,    -4147.8 ,     3604.4 , 
  ;             0   ,    -7626.3 ,   -10948.4 ,     2633.3 , 
  ;          5394.2 ,    -6050.8 ,    -1526.0 ,    -5684.2 , 
  ;             0   ,    -2409.7 ,   -11614.6 ,     9287.0 , 
  ;          -911.9 ,     2754.5 ,    -2446.1 ,     1207.3 , 
  ;         -2887.3 ,        0   ,    -8467.4 ,    -1383.8 , 
  ;          5697.7 ,    -2056.3 ,     3081.5 ,     -721.2 , 
  ;          1352.5 ,     -210.1 ,     1567.6 ,        0   , 
  ;         -4692.6 ,     4445.8 ,    -2378.6 ,    -2204.3 , 
  ;           164.1 ,    -1361.6 ,    -2031.5 ,     1411.8 , 
  ;          -714.3 ,     1676.5 ]
  ; These arrays are then extended and manipulated to make larger g and h arrays, and a rec array.
; ######################################################################
; The following is the Python code that was used to expand and process the
; g and h arrays, and create the rec array for pasting the numbers in to
; this source code:
;
; degree = 10      # = order
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
              0d                        ,             0.25174825174825177231952d,             0.24475524475524476630817d,             0.22377622377622377602968d, $
              0.18881118881118880148406d,             0.13986013986013987042689d,             0.07692307692307692734701d,             0d                        , $
              0.25128205128205127749652d,             0.24615384615384616751044d,             0.23076923076923078204103d,             0.20512820512820512108831d, $
              0.16923076923076924016343d,             0.12307692307692308375522d,             0.06666666666666666574148d,             0d                        , $
              0.25098039215686274161499d,             0.24705882352941177515504d,             0.23529411764705882026405d,             0.21568627450980393245317d, $
              0.18823529411764705621124d,             0.15294117647058824704942d,             0.10980392156862744945656d,             0.05882352941176470506601d, $
              0d                        ,             0.25077399380804954454049d,             0.24767801857585139413409d,             0.23839009287925697067045d, $
              0.22291021671826624639401d,             0.20123839009287924906033d,             0.17337461300309597866942d,             0.13931888544891640746570d, $
              0.09907120743034056320475d,             0.05263157894736841813099d,             0d                        ,             0.25062656641604008633806d, $
              0.24812030075187968547468d,             0.24060150375939848288454d,             0.22807017543859647856763d,             0.21052631578947367252397d, $
              0.18796992481203006475354d,             0.16040100250626565525636d,             0.12781954887218044403241d,             0.09022556390977443108170d, $
              0.04761904761904761640423d,             0d                        ]

  ; This is the modified g array, not the original g coefficients, and will be further modified.
  g = [       0d                        , $
              0d                        ,        410244.70000000001164153218269d,        -71498.30000000000291038304567d,         17505.59999999999854480847716d, $
         -98442.49328882319969125092030d,         42166.34389756242308067157865d,         10046.50000000000000000000000d,       -115711.13977311669441405683756d, $
          30841.14733335159326088614762d,         -2142.83839947159822258981876d,       -151573.62500000000000000000000d,        -45642.10215250826877309009433d, $
          -9415.35553115892798814456910d,        -23183.43100524596593459136784d,        -13190.78728838805909617803991d,       -141935.84999999997671693563461d, $
          47619.25007516491314163431525d,        124193.04328343032102566212416d,        -77191.29987306696421001106501d,         -5769.73075946518656564876437d, $
          -2568.20347420563030027551576d,       -300582.97499999997671693563461d,        188897.03523126235813833773136d,        176219.39807558275060728192329d, $
        -125279.49167958697944413870573d,         14568.18469249509325891267508d,          2590.20913175944133399752900d,          5094.72643062895076582208276d, $
          16044.59999999999854480847716d,        165497.62303578440332785248756d,       -188120.73349108296679332852364d,        -51533.85619684254197636619210d, $
         -79632.08074599411338567733765d,         11455.48572598611099238041788d,         -7006.09637449074671167181805d,          1921.06723268604082477395423d, $
         505710.56250000005820766091347d,        129665.25000000000000000000000d,       -375914.50046967255184426903725d,          6366.18842053944717918056995d, $
        -110265.51694787663291208446026d,        -12861.08441189933364512398839d,        -25680.06972210412277490831912d,         -1835.99981426163253672712017d, $
          -1524.90263109687475662212819d,        918443.19531249988358467817307d,       -388096.44079238711856305599213d,         28346.79585738653622684068978d, $
         171882.35647868152591399848461d,        187708.65848237561294808983803d,        -85005.63010931370081380009651d,         31091.84173000367445638403296d, $
          -8649.99331942896606051363051d,          3298.44757452672547515248880d,         -1203.96883845257525535998866d,       -414889.08398437500000000000000d, $
         488932.02253022126387804746628d,        448310.26802277535898610949516d,        578200.21502919984050095081329d,        346825.93308208207599818706512d, $
           1204.82452975374712877965067d,         74644.73375998696428723633289d,          -932.00811829940028019336751d,         23711.52712466476441477425396d, $
           1525.17384009177953885227907d,           771.06330156869501024630154d]

  ; This is the modified h array, not the original h coefficients, and will be further modified.
  h = [       0d                        , $
              0d                        ,             0d                        ,         21330.50000000000000000000000d,             0d                        , $
         -72793.41890493947721552103758d,         16760.36284452099789632484317d,             0d                        ,       -100910.71037478504877071827650d, $
          81496.28381673301919363439083d,        -21775.60208180246991105377674d,             0d                        ,        177057.11204695011838339269161d, $
         108828.53403772377350833266973d,         -1937.07712643301852040167432d,           271.47411104698363715215237d,             0d                        , $
         461033.11139938322594389319420d,         -5756.22459277780399133916944d,         29500.89399184978901757858694d,         24092.34750470571452751755714d, $
           6740.87640657726115023251623d,             0d                        ,        274721.00218349619535729289055d,       -158302.75800766979227773845196d, $
           5664.86031229477885062806308d,         70239.09162317455047741532326d,         -9651.15831540766339458059520d,          2421.05129224630400130990893d, $
              0d                        ,       -270501.83727851061848923563957d,       -317074.53216031729243695735931d,         53925.73158082475129049271345d, $
          66612.60292471760476473718882d,        -37360.45546854781423462554812d,         -3695.70433387703633343335241d,         -3679.15443532142626281711273d, $
              0d                        ,       -161525.20312500000000000000000d,       -651374.26444599486421793699265d,        384663.57749869779217988252640d, $
         -24380.75866950407362310215831d,         40850.84987612628174247220159d,        -16793.03281485293700825423002d,          3026.49177465601997027988546d, $
          -1809.49012278727877855999395d,             0d                        ,      -1078776.11541115446016192436218d,       -150349.92758701223647221922874d, $
         472811.32743136369390413165092d,       -115925.43081370404979679733515d,        103818.65529778851487208157778d,        -12547.38753045641169592272490d, $
          10189.09246170325241109821945d,          -542.89372143209163823485142d,           954.74582717435089307400631d,             0d                        , $
       -1141644.22994741331785917282104d,        936694.13928736466914415359497d,       -393135.81781678373226895928383d,       -257618.41363149805692955851555d, $
          12129.55247439201775705441833d,        -56261.42789238762634340673685d,        -40717.73101774691895116120577d,         11552.18924515208345837891102d, $
          -1896.31274843787309691833798d,           995.21720307946509365137899d]

  ; ============
  ; End parts that are hard-coded for JRM09_ORDER10
  ; ============

  IF scalar_input THEN BEGIN
    a         = [ 0d, 0d, 0d, 0d, 0d, 0d, 0d, 0d, 0d, 0d, 0d, 0d]  ; = DBLARR(k+1)
    DINDGEN_k = [ 0d, 1d, 2d, 3d, 4d, 5d, 6d, 7d, 8d, 9d,10d,11d]  ; = DINDGEN(k+1), done manually for speed
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
