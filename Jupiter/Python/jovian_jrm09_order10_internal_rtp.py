import numpy as np

def jovian_jrm09_order10_internal_rtp( r_rj, colat_rads, elong_rads):
    # Code to calculate the JRM09_ORDER10 model of Jupiter's internal magnetic field model
    # with Degree 10 and Order 10.
    # Reference: Connerney et al. (2018), https://doi.org/10.1002/2018GL077312
    #
    # Required inputs (System III (1965) Spherical, right handed, and assuming 1 Rj = 71492 km):
    #  r_rj       - radial distance, in Rj.
    #  colat_rads - colatitude, in radians.                    Value(s) should be 0 <= colat_rads <=  pi.
    #  elong_rads - East longitude, right handed, in radians.  Value(s) should be 0 <= elong_rads <= 2pi.
    #
    # Outputs:
    #  B - Spherical Magnetic field vector the JRM09_ORDER10 internal magnetic field model, [Br, Btheta, Bphi], units of nT.
    #
    # Usage:
    # For internal field only: B = jovian_jrm09_order10_internal_rtp(r_rj, colat_rads, elong_rads)
    #
    # This code was written by Marissa Vogt (mvogt@bu.edu) and Rob Wilson (rob.wilson@lasp.colorado.edu).
    # It is based on a routine originally written by K. Khurana, translated into IDL by Marissa Vogt in 2009.
    # Thanks to Masafumi Imai for providing code for his version of the JRM09 model, which was used to test and validate this code.
    #
    # Citation Info:
    #  DOI: 10.5281/zenodo.6814109     This DOI links to all versions of code at the Github.
    #  Github: https://github.com/rjwilson-LASP/PSH
    #  Individual versions released on the Github repository can have a different DOI,
    #  See the DOI above for a list of DOIs for each specific Github released version.
    #
    # Version Info:
    #  Last update of this file: 2024-03-13 11:48:53.848658 by user wilsonr. 
    #  This code was re-written/re-formatted by the Mother_Source python code:
    #   /Users/wilsonr/Documents/Publications_Presentations/GitHub_Desktop/PSH/Mother_Source/MOP_spherical.py
    #   which itself was last updated at UTC 2024-03-13T17:48:37.
    #
    #  The Spherical Harmonic g and h values used for this order 10 code are below: 
    #  
    #  g[i,j] values (nT) used are:
    # g[ 1, 0] =   410244.7, g[ 1, 1] =   -71498.3, 
    # g[ 2, 0] =    11670.4, g[ 2, 1] =   -56835.8, g[ 2, 2] =    48689.5, 
    # g[ 3, 0] =     4018.6, g[ 3, 1] =   -37791.1, g[ 3, 2] =    15926.3, g[ 3, 3] =    -2710.5, 
    # g[ 4, 0] =   -34645.4, g[ 4, 1] =    -8247.6, g[ 4, 2] =    -2406.1, g[ 4, 3] =   -11083.8, g[ 4, 4] =   -17837.2, 
    # g[ 5, 0] =   -18023.6, g[ 5, 1] =     4683.9, g[ 5, 2] =    16160.0, g[ 5, 3] =   -16402.0, g[ 5, 4] =    -2600.7, g[ 5, 5] =    -3660.7, 
    # g[ 6, 0] =   -20819.6, g[ 6, 1] =     9992.9, g[ 6, 2] =    11791.8, g[ 6, 3] =   -12574.7, g[ 6, 4] =     2669.7, g[ 6, 5] =     1113.2, g[ 6, 6] =     7584.9, 
    # g[ 7, 0] =      598.4, g[ 7, 1] =     4665.9, g[ 7, 2] =    -6495.7, g[ 7, 3] =    -2516.5, g[ 7, 4] =    -6448.5, g[ 7, 5] =     1855.3, g[ 7, 6] =    -2892.9, g[ 7, 7] =     2968.0, 
    # g[ 8, 0] =    10059.2, g[ 8, 1] =     1934.4, g[ 8, 2] =    -6702.9, g[ 8, 3] =      153.7, g[ 8, 4] =    -4124.2, g[ 8, 5] =     -867.2, g[ 8, 6] =    -3740.6, g[ 8, 7] =     -732.4, g[ 8, 8] =    -2433.2, 
    # g[ 9, 0] =     9671.8, g[ 9, 1] =    -3046.2, g[ 9, 2] =      260.9, g[ 9, 3] =     2071.3, g[ 9, 4] =     3329.6, g[ 9, 5] =    -2523.1, g[ 9, 6] =     1787.1, g[ 9, 7] =    -1148.2, g[ 9, 8] =     1276.5, g[ 9, 9] =    -1976.8, 
    # g[10, 0] =    -2299.5, g[10, 1] =     2009.7, g[10, 2] =     2127.8, g[10, 3] =     3498.3, g[10, 4] =     2967.6, g[10, 5] =       16.3, g[10, 6] =     1806.5, g[10, 7] =      -46.5, g[10, 8] =     2897.8, g[10, 9] =      574.5, g[10,10] =     1298.9, 
    #
    #  h[i,j] values (nT) used are:
    #                        h[ 1, 1] =    21330.5, 
    #                        h[ 2, 1] =   -42027.3, h[ 2, 2] =    19353.2, 
    #                        h[ 3, 1] =   -32957.3, h[ 3, 2] =    42084.5, h[ 3, 3] =   -27544.2, 
    #                        h[ 4, 1] =    31994.5, h[ 4, 2] =    27811.2, h[ 4, 3] =     -926.1, h[ 4, 4] =      367.1, 
    #                        h[ 5, 1] =    45347.9, h[ 5, 2] =     -749.0, h[ 5, 3] =     6268.5, h[ 5, 4] =    10859.6, h[ 5, 5] =     9608.4, 
    #                        h[ 6, 1] =    14533.1, h[ 6, 2] =   -10592.9, h[ 6, 3] =      568.6, h[ 6, 4] =    12871.7, h[ 6, 5] =    -4147.8, h[ 6, 6] =     3604.4, 
    #                        h[ 7, 1] =    -7626.3, h[ 7, 2] =   -10948.4, h[ 7, 3] =     2633.3, h[ 7, 4] =     5394.2, h[ 7, 5] =    -6050.8, h[ 7, 6] =    -1526.0, h[ 7, 7] =    -5684.2, 
    #                        h[ 8, 1] =    -2409.7, h[ 8, 2] =   -11614.6, h[ 8, 3] =     9287.0, h[ 8, 4] =     -911.9, h[ 8, 5] =     2754.5, h[ 8, 6] =    -2446.1, h[ 8, 7] =     1207.3, h[ 8, 8] =    -2887.3, 
    #                        h[ 9, 1] =    -8467.4, h[ 9, 2] =    -1383.8, h[ 9, 3] =     5697.7, h[ 9, 4] =    -2056.3, h[ 9, 5] =     3081.5, h[ 9, 6] =     -721.2, h[ 9, 7] =     1352.5, h[ 9, 8] =     -210.1, h[ 9, 9] =     1567.6, 
    #                        h[10, 1] =    -4692.6, h[10, 2] =     4445.8, h[10, 3] =    -2378.6, h[10, 4] =    -2204.3, h[10, 5] =      164.1, h[10, 6] =    -1361.6, h[10, 7] =    -2031.5, h[10, 8] =     1411.8, h[10, 9] =     -714.3, h[10,10] =     1676.5, 
    
    # Check inputs r_rj, colat_rads and elong_rads are all numbers, and convert to numpy doubles here.
    try:
        r_rj       = np.float64(    r_rj  )
        colat_rads = np.float64(colat_rads)
        elong_rads = np.float64(elong_rads)
    except Exception as e:
        print('ERROR: Inputs must be numeric.')
        raise SystemExit
    
    # Check inputs are same size.
    N_input = r_rj.size
    scalar_input = (N_input == 1)  # scalar or not
    
    # Check inputs r_rj, colat_rads and elong_rads are all arrays of the same size (scalar or 1D only)
    if (r_rj.ndim > 1):
        print('ERROR: First  argument    r_rj    must be a scalar number or 1D array of numbers')
        raise SystemExit
    if (colat_rads.ndim > 1):
        print('ERROR: Second argument colat_rads must be a scalar number or 1D array of numbers')
        raise SystemExit
    if (elong_rads.ndim > 1):
        print('ERROR: Third  argument elong_rads must be a scalar number or 1D array of numbers')
        raise SystemExit
    if (N_input != colat_rads.size):
        print('ERROR: First argument r_rj must be the same size as 2nd argument colat_rads')
        raise SystemExit
    if (N_input != elong_rads.size):
        print('ERROR: First argument r_rj must be the same size as 3rd argument elong_rads')
        raise SystemExit
    
    
    # Changing inputs to Doubles, and not using input names (so as not to alter inputs, an IDL issue)
    r_rj_dbl       =              r_rj
    colat_rads_dbl =          colat_rads
    elong_rads_dbl =          elong_rads
    
    # ============
    # Begin hard-coding for JRM09_ORDER10
    # Values from Connerney et al. (2018), https://doi.org/10.1002/2018GL077312
    # See supplemental online information Table S1, https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2F2018GL077312&file=grl57087-sup-0005-2018GL077312-ds01.txt
    # ============
    
    # order = 10; # degree = order for this code 
    # k     = order + 1
    # k     = 11  # order + 1 
    k_plus1 = 12  # k+1, used in for loops when I want to go up to k 
    
    # Arrays rec, g and h are processed (depending on degree) but otherwise do not
    # change. So we calculate them once and use in the code. The initial g and h 
    # values are given in the comments at the top of this code, and are reformatted
    # here in to 1D arrays.
    # g = [       0   ,   410244.7 ,   -71498.3 ,    11670.4 , 
    #        -56835.8 ,    48689.5 ,     4018.6 ,   -37791.1 , 
    #         15926.3 ,    -2710.5 ,   -34645.4 ,    -8247.6 , 
    #         -2406.1 ,   -11083.8 ,   -17837.2 ,   -18023.6 , 
    #          4683.9 ,    16160.0 ,   -16402.0 ,    -2600.7 , 
    #         -3660.7 ,   -20819.6 ,     9992.9 ,    11791.8 , 
    #        -12574.7 ,     2669.7 ,     1113.2 ,     7584.9 , 
    #           598.4 ,     4665.9 ,    -6495.7 ,    -2516.5 , 
    #         -6448.5 ,     1855.3 ,    -2892.9 ,     2968.0 , 
    #         10059.2 ,     1934.4 ,    -6702.9 ,      153.7 , 
    #         -4124.2 ,     -867.2 ,    -3740.6 ,     -732.4 , 
    #         -2433.2 ,     9671.8 ,    -3046.2 ,      260.9 , 
    #          2071.3 ,     3329.6 ,    -2523.1 ,     1787.1 , 
    #         -1148.2 ,     1276.5 ,    -1976.8 ,    -2299.5 , 
    #          2009.7 ,     2127.8 ,     3498.3 ,     2967.6 , 
    #            16.3 ,     1806.5 ,      -46.5 ,     2897.8 , 
    #           574.5 ,     1298.9 ]
    # h = [       0   ,        0   ,    21330.5 ,        0   , 
    #        -42027.3 ,    19353.2 ,        0   ,   -32957.3 , 
    #         42084.5 ,   -27544.2 ,        0   ,    31994.5 , 
    #         27811.2 ,     -926.1 ,      367.1 ,        0   , 
    #         45347.9 ,     -749.0 ,     6268.5 ,    10859.6 , 
    #          9608.4 ,        0   ,    14533.1 ,   -10592.9 , 
    #           568.6 ,    12871.7 ,    -4147.8 ,     3604.4 , 
    #             0   ,    -7626.3 ,   -10948.4 ,     2633.3 , 
    #          5394.2 ,    -6050.8 ,    -1526.0 ,    -5684.2 , 
    #             0   ,    -2409.7 ,   -11614.6 ,     9287.0 , 
    #          -911.9 ,     2754.5 ,    -2446.1 ,     1207.3 , 
    #         -2887.3 ,        0   ,    -8467.4 ,    -1383.8 , 
    #          5697.7 ,    -2056.3 ,     3081.5 ,     -721.2 , 
    #          1352.5 ,     -210.1 ,     1567.6 ,        0   , 
    #         -4692.6 ,     4445.8 ,    -2378.6 ,    -2204.3 , 
    #           164.1 ,    -1361.6 ,    -2031.5 ,     1411.8 , 
    #          -714.3 ,     1676.5 ]
    # These arrays are then extended and manipulated to make larger g and h arrays, and a rec array.
    # ######################################################################
    # The following is the Python code that was used to expand and process the
    # g and h arrays, and create the rec array for pasting the numbers in to
    # this source code:
    #
    # degree = 10      # = order
    # g, h, rec = expand_out_g_and_h(degree,order,g,h)
    # ----------------------------------------------------------------------
    # import numpy as np
    # def expand_out_g_and_h(degree,sh_order,g,h):
    #
    #     # Expand out g and h for later use. i.e. want length = 232 if degree is 20
    #     max_gh_len = int( (degree +1)*(degree)/2+1 + degree + 1 )
    #     # if g and h arrays aren't long enough, pad them to correct size with zeros
    #     if (max_gh_len > len(g)):
    #         g = np.append(g,np.zeros(max_gh_len - len(g),dtype='float64'))
    #     if (max_gh_len > len(h)):
    #         h = np.append(h,np.zeros(max_gh_len - len(h),dtype='float64'))
    #
    #     one_float = np.float64(1)  # = 1.0
    #     two_float = np.float64(2)  # = 2.0
    #     rec = np.zeros(max_gh_len,dtype='float64')
    #
    #     for n in range(1, degree +1 +1):
    #         n2 = np.float64( 2*n-1 )
    #         n2 = n2 * (n2 - two_float)
    #         for m in range(1, n +1):
    #             mn = int( n*(n-1)/2 + m )
    #             rec[mn] = np.float64( (n-m)*(n+m-2) )/n2
    #
    #     s = one_float.copy() # = 1.0
    #     for n in range(2, degree+1 +1):
    #         mn = int( n*(n-1)/2 + 1 )
    #         s = s * np.float64( 2*n - 3 )/np.float64( n - 1 )
    #         p = s.copy() # = a copy of s, not a pointer to s
    #         g[mn] = g[mn] * s
    #         h[mn] = h[mn] * s
    #         for m in range (2, n +1):
    #             if (m == 2):
    #                 aa = two_float.copy() # = 2.0
    #             else:
    #                 aa = one_float.copy() # = 1.0
    #             p = p * np.sqrt( aa*np.float64( n-m+1 )/np.float64( n+m-2 ) )
    #             mnn = int( mn+m-1 )
    #             g[mnn] = g[mnn] * p;
    #             h[mnn] = h[mnn] * p;
    #
    #     # In use, max index called is k*(k-1)/2 + k , where k = order + 1.
    #     # so for k = 11, that's index 66, so size 67 (as indexes start at 0 in Python)
    #     k = sh_order + 1
    #     max_index = int( k*(k-1)/2 + k )
    #     if (len(g) > max_index +1 ):  # +1 for index 0
    #         g   =   g[0:(max_index +1)]
    #     if (len(h) > max_index +1 ):  # +1 for index 0
    #         h   =   h[0:(max_index +1)]
    #     if (len(rec) > max_index +1 ):  # +1 for index 0
    #         rec = rec[0:(max_index +1)]
    #
    #     # Done, return arrays back to main code
    #     return g, h, rec
    # ----------------------------------------------------------------------
    # ######################################################################
    
    rec = np.array([0                         , \
                    0                         ,             0.33333333333333331482962 ,             0                         ,             0.26666666666666666296592 , \
                    0.20000000000000001110223 ,             0                         ,             0.25714285714285711748062 ,             0.22857142857142856429142 , \
                    0.14285714285714284921269 ,             0                         ,             0.25396825396825395415590 ,             0.23809523809523808202115 , \
                    0.19047619047619046561692 ,             0.11111111111111110494321 ,             0                         ,             0.25252525252525254151337 , \
                    0.24242424242424243097105 ,             0.21212121212121212709967 ,             0.16161616161616162989922 ,             0.09090909090909091161414 , \
                    0                         ,             0.25174825174825177231952 ,             0.24475524475524476630817 ,             0.22377622377622377602968 , \
                    0.18881118881118880148406 ,             0.13986013986013987042689 ,             0.07692307692307692734701 ,             0                         , \
                    0.25128205128205127749652 ,             0.24615384615384616751044 ,             0.23076923076923078204103 ,             0.20512820512820512108831 , \
                    0.16923076923076924016343 ,             0.12307692307692308375522 ,             0.06666666666666666574148 ,             0                         , \
                    0.25098039215686274161499 ,             0.24705882352941177515504 ,             0.23529411764705882026405 ,             0.21568627450980393245317 , \
                    0.18823529411764705621124 ,             0.15294117647058824704942 ,             0.10980392156862744945656 ,             0.05882352941176470506601 , \
                    0                         ,             0.25077399380804954454049 ,             0.24767801857585139413409 ,             0.23839009287925697067045 , \
                    0.22291021671826624639401 ,             0.20123839009287924906033 ,             0.17337461300309597866942 ,             0.13931888544891640746570 , \
                    0.09907120743034056320475 ,             0.05263157894736841813099 ,             0                         ,             0.25062656641604008633806 , \
                    0.24812030075187968547468 ,             0.24060150375939848288454 ,             0.22807017543859647856763 ,             0.21052631578947367252397 , \
                    0.18796992481203006475354 ,             0.16040100250626565525636 ,             0.12781954887218044403241 ,             0.09022556390977443108170 , \
                    0.04761904761904761640423 ,             0                         ], dtype='float64')
    
    # This is the modified g array, not the original g coefficients, and will be further modified.
    g = np.array([  0                         , \
                    0                         ,        410244.70000000001164153218269 ,        -71498.30000000000291038304567 ,         17505.59999999999854480847716 , \
               -98442.49328882319969125092030 ,         42166.34389756242308067157865 ,         10046.50000000000000000000000 ,       -115711.13977311669441405683756 , \
                30841.14733335159326088614762 ,         -2142.83839947159822258981876 ,       -151573.62500000000000000000000 ,        -45642.10215250826877309009433 , \
                -9415.35553115892798814456910 ,        -23183.43100524596593459136784 ,        -13190.78728838805909617803991 ,       -141935.84999999997671693563461 , \
                47619.25007516491314163431525 ,        124193.04328343032102566212416 ,        -77191.29987306696421001106501 ,         -5769.73075946518656564876437 , \
                -2568.20347420563030027551576 ,       -300582.97499999997671693563461 ,        188897.03523126235813833773136 ,        176219.39807558275060728192329 , \
              -125279.49167958697944413870573 ,         14568.18469249509325891267508 ,          2590.20913175944133399752900 ,          5094.72643062895076582208276 , \
                16044.59999999999854480847716 ,        165497.62303578440332785248756 ,       -188120.73349108296679332852364 ,        -51533.85619684254197636619210 , \
               -79632.08074599411338567733765 ,         11455.48572598611099238041788 ,         -7006.09637449074671167181805 ,          1921.06723268604082477395423 , \
               505710.56250000005820766091347 ,        129665.25000000000000000000000 ,       -375914.50046967255184426903725 ,          6366.18842053944717918056995 , \
              -110265.51694787663291208446026 ,        -12861.08441189933364512398839 ,        -25680.06972210412277490831912 ,         -1835.99981426163253672712017 , \
                -1524.90263109687475662212819 ,        918443.19531249988358467817307 ,       -388096.44079238711856305599213 ,         28346.79585738653622684068978 , \
               171882.35647868152591399848461 ,        187708.65848237561294808983803 ,        -85005.63010931370081380009651 ,         31091.84173000367445638403296 , \
                -8649.99331942896606051363051 ,          3298.44757452672547515248880 ,         -1203.96883845257525535998866 ,       -414889.08398437500000000000000 , \
               488932.02253022126387804746628 ,        448310.26802277535898610949516 ,        578200.21502919984050095081329 ,        346825.93308208207599818706512 , \
                 1204.82452975374712877965067 ,         74644.73375998696428723633289 ,          -932.00811829940028019336751 ,         23711.52712466476441477425396 , \
                 1525.17384009177953885227907 ,           771.06330156869501024630154 ], dtype='float64')
    
    # This is the modified h array, not the original h coefficients, and will be further modified.
    h = np.array([  0                         , \
                    0                         ,             0                         ,         21330.50000000000000000000000 ,             0                         , \
               -72793.41890493947721552103758 ,         16760.36284452099789632484317 ,             0                         ,       -100910.71037478504877071827650 , \
                81496.28381673301919363439083 ,        -21775.60208180246991105377674 ,             0                         ,        177057.11204695011838339269161 , \
               108828.53403772377350833266973 ,         -1937.07712643301852040167432 ,           271.47411104698363715215237 ,             0                         , \
               461033.11139938322594389319420 ,         -5756.22459277780399133916944 ,         29500.89399184978901757858694 ,         24092.34750470571452751755714 , \
                 6740.87640657726115023251623 ,             0                         ,        274721.00218349619535729289055 ,       -158302.75800766979227773845196 , \
                 5664.86031229477885062806308 ,         70239.09162317455047741532326 ,         -9651.15831540766339458059520 ,          2421.05129224630400130990893 , \
                    0                         ,       -270501.83727851061848923563957 ,       -317074.53216031729243695735931 ,         53925.73158082475129049271345 , \
                66612.60292471760476473718882 ,        -37360.45546854781423462554812 ,         -3695.70433387703633343335241 ,         -3679.15443532142626281711273 , \
                    0                         ,       -161525.20312500000000000000000 ,       -651374.26444599486421793699265 ,        384663.57749869779217988252640 , \
               -24380.75866950407362310215831 ,         40850.84987612628174247220159 ,        -16793.03281485293700825423002 ,          3026.49177465601997027988546 , \
                -1809.49012278727877855999395 ,             0                         ,      -1078776.11541115446016192436218 ,       -150349.92758701223647221922874 , \
               472811.32743136369390413165092 ,       -115925.43081370404979679733515 ,        103818.65529778851487208157778 ,        -12547.38753045641169592272490 , \
                10189.09246170325241109821945 ,          -542.89372143209163823485142 ,           954.74582717435089307400631 ,             0                         , \
             -1141644.22994741331785917282104 ,        936694.13928736466914415359497 ,       -393135.81781678373226895928383 ,       -257618.41363149805692955851555 , \
                12129.55247439201775705441833 ,        -56261.42789238762634340673685 ,        -40717.73101774691895116120577 ,         11552.18924515208345837891102 , \
                -1896.31274843787309691833798 ,           995.21720307946509365137899 ], dtype='float64')
    
    # ============
    # End parts that are hard-coded for JRM09_ORDER10
    # ============
    
    if scalar_input:
        a         = np.array([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],dtype='float64') # = np.zeros(k_plus1,dtype='float64')
        DINDGEN_k = np.array([ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11],dtype='float64') # = 0:k, done manually for speed
    else:
        a         = np.zeros((N_input,k_plus1),dtype='float64')
        DINDGEN_k = a.copy()
        for i in range(k_plus1):
            DINDGEN_k[:,i] = i
    
    da = np.float64(1)/r_rj_dbl
    if scalar_input:
        a[0] = da
        for i in range(1,k_plus1):
            a[i] = a[i-1]*da
    else:
        a[:,0] = da
        for i in range(1,k_plus1):
            a[:,i] = a[:,i-1]*da
    
    b = a  * DINDGEN_k
    
    cos_phi   = np.cos(elong_rads_dbl,dtype='float64')
    sin_phi   = np.sin(elong_rads_dbl,dtype='float64')
    cos_theta = np.cos(colat_rads_dbl,dtype='float64')
    sin_theta = np.sin(colat_rads_dbl,dtype='float64')
    not_bk = (sin_theta >= 0.00001 )  # = 1d-5 - also see bk both times below
    if scalar_input:
        # bk = (sin_theta <  0.00001 )  # bk not needed for scalar
        zero_array = np.float64(0)
        p   = np.float64(1)
        d   = zero_array.copy()
        bbr = zero_array.copy()
        bbt = zero_array.copy()
        bbf = zero_array.copy()
        x = zero_array.copy()
        y = p.copy()
    else:
        bk = (sin_theta <  0.00001 )
        zero_array = np.zeros(N_input,dtype='float64')
        p   = zero_array + np.float64(1)
        d   = zero_array.copy()
        bbr = zero_array.copy()
        bbt = zero_array.copy()
        bbf = zero_array.copy()
        x = zero_array.copy()
        y = p.copy() # 1s
    
    for m in range(1, k_plus1):
        bm  = (m != 1)
        if bm:
            m_minus_1 = np.float64(m - 1)
            w = x.copy()
            x = w *cos_phi + y *sin_phi
            y = y *cos_phi - w *sin_phi
        q = p.copy()
        z = d.copy()
        bi = zero_array.copy()
        p2 = zero_array.copy()
        d2 = zero_array.copy()
        for n in range(m, k_plus1):
            mn = int( n*(n-1)/2 + m )
            w  = g[mn]*y + h[mn]*x
            if scalar_input:
                bbr += b[  n]*w*q
                bbt -= a[  n]*w*z
                if bm:
                    if not_bk:
                        bi += a[n] * (g[mn]*x-h[mn]*y) * q
                    else:
                        bi += a[n] * (g[mn]*x-h[mn]*y) * z
            else:
                bbr += b[:,n]*w*q
                bbt -= a[:,n]*w*z
                if bm:
                    qq = q.copy()
                    ind = np.where(bk)[0]
                    if (len(ind) != 0):
                        qq[ind] = z[ind]
                    bi += a[:,n] * (g[mn]*x-h[mn]*y) * qq
            xk = rec[mn] # faster to write this to xk, to use below twice
            dp = cos_theta *z - sin_theta *q - d2*xk
            pm = cos_theta *q                - p2*xk
            d2 = z.copy()
            p2 = q.copy()
            z = dp.copy()
            q = pm.copy()
        d = sin_theta *d + cos_theta *p
        p = sin_theta *p
        if bm:
            bi  *= m_minus_1
            bbf += bi
    
    # br = bbr  # This doesn't change again
    # bt = bbt  # This doesn't change again
    if scalar_input:
        if not_bk:
            bf = bbf/sin_theta
        else:
            if (cos_theta >= 0):
                bf =  bbf.copy()
            else:
                bf = np.float64(-1)*bbf
    else:
        bf = bbf.copy() # set size of array and do the 3rd case
        ind = np.where((bk == 1) & (cos_theta < 0))[0]
        if (len(ind) != 0):
            bf[ind] = -bbf[ind]
        ind = np.where(bk == 0)[0]
        if (len(ind) != 0):
            bf[ind] =  bbf[ind]/sin_theta[ind]
    
    if scalar_input:
        return             np.array([[bbr,bbt,bf]])
    else:
        return np.transpose(np.array([bbr,bbt,bf]))
