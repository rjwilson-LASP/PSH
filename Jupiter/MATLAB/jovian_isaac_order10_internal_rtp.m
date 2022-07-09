function Brtp = jovian_isaac_order10_internal_rtp( r_rj, colat_rads, elong_rads)
% Code to calculate the ISAAC_ORDER10 model of Jupiter's internal magnetic field model
% with Degree 10 and Order 10.
% Reference: Hess et al. (2017), https://doi.org/10.1553/PRE8s157
%
% Required inputs (System III (1965) Spherical, right handed, and assuming 1 Rj = 71492 km):
%  r_rj       - radial distance, in Rj.
%  colat_rads - colatitude, in radians.                    Value(s) should be 0 <= colat_rads <=  pi.
%  elong_rads - East longitude, right handed, in radians.  Value(s) should be 0 <= elong_rads <= 2pi.
%
% Outputs:
%  B - Spherical Magnetic field vector the ISAAC_ORDER10 internal magnetic field model, [Br, Btheta, Bphi], units of nT.
%
% Usage:
% For internal field only: B = jovian_isaac_order10_internal_rtp(r_rj, colat_rads, elong_rads)
%
% This code was written by Marissa Vogt (mvogt@bu.edu) and Rob Wilson (rob.wilson@lasp.colorado.edu).
% It is based on a routine originally written by K. Khurana, translated into IDL by Marissa Vogt in 2009.
%
% Citation Info:
%  DOI: 10.5281/zenodo.6814109     This DOI links to all versions of code at the Github.
%  Github: https://github.com/rjwilson-LASP/PSH
%  Individual versions released on the Github repository can have a different DOI,
%  See the DOI above for a list of DOIs for each specific Github released version.
%
% Version Info:
%  Last update of this file: 2022-07-09 11:51:05.688642 by user wilsonr. 
%  This code was re-written/re-formatted by the Mother_Source python code:
%   /Users/wilsonr/Documents/JADE/Level2_Processing_Code/IDL/Field_Model/2022/Git_initial/Mother_Source/MOP_spherical.py
%   which itself was last updated at UTC 2022-07-09T17:50:57.
%
%  The Spherical Harmonic g and h values used for this order 10 code are below: 
%  
%  g[i,j] values (nT) used are:
% g[ 1, 0] =   406650.0, g[ 1, 1] =   -71420.0, 
% g[ 2, 0] =   -12860.0, g[ 2, 1] =   -69810.0, g[ 2, 2] =    38520.0, 
% g[ 3, 0] =    -4790.0, g[ 3, 1] =   -46420.0, g[ 3, 2] =    28670.0, g[ 3, 3] =    -9340.0, 
% g[ 4, 0] =   -22300.0, g[ 4, 1] =    18930.0, g[ 4, 2] =     2760.0, g[ 4, 3] =   -13170.0, g[ 4, 4] =     1110.0, 
% g[ 5, 0] =    -1650.0, g[ 5, 1] =     7550.0, g[ 5, 2] =     6230.0, g[ 5, 3] =    -1500.0, g[ 5, 4] =    -4000.0, g[ 5, 5] =     1420.0, 
% g[ 6, 0] =    -6370.0, g[ 6, 1] =     3240.0, g[ 6, 2] =     8020.0, g[ 6, 3] =      270.0, g[ 6, 4] =    -4070.0, g[ 6, 5] =     2790.0, g[ 6, 6] =     -680.0, 
% g[ 7, 0] =      640.0, g[ 7, 1] =     8920.0, g[ 7, 2] =     3660.0, g[ 7, 3] =    -6980.0, g[ 7, 4] =    -1140.0, g[ 7, 5] =     2390.0, g[ 7, 6] =    -1230.0, g[ 7, 7] =      310.0, 
% g[ 8, 0] =    -5110.0, g[ 8, 1] =     7200.0, g[ 8, 2] =     4530.0, g[ 8, 3] =    -3710.0, g[ 8, 4] =     -200.0, g[ 8, 5] =     3120.0, g[ 8, 6] =    -2100.0, g[ 8, 7] =      710.0, g[ 8, 8] =     -133.0, 
% g[ 9, 0] =      100.0, g[ 9, 1] =      330.0, g[ 9, 2] =     -960.0, g[ 9, 3] =      370.0, g[ 9, 4] =     3390.0, g[ 9, 5] =     -740.0, g[ 9, 6] =     -680.0, g[ 9, 7] =      400.0, g[ 9, 8] =      -92.0, g[ 9, 9] =        9.0, 
% g[10, 0] =      280.0, g[10, 1] =      490.0, g[10, 2] =     -630.0, g[10, 3] =      910.0, g[10, 4] =     2850.0, g[10, 5] =     -380.0, g[10, 6] =     -620.0, g[10, 7] =      440.0, g[10, 8] =     -169.0, g[10, 9] =       31.0, g[10,10] =       -0.8, 
%
%  h[i,j] values (nT) used are:
%                        h[ 1, 1] =    23530.0, 
%                        h[ 2, 1] =   -31700.0, h[ 2, 2] =     7950.0, 
%                        h[ 3, 1] =    -7503.0, h[ 3, 2] =    40310.0, h[ 3, 3] =   -36860.0, 
%                        h[ 4, 1] =    20780.0, h[ 4, 2] =    32930.0, h[ 4, 3] =   -16950.0, h[ 4, 4] =     5960.0, 
%                        h[ 5, 1] =       30.0, h[ 5, 2] =    -3340.0, h[ 5, 3] =     2540.0, h[ 5, 4] =     3490.0, h[ 5, 5] =     -840.0, 
%                        h[ 6, 1] =    10510.0, h[ 6, 2] =    -9600.0, h[ 6, 3] =     5030.0, h[ 6, 4] =     1270.0, h[ 6, 5] =    -1060.0, h[ 6, 6] =      180.0, 
%                        h[ 7, 1] =     -270.0, h[ 7, 2] =    -2970.0, h[ 7, 3] =     7390.0, h[ 7, 4] =     1980.0, h[ 7, 5] =    -2120.0, h[ 7, 6] =      350.0, h[ 7, 7] =       19.0, 
%                        h[ 8, 1] =     7810.0, h[ 8, 2] =    -7380.0, h[ 8, 3] =     8520.0, h[ 8, 4] =     -870.0, h[ 8, 5] =    -1680.0, h[ 8, 6] =      720.0, h[ 8, 7] =       -8.0, h[ 8, 8] =      -46.0, 
%                        h[ 9, 1] =     -290.0, h[ 9, 2] =     1890.0, h[ 9, 3] =     1850.0, h[ 9, 4] =      690.0, h[ 9, 5] =      320.0, h[ 9, 6] =     -430.0, h[ 9, 7] =      270.0, h[ 9, 8] =     -120.0, h[ 9, 9] =       27.0, 
%                        h[10, 1] =      230.0, h[10, 2] =     1720.0, h[10, 3] =     1250.0, h[10, 4] =     -110.0, h[10, 5] =     -110.0, h[10, 6] =      -15.0, h[10, 7] =      360.0, h[10, 8] =     -250.0, h[10, 9] =       80.0, h[10,10] =      -14.0, 

%%
%%
% Check inputs are same size.
N_input = numel(r_rj);
scalar_input = (N_input == 1); % scalar or not

% Check inputs r_rj, colat_rads and elong_rads are all numbers,  and same size (scalar or 1D only)
if (~isnumeric(    r_rj  )) || (size(    r_rj  ,2) ~= 1), error('ERROR: First  argument    r_rj    must be a scalar number or 1D column array of numbers'); end
if (~isnumeric(colat_rads)) || (size(colat_rads,2) ~= 1), error('ERROR: Second argument colat_rads must be a scalar number or 1D column array of numbers'); end
if (~isnumeric(elong_rads)) || (size(elong_rads,2) ~= 1), error('ERROR: Third  argument elong_rads must be a scalar number or 1D column array of numbers'); end
if (N_input ~= length(colat_rads)), error('ERROR: First argument r_rj must be the same size as 2nd argument colat_rads'); end
if (N_input ~= length(elong_rads)), error('ERROR: First argument r_rj must be the same size as 3rd argument elong_rads'); end

% Do this check to be sure that user hasn't got position in km, must be in planetary radii.
if scalar_input
    if (    r_rj   <= 0 ) || (    r_rj   >= 200 ), error('ERROR: First  argument, Position    r_rj   , must be in units of Rj and >0 and <200 only, and not outside that range (did you use km instead?)'),end
    if (colat_rads <  0 ) || (colat_rads >  pi  ), error('ERROR: Second argument, Position colat_rads, must be in units of radians and >= 0 and <=   Pi only, and not outside that range (did you use degrees instead?)'),end
    if (elong_rads <  0 ) || (elong_rads >  pi*2), error('ERROR: Third  argument, Position elong_rads, must be in units of radians and >= 0 and <= 2*Pi only, and not outside that range (did you use degress instead?)'),end
else
    if (min(r_rj )      <=    0) || (max(r_rj ) >=      200), error('ERROR: First  argument, Position    r_rj   , must be in units of Rj and >0 and <200 only, and not outside that range (did you use km instead?)'); end
    if (min(colat_rads) <     0) || (max(colat_rads) > pi  ), error('ERROR: Second argument, Position colat_rads, must be in units of radians and >= 0 and <=   Pi only, and not outside that range (did you use degrees instead?)'); end
    if (min(elong_rads) <     0) || (max(elong_rads) > pi*2), error('ERROR: Third  argument, Position elong_rads, must be in units of radians and >= 0 and <= 2*Pi only, and not outside that range (did you use degrees instead?)'); end
end

%%
% Changing inputs to Doubles, and not using input names (so as not to alter inputs, an IDL issue)
r_rj_dbl       = double(      r_rj  );
colat_rads_dbl = double(  colat_rads);
elong_rads_dbl = double(  elong_rads);

%%
% ============
% Begin hard-coding for ISAAC_ORDER10
% Values from Hess et al. (2017), https://doi.org/10.1553/PRE8s157
% Values are give in nT here, but in the original paper were given to 4 (mostly) to 6 decimal places in units of G.
% ============

% order = 10; % degree = order for this code 
% k     = order + 1;
k       = 11; % order + 1 

%%
% Arrays rec, g and h are processed (depending on degree) but otherwise do not
% change. So we calculate them once and use in the code. The initial g and h 
% values are given in the comments at the top of this code, and are reformatted
% here in to 1D arrays.
% g = [       0   ,   406650.0 ,   -71420.0 ,   -12860.0 , 
%        -69810.0 ,    38520.0 ,    -4790.0 ,   -46420.0 , 
%         28670.0 ,    -9340.0 ,   -22300.0 ,    18930.0 , 
%          2760.0 ,   -13170.0 ,     1110.0 ,    -1650.0 , 
%          7550.0 ,     6230.0 ,    -1500.0 ,    -4000.0 , 
%          1420.0 ,    -6370.0 ,     3240.0 ,     8020.0 , 
%           270.0 ,    -4070.0 ,     2790.0 ,     -680.0 , 
%           640.0 ,     8920.0 ,     3660.0 ,    -6980.0 , 
%         -1140.0 ,     2390.0 ,    -1230.0 ,      310.0 , 
%         -5110.0 ,     7200.0 ,     4530.0 ,    -3710.0 , 
%          -200.0 ,     3120.0 ,    -2100.0 ,      710.0 , 
%          -133.0 ,      100.0 ,      330.0 ,     -960.0 , 
%           370.0 ,     3390.0 ,     -740.0 ,     -680.0 , 
%           400.0 ,      -92.0 ,        9.0 ,      280.0 , 
%           490.0 ,     -630.0 ,      910.0 ,     2850.0 , 
%          -380.0 ,     -620.0 ,      440.0 ,     -169.0 , 
%            31.0 ,       -0.8 ]
% h = [       0   ,        0   ,    23530.0 ,        0   , 
%        -31700.0 ,     7950.0 ,        0   ,    -7503.0 , 
%         40310.0 ,   -36860.0 ,        0   ,    20780.0 , 
%         32930.0 ,   -16950.0 ,     5960.0 ,        0   , 
%            30.0 ,    -3340.0 ,     2540.0 ,     3490.0 , 
%          -840.0 ,        0   ,    10510.0 ,    -9600.0 , 
%          5030.0 ,     1270.0 ,    -1060.0 ,      180.0 , 
%             0   ,     -270.0 ,    -2970.0 ,     7390.0 , 
%          1980.0 ,    -2120.0 ,      350.0 ,       19.0 , 
%             0   ,     7810.0 ,    -7380.0 ,     8520.0 , 
%          -870.0 ,    -1680.0 ,      720.0 ,       -8.0 , 
%           -46.0 ,        0   ,     -290.0 ,     1890.0 , 
%          1850.0 ,      690.0 ,      320.0 ,     -430.0 , 
%           270.0 ,     -120.0 ,       27.0 ,        0   , 
%           230.0 ,     1720.0 ,     1250.0 ,     -110.0 , 
%          -110.0 ,      -15.0 ,      360.0 ,     -250.0 , 
%            80.0 ,      -14.0 ]
% These arrays are then extended and manipulated to make larger g and h arrays, and a rec array.
% ######################################################################
% The following is the Python code that was used to expand and process the
% g and h arrays, and create the rec array for pasting the numbers in to
% this source code:
%
% degree = 10      # = order
% g, h, rec = expand_out_g_and_h(degree,order,g,h)
% ----------------------------------------------------------------------
% import numpy as np
% def expand_out_g_and_h(degree,sh_order,g,h):
%
%     # Expand out g and h for later use. i.e. want length = 232 if degree is 20
%     max_gh_len = int( (degree +1)*(degree)/2+1 + degree + 1 )
%     # if g and h arrays aren't long enough, pad them to correct size with zeros
%     if (max_gh_len > len(g)):
%         g = np.append(g,np.zeros(max_gh_len - len(g),dtype='float64'))
%     if (max_gh_len > len(h)):
%         h = np.append(h,np.zeros(max_gh_len - len(h),dtype='float64'))
%
%     one_float = np.float64(1)  # = 1.0
%     two_float = np.float64(2)  # = 2.0
%     rec = np.zeros(max_gh_len,dtype='float64')
%
%     for n in range(1, degree +1 +1):
%         n2 = np.float64( 2*n-1 )
%         n2 = n2 * (n2 - two_float)
%         for m in range(1, n +1):
%             mn = int( n*(n-1)/2 + m )
%             rec[mn] = np.float64( (n-m)*(n+m-2) )/n2
%
%     s = one_float.copy() # = 1.0
%     for n in range(2, degree+1 +1):
%         mn = int( n*(n-1)/2 + 1 )
%         s = s * np.float64( 2*n - 3 )/np.float64( n - 1 )
%         p = s.copy() # = a copy of s, not a pointer to s
%         g[mn] = g[mn] * s
%         h[mn] = h[mn] * s
%         for m in range (2, n +1):
%             if (m == 2):
%                 aa = two_float.copy() # = 2.0
%             else:
%                 aa = one_float.copy() # = 1.0
%             p = p * np.sqrt( aa*np.float64( n-m+1 )/np.float64( n+m-2 ) )
%             mnn = int( mn+m-1 )
%             g[mnn] = g[mnn] * p;
%             h[mnn] = h[mnn] * p;
%
%     # In use, max index called is k*(k-1)/2 + k , where k = order + 1.
%     # so for k = 11, that's index 66, so size 67 (as indexes start at 0 in Python)
%     k = sh_order + 1
%     max_index = int( k*(k-1)/2 + k )
%     if (len(g) > max_index +1 ):  # +1 for index 0
%         g   =   g[0:(max_index +1)]
%     if (len(h) > max_index +1 ):  # +1 for index 0
%         h   =   h[0:(max_index +1)]
%     if (len(rec) > max_index +1 ):  # +1 for index 0
%         rec = rec[0:(max_index +1)]
%
%     # Done, return arrays back to main code
%     return g, h, rec
% ----------------------------------------------------------------------
% ######################################################################

rec = [... % MATLAB starts at index 1, not 0, so no value on this first line compared to IDL & Python output
                0                         ,             0.33333333333333331482962 ,             0                         ,             0.26666666666666666296592 , ...
                0.20000000000000001110223 ,             0                         ,             0.25714285714285711748062 ,             0.22857142857142856429142 , ...
                0.14285714285714284921269 ,             0                         ,             0.25396825396825395415590 ,             0.23809523809523808202115 , ...
                0.19047619047619046561692 ,             0.11111111111111110494321 ,             0                         ,             0.25252525252525254151337 , ...
                0.24242424242424243097105 ,             0.21212121212121212709967 ,             0.16161616161616162989922 ,             0.09090909090909091161414 , ...
                0                         ,             0.25174825174825177231952 ,             0.24475524475524476630817 ,             0.22377622377622377602968 , ...
                0.18881118881118880148406 ,             0.13986013986013987042689 ,             0.07692307692307692734701 ,             0                         , ...
                0.25128205128205127749652 ,             0.24615384615384616751044 ,             0.23076923076923078204103 ,             0.20512820512820512108831 , ...
                0.16923076923076924016343 ,             0.12307692307692308375522 ,             0.06666666666666666574148 ,             0                         , ...
                0.25098039215686274161499 ,             0.24705882352941177515504 ,             0.23529411764705882026405 ,             0.21568627450980393245317 , ...
                0.18823529411764705621124 ,             0.15294117647058824704942 ,             0.10980392156862744945656 ,             0.05882352941176470506601 , ...
                0                         ,             0.25077399380804954454049 ,             0.24767801857585139413409 ,             0.23839009287925697067045 , ...
                0.22291021671826624639401 ,             0.20123839009287924906033 ,             0.17337461300309597866942 ,             0.13931888544891640746570 , ...
                0.09907120743034056320475 ,             0.05263157894736841813099 ,             0                         ,             0.25062656641604008633806 , ...
                0.24812030075187968547468 ,             0.24060150375939848288454 ,             0.22807017543859647856763 ,             0.21052631578947367252397 , ...
                0.18796992481203006475354 ,             0.16040100250626565525636 ,             0.12781954887218044403241 ,             0.09022556390977443108170 , ...
                0.04761904761904761640423 ,             0                         ];

% This is the modified g array, not the original g coefficients, and will be further modified.
g = [  ... % MATLAB starts at index 1, not 0, so no value on this first line compared to IDL & Python output
                0                         ,        406650.00000000000000000000000 ,        -71420.00000000000000000000000 ,        -19290.00000000000000000000000 , ...
          -120914.46687638331786729395390 ,         33359.29855377657804638147354 ,        -11975.00000000000000000000000 ,       -142131.64232499391073361039162 , ...
            55519.21626788332650903612375 ,         -7383.91833649316595256095752 ,        -97562.50000000000000000000000 ,        104758.35318722798547241836786 , ...
            10800.20833132398547604680061 ,        -27547.03137363443966023623943 ,           820.85606990507176305982284 ,        -12993.75000000000000000000000 , ...
            76757.68869264824024867266417 ,         47878.87745394621742889285088 ,         -7059.31897388126071746228263 ,         -8874.11967464942063088528812 , ...
              996.21627922856157510977937 ,        -91966.87500000000000000000000 ,         61246.12416308479441795498133 ,        119852.74280145303055178374052 , ...
             2689.96180851141434686724097 ,        -22209.42866181781937484629452 ,          6491.81052605896547902375460 ,          -456.75143677934931929485174 , ...
            17160.00000000000000000000000 ,        316388.86334452027222141623497 ,        105996.56458539706363808363676 ,       -142939.12825510071706958115101 , ...
           -14077.78119724483076424803585 ,         14756.97239535751941730268300 ,         -2978.84425338712617303826846 ,           200.65055327920237004946102 , ...
          -256897.26562500000000000000000 ,        482625.00000000000000000000000 ,        254053.12433836350101046264172 ,       -153666.61704750391072593629360 , ...
            -5347.24392356707448925590143 ,         46271.42915720240125665441155 ,        -14416.97760156623553484678268 ,          1779.84689804172467120224610 , ...
              -83.35198501392584091718163 ,          9496.09375000000000000000000 ,         42043.14406850757222855463624 ,       -104304.03995052157551981508732 , ...
            30703.65079762089226278476417 ,        191113.75308002563542686402798 ,        -24931.30128845156650641001761 ,        -11830.59279077975406835321337 , ...
             3013.40997018950156416394748 ,          -237.72595131724148131979746 ,             5.48144452957971317630381 ,         50519.21875000000000000000000 , ...
           119210.17616550152888521552086 ,       -132735.90979149754275567829609 ,        150405.10981807499774731695652 ,        333081.92117668618448078632355 , ...
           -28087.93382247999397804960608 ,        -25618.45277121058461489155889 ,          8819.00154949970055895391852 ,         -1382.85874941967858831048943 , ...
               82.29832731565738868084736 ,            -0.47490233370925860612033 ];

% This is the modified h array, not the original h coefficients, and will be further modified.
h = [  ... % MATLAB starts at index 1, not 0, so no value on this first line compared to IDL & Python output
                0                         ,             0                         ,         23530.00000000000000000000000 ,             0                         , ...
           -54906.01059993340459186583757 ,          6884.90196008628663548734039 ,             0                         ,        -22973.15192512773137423209846 , ...
            78059.97934281048947013914585 ,        -29140.38863845161904464475811 ,             0                         ,        114996.22711202311620581895113 , ...
           128859.00737336913880426436663 ,        -35453.46862438145035412162542 ,          4407.47943840921379887731746 ,             0                         , ...
              304.99743851383408355104621 ,        -25668.61166872879039146937430 ,         11953.78012910560210002586246 ,          7742.66941613161998247960582 , ...
             -589.31103841689559885708150 ,             0                         ,        198671.84103519172640517354012 ,       -143464.62978727545123547315598 , ...
            50112.99221041634882567450404 ,          6930.21484041981057089287788 ,         -2466.42263714068212721031159 ,           120.90479208865129123751103 , ...
                0                         ,         -9576.79294876911080791614950 ,        -86013.60568815008446108549833 ,        151335.26616120262769982218742 , ...
            24450.88313205680969986133277 ,        -13089.86672726273718581069261 ,           847.63860868739368470414774 ,            12.29793713646724206967065 , ...
                0                         ,        523514.06250000000000000000000 ,       -413887.87143865844700485467911 ,        352894.76475599280092865228653 , ...
           -23260.51106751677070860750973 ,        -24915.38493080129046575166285 ,          4942.96374910842314420733601 ,           -20.05461293568140490606311 , ...
              -28.82850609504201955246572 ,             0                         ,        -36947.00539353695785393938422 ,        205348.57865258934907615184784 , ...
           153518.25398810446495190262794 ,         38899.25947646539134439080954 ,         10781.10325987094620359130204 ,         -7481.11014711072766658617184 , ...
             2034.05172987791365812881850 ,          -310.07732780509758185871760 ,            16.44433358873913775255460 ,             0                         , ...
            55955.79697564357775263488293 ,        362390.10292281868169084191322 ,        206600.42557427886640653014183 ,        -12855.79344892473091022111475 , ...
            -8130.71768545473514677723870 ,          -619.80127672283663287089439 ,          7215.54672231793756509432569 ,         -2045.64903760307493030268233 , ...
              212.38278016943840498242935 ,            -8.31079083991202516301655 ];

% ============
% End parts that are hard-coded for ISAAC_ORDER10
% ============

%%
if scalar_input
    a         = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; % = zeros(1, k)
    DINDGEN_k = [ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11]; % = 1:k, done manually for speed
else
    a         = zeros( N_input,k  );
    DINDGEN_k = a;
    for i = 1:k
        DINDGEN_k(:,i) = i;
    end
end

%%
da = 1./r_rj_dbl;
a(:,1) = da.*da;
for i=2:k
    a(:,i) = a(:,i-1).*da;
end

b = a .* DINDGEN_k;

cos_phi   = cos(elong_rads_dbl);
sin_phi   = sin(elong_rads_dbl);
cos_theta = cos(colat_rads_dbl);
sin_theta = sin(colat_rads_dbl);
not_bk = (sin_theta >= 0.00001 ); % = 1e-5 - also see bk both times below
if scalar_input
    % bk = (sin_theta <  0.00001 ); % bk not needed for scalar
    zero_array = 0;
    p   = 1;
    d   = 0;
    bbr = 0;
    bbt = 0;
    bbf = 0;
    x = 0;
    y = 1;
else
    bk = (sin_theta <  0.00001 );
    zero_array = zeros( N_input,1);
    p   =         ones( N_input,1);
    d   = zero_array;
    bbr = zero_array;
    bbt = zero_array;
    bbf = zero_array;
    x = zero_array;
    y = p; % 1s
end

for m = 1:k
    bm  = (m ~= 1);
    if bm
        m_minus_1 =        m - 1;
        w = x;
        x = w.*cos_phi + y.*sin_phi;
        y = y.*cos_phi - w.*sin_phi;
    end
    q = p;
    z = d;
    bi = zero_array;
    p2 = zero_array;
    d2 = zero_array;
    for n = m:k
        mn = n*(n-1)/2 + m;
        w  = g(mn)*y + h(mn)*x;
        bbr = bbr + b(:,n).*w.*q;
        bbt = bbt - a(:,n).*w.*z;
        if bm
            if scalar_input
                if not_bk
                    bi = bi + a(  n)  * (g(mn)*x-h(mn)*y)  * q;
                else
                    bi = bi + a(  n)  * (g(mn)*x-h(mn)*y)  * z;
                end
            else
                qq = q;
                ind = find(bk);
                if numel(ind)
                    qq(ind) = z(ind);
                end
                bi = bi + a(:,n) .* (g(mn)*x-h(mn)*y) .* qq;
            end
        end
        dp = cos_theta.*z - sin_theta.*q - d2*rec(mn);
        pm = cos_theta.*q                - p2*rec(mn);
        d2 = z;
        p2 = q;
        z = dp;
        q = pm;
    end
    d = sin_theta.*d + cos_theta.*p;
    p = sin_theta.*p;
    if bm
        bi  = bi  * m_minus_1;
        bbf = bbf + bi;
    end
end

% br = bbr; % This doesn't change again
% bt = bbt; % This doesn't change again
if scalar_input
    if not_bk
        bf = bbf/sin_theta;
    else
        if (cos_theta >= 0)
            bf =  bbf;
        else
            bf = -bbf;
        end
    end
else
    bf = bbf; % set size of array and do the 3rd case
    ind = find(bk & (cos_theta <  0));
    if numel(ind)
        bf(ind) = -bbf(ind);
    end
    ind = find(~bk); % find(bk == 0)
    if numel(ind)
        bf(ind) = bbf(ind)./sin_theta(ind);
    end
end

%%
Brtp = [ bbr , bbt , bf ];
