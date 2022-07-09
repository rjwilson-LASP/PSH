function Brtp = jovian_vip4_order04_internal_rtp( r_rj, colat_rads, elong_rads)
% Code to calculate the VIP4_ORDER04 model of Jupiter's internal magnetic field model
% with Degree 4 and Order 4.
% Reference: Connerney et al. (1998), https://doi.org/10.1029/97JA03726
%
% Required inputs (System III (1965) Spherical, right handed, and assuming 1 Rj = 71492 km):
%  r_rj       - radial distance, in Rj.
%  colat_rads - colatitude, in radians.                    Value(s) should be 0 <= colat_rads <=  pi.
%  elong_rads - East longitude, right handed, in radians.  Value(s) should be 0 <= elong_rads <= 2pi.
%
% Outputs:
%  B - Spherical Magnetic field vector the VIP4_ORDER04 internal magnetic field model, [Br, Btheta, Bphi], units of nT.
%
% Usage:
% For internal field only: B = jovian_vip4_order04_internal_rtp(r_rj, colat_rads, elong_rads)
%
% This code was written by Marissa Vogt (mvogt@bu.edu) and Rob Wilson (rob.wilson@lasp.colorado.edu).
% It is based on a routine originally written by K. Khurana, translated into IDL by Marissa Vogt in 2009.
%
% Citation Info:
%  DOI: 10.5281/zenodo.6814109     This doi links to all versions of code at the Github.
%  Github: https://github.com/rjwilson-LASP/PSH
%  Individual versions released on the Github repository can have a different DOI,
%  See the DOI above for a list of DOIs for each specific Github released version.
% Version Info:
%  Last update of this file: 2022-07-09 11:45:40.165977 by user wilsonr. 
%  This code was re-written/re-formatted by Rob's python code:
%   /Users/wilsonr/Documents/JADE/Level2_Processing_Code/IDL/Field_Model/2022/Git_initial/Mother_Source/MOP_spherical.py
%   which itself was last updated at UTC 2022-07-09T17:45:29.
%
%  The Spherical Harmonic g and h values used for this order 4 code are below: 
%  
%  g[i,j] values (nT) used are:
% g[ 1, 0] =     420543, g[ 1, 1] =     -65920, 
% g[ 2, 0] =      -5118, g[ 2, 1] =     -61904, g[ 2, 2] =      49690, 
% g[ 3, 0] =      -1576, g[ 3, 1] =     -52036, g[ 3, 2] =      24386, g[ 3, 3] =     -17597, 
% g[ 4, 0] =     -16758, g[ 4, 1] =      22210, g[ 4, 2] =      -6074, g[ 4, 3] =     -20243, g[ 4, 4] =       6643, 
%
%  h[i,j] values (nT) used are:
%                        h[ 1, 1] =      24992, 
%                        h[ 2, 1] =     -36052, h[ 2, 2] =       5250, 
%                        h[ 3, 1] =      -8804, h[ 3, 2] =      40829, h[ 3, 3] =     -31586, 
%                        h[ 4, 1] =       7557, h[ 4, 2] =      40411, h[ 4, 3] =     -16597, h[ 4, 4] =       3866, 

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

% Scaling distances since vip4_order04 expects 1Rj to be 71323 km (not the 71492 km that the inputs expect)
r_rj_dbl = r_rj_dbl * double(71492.000000)  /  double(71323.000000);

%%
% ============
% Begin hard-coding for VIP4_ORDER04
% Values from Connerney et al. (1998), https://doi.org/10.1029/97JA03726
% Original paper is Connerney et al (1998) [https://doi.org/10.1029/97JA03726], however table 3 of Connerney (2007) [https://doi.org/10.1016/B978-044452748-6.00159-0] provides the g and h values to more significant figures, which are used here. i.e. 4.205 G (1998) -> 420543 nT (2007)
% ============

% order = 4; % degree = order for this code 
% k     = order + 1;
k       = 5; % order + 1 

%%
% Arrays rec, g and h are processed (depending on degree) but otherwise do not
% change. So we calculate them once and use in the code. The initial g and h 
% values are given in the comments at the top of this code, and are reformatted
% here in to 1D arrays.
% g = [            0                         ,        420543.00000000000000000000000 ,        -65920.00000000000000000000000 ,         -5118.00000000000000000000000 , 
%             -61904.00000000000000000000000 ,         49690.00000000000000000000000 ,         -1576.00000000000000000000000 ,        -52036.00000000000000000000000 , 
%              24386.00000000000000000000000 ,        -17597.00000000000000000000000 ,        -16758.00000000000000000000000 ,         22210.00000000000000000000000 , 
%              -6074.00000000000000000000000 ,        -20243.00000000000000000000000 ,          6643.00000000000000000000000 ]
% h = [            0                         ,             0                         ,         24992.00000000000000000000000 ,             0                         , 
%             -36052.00000000000000000000000 ,          5250.00000000000000000000000 ,             0                         ,         -8804.00000000000000000000000 , 
%              40829.00000000000000000000000 ,        -31586.00000000000000000000000 ,             0                         ,          7557.00000000000000000000000 , 
%              40411.00000000000000000000000 ,        -16597.00000000000000000000000 ,          3866.00000000000000000000000 ]
% These arrays are then extended and manipulated to make larger g and h arrays, and a rec array.
% ######################################################################
% The following is the Python code that was used to expand and process the
% g and h arrays, and create the rec array for pasting the numbers in to
% this source code:
%
% degree = 4      # = order
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
                0.19047619047619046561692 ,             0.11111111111111110494321 ,             0                         ];

% This is the modified g array, not the original g coefficients, and will be further modified.
g = [  ... % MATLAB starts at index 1, not 0, so no value on this first line compared to IDL & Python output
                0                         ,        420543.00000000000000000000000 ,        -65920.00000000000000000000000 ,         -7677.00000000000000000000000 , ...
          -107220.87319174376898445188999 ,         43032.80231404875667067244649 ,         -3940.00000000000000000000000 ,       -159327.06031933182384818792343 , ...
            47223.28594030703243333846331 ,        -13911.64999649574383511207998 ,        -73316.25000000000000000000000 ,        122909.82695659447927027940750 , ...
           -23768.28456683401600457727909 ,        -42341.27229282323241932317615 ,          4912.56474989134403585921973 ];

% This is the modified h array, not the original h coefficients, and will be further modified.
h = [  ... % MATLAB starts at index 1, not 0, so no value on this first line compared to IDL & Python output
                0                         ,             0                         ,         24992.00000000000000000000000 ,             0                         , ...
           -62443.89571447316120611503720 ,          4546.63336986830290697980672 ,             0                         ,        -26956.63461932887366856448352 , ...
            79065.01852115131623577326536 ,        -24970.92554351961007341742516 ,             0                         ,         41820.33148631177027709782124 , ...
           158133.05031780200079083442688 ,        -34715.11615096512832678854465 ,          2858.94555518288962048245594 ];

% ============
% End parts that are hard-coded for VIP4_ORDER04
% ============

%%
if scalar_input
    a         = [ 0, 0, 0, 0, 0]; % = zeros(1, k)
    DINDGEN_k = [ 1, 2, 3, 4, 5]; % = 1:k, done manually for speed
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
