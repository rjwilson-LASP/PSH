function Bxyz = jovian_vipal_order05_internal_xyz( x_rj, y_rj, z_rj)
% Code to calculate the VIPAL_ORDER05 model of Jupiter's internal magnetic field model
% with Degree 5 and Order 5.
% Reference: Hess et al. (2011), https://doi.org/10.1029/2010JA016262
%
% Required inputs (System III (1965) Cartesian, right handed, and assuming 1 Rj = 71492 km):
%  x_rj       - Jupiter SYSIII right-handed position in x, in Rj.
%  y_rj       - Jupiter SYSIII right-handed position in y, in Rj.
%  z_rj       - Jupiter SYSIII right-handed position in z, in Rj.
%
% Outputs:
%  B - Cartesian Magnetic field vector the VIPAL_ORDER05 internal magnetic field model, [Bx, By, Bz], units of nT.
%
% Usage:
% For internal field only: B = jovian_vipal_order05_internal_xyz(x_rj, y_rj, z_rj)
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
%  Last update of this file: 2022-08-31 11:48:40.505033 by user wilsonr. 
%  This code was re-written/re-formatted by the Mother_Source python code:
%   /Volumes/wilsonr/Documents/JADE/Level2_Processing_Code/IDL/Field_Model/2022/Git_initial/Mother_Source/MOP_spherical.py
%   which itself was last updated at UTC 2022-08-31T17:46:34.
%
%  The Spherical Harmonic g and h values used for this order 5 code are below: 
%  
%  g[i,j] values (nT) used are:
% g[ 1, 0] =     420000, g[ 1, 1] =     -69750, 
% g[ 2, 0] =      64410, g[ 2, 1] =     -86720, g[ 2, 2] =      95980, 
% g[ 3, 0] =     -10580, g[ 3, 1] =     -59000, g[ 3, 2] =      63220, g[ 3, 3] =      46710, 
% g[ 4, 0] =     -74660, g[ 4, 1] =      32820, g[ 4, 2] =     -33800, g[ 4, 3] =      18260, g[ 4, 4] =     -14290, 
% g[ 5, 0] =      -6600, g[ 5, 1] =       7370, g[ 5, 2] =     -17110, g[ 5, 3] =     -17930, g[ 5, 4] =       -770, g[ 5, 5] =      -7400, 
%
%  h[i,j] values (nT) used are:
%                        h[ 1, 1] =      19730, 
%                        h[ 2, 1] =     -40410, h[ 2, 2] =      60300, 
%                        h[ 3, 1] =     -23100, h[ 3, 2] =      51600, h[ 3, 3] =     -11310, 
%                        h[ 4, 1] =      32830, h[ 4, 2] =     -21310, h[ 4, 3] =      -6060, h[ 4, 4] =      -4860, 
%                        h[ 5, 1] =      20650, h[ 5, 2] =     -11670, h[ 5, 3] =      -2880, h[ 5, 4] =       -500, h[ 5, 5] =     -22790, 

%%
% Check inputs are same size.
N_input = numel(x_rj);
scalar_input = (N_input == 1); % scalar or not

% Check inputs x_rj, y_rj and z_rj are all numbers,  and same size (scalar or 1D only)
if (N_input ~= length(y_rj)), error('ERROR: First argument x_rj must be the same size as 2nd argument y_rj'); end
if (N_input ~= length(z_rj)), error('ERROR: First argument x_rj must be the same size as 3rd argument z_rj'); end
if (~isnumeric(x_rj)) || (size(x_rj,2) ~= 1), error('ERROR: First  argument x_rj must be a scalar number or 1D column array of numbers'); end
if (~isnumeric(y_rj)) || (size(y_rj,2) ~= 1), error('ERROR: Second argument y_rj must be a scalar number or 1D column array of numbers'); end
if (~isnumeric(z_rj)) || (size(z_rj,2) ~= 1), error('ERROR: Third  argument z_rj must be a scalar number or 1D column array of numbers'); end

% Changing inputs to Doubles, and not using input names (so as not to alter inputs, an IDL issue)
x_in = double(x_rj); % X in SYSIII, units Rj
y_in = double(y_rj); % Y in SYSIII, units Rj
z_in = double(z_rj); % Z in SYSIII, units Rj

rho_rj_sq = x_in.*x_in + y_in.*y_in;
r_rj = sqrt(rho_rj_sq + z_in.*z_in);

colat_rads = acos(z_in./r_rj);
elong_rads = atan2(y_in,x_in);

%%
% ######################################################################
% Start of RTP code.
% ######################################################################
%%
% Do this check to be sure that user hasn't got position in km, must be in planetary radii.
if scalar_input
    if (    r_rj   <= 0 ) || (    r_rj   >= 200), error('ERROR: First  argument, Position    r_rj   , must be in units of Rj and >0 and <200 only, and not outside that range (did you use km instead?)'),end
else
    if (min(r_rj ) <= 0 ) || (max(r_rj ) >= 200), error('ERROR: First  argument, Position    r_rj   , must be in units of Rj and >0 and <200 only, and not outside that range (did you use km instead?)'); end
end

%%
% Code is not using input names (so as not to alter inputs, an IDL issue)
r_rj_dbl       =  r_rj;
colat_rads_dbl = colat_rads;
elong_rads_dbl = elong_rads;

%%
% ============
% Begin hard-coding for VIPAL_ORDER05
% Values from Hess et al. (2011), https://doi.org/10.1029/2010JA016262
% Values are give to whole nT here, but in the original paper were given to 1 to 4 decimal places in units of G.
% ============

% order = 5; % degree = order for this code 
% k     = order + 1;
k       = 6; % order + 1 

%%
% Arrays rec, g and h are processed (depending on degree) but otherwise do not
% change. So we calculate them once and use in the code. The initial g and h 
% values are given in the comments at the top of this code, and are reformatted
% here in to 1D arrays.
% g = [            0                         ,        420000.00000000000000000000000 ,        -69750.00000000000000000000000 ,         64410.00000000000000000000000 , 
%             -86720.00000000000000000000000 ,         95980.00000000000000000000000 ,        -10580.00000000000000000000000 ,        -59000.00000000000000000000000 , 
%              63220.00000000000000000000000 ,         46710.00000000000000000000000 ,        -74660.00000000000000000000000 ,         32820.00000000000000000000000 , 
%             -33800.00000000000000000000000 ,         18260.00000000000000000000000 ,        -14290.00000000000000000000000 ,         -6600.00000000000000000000000 , 
%               7370.00000000000000000000000 ,        -17110.00000000000000000000000 ,        -17930.00000000000000000000000 ,          -770.00000000000000000000000 , 
%              -7400.00000000000000000000000 ]
% h = [            0                         ,             0                         ,         19730.00000000000000000000000 ,             0                         , 
%             -40410.00000000000000000000000 ,         60300.00000000000000000000000 ,             0                         ,        -23100.00000000000000000000000 , 
%              51600.00000000000000000000000 ,        -11310.00000000000000000000000 ,             0                         ,         32830.00000000000000000000000 , 
%             -21310.00000000000000000000000 ,         -6060.00000000000000000000000 ,         -4860.00000000000000000000000 ,             0                         , 
%              20650.00000000000000000000000 ,        -11670.00000000000000000000000 ,         -2880.00000000000000000000000 ,          -500.00000000000000000000000 , 
%             -22790.00000000000000000000000 ]
% These arrays are then extended and manipulated to make larger g and h arrays, and a rec array.
% ######################################################################
% The following is the Python code that was used to expand and process the
% g and h arrays, and create the rec array for pasting the numbers in to
% this source code:
%
% degree = 5      # = order
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
                0                         ];

% This is the modified g array, not the original g coefficients, and will be further modified.
g = [  ... % MATLAB starts at index 1, not 0, so no value on this first line compared to IDL & Python output
                0                         ,        420000.00000000000000000000000 ,        -69750.00000000000000000000000 ,         96615.00000000000000000000000 , ...
          -150203.44603237303090281784534 ,         83121.11825523042352870106697 ,        -26450.00000000000000000000000 ,       -180649.86853025938034988939762 , ...
           122425.00357361645728815346956 ,         36927.49737661625113105401397 ,       -326637.50000000000000000000000 ,        181625.41741177084622904658318 , ...
          -132263.42086911256774328649044 ,         38193.53021128055115696042776 ,        -10567.59751256168965483084321 ,        -51975.00000000000000000000000 , ...
            74927.70406156523677054792643 ,       -131493.99570417654467746615410 ,        -84382.39280112733831629157066 ,         -1708.26803737001364424941130 , ...
            -5191.54962414884175814222544 ];

% This is the modified h array, not the original h coefficients, and will be further modified.
h = [  ... % MATLAB starts at index 1, not 0, so no value on this first line compared to IDL & Python output
                0                         ,             0                         ,         19730.00000000000000000000000 ,             0                         , ...
           -69992.17313385833404026925564 ,         52221.33184820164751727133989 ,             0                         ,        -70729.01632286426320206373930 , ...
            99922.97033215136616490781307 ,         -8941.34008412609364313539118 ,             0                         ,        181680.75727082381490617990494 , ...
           -83388.56505091091094072908163 ,        -12675.39940199124430364463478 ,         -3594.01846823301684707985260 ,             0                         , ...
           209939.90351035579806193709373 ,        -89686.43657906138105317950249 ,        -13553.89242985202145064249635 ,         -1109.26495933117757886066101 , ...
           -15988.56972085839333885814995 ];

% ============
% End parts that are hard-coded for VIPAL_ORDER05
% ============

%%
if scalar_input
    a         = [ 0, 0, 0, 0, 0, 0]; % = zeros(1, k)
    DINDGEN_k = [ 1, 2, 3, 4, 5, 6]; % = 1:k, done manually for speed
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
% ######################################################################
% End of RTP code.
% ######################################################################
% Brtp = [ bbr , bbt , bf ];

%%
% Convert to cartesian coordinates
Bxyz = [ ...
    bbr.*sin_theta.*cos_phi + bbt.*cos_theta.*cos_phi - bf.*sin_phi   ... % Bx
    bbr.*sin_theta.*sin_phi + bbt.*cos_theta.*sin_phi + bf.*cos_phi   ... % By
    bbr.*cos_theta          - bbt.*sin_theta                          ... % Bz
    ];  % size n x 3
