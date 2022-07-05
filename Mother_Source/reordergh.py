import numpy as np
def expand_out_g_and_h(degree,sh_order,g,h):

    # Expand out g and h for later use. i.e. want length = 232 if degree is 20
    max_gh_len = int( (degree +1)*(degree)/2+1 + degree + 1 )
    # if g and h arrays aren't long enough, pad them to correct size with zeros
    if (max_gh_len > len(g)):
        g = np.append(g,np.zeros(max_gh_len - len(g),dtype='float64'))
    if (max_gh_len > len(h)):
        h = np.append(h,np.zeros(max_gh_len - len(h),dtype='float64'))

    one_float = np.float64(1)  # = 1.0
    two_float = np.float64(2)  # = 2.0
    rec = np.zeros(max_gh_len,dtype='float64')

    for n in range(1, degree +1 +1):
        n2 = np.float64( 2*n-1 )
        n2 = n2 * (n2 - two_float)
        for m in range(1, n +1):
            mn = int( n*(n-1)/2 + m )
            rec[mn] = np.float64( (n-m)*(n+m-2) )/n2

    s = one_float.copy() # = 1.0
    for n in range(2, degree+1 +1):
        mn = int( n*(n-1)/2 + 1 )
        s = s * np.float64( 2*n - 3 )/np.float64( n - 1 )
        p = s.copy() # = a copy of s, not a pointer to s
        g[mn] = g[mn] * s
        h[mn] = h[mn] * s
        for m in range (2, n +1):
            if (m == 2):
                aa = two_float.copy() # = 2.0
            else:
                aa = one_float.copy() # = 1.0
            p = p * np.sqrt( aa*np.float64( n-m+1 )/np.float64( n+m-2 ) )
            mnn = int( mn+m-1 )
            g[mnn] = g[mnn] * p;
            h[mnn] = h[mnn] * p;

    # In use, max index called is k*(k-1)/2 + k , where k = order + 1.
    # so for k = 11, that's index 66, so size 67 (as indexes start at 0 in Python)
    k = sh_order + 1
    max_index = int( k*(k-1)/2 + k )
    if (len(g) > max_index +1 ):  # +1 for index 0
        g   =   g[0:(max_index +1)]
    if (len(h) > max_index +1 ):  # +1 for index 0
        h   =   h[0:(max_index +1)]
    if (len(rec) > max_index +1 ):  # +1 for index 0
        rec = rec[0:(max_index +1)]

    # Done, return arrays back to main code
    return g, h, rec