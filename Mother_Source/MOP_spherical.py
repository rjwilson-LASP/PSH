# RJW's Python code to make IDL and Matlab codes... 
# Use "" for strings in this Python code, but '' for strings in the code I'm writing out
import numpy as np
from datetime import datetime, timezone # for today's date
import pathlib # for full path of this files
import os
from reordergh import expand_out_g_and_h # I need this as a separate file

##############################################################################
def add_commented_line(IDLpro,MATLAB,PYTHON,text_list,standards,IDL_indent):
    for i in range(len(text_list)):
        if len(text_list[i])==0:
            IDLpro.extend(["%s%s"%(IDL_indent,standards['comment_IDLpro'])])
            MATLAB.extend(["%s"%standards['comment_Matlab']])
            PYTHON.extend(["%s"%standards['comment_Python']])
        else:
            IDLpro.extend(["%s%s %s"%(IDL_indent,standards['comment_IDLpro'],text_list[i])])
            MATLAB.extend(["%s %s"%(standards['comment_Matlab'],text_list[i])])
            PYTHON.extend(["%s %s"%(standards['comment_Python'],text_list[i])])
    return IDLpro,MATLAB,PYTHON
##############################################################################
def add_blank_line(IDLpro,MATLAB,PYTHON):
    IDLpro.extend([""])
    MATLAB.extend([""])
    PYTHON.extend([""])
    return IDLpro,MATLAB,PYTHON
##############################################################################
def write_out_original_g_h_by_index(readme,sh_order,g,h):
    readme.append(" The Spherical Harmonic g and h values used for this order %d code are below: "%sh_order)
    readme.append(" ")
    readme.append(" g[i,j] values (nT) used are:")

    cnt = 1
    for i in range(sh_order):
        line = ""
        for j in range(i+2):
            cnt += 1
            if (cnt < len(g)):
                if gh_dp == 1:
                    line += 'g[%2d,%2d] = %10.1f, '%(i+1,j,g[cnt])
                elif gh_dp == 0:
                    line += 'g[%2d,%2d] = %10d, '%(  i+1,j,g[cnt])
                else:
                    print("Should not get here (A)")
                    raise SystemExit 
        readme.append(line)
    readme.append("")
    readme.append(" h[i,j] values (nT) used are:")
    cnt = 1
    for i in range(sh_order):
        line = ""
        for j in range(i+2):
            cnt += 1
            if (j==0):
                if (h[cnt] != 0):
                    print("Error: h[%d,0] not zero"%(i+1))
                    raise SystemExit
                else:
                    if gh_dp == 1:
                        line += '                       '
                    elif gh_dp == 0:
                        line += '                       '
                    else:
                        print("Should not get here (B)")
                        raise SystemExit 
            else:
                if (cnt < len(g)):
                    if gh_dp == 1:
                        line += 'h[%2d,%2d] = %10.1f, '%(i+1,j,h[cnt])
                    elif gh_dp == 0:
                        line += 'h[%2d,%2d] = %10d, '%(  i+1,j,h[cnt])
                    else:
                        print("Should not get here (C)")
                        raise SystemExit 
        readme.append(line)

    return readme
##############################################################################
def add_original_g_h(readme,start_txt,x,gh_dp,yIDL,yMatlab,yPython,lineend,standards,line_start):
    if (yIDL):
        IDLstr = 'd'
    else:
        IDLstr = ''

    if (start_txt != ''):
        if (yPython):
            line_to_add = '%s = np.array(['%start_txt
        else:
            line_to_add = '%s = ['%start_txt
        if (lineend !=''):
            if (len(start_txt)==1):
                line_to_add += "  "
            
            if (yIDL):        # add IDL indent:
                line_to_add = "  "+line_to_add
            if (yMatlab == 1):
                # Uses Matlab comment, but hardcoded
                line_to_add += lineend
                line_to_add += ' %s MATLAB starts at index 1, not 0, so no value on this first line compared to IDL & Python output'%standards['comment_Matlab']
            else:
                if (x[0]!=0):
                    print("Error.... x[0] not 0")
                    raise SystemExit
                if (yIDL):
                    line_to_add += '  '
                    line_to_add += '   0%s                        '%IDLstr +', '+lineend
                elif (yPython):
                    line_to_add += '0                         , '+lineend
                else:
                    line_to_add += '0                         , '+lineend
            readme.extend([line_to_add])
            line_to_add = "%s"%line_start
    else:
        line_to_add = "%s"%line_start


    for i in range(len(x)):
        if (i==0):
            continue
        if (x[i]==0):
            if gh_dp == 1:
                line_to_add += '       0%s  '%IDLstr
            else:
                #line_to_add += '       0%s                        '%IDLstr
                line_to_add += '            0%s                        '%IDLstr
            if (yIDL == 0):
                line_to_add += ' '
        else:
            if gh_dp == 1:
                line_to_add += '%10.1f'%x[i]
            else:
                #line_to_add += '%32.23f'%x[i]
                line_to_add += '%37.23f'%x[i]
            if (yIDL):
                line_to_add +=IDLstr
            else:
                line_to_add += ' '

        if (i < len(x)-1):
            line_to_add += ', '

#        # Exception for Matlab and index 0
#        if (yMatlab == 1) and (i == 0):
#            line_to_add = ' ' * len(line_to_add) # replaces string with spaces, same number of chars

        if (i == len(x)-1):
            line_to_add += ']'
            if (yMatlab):
                line_to_add += ";"
            elif (yPython):
                line_to_add += ", dtype='float64')"
            readme.extend([line_to_add])
            break
        elif (i%4 == 0 ):
            readme.extend([line_to_add+lineend])
            line_to_add = "%s"%line_start
    
    return readme
##############################################################################
def write_out_g_h_rec(degree, sh_order, g, h, IDLpro, MATLAB, PYTHON, standards, IDL_indent):

    g, h, rec = expand_out_g_and_h(degree, sh_order , g, h)

    IDLpro = add_original_g_h(IDLpro,'rec',rec,23,1,0,0,standards['IDLpro_contline'],standards,'  '  )
    MATLAB = add_original_g_h(MATLAB,'rec',rec,23,0,1,0,standards['Matlab_contline'],standards,'    ')
    PYTHON = add_original_g_h(PYTHON,'rec',rec,23,0,0,1,standards['Python_contline'],standards,'    ')

    # Add Blank Line
    (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)

    IDL_indent = '  '
    readme = ['This is the modified g array, not the original g coefficients, and will be further modified.']
    (IDLpro,MATLAB,PYTHON) = add_commented_line(IDLpro,MATLAB,PYTHON,readme,standards,IDL_indent)
    IDLpro = add_original_g_h(IDLpro,'g',g,23,1,0,0,standards['IDLpro_contline'],standards,'  '  )
    MATLAB = add_original_g_h(MATLAB,'g',g,23,0,1,0,standards['Matlab_contline'],standards,'    ')
    PYTHON = add_original_g_h(PYTHON,'g',g,23,0,0,1,standards['Python_contline'],standards,'    ')

    # Add Blank Line
    (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)

    readme = ['This is the modified h array, not the original h coefficients, and will be further modified.']
    (IDLpro,MATLAB,PYTHON) = add_commented_line(IDLpro,MATLAB,PYTHON,readme,standards,IDL_indent)
    IDLpro = add_original_g_h(IDLpro,'h',h,23,1,0,0,standards['IDLpro_contline'],standards,'  '  )
    MATLAB = add_original_g_h(MATLAB,'h',h,23,0,1,0,standards['Matlab_contline'],standards,'    ')
    PYTHON = add_original_g_h(PYTHON,'h',h,23,0,0,1,standards['Python_contline'],standards,'    ')

    # Add Blank Line
    (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)

    return IDLpro,MATLAB,PYTHON
##############################################################################
line_end = '\n'

standards = {\
'comment_IDLpro' : ';', \
'comment_Matlab' : '%', \
'comment_Python' : '#', \
'IDLpro_contline': '$'  , \
'Matlab_contline': '...', \
'Python_contline': '\\' , \
'IDLpro_indent'  : '', \
'Matlab_indent'  : '', \
'Python_indent'  : ''}

#models_to_do = (['jrm09_order10','jrm33_order13','jrm33_order18','jrm33_order30','vip4_order04','vit4_order04','vipal_order05','isaac_order10','o6_order03'])
models_to_do = (['jrm09_order10','jrm33_order13','jrm33_order18','vip4_order04','vit4_order04','vipal_order05','isaac_order10','o6_order03'])

#root_output_dir = '/Users/wilsonr/Documents/JADE/Level2_Processing_Code/IDL/Field_Model/2022'
# Moving to relative paths: Go back one directory.
root_output_dir = os.path.dirname(os.getcwd())
print("Writing out language sub-directories and files to: %s"%root_output_dir)


for model in (models_to_do): # loop for different models?

    for coord in ['rtp','xyz']: # loop for rtp or xyz

        if (model == 'jrm09_order10'):
            planet = 'Jupiter'
            ref = 'Connerney et al. (2018), https://doi.org/10.1002/2018GL077312'
            gh_notes = 'See supplemental online information Table S1, https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2F2018GL077312&file=grl57087-sup-0005-2018GL077312-ds01.txt'
# The above link goes to 10 degrees

            #degree = 20 # this was originally 18, but we only go up to index 66 (size 67), so didn't matter
            degree = 10 # this was originally 18, but we only go up to index 66 (size 67), so didn't matter
            sh_order  = 10
            gh_dp = 1 # number of decimal places given for original g and h
            r_ref =  71492 # 1 Rj in km
            g   = np.array([\
    0., 0.,\
    410244.7, -71498.3,\
    11670.4, -56835.8, 48689.5, \
    4018.6, -37791.1, 15926.3, -2710.5, \
    -34645.4, -8247.6, -2406.1, -11083.8, -17837.2, \
    -18023.6, 4683.9, 16160.0, -16402.0, -2600.7, -3660.7, \
    -20819.6, 9992.9, 11791.8, -12574.7, 2669.7, 1113.2, 7584.9, \
    598.4, 4665.9, -6495.7, -2516.5, -6448.5, 1855.3, -2892.9, 2968., \
    10059.2, 1934.4, -6702.9, 153.7, -4124.2, -867.2, -3740.6, -732.4, -2433.2, \
    9671.8, -3046.2, 260.9, 2071.3, 3329.6, -2523.1, 1787.1, -1148.2, 1276.5, -1976.8, \
    -2299.5, 2009.7, 2127.8, 3498.3, 2967.6, 16.3, 1806.5, -46.5, 2897.8, 574.5, 1298.9 ],dtype='float64')
            h = np.array([\
    0., 0.,\
    0., 21330.5,\
    0., -42027.3, 19353.2,\
    0., -32957.3, 42084.5, -27544.2, \
    0., 31994.5, 27811.2, -926.1, 367.1,\
    0., 45347.9, -749.0, 6268.5, 10859.6, 9608.4,\
    0., 14533.1, -10592.9, 568.6, 12871.7, -4147.8, 3604.4, \
    0., -7626.3, -10948.4, 2633.3, 5394.2, -6050.8, -1526., -5684.2, \
    0., -2409.7, -11614.6, 9287.0, -911.9, 2754.5, -2446.1, 1207.3, -2887.3, \
    0., -8467.4, -1383.8, 5697.7, -2056.3, 3081.5, -721.2, 1352.5, -210.1, 1567.6, \
    0., -4692.6, 4445.8, -2378.6, -2204.3, 164.1, -1361.6, -2031.5, 1411.8, -714.3, 1676.5],dtype='float64')


        elif (model == 'jrm33_order13') or (model == 'jrm33_order18') or (model == 'jrm33_order30'):
            planet = 'Jupiter'
            ref = 'Connerney et al. (2021),  https://doi.org/10.1029/2021JE007055'
            gh_notes = 'See supplemental online information, file 2021JE007055-sup-0002-Supporting Information SI-S02.mod'
            if (model == 'jrm33_order13'):
                degree   = 13
                sh_order = 13
            elif (model == 'jrm33_order18'):
                degree   = 18
                sh_order = 18
            elif (model == 'jrm33_order30'):
                degree   = 30
                sh_order = 30
            gh_dp = 1 # number of decimal places given for original g and h
            r_ref =  71492 # 1 Rj in km
            g   = np.array([\
     0.0,     0.0, \
410993.4,-71305.9, \
 11796.7,-56972.4, 48250.2, \
  2799.3,-37488.4, 15396.8, -1489.8, \
-34402.0, -8080.8, -2440.5,-10848.3,-17919.1, \
-18265.7,  4221.8, 16599.5,-17345.8, -2544.5, -4987.7, \
-20968.0,  9887.6, 12192.4,-12548.7,  2742.2,  1557.6,  8018.2, \
    59.9,  5366.1, -7099.5, -1533.4, -7055.7,  3060.6, -2488.3,  3700.8, \
 10849.5,  1323.8, -6952.2,   -95.0, -4746.6, -1301.7, -4284.6, -1436.7, -3024.6, \
  8914.4, -3506.7,   288.1,   773.6,  3592.7, -3170.8,  1406.3, -1526.6,  1313.3, -2314.6, \
 -2516.5,  1883.2,  2836.1,  4259.4,  3776.7,   764.4,  2112.3,   724.8,  2496.3,   840.1,  1179.4, \
  1311.6,  3056.9, -2060.6,  3550.1,   589.7,   194.8, -1073.3,  -521.7,   685.3,   350.2,  -326.6,  1088.8, \
  2300.5,  1688.4,  -291.7,   483.2,  -581.9,  -818.4, -1095.6, -1987.8,   484.5, -1473.1,  -666.5,  -755.4,  -220.8, \
   751.5,  -456.4,  1160.4, -2307.8,   -77.1,   -16.9,   594.8,  -950.7,  1042.3,  -641.8,    41.6,  -229.1,   -40.3,  -287.8, \
  1164.6, -2478.0,   907.7, -1280.7,   148.5,  1468.8,    -8.1,   604.5,   787.2,   327.7,   111.6,   245.2,   127.3,   425.9,     6.3, \
   469.4, -1929.2,   -38.7,  -815.9,  -775.7,    31.0,  -149.9,   318.2,  -298.3,  -149.3,   241.4,   -69.4,   367.0,   -22.2,    91.8,   -43.7, \
  -809.2,  -888.9,  -604.0,  -339.1,  -350.6,  -503.8,  -235.5,     9.0,  -498.7,  -283.2,   -50.9,  -515.2,  -255.3,    15.6,   -24.7,  -176.5,    49.2, \
 -1238.7,  -190.2,  -101.5,   -58.0,    10.0,  -525.9,   468.4,    88.1,  -385.1,    77.6,    75.7,  -148.9,   201.8,    61.7,  -243.6,   168.2,   -42.2,     7.3, \
  -796.0,   443.6,   266.9,   391.2,   366.2,   175.8,   512.5,   151.0,   214.9,   243.4,   337.5,   -10.2,   -34.9,   161.2,   222.3,  -119.6,    28.0,   132.6,  -180.6, \
  -493.7,   486.7,   -62.0,   539.2,   187.4,   136.8,    -9.4,  -581.2,   291.8,  -198.1,    95.5,  -372.6,   190.7,  -119.6,   -11.9,    74.4,  -391.4,   -33.1,  -118.0,   -57.2, \
  -319.9,   113.6,  -262.2,    94.2,  -306.7,  -139.2,  -291.5,  -716.6,    55.5,  -219.0,  -165.1,  -206.6,  -115.7,  -117.5,   -56.9,    -6.9,    31.4,  -146.9,   453.7,  -256.5,   424.2, \
   -71.5,  -157.0,  -127.9,   -73.9,   -98.0,   222.1,  -299.3,  -124.8,   155.6,    61.0,  -137.5,   -52.1,  -177.1,   -51.6,    44.3,    35.7,  -365.8,    91.4,  -200.6,    28.0,  -113.8,    16.7, \
   -48.6,  -197.9,   -16.4,    84.0,   155.6,   194.9,   -91.6,    55.5,   183.6,    43.1,   -16.4,  -149.6,    33.1,   -39.4,   123.4,   -62.2,   -74.9,   -28.0,   142.1,   147.2,  -307.5,   228.2,  -365.4, \
   -53.7,  -120.7,    23.0,    53.9,   -25.8,  -185.5,    65.2,   -45.7,   -22.3,   -34.5,    10.4,  -158.9,   -15.6,   -59.4,   122.2,  -119.5,   -18.3,   108.2,   -35.6,   -94.6,   147.0,   -88.5,   199.1,  -121.6, \
   128.5,   -31.9,    38.5,   -80.7,  -118.6,  -177.4,    44.8,   -26.8,  -153.3,     7.6,   -49.6,   -48.7,  -122.1,    -1.7,    40.1,   -57.1,   -82.7,   170.5,  -109.1,    28.0,   -45.5,   -42.1,    16.4,   -42.7,   174.8, \
   203.0,    42.3,    17.1,   -79.5,    -1.5,    31.4,    18.8,    21.5,   -94.6,    35.0,   -60.6,   -11.7,   -84.7,    60.1,     5.8,    -4.9,    25.9,   117.3,     9.2,    99.6,  -126.2,    34.2,    25.9,    23.7,    89.8,    63.2, \
    42.0,    69.7,   -12.6,    -7.0,    58.4,    71.8,    30.1,    24.8,     0.7,     8.7,    10.8,   -21.4,    -0.1,    27.6,    17.1,    23.1,    59.9,    89.4,    40.0,   -30.5,   -42.5,    29.0,   -21.3,     5.1,  -205.5,    60.1,  -142.3, \
   -97.3,    29.4,   -17.8,    30.4,    23.1,    16.1,    13.1,     8.6,    28.3,   -14.0,    52.4,    10.3,    17.4,   -33.2,    12.0,    22.1,    15.7,    42.5,    18.8,   -62.5,   -23.7,    15.4,   -38.5,   -10.2,  -106.4,   -25.5,    -2.6,   -86.8, \
   -66.8,   -19.6,    -3.9,    27.8,    -9.3,    -6.2,   -12.2,    -6.2,    16.7,    -9.1,    17.8,    34.2,    -0.4,   -30.6,    -3.2,     2.9,     1.7,   -19.7,     2.2,    16.7,   -20.5,    -3.5,     9.8,     7.6,    12.8,   -35.4,    43.7,    10.5,   -23.6, \
    15.3,   -22.2,     8.2,     5.4,   -10.1,    -4.4,   -13.1,    -6.5,     3.0,     4.2,   -19.7,    13.2,     1.2,     6.9,    -7.7,    -1.8,     6.1,   -27.1,   -11.2,    40.6,     8.1,   -14.3,    23.9,    11.5,     2.0,    -9.8,     5.2,     4.1,    42.8,     5.3, \
    32.2,    -0.1,     8.0,   -12.8,    -0.8,    -4.6,    -3.2,     0.8,    -5.9,     4.4,   -14.3,    -9.4,     8.0,    16.1,    -3.3,     3.6,     3.7,     1.0,   -10.9,     4.3,    21.2,    -5.4,     0.7,    -5.2,   -13.0,     5.0,   -19.2,   -14.3,    32.1,    28.0,    26.9, \
],dtype='float64')
            h = np.array([\
0.0,      0.0, \
0.0,  20958.4, \
0.0, -42549.0, 20221.5, \
0.0, -32890.6, 42518.4,-27397.7, \
0.0,  32452.4, 27438.6,  -501.4, -1325.1, \
0.0,  45363.1,  -826.2,  6000.6, 10568.8, 10091.5, \
0.0,  14016.9,-10119.1,  -294.9, 13948.3, -3686.9,  4783.1, \
0.0,  -7654.8,-11398.6,  2171.0,  5301.6, -6618.1, -1933.8, -5802.9, \
0.0,  -2297.4,-12833.5, 10019.6, -1725.6,  2387.1, -3237.1,   906.1, -3178.3, \
0.0,  -7899.6, -1328.3,  6566.5, -1275.3,  3617.1,  -194.8,  1960.3,   956.9,  1831.7, \
0.0,  -5689.8,  5570.3, -2541.6, -1795.2,  -217.8,  -628.3, -1619.6,   454.6,  -840.0,  1581.4, \
0.0,    549.1,  4527.6, -4223.0, -2761.9,  -476.0, -2714.4, -1349.2, -1649.0, -1471.5,  -583.0,  -499.8, \
0.0,   4204.2,  2228.0, -1901.1, -1271.2,  1073.4, -1111.9,   872.9,  -837.9,  -462.8,   -22.0,   110.3,  -712.7, \
0.0,   4125.6,   -44.0,  -455.4,  2160.6,   255.4,  1105.1,  1214.2,   196.1,  -207.7,  1195.7,   472.4,   721.6,    51.3, \
0.0,   -847.8,  -124.9,   108.8,  1456.2, -1284.6,   896.5,   237.0,   225.8,   -80.9,   696.9,   326.5,   225.9,    71.0,   124.4, \
0.0,  -1146.4,  -580.4,   527.3,   120.0,  -478.9,   148.4,    45.7,  -882.5,  -206.6,   -74.9,    51.0,    36.4,   -44.4,  -285.9,   135.1, \
0.0,   -932.8,  -608.1,   510.2,  -587.4,   -39.4,   358.1,   279.8,   244.5,   487.3,   366.7,    93.9,   282.8,  -267.7,    39.8,  -213.1,    32.4, \
0.0,   -209.2,   -15.8,    54.4,   308.7,   722.0,   233.3,   -82.1,   221.5,   488.6,   206.3,   176.7,   353.4,   204.8,   -29.5,    12.8,   -89.6,  -166.6, \
0.0,    670.9,  -176.1,  -340.3,    37.0,   304.3,  -348.8,  -291.9,   165.6,   360.9,  -119.0,   100.1,    26.9,     1.0,   -60.2,    66.5,   277.8,    29.1,    15.3, \
0.0,     85.8,  -412.1,  -502.0,  -425.1,  -232.2,  -477.0,  -463.5,    73.6,  -228.3,   -96.8,  -256.0,  -228.0,   119.3,   -76.9,  -319.8,  -284.0,   197.9,  -239.6,   302.6, \
0.0,   -297.6,   -13.7,  -118.4,   105.6,   101.9,    45.3,  -213.0,   372.5,  -382.8,   154.2,  -290.3,    17.1,     6.7,    13.2,   -49.3,   -32.2,  -152.8,    62.2,     8.5,    35.0, \
0.0,   -155.8,   312.3,   291.1,   207.9,   165.8,   291.2,   198.8,   442.1,  -104.0,   209.6,   -45.0,  -139.6,   229.4,   -59.8,   -30.0,  -182.0,   162.1,   131.7,  -356.9,   379.5,  -318.6, \
0.0,    -52.5,   180.6,   182.5,   -69.5,   -24.5,    93.2,   188.1,    54.5,   -68.5,    93.8,   -87.3,  -175.3,   202.4,   -49.8,  -129.2,  -241.1,   -63.6,    68.4,   168.1,  -160.0,    62.6,  -117.5, \
0.0,     52.8,    -2.1,   -51.9,   -92.9,    -7.9,    21.1,    16.4,   -51.3,  -135.8,    56.1,  -172.6,    40.8,    21.3,    19.6,   -63.1,  -135.9,    -3.1,    55.9,   -63.5,    -6.3,   188.1,  -119.3,   313.2, \
0.0,     59.4,   -45.4,   -76.0,    13.9,    26.5,    63.0,     9.1,    29.6,    -7.6,    38.3,   -81.6,    63.3,    40.0,    50.3,   -52.5,   -88.2,    99.9,    59.3,   -52.2,   106.4,   -62.1,    34.6,  -171.1,    43.5, \
0.0,     -6.2,   -62.4,   -17.5,    72.8,    -6.1,   -12.2,    28.2,   -18.0,   122.4,   -33.4,    24.0,    19.7,    70.6,    62.2,   -60.1,   -54.4,    54.0,   -48.5,    45.4,    24.9,   -77.3,    36.9,     4.8,    85.7,  -262.0, \
0.0,      5.1,   -50.8,     9.5,    33.9,    -7.2,   -73.0,   -12.9,   -48.1,    55.4,   -60.8,    39.0,    46.9,    43.8,    29.2,    -1.9,     5.8,    36.0,  -113.8,    29.2,   -22.9,   -10.0,   -23.4,   -30.5,    66.7,   126.9,   -43.3, \
0.0,     23.2,     8.1,    16.8,   -44.1,     6.2,   -17.0,   -30.2,     0.7,   -45.5,   -18.1,    13.2,    31.8,    20.7,   -14.4,    30.7,    41.0,    23.0,   -56.4,     1.6,   -14.7,   -16.3,   -31.1,   -32.5,   -49.2,    79.7,   -45.3,    35.8, \
0.0,     -9.5,    37.4,    14.0,   -46.1,    -2.3,    34.7,    -6.7,    21.6,   -33.1,    16.6,    -0.5,   -20.5,    -3.9,   -16.5,    14.5,    20.4,   -14.2,    15.7,   -12.9,   -11.8,   -42.8,    -8.6,    27.4,   -67.5,    -2.2,   -19.2,    28.9,    79.4, \
0.0,    -23.6,    13.4,     0.8,     9.2,    -9.7,    17.2,    11.0,     4.0,    14.7,    13.9,    -2.4,   -28.1,   -21.7,     1.6,     1.5,    -9.8,   -21.9,    26.9,   -15.6,   -12.6,   -28.6,     8.0,    25.6,   -11.0,   -10.3,   -17.9,     2.0,     5.8,    40.6, \
0.0,     -1.0,   -10.2,    -8.8,    29.1,    -2.4,   -10.4,     6.7,    -4.9,    16.9,     2.3,    -1.6,    -3.2,   -14.0,     8.2,     0.5,    -9.8,    -4.4,     6.2,    -5.8,    -4.7,     0.1,    10.2,    -3.4,     9.9,    -2.5,   -11.7,   -10.5,   -21.5,    14.6,    38.3, \
],dtype='float64')


        elif (model == 'isaac_order10'):
            planet = 'Jupiter'
            ref = 'Hess et al. (2017), https://doi.org/10.1553/PRE8s157'
            gh_notes = 'Values are give in nT here, but in the original paper were given to 4 (mostly) to 6 decimal places in units of G.'
            degree = 10
            sh_order  = 10
            gh_dp = 1 # number of decimal places given for original g and h - actually it's not for this, see gh_notes
            r_ref =  71492 # 1 Rj in km
            g   = np.array([\
    0., 0., \
    406650, -71420, \
    -12860, -69810, 38520, \
     -4790, -46420, 28670,  -9340, \
    -22300,  18930,  2760, -13170,  1110, \
     -1650,   7550,  6230,  -1500, -4000, 1420, \
     -6370,   3240,  8020,    270, -4070, 2790,  -680, \
       640,   8920,  3660,  -6980, -1140, 2390, -1230, 310, \
     -5110,   7200,  4530,  -3710,  -200, 3120, -2100, 710, -133, \
       100,    330,  -960,    370,  3390, -740,  -680, 400,  -92,  9, \
       280,    490,  -630,    910,  2850, -380,  -620, 440, -169, 31, -0.8],dtype='float64')
            h = np.array([\
    0., 0.,\
    0,  23530, \
    0, -31700,  7950, \
    0,  -7503, 40310, -36860, \
    0,  20780, 32930, -16950, 5960, \
    0,     30, -3340,   2540, 3490,  -840, \
    0,  10510, -9600,   5030, 1270, -1060,  180, \
    0,   -270, -2970,   7390, 1980, -2120,  350,  19, \
    0,   7810, -7380,   8520, -870, -1680,  720,  -8,  -46, \
    0,   -290,  1890,   1850,  690,   320, -430, 270, -120, 27, \
    0,    230,  1720,   1250, -110,  -110,  -15, 360, -250, 80, -14],dtype='float64')

        elif (model == 'vip4_order04'):
            planet = 'Jupiter'
            ref = 'Connerney et al. (1998), https://doi.org/10.1029/97JA03726'
            gh_notes = 'Original paper is Connerney et al (1998) [https://doi.org/10.1029/97JA03726], however table 3 of Connerney (2007) [https://doi.org/10.1016/B978-044452748-6.00159-0] provides the g and h values to more significant figures, which are used here. i.e. 4.205 G (1998) -> 420543 nT (2007)'
            degree = 4
            sh_order  = 4
            gh_dp = 0 # number of decimal places given for original g and h
            r_ref =  71323 # 1 Rj in km
            g   = np.array([\
    0., 0.,\
    420543, -65920, \
     -5118, -61904, 49690, \
     -1576, -52036, 24386, -17597, \
    -16758,  22210, -6074, -20243, 6643],dtype='float64')
            h = np.array([\
    0., 0.,\
    0,  24992, \
    0, -36052,  5250, \
    0,  -8804, 40829, -31586, \
    0,   7557, 40411, -16597, 3866],dtype='float64')


        elif (model == 'vit4_order04'):
            planet = 'Jupiter'
            ref = 'Connerney (2007), https://doi.org/10.1016/B978-044452748-6.00159-0'
            gh_notes = 'Original paper is Connerney et al (1998) [https://doi.org/10.1029/97JA03726], however table 3 of Connerney (2007) [https://doi.org/10.1016/B978-044452748-6.00159-0] provides the g and h values to more significant figures, which are used here. i.e. 4.205 G (1998) -> 420543 nT (2007)'
            degree = 4
            sh_order  = 4
            gh_dp = 0 # number of decimal places given for original g and h
            r_ref =  71323 # 1 Rj in km
            g   = np.array([\
    0., 0.,\
    428077, -75306, \
     -4283, -59426, 44386, \
      8906, -21447, 21130, -1190, \
    -22925,  18940, -3851,  9926, 1271],dtype='float64')
            h = np.array([\
    0., 0.,\
    0,  24616, \
    0, -50154, 38452, \
    0, -17187, 40667, -35263, \
    0,  16088, 11807,   6195, 12641],dtype='float64')

        elif (model == 'vipal_order05'):
            planet = 'Jupiter'
            ref = 'Hess et al. (2011), https://doi.org/10.1029/2010JA016262'
            gh_notes = 'Values are give to whole nT here, but in the original paper were given to 1 to 4 decimal places in units of G.'
            degree = 5
            sh_order  = 5
            gh_dp = 0 # number of decimal places given for original g and h - actually it's not for this, see gh_notes
            r_ref =  71492 # 1 Rj in km
            g   = np.array([\
    0., 0.,\
    420000 ,-69750, \
     64410, -86720,  95980, \
    -10580, -59000,  63220,  46710, \
    -74660,  32820, -33800,  18260, -14290, \
     -6600,   7370, -17110, -17930,   -770, -7400],dtype='float64')
            h = np.array([\
    0., 0.,\
    0,  19730, \
    0, -40410,  60300, \
    0, -23100,  51600, -11310, \
    0,  32830, -21310,  -6060, -4860, \
    0,  20650, -11670,  -2880,  -500, -22790],dtype='float64')

        elif (model == 'o6_order03'):
            planet = 'Jupiter'
            ref = 'Connerney (1992) (No known DOI)'
            gh_notes = 'This reference does not have a DOI, but we found a NASA ADS page: https://ui.adsabs.harvard.edu/abs/1992pre3.conf...13C/abstract'
            degree = 3
            sh_order  = 3
            gh_dp = 0
            r_ref =  71372 # 1 Rj in km
            g   = np.array([\
    0., 0.,\
    424202, -65929, \
     -2181, -71106, 48714, \
      7565, -15493, 19775, -17958],dtype='float64')
            h = np.array([\
    0., 0.,\
    0,  24116, \
    0, -40304,  7179, \
    0, -38824, 34243, -22439],dtype='float64')


        else:
            print('Did not recognize model %s'%model)
            raise SystemExit

        if planet=='Jupiter':
            OneRExpected = 71492 # km expected by this code.
        else:
            print('Did not recognize planet: %s'%planet)
            raise SystemExit
            



        #outfile_IDLpro = '%s/python_out/jovian_%s_internal_%s.pro'%(root_output_dir,model,coord)
        #outfile_MATLAB = '%s/python_out/jovian_%s_internal_%s.m'%(  root_output_dir,model,coord)
        #outfile_PYTHON = '%s/python_out/jovian_%s_internal_%s.py'%( root_output_dir,model,coord)

        # Make relative paths for system             planet
        outfile_IDLpro = os.path.join(root_output_dir,planet,'IDL'   )
        outfile_MATLAB = os.path.join(root_output_dir,planet,'MATLAB')
        outfile_PYTHON = os.path.join(root_output_dir,planet,'Python')

        if not os.path.exists(outfile_IDLpro):
            os.makedirs(outfile_IDLpro) # Create a new directory because it does not exist
        if not os.path.exists(outfile_MATLAB):
            os.makedirs(outfile_MATLAB) # Create a new directory because it does not exist
        if not os.path.exists(outfile_PYTHON):
            os.makedirs(outfile_PYTHON) # Create a new directory because it does not exist

        outfile_IDLpro = os.path.join(outfile_IDLpro,'jovian_%s_internal_%s.pro'%(model,coord)) 
        outfile_MATLAB = os.path.join(outfile_MATLAB,'jovian_%s_internal_%s.m'%(  model,coord))
        outfile_PYTHON = os.path.join(outfile_PYTHON,'jovian_%s_internal_%s.py'%( model,coord))

        IDLpro = ([])
        MATLAB = ([])
        PYTHON = ([])
        PYTHON.extend(["import numpy as np",""])
#        PYTHON.extend(["from numba import jit","import numpy as np","","@jit(nopython=True)",""])
        
        if (coord == 'rtp'):
            IDLpro.extend(["FUNCTION "   +"jovian_%s_internal_%s, r_rj, colat_rads, elong_rads"%(        model,coord)])
            MATLAB.extend(["function B%s = jovian_%s_internal_%s( r_rj, colat_rads, elong_rads)"%( coord,model,coord)])
            PYTHON.extend(["def "        +"jovian_%s_internal_%s( r_rj, colat_rads, elong_rads):"%(      model,coord)])
        elif (coord == 'xyz'):
            IDLpro.extend(["FUNCTION "   +"jovian_%s_internal_%s, x_rj, y_rj, z_rj"%(        model,coord)])
            MATLAB.extend(["function B%s = jovian_%s_internal_%s( x_rj, y_rj, z_rj)"%( coord,model,coord)])
            PYTHON.extend(["def "        +"jovian_%s_internal_%s( x_rj, y_rj, z_rj):"%(      model,coord)])
        else:
            print("Error: Should not get to this part of code!")
            raise SystemExit

        readme = [
"Code to calculate the %s model of Jupiter's internal magnetic field model"%model.upper(),
"with Degree %d and Order %d."%(degree,sh_order),
"Reference: %s"%ref,
""]

        if (coord == 'rtp'):
            readme.extend([
"Required inputs (System III (1965) Spherical, right handed, and assuming 1 Rj = %d km):"%OneRExpected,
" r_rj       - radial distance, in Rj.",
" colat_rads - colatitude, in radians.                    Value(s) should be 0 <= colat_rads <=  pi.",
" elong_rads - East longitude, right handed, in radians.  Value(s) should be 0 <= elong_rads <= 2pi.",
"",
"Outputs:",
" B - Spherical Magnetic field vector the %s internal magnetic field model, [Br, Btheta, Bphi], units of nT."%model.upper(),
"",
"Usage:",
"For internal field only: B = jovian_%s_internal_rtp(r_rj, colat_rads, elong_rads)"%model])
        elif (coord == 'xyz'):
            readme.extend([
"Required inputs (System III (1965) Cartesian, right handed, and assuming 1 Rj = %d km):"%OneRExpected,
" x_rj       - Jupiter SYSIII right-handed position in x, in Rj.",
" y_rj       - Jupiter SYSIII right-handed position in y, in Rj.",
" z_rj       - Jupiter SYSIII right-handed position in z, in Rj.",
"",
"Outputs:",
" B - Cartesian Magnetic field vector the %s internal magnetic field model, [Bx, By, Bz], units of nT."%model.upper(),
"",
"Usage:",
"For internal field only: B = jovian_%s_internal_xyz(x_rj, y_rj, z_rj)"%model])
        else:
            print("Error: Should not get to this part of code!")
            raise SystemExit

        readme.extend([
"",
"This code was written by Marissa Vogt (mvogt@bu.edu) and Rob Wilson (rob.wilson@lasp.colorado.edu).",
"It is based on a routine originally written by K. Khurana, translated into IDL by Marissa Vogt in 2009."])

        # Special case for jrm09_order10
        if model=='jrm09_order10':
            readme.append("Thanks to Masafumi Imai for providing code for his version of the JRM09 model, which was used to test and validate this code.")
        
        this_file_path = pathlib.Path(__file__).resolve()
        this_file_time = datetime.fromtimestamp(os.path.getmtime(this_file_path), tz=timezone.utc)
         
         
        readme.extend([
"",
"Citation Info:",
" DOI: 10.5281/zenodo.6814109     This doi links to all versions of code at the Github.",
" Github: https://github.com/rjwilson-LASP/PSH",
" Individual versions released on the Github repository can have a different DOI,",
" See the DOI above for a list of DOIs for each specific Github released version."
"",
"Version Info:",
" Last update of this file: %s by user %s. "%(datetime.today(),os.getlogin()),
" This code was re-written/re-formatted by Rob's python code:",
"  %s"%this_file_path,
"  which itself was last updated at UTC %s."%this_file_time.strftime("%Y-%m-%dT%H:%M:%S"),
""])

        readme = write_out_original_g_h_by_index(readme,sh_order,g,h)
        
        # now add to output list
        IDL_indent = '  '
        (IDLpro,MATLAB,PYTHON) = add_commented_line(IDLpro,MATLAB,PYTHON,readme,standards,IDL_indent)

        # Add Blank Line
        (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)

        # IDL only lines for errors
        IDLpro.extend(["  ON_ERROR, 2 %s Exit code if an error in main, don't stop in code - no MATLAB equivalent, just delete line in MATLAB"%standards['comment_IDLpro']])
        IDLpro.extend([""])

        # Start new section
        MATLAB.extend(["%%"]) # New Cell


        if (coord == 'xyz'):
            #For Python, need to convert all to numpy immediately
            PYTHON.extend([
"%s Check inputs x_rj, y_rj and z_rj are all numbers, and convert to numpy doubles here."%standards['comment_Python'],
"try:",
"    x_rj = np.float64(x_rj)",
"    y_rj = np.float64(y_rj)",
"    z_rj = np.float64(z_rj)",
"except Exception as e:",
"    print('ERROR: Inputs must be numeric.')",
"    raise SystemExit",
""])

            #Check Inputs
            IDLpro.extend(["  %s Check inputs are same size."%standards['comment_IDLpro']])
            MATLAB.extend([  "%s Check inputs are same size."%standards['comment_Matlab']])
            PYTHON.extend([  "%s Check inputs are same size."%standards['comment_Python']])

            IDLpro.extend(["  N_input = N_ELEMENTS(x_rj)","  scalar_input = (N_input EQ 1)  %s scalar or not"%standards['comment_IDLpro']])
            MATLAB.extend([  "N_input = numel(x_rj);",      "scalar_input = (N_input == 1); %s scalar or not"%standards['comment_Matlab']])
            PYTHON.extend([  "N_input = x_rj.size",         "scalar_input = (N_input == 1)  %s scalar or not"%standards['comment_Python']])

            # Add Blank Line
            (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)

            IDLpro.extend(["  %s Check inputs x_rj, y_rj and z_rj are all numbers,  and same size (scalar or 1D only)"%standards['comment_IDLpro']])
            MATLAB.extend([  "%s Check inputs x_rj, y_rj and z_rj are all numbers,  and same size (scalar or 1D only)"%standards['comment_Matlab']])
            PYTHON.extend([  "%s Check inputs x_rj, y_rj and z_rj are all arrays of the same size (scalar or 1D only)"%standards['comment_Python']])

            IDLpro.extend(["  IF (N_input NE N_ELEMENTS(y_rj)) THEN MESSAGE,'ERROR: First argument x_rj must be the same size as 2nd argument y_rj'"])
            MATLAB.extend([  "if (N_input ~= "+ "length(y_rj)), "+   "error('ERROR: First argument x_rj must be the same size as 2nd argument y_rj'); end"])
            PYTHON.extend([  "if (N_input != "+ "y_rj.size):",   "    print('ERROR: First argument x_rj must be the same size as 2nd argument y_rj')","    raise SystemExit"])

            IDLpro.extend(["  IF (N_input NE N_ELEMENTS(z_rj)) THEN MESSAGE,'ERROR: First argument x_rj must be the same size as 3rd argument z_rj'"])
            MATLAB.extend([  "if (N_input ~= "+ "length(z_rj)), "+   "error('ERROR: First argument x_rj must be the same size as 3rd argument z_rj'); end"])
            PYTHON.extend([  "if (N_input != "+ "z_rj.size):",   "    print('ERROR: First argument x_rj must be the same size as 3rd argument z_rj')","    raise SystemExit"])

            IDLpro.extend(["  IF (ISA(x_rj, NUMBER=1) EQ 0) OR (SIZE(x_rj, N_DIMENSIONS=1) GT 1) THEN MESSAGE,'ERROR: First  argument x_rj must be a scalar number or 1D "+    "array of numbers'"])
            MATLAB.extend([  "if (~isnumeric(x_rj)) "+     "|| (size(x_rj,2) ~= 1), "+                 "error('ERROR: First  argument x_rj must be a scalar number or 1D column array of numbers'); end"])
            PYTHON.extend([  "if "+                                "(x_rj.ndim "+          "> 1):","    print('ERROR: First  argument x_rj must be a scalar number or 1D "+    "array of numbers')","    raise SystemExit"])

            IDLpro.extend(["  IF (ISA(y_rj, NUMBER=1) EQ 0) OR (SIZE(y_rj, N_DIMENSIONS=1) GT 1) THEN MESSAGE,'ERROR: Second argument y_rj must be a scalar number or 1D "+    "array of numbers'"])
            MATLAB.extend([  "if (~isnumeric(y_rj)) "+     "|| (size(y_rj,2) ~= 1), "+                 "error('ERROR: Second argument y_rj must be a scalar number or 1D column array of numbers'); end"])
            PYTHON.extend([  "if "+                                "(y_rj.ndim "+          "> 1):","    print('ERROR: Second argument y_rj must be a scalar number or 1D "+    "array of numbers')","    raise SystemExit"])

            IDLpro.extend(["  IF (ISA(z_rj, NUMBER=1) EQ 0) OR (SIZE(z_rj, N_DIMENSIONS=1) GT 1) THEN MESSAGE,'ERROR: Third  argument z_rj must be a scalar number or 1D "+    "array of numbers'"])
            MATLAB.extend([  "if (~isnumeric(z_rj)) "+     "|| (size(z_rj,2) ~= 1), "+                 "error('ERROR: Third  argument z_rj must be a scalar number or 1D column array of numbers'); end"])
            PYTHON.extend([  "if "+                                "(z_rj.ndim "+          "> 1):","    print('ERROR: Third  argument z_rj must be a scalar number or 1D "+    "array of numbers')","    raise SystemExit"])

            # Add Blank Line
            (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)

            IDLpro.extend(["  %s Changing inputs to Doubles, and not using input names (so as not to alter inputs, an IDL issue)"%standards['comment_IDLpro']])
            MATLAB.extend([  "%s Changing inputs to Doubles, and not using input names (so as not to alter inputs, an IDL issue)"%standards['comment_Matlab']])
            PYTHON.extend([  "%s Changing inputs to Doubles, and not using input names (so as not to alter inputs, an IDL issue)"%standards['comment_Python']])

            IDLpro.extend(["  x_in = double(x_rj)  %s X in SYSIII, units Rj"%standards['comment_IDLpro']])
            MATLAB.extend([  "x_in = double(x_rj); %s X in SYSIII, units Rj"%standards['comment_Matlab']])
            PYTHON.extend([  "x_in =        x_rj   %s X in SYSIII, units Rj"%standards['comment_Python']]) # Already float64 from above

            IDLpro.extend(["  y_in = double(y_rj)  %s Y in SYSIII, units Rj"%standards['comment_IDLpro']])
            MATLAB.extend([  "y_in = double(y_rj); %s Y in SYSIII, units Rj"%standards['comment_Matlab']])
            PYTHON.extend([  "y_in =        y_rj   %s Y in SYSIII, units Rj"%standards['comment_Python']]) # Already float64 from above

            IDLpro.extend(["  z_in = double(z_rj)  %s Z in SYSIII, units Rj"%standards['comment_IDLpro']])
            MATLAB.extend([  "z_in = double(z_rj); %s Z in SYSIII, units Rj"%standards['comment_Matlab']])
            PYTHON.extend([  "z_in =        z_rj   %s Z in SYSIII, units Rj"%standards['comment_Python']]) # Already float64 from above

            # Adjust distances if needed
            if (OneRExpected != r_ref):
                # Add Blank Line
                (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)

                IDLpro.extend(["  %s Scaling distances since %s expects 1Rj to be %s km (not the %d km that the inputs expect)"%(standards['comment_IDLpro'],model,r_ref,OneRExpected)])
                MATLAB.extend([  "%s Scaling distances since %s expects 1Rj to be %s km (not the %d km that the inputs expect)"%(standards['comment_Matlab'],model,r_ref,OneRExpected)])
                PYTHON.extend([  "%s Scaling distances since %s expects 1Rj to be %s km (not the %d km that the inputs expect)"%(standards['comment_Python'],model,r_ref,OneRExpected)])

                IDLpro.extend(["  r_scale = "+ "double(%f)  /  double(%f)"%( OneRExpected,r_ref),"  x_in = x_in * r_scale" ,"  y_in = y_in * r_scale" ,"  z_in = z_in * r_scale" ])
                MATLAB.extend([  "r_scale = "+ "double(%f)  /  double(%f);"%(OneRExpected,r_ref),  "x_in = x_in * r_scale;",  "y_in = y_in * r_scale;",  "z_in = z_in * r_scale;"])
                PYTHON.extend([  "r_scale = np.float64(%f)/np.float64(%f)"%( OneRExpected,r_ref),  "x_in = x_in * r_scale" ,  "y_in = y_in * r_scale" ,  "z_in = z_in * r_scale" ])

            # Add Blank Line
            (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)

            IDLpro.extend(["  rho_rj_sq = x_in *x_in + y_in *y_in" ])
            MATLAB.extend([  "rho_rj_sq = x_in.*x_in + y_in.*y_in;"])
            PYTHON.extend([  "rho_rj_sq = x_in *x_in + y_in *y_in" ])

            IDLpro.extend(["  r_rj = "+"sqrt(rho_rj_sq + z_in *z_in)" ])
            MATLAB.extend([  "r_rj = "+"sqrt(rho_rj_sq + z_in.*z_in);"])
            PYTHON.extend([  "r_rj = np.sqrt(rho_rj_sq + z_in *z_in)" ])

            # Add Blank Line
            (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)

            IDLpro.extend(["  colat_rads = "  +"acos(z_in /r_rj)" ])
            MATLAB.extend([  "colat_rads = "  +"acos(z_in./r_rj);"])
            PYTHON.extend([  "colat_rads = np.arccos(z_in /r_rj)" ])

            IDLpro.extend(["  elong_rads = "  +"atan( y_in,x_in)" ])
            MATLAB.extend([  "elong_rads = "  +"atan2(y_in,x_in);"])
            PYTHON.extend([  "elong_rads = np.arctan2(y_in,x_in)" ])


#            IDLpro.extend(["  IF scalar_input THEN BEGIN"])
#            MATLAB.extend([  "if scalar_input" ])
#            PYTHON.extend([  "if scalar_input:"])

#            IDLpro.extend(["    IF (elong_rads LT 0)      THEN elong_rads      += 2d*!DPI"])
#            MATLAB.extend([  "    if (elong_rads<0)", "        elong_rads = elong_rads + 2*pi;","    end"])
#            PYTHON.extend([  "    if (elong_rads<0):","        elong_rads += np.float64(2)*np.pi"])

#            IDLpro.extend(["  ENDIF ELSE BEGIN"])
#            MATLAB.extend([  "else" ])
#            PYTHON.extend([  "else:"])

#            IDLpro.extend(["    ind = WHERE(   elong_rads LT 0, NULL = 1)"])
#            MATLAB.extend(["    ind = find(    elong_rads <  0);"])
#            PYTHON.extend(["    ind = np.where(elong_rads <  0)[0]"])

#            IDLpro.extend(["    IF (N_ELEMENTS(ind) NE 0) THEN  elong_rads[ind] += 2d*!DPI"])
#            MATLAB.extend([  "    if numel(ind)",      "        elong_rads(ind) = elong_rads(ind) + 2*pi;","    end"])
#            PYTHON.extend([  "    if (len(ind) != 0):","        elong_rads[ind] += np.float64(2)*np.pi"])

#            IDLpro.extend(["  ENDELSE"])
#            MATLAB.extend([  "end"])

            # Add Blank Line
            (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)
            MATLAB.extend(["%%"]) # New Cell

            readme = [
"######################################################################",
"Start of RTP code.",
"######################################################################"]
            IDL_indent = '  '
            (IDLpro,MATLAB,PYTHON) = add_commented_line(IDLpro,MATLAB,PYTHON,readme,standards, IDL_indent)

        ######################################################################
        # Main RTP Code starts here
        ######################################################################
        # Start new section
        MATLAB.extend(["%%"]) # New Cell

        if (coord == 'rtp'):
            #For Python, need to convert all to numpy immediately
            PYTHON.extend([
"%s Check inputs r_rj, colat_rads and elong_rads are all numbers, and convert to numpy doubles here."%standards['comment_Python'],
"try:",
"    r_rj       = np.float64(    r_rj  )",
"    colat_rads = np.float64(colat_rads)",
"    elong_rads = np.float64(elong_rads)",
"except Exception as e:",
"    print('ERROR: Inputs must be numeric.')",
"    raise SystemExit",
""])


            #Check Inputs
            IDLpro.extend(["  %s Check inputs are same size."%standards['comment_IDLpro']])
            MATLAB.extend([  "%s Check inputs are same size."%standards['comment_Matlab']])
            PYTHON.extend([  "%s Check inputs are same size."%standards['comment_Python']])


            IDLpro.extend(["  N_input = N_ELEMENTS(r_rj)","  scalar_input = (N_input EQ 1)  %s scalar or not"%standards['comment_IDLpro']])
            MATLAB.extend([  "N_input = numel(r_rj);",      "scalar_input = (N_input == 1); %s scalar or not"%standards['comment_Matlab']])
            PYTHON.extend([  "N_input = r_rj.size",         "scalar_input = (N_input == 1)  %s scalar or not"%standards['comment_Python']])

            # Add Blank Line
            (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)

            IDLpro.extend(["  %s Check inputs r_rj, colat_rads and elong_rads are all numbers,  and same size (scalar or 1D only)"%standards['comment_IDLpro']])
            MATLAB.extend([  "%s Check inputs r_rj, colat_rads and elong_rads are all numbers,  and same size (scalar or 1D only)"%standards['comment_Matlab']])
            PYTHON.extend([  "%s Check inputs r_rj, colat_rads and elong_rads are all arrays of the same size (scalar or 1D only)"%standards['comment_Python']])


            IDLpro.extend(["  IF (ISA(r_rj      , NUMBER=1) EQ 0) OR (SIZE(r_rj      , N_DIMENSIONS=1) GT 1) THEN MESSAGE,'ERROR: First  argument    r_rj    must be a scalar number or 1D "+    "array of numbers'"])
            MATLAB.extend(["if (~isnumeric(    r_rj  )) "+       "|| (size(    r_rj  ,2) ~= 1), "+                 "error('ERROR: First  argument    r_rj    must be a scalar number or 1D column array of numbers'); end"])
            PYTHON.extend([  "if "+                                          "(r_rj.ndim "+            "> 1):","    print('ERROR: First  argument    r_rj    must be a scalar number or 1D "+    "array of numbers')","    raise SystemExit"])

            IDLpro.extend(["  IF (ISA(colat_rads, NUMBER=1) EQ 0) OR (SIZE(colat_rads, N_DIMENSIONS=1) GT 1) THEN MESSAGE,'ERROR: Second argument colat_rads must be a scalar number or 1D "+    "array of numbers'"])
            MATLAB.extend(["if (~isnumeric(colat_rads)) "+       "|| (size(colat_rads,2) ~= 1),"+                 " error('ERROR: Second argument colat_rads must be a scalar number or 1D column array of numbers'); end"])
            PYTHON.extend([  "if "+                                      "(colat_rads.ndim "+          "> 1):","    print('ERROR: Second argument colat_rads must be a scalar number or 1D "+    "array of numbers')","    raise SystemExit"])

            IDLpro.extend(["  IF (ISA(elong_rads, NUMBER=1) EQ 0) OR (SIZE(elong_rads, N_DIMENSIONS=1) GT 1) THEN MESSAGE,'ERROR: Third  argument elong_rads must be a scalar number or 1D "+    "array of numbers'"])
            MATLAB.extend(["if (~isnumeric(elong_rads)) "+       "|| (size(elong_rads,2) ~= 1),"+                 " error('ERROR: Third  argument elong_rads must be a scalar number or 1D column array of numbers'); end"])
            PYTHON.extend([  "if "+                                      "(elong_rads.ndim "+          "> 1):","    print('ERROR: Third  argument elong_rads must be a scalar number or 1D "+    "array of numbers')","    raise SystemExit"])

            IDLpro.extend(["  IF (N_input NE N_ELEMENTS(colat_rads)) THEN MESSAGE,'ERROR: First argument r_rj must be the same size as 2nd argument colat_rads'"])
            MATLAB.extend([  "if (N_input ~= length(colat_rads)),"+       " error('ERROR: First argument r_rj must be the same size as 2nd argument colat_rads'); end"])
            PYTHON.extend([  "if (N_input != "+    "colat_rads.size):","    print('ERROR: First argument r_rj must be the same size as 2nd argument colat_rads')","    raise SystemExit"])

            IDLpro.extend(["  IF (N_input NE N_ELEMENTS(elong_rads)) THEN MESSAGE,'ERROR: First argument r_rj must be the same size as 3rd argument elong_rads'"])
            MATLAB.extend([  "if (N_input ~= length(elong_rads)),"+       " error('ERROR: First argument r_rj must be the same size as 3rd argument elong_rads'); end"])
            PYTHON.extend([  "if (N_input != "+    "elong_rads.size):","    print('ERROR: First argument r_rj must be the same size as 3rd argument elong_rads')","    raise SystemExit"])

            # Add Blank Line
            (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)

        IDLpro.extend(["  %s Do this check to be sure that user hasn't got position in km, must be in planetary radii."%standards['comment_IDLpro']])
        MATLAB.extend([  "%s Do this check to be sure that user hasn't got position in km, must be in planetary radii."%standards['comment_Matlab']])

        IDLpro.extend(["  IF scalar_input THEN BEGIN"])
        MATLAB.extend([  "if scalar_input"])

        if (coord == 'xyz'):
            IDLpro.extend(["    IF (    r_rj   LE 0d) OR (    r_rj   GE 200d   ) THEN MESSAGE,'ERROR: First  argument, Position    r_rj   , must be in units of Rj and >0 and <200 only, and not outside that range (did you use km instead?)'"])
            MATLAB.extend(["    if (    r_rj   <= 0 ) || (    r_rj   >= 200),"+       " error('ERROR: First  argument, Position    r_rj   , must be in units of Rj and >0 and <200 only, and not outside that range (did you use km instead?)'),end"])
        elif (coord == 'rtp'):
            IDLpro.extend(["    IF (    r_rj   LE 0d) OR (    r_rj   GE 200d   ) THEN MESSAGE,'ERROR: First  argument, Position    r_rj   , must be in units of Rj and >0 and <200 only, and not outside that range (did you use km instead?)'"])
            MATLAB.extend(["    if (    r_rj   <= 0 ) || (    r_rj   >= 200 ),"+      " error('ERROR: First  argument, Position    r_rj   , must be in units of Rj and >0 and <200 only, and not outside that range (did you use km instead?)'),end"])
            IDLpro.extend(["    IF (colat_rads LT 0d) OR (colat_rads GT !DPI   ) THEN MESSAGE,'ERROR: Second argument, Position colat_rads, must be in units of radians and >= 0 and <=   Pi only, and not outside that range (did you use degrees instead?)'"])
            MATLAB.extend(["    if (colat_rads <  0 ) || (colat_rads >  pi  ),"+      " error('ERROR: Second argument, Position colat_rads, must be in units of radians and >= 0 and <=   Pi only, and not outside that range (did you use degrees instead?)'),end"])
            IDLpro.extend(["    IF (elong_rads LT 0d) OR (elong_rads GT !DPI*2d) THEN MESSAGE,'ERROR: Third  argument, Position elong_rads, must be in units of radians and >= 0 and <= 2*Pi only, and not outside that range (did you use degress instead?)'"])
            MATLAB.extend(["    if (elong_rads <  0 ) || (elong_rads >  pi*2),"+      " error('ERROR: Third  argument, Position elong_rads, must be in units of radians and >= 0 and <= 2*Pi only, and not outside that range (did you use degress instead?)'),end"])

        IDLpro.extend(["  ENDIF ELSE BEGIN"])
        MATLAB.extend(["else"])

        IDLpro.extend(["    min_x = MIN(r_rj, MAX=max_x)"])
        if (coord == 'xyz'):
            IDLpro.extend(["    IF (   min_x   LE 0d) "+       "OR (   max_x   GE 200d   ) THEN MESSAGE,'ERROR: First  argument, Position    r_rj   , must be in units of Rj and >0 and <200 only, and not outside that range (did you use km instead?)'"])
            MATLAB.extend(["    if (min(r_rj ) <= 0 ) || (max(r_rj ) >= 200), "+       "error('ERROR: First  argument, Position    r_rj   , must be in units of Rj and >0 and <200 only, and not outside that range (did you use km instead?)'); end"])
        elif (coord == 'rtp'):
            IDLpro.extend(["    IF (min_x "+  "LE 0d) "+        "OR (max_x "+ "GE 200d   ) THEN "+ "MESSAGE,'ERROR: First  argument, Position    r_rj   , must be in units of Rj and >0 and <200 only, and not outside that range (did you use km instead?)'"])
            MATLAB.extend(["    if (min(r_rj )      <=    0) || (max(r_rj ) >=      200), "+"error('ERROR: First  argument, Position    r_rj   , must be in units of Rj and >0 and <200 only, and not outside that range (did you use km instead?)'); end"])
            IDLpro.extend(["    min_y = MIN(colat_rads, MAX=max_y)"])
            IDLpro.extend(["    IF (min_y LT 0d) "+         "OR (max_y GT !DPI   ) THEN MESSAGE,'ERROR: Second argument, Position colat_rads, must be in units of radians and >= 0 and <=   Pi only, and not outside that range (did you use degrees instead?)'"])
            MATLAB.extend(["    if (min(colat_rads) <     0) || (max(colat_rads) > pi  ), error('ERROR: Second argument, Position colat_rads, must be in units of radians and >= 0 and <=   Pi only, and not outside that range (did you use degrees instead?)'); end"])

            IDLpro.extend(["    min_z = MIN(elong_rads, MAX=max_z)"])
            IDLpro.extend(["    IF (min_z LT 0d) "+         "OR (max_z GT !DPI*2d) THEN MESSAGE,'ERROR: Third  argument, Position elong_rads, must be in units of radians and >= 0 and <= 2*Pi only, and not outside that range (did you use degress instead?)'"])
            MATLAB.extend(["    if (min(elong_rads) <     0) || (max(elong_rads) > pi*2), error('ERROR: Third  argument, Position elong_rads, must be in units of radians and >= 0 and <= 2*Pi only, and not outside that range (did you use degrees instead?)'); end"])

        IDLpro.extend(["  ENDELSE"])
        MATLAB.extend(["end"])

        # Add Blank Line
        (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)
        MATLAB.extend(["%%"]) # New Matlab Cell


        if (coord == 'rtp'):
            IDLpro.extend(["  %s Changing inputs to Doubles, and not using input names (so as not to alter inputs, an IDL issue)"%standards['comment_IDLpro']])
            MATLAB.extend([  "%s Changing inputs to Doubles, and not using input names (so as not to alter inputs, an IDL issue)"%standards['comment_Matlab']])
            PYTHON.extend([  "%s Changing inputs to Doubles, and not using input names (so as not to alter inputs, an IDL issue)"%standards['comment_Python']])

            IDLpro.extend(["  r_rj_dbl       = DOUBLE(      r_rj  )" ,"  colat_rads_dbl = DOUBLE(  colat_rads)" ,"  elong_rads_dbl = DOUBLE(  elong_rads)" ])
            MATLAB.extend([  "r_rj_dbl       = double(      r_rj  );",  "colat_rads_dbl = double(  colat_rads);",  "elong_rads_dbl = double(  elong_rads);"])
            PYTHON.extend([  "r_rj_dbl       =              r_rj"    ,  "colat_rads_dbl =          colat_rads"  ,  "elong_rads_dbl =          elong_rads"  ]) # Already float64 from above # Pointer, rather deep copy, is fine here.

            # Adjust distances if needed - PUT THIS RTP AND XYZ CASE NEXT TO EACH OTHER IN THIS CODE!!!!
            if (OneRExpected != r_ref):
                # Add Blank Line
                (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)

                IDLpro.extend(["  %s Scaling distances since %s expects 1Rj to be %s km (not the %d km that the inputs expect)"%(standards['comment_IDLpro'],model,r_ref,OneRExpected)])
                MATLAB.extend([  "%s Scaling distances since %s expects 1Rj to be %s km (not the %d km that the inputs expect)"%(standards['comment_Matlab'],model,r_ref,OneRExpected)])
                PYTHON.extend([  "%s Scaling distances since %s expects 1Rj to be %s km (not the %d km that the inputs expect)"%(standards['comment_Python'],model,r_ref,OneRExpected)])

                IDLpro.extend(["  r_rj_dbl = r_rj_dbl * "+ "double(%f)  /  double(%f)"%( OneRExpected,r_ref)])
                MATLAB.extend([  "r_rj_dbl = r_rj_dbl * "+ "double(%f)  /  double(%f);"%(OneRExpected,r_ref)])
                PYTHON.extend([  "r_rj_dbl = r_rj_dbl * np.float64(%f)/np.float64(%f)"%( OneRExpected,r_ref)])

        elif (coord == 'xyz'):
            IDLpro.extend(["  %s Code is not using input names (so as not to alter inputs, an IDL issue)"%standards['comment_IDLpro']])
            MATLAB.extend([  "%s Code is not using input names (so as not to alter inputs, an IDL issue)"%standards['comment_Matlab']])
            PYTHON.extend([  "%s Code is not using input names (so as not to alter inputs, an IDL issue)"%standards['comment_Python']])

            IDLpro.extend(["  r_rj_dbl       =  r_rj" ,"  colat_rads_dbl = colat_rads" ,"  elong_rads_dbl = elong_rads" ])
            MATLAB.extend([  "r_rj_dbl       =  r_rj;",  "colat_rads_dbl = colat_rads;",  "elong_rads_dbl = elong_rads;"])
            PYTHON.extend([  "r_rj_dbl       =  r_rj",   "colat_rads_dbl = colat_rads",   "elong_rads_dbl = elong_rads" ]) # Already float64 from above # Pointer, rather deep copy, is fine here.

        # Add Blank Line
        (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)
        MATLAB.extend(["%%"]) # New Matlab Cell


        readme = [
"============",
"Begin hard-coding for %s"%model.upper(),
"Values from %s"%ref]
        if (gh_notes != ''):
            readme.extend([gh_notes])
        readme.extend([
"============"])
        IDL_indent = '  '
        (IDLpro,MATLAB,PYTHON) = add_commented_line(IDLpro,MATLAB,PYTHON,readme,standards,IDL_indent)

        # Add Blank Line
        (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)

        IDLpro.extend(["  %s order = %d  %s degree = order for this code "%(standards['comment_IDLpro'],sh_order,standards['comment_IDLpro'])])
        MATLAB.extend([  "%s order = %d; %s degree = order for this code "%(standards['comment_Matlab'],sh_order,standards['comment_Matlab'])])
        PYTHON.extend([  "%s order = %d; %s degree = order for this code "%(standards['comment_Python'],sh_order,standards['comment_Python'])])

        IDLpro.extend(["  %s k     = order + 1"%( standards['comment_IDLpro'])])
        MATLAB.extend([  "%s k     = order + 1;"%(standards['comment_Matlab'])])
        PYTHON.extend([  "%s k     = order + 1"%( standards['comment_Python'])])
        
        k = sh_order + 1
        IDLpro.extend([ "  k       = %d  %s order + 1 "%(k,standards['comment_IDLpro'])])
        MATLAB.extend([   "k       = %d; %s order + 1 "%(k,standards['comment_Matlab'])])
        #PYTHON.extend([   "k       = %d  %s order + 1 "%(k,standards['comment_Python'])])
        PYTHON.extend([ "%s k     = %d  %s order + 1 "%(standards['comment_Python'],k,standards['comment_Python'])])
        PYTHON.extend([   "k_plus1 = %d  %s k+1, used in for loops when I want to go up to k "%(k+1,standards['comment_Python'])])

        # Add Blank Line
        (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)
        MATLAB.extend(["%%"])


#        readme = [
#"Used IDL function _jovian_jrm09_internal_rjw_setup (below in this file) to copy/paste rec, g and h",
#"rec, g and h never change - so why calculate them each time?"]
        readme = [
"Arrays rec, g and h are processed (depending on degree) but otherwise do not",
"change. So we calculate them once and use in the code. The initial g and h ",
"values are given in the comments at the top of this code, and are reformatted",
"here in to 1D arrays."]
        readme = add_original_g_h(readme,'g',g,gh_dp,0,0,0,'',standards,'     ')
        readme = add_original_g_h(readme,'h',h,gh_dp,0,0,0,'',standards,'     ')
        readme.extend([
"These arrays are then extended and manipulated to make larger g and h arrays, and a rec array."])
        IDL_indent = '  '
        (IDLpro,MATLAB,PYTHON) = add_commented_line(IDLpro,MATLAB,PYTHON,readme,standards, IDL_indent)

        # Copy in the python code used to expand out g, h and rec.
        readme = [
"######################################################################",
"The following is the Python code that was used to expand and process the",
"g and h arrays, and create the rec array for pasting the numbers in to",
"this source code:",
"",
"degree = %d      # = order"%degree,
"g, h, rec = expand_out_g_and_h(degree,order,g,h)"]

        readme.append("----------------------------------------------------------------------")
        file = open("reordergh.py", "r")
        content_list = file.readlines()
        file.close()
        for i in range(len(content_list)):
            readme.append(content_list[i].rstrip())
        readme.append("----------------------------------------------------------------------")

        readme.append("######################################################################")
        IDL_indent = '' # No indent here!
        (IDLpro,MATLAB,PYTHON) = add_commented_line(IDLpro,MATLAB,PYTHON,readme,standards, IDL_indent)

        # Add Blank Line
        (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)

        (IDLpro,MATLAB,PYTHON) = write_out_g_h_rec(degree, sh_order, g, h, IDLpro, MATLAB, PYTHON, standards, IDL_indent)

        readme = [
"============",
"End parts that are hard-coded for %s"%model.upper(),
"============"]
        IDL_indent = '  '
        (IDLpro,MATLAB,PYTHON) = add_commented_line(IDLpro,MATLAB,PYTHON,readme,standards,IDL_indent)


        # Add Blank Line
        (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)
        MATLAB.extend(["%%"])

        IDLpro.extend(["  IF scalar_input THEN BEGIN"])
        MATLAB.extend([  "if scalar_input"])
        PYTHON.extend([  "if scalar_input:"])

        # write out array of zeros and increments up to k
        part = "    a         = [ 0d"
        for i in range(k):
            part += ", 0d"
        part += "]  %s = DBLARR(k+1)"%standards['comment_IDLpro']
        IDLpro.append(part)
        part = "    a         = [ 0"
        for i in range(k-1): # -1 as one less for Matlab
            part += ", 0"
        part += "]; %s = zeros(1, k)"%standards['comment_Matlab']
        MATLAB.append(part)
        part = "    a         = np.array([ 0"
        for i in range(k):
            part += ", 0"
        part += "],dtype='float64') %s = np.zeros(k_plus1,dtype='float64')"%standards['comment_Python']
        PYTHON.append(part)

        part = "    DINDGEN_k = [ 0d"
        for i in range(k):
            part += ",%2dd"%(i+1)
        part += "]  %s = DINDGEN(k+1)"%(standards['comment_IDLpro'])+", done manually for speed"
        IDLpro.append(part)
        part = "    DINDGEN_k = [ 1"
        for i in range(k-1):  # -1 as one less for Matlab
            part += ",%2d"%(i+2)
        part += "]; %s = 1:k"%standards['comment_Matlab']+", done manually for speed"
        MATLAB.append(part)
        part = "    DINDGEN_k = np.array([ 0"
        for i in range(k):
            part += ",%2d"%(i+1)
        part += "],dtype='float64') %s = 0:k"%(standards['comment_Python'])+", done manually for speed"
        PYTHON.append(part)

        IDLpro.extend(["  ENDIF ELSE BEGIN"])
        MATLAB.extend([  "else"])
        PYTHON.extend([  "else:"])

        IDLpro.extend(["    a         = DBLARR(N_input,k+1)"])
        MATLAB.extend(["    a         = zeros( N_input,k  );"])
        PYTHON.extend(["    a         = np.zeros((N_input,k_plus1),dtype='float64')"])

        IDLpro.extend(["    DINDGEN_k = a"])
        MATLAB.extend(["    DINDGEN_k = a;"])
        PYTHON.extend(["    DINDGEN_k = a.copy()"])

        IDLpro.extend(["    FOR i = 0,k DO "+                 "DINDGEN_k[*,i] = i"])
        MATLAB.extend(["    for i = 1:k",             "        DINDGEN_k(:,i) = i;","    end"])
        PYTHON.extend(["    for i in range(k_plus1):","        DINDGEN_k[:,i] = i"])

        IDLpro.extend(["  ENDELSE"])
        MATLAB.extend([  "end"])

        # Add Blank Line
        (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)
        MATLAB.extend(["%%"])

#        IDLpro.extend(["  %s Instead of: a = (1d/r_rj_dbl)^DINDGEN_kp1,  do the da to end of for loop"%standards['comment_IDLpro']])
#        MATLAB.extend([  "%s Instead of: a = (1/r_rj_dbl).^DINDGEN_kp1;, do the da to end of for loop"%standards['comment_Matlab']])
#        PYTHON.extend([  "%s Instead of: a = (1d/r_rj_dbl)^DINDGEN_kp1,  do the da to end of for loop"%standards['comment_Python']])
        IDLpro.extend(["  da = 1d"       + "/r_rj_dbl" ])
        MATLAB.extend([  "da = 1"        +"./r_rj_dbl;"])
        PYTHON.extend([  "da = np.float64(1)/r_rj_dbl" ])

        IDLpro.extend([
"  IF scalar_input THEN BEGIN",
"    a[0] = da",
"    FOR i=1,k DO a[i] = a[i-1]*da",
"  ENDIF ELSE BEGIN %s the following vectorized 2 lines works for scalars, but is slower."%standards['comment_IDLpro'],
"    a[*,0] = da",
"    FOR i=1,k DO a[*,i] = a[*,i-1]*da",
"  ENDELSE"])
        # Different code format
        MATLAB.extend([
"a(:,1) = da.*da;", # MATLAB.extend(["a(:,1) = da.*da; %s in IDL a[0] = da, and for i = 1:k"%standards['comment_Matlab']])
"for i=2:k",
"    a(:,i) = a(:,i-1).*da;",
"end"])
        PYTHON.extend([
"if scalar_input:",
"    a[0] = da",
"    for i in range(1,k_plus1):",
"        a[i] = a[i-1]*da",
"else:",
"    a[:,0] = da",
"    for i in range(1,k_plus1):",
"        a[:,i] = a[:,i-1]*da"])

        # Add Blank Line
        (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)

        IDLpro.extend(["  b = a  * DINDGEN_k" ])
        MATLAB.extend([  "b = a .* DINDGEN_k;"])
        PYTHON.extend([  "b = a  * DINDGEN_k" ])


        # Add Blank Line
        (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)

        IDLpro.extend(["  cos_phi   = "+"cos(elong_rads_dbl)"])
        MATLAB.extend([  "cos_phi   = "+"cos(elong_rads_dbl);"])
        PYTHON.extend([  "cos_phi   = np.cos(elong_rads_dbl,dtype='float64')"])

        IDLpro.extend(["  sin_phi   = "+"sin(elong_rads_dbl)"])
        MATLAB.extend([  "sin_phi   = "+"sin(elong_rads_dbl);"])
        PYTHON.extend([  "sin_phi   = np.sin(elong_rads_dbl,dtype='float64')"])

        IDLpro.extend(["  cos_theta = "+"cos(colat_rads_dbl)"])
        MATLAB.extend([  "cos_theta = "+"cos(colat_rads_dbl);"])
        PYTHON.extend([  "cos_theta = np.cos(colat_rads_dbl,dtype='float64')"])

        IDLpro.extend(["  sin_theta = "+"sin(colat_rads_dbl)"])
        MATLAB.extend([  "sin_theta = "+"sin(colat_rads_dbl);"])
        PYTHON.extend([  "sin_theta = np.sin(colat_rads_dbl,dtype='float64')"])

        IDLpro.extend(["  not_bk = (sin_theta GE 0.00001d)  %s = 1d-5 - also see bk both times below"%standards['comment_IDLpro']])
        MATLAB.extend(  ["not_bk = (sin_theta >= 0.00001 ); %s = 1e-5 - also see bk both times below"%standards['comment_Matlab']])
        PYTHON.extend(  ["not_bk = (sin_theta >= 0.00001 )  %s = 1d-5 - also see bk both times below"%standards['comment_Python']])

        IDLpro.extend(["  IF scalar_input THEN BEGIN"])
        MATLAB.extend([  "if scalar_input"])
        PYTHON.extend([  "if scalar_input:"])

        IDLpro.extend(["    %s bk = (sin_theta LT 0.00001d)  %s bk not needed for scalar"%(standards['comment_IDLpro'],standards['comment_IDLpro'])])
        MATLAB.extend(["    %s bk = (sin_theta <  0.00001 ); %s bk not needed for scalar"%(standards['comment_Matlab'],standards['comment_Matlab'])])
        PYTHON.extend(["    %s bk = (sin_theta <  0.00001 )  %s bk not needed for scalar"%(standards['comment_Python'],standards['comment_Python'])])

        IDLpro.extend(["    zero_array = "+        "0d","    p   = "+        "1d","    d   = 0d"               ,"    bbr = 0d"               ,"    bbt = 0d"               ,"    bbf = 0d"               ,"    x = 0d"               ,"    y = 1d"      ])
        MATLAB.extend(["    zero_array = "+        "0;","    p   = "+        "1;","    d   = 0;"               ,"    bbr = 0;"               ,"    bbt = 0;"               ,"    bbf = 0;"               ,"    x = 0;"               ,"    y = 1;"      ])
        PYTHON.extend(["    zero_array = np.float64(0)","    p   = np.float64(1)","    d   = zero_array.copy()","    bbr = zero_array.copy()","    bbt = zero_array.copy()","    bbf = zero_array.copy()","    x = zero_array.copy()","    y = p.copy()"])



        IDLpro.extend(["  ENDIF ELSE BEGIN"])
        MATLAB.extend(["else"])
        PYTHON.extend(["else:"])

        IDLpro.extend(["    bk = (sin_theta LT 0.00001d)" ])
        MATLAB.extend(["    bk = (sin_theta <  0.00001 );"])
        PYTHON.extend(["    bk = (sin_theta <  0.00001 )" ]) 

        IDLpro.extend(["    zero_array = DBLARR(N_input)"                  ,"    p   = zero_array + 1d"])
        MATLAB.extend(["    zero_array = zeros( N_input,1);"               ,"    p   =         ones( N_input,1);"])
        PYTHON.extend(["    zero_array = np.zeros(N_input,dtype='float64')","    p   = zero_array + np.float64(1)"])
        
        IDLpro.extend(["    d   = zero_array"       ,"    bbr = zero_array"       ,"    bbt = zero_array"       ,"    bbf = zero_array"       ])
        MATLAB.extend(["    d   = zero_array;"      ,"    bbr = zero_array;"      ,"    bbt = zero_array;"      ,"    bbf = zero_array;"      ])
        PYTHON.extend(["    d   = zero_array.copy()","    bbr = zero_array.copy()","    bbt = zero_array.copy()","    bbf = zero_array.copy()"])

        IDLpro.extend(["    x = zero_array"       ,"    y = p  "   +"%s 1s"%standards['comment_IDLpro']])
        MATLAB.extend(["    x = zero_array;"      ,"    y = p; "   +"%s 1s"%standards['comment_Matlab']])
        PYTHON.extend(["    x = zero_array.copy()","    y = p.copy() %s 1s"%standards['comment_Python']])

        IDLpro.extend(["  ENDELSE"])
        MATLAB.extend([  "end"])

        # Add Blank Line
        (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)

        IDLpro.extend(["  FOR m = 1,k DO BEGIN"])
        MATLAB.extend([  "for m = 1:k"])
        PYTHON.extend([  "for m in range(1, k_plus1):"])

        IDLpro.extend(["    bm  = (m NE 1)" ])
        MATLAB.extend(["    bm  = (m ~= 1);"])
        PYTHON.extend(["    bm  = (m != 1)" ])

        IDLpro.extend(["    IF bm THEN BEGIN"])
        MATLAB.extend(["    if bm" ])
        PYTHON.extend(["    if bm:"])

        IDLpro.extend([  "      m_minus_1 = DOUBLE(m - 1)"])
        MATLAB.extend(["        m_minus_1 =        m - 1;"])
        PYTHON.extend(["        m_minus_1 = np.float64(m - 1)"])

        IDLpro.extend([  "      w = x"       ,  "      x = w *cos_phi + y *sin_phi" ,  "      y = y *cos_phi - w *sin_phi" ])
        MATLAB.extend(["        w = x;"      ,"        x = w.*cos_phi + y.*sin_phi;","        y = y.*cos_phi - w.*sin_phi;"])
        PYTHON.extend(["        w = x.copy()","        x = w *cos_phi + y *sin_phi" ,"        y = y *cos_phi - w *sin_phi" ])

        IDLpro.extend(["    ENDIF"])
        MATLAB.extend(["    end"])

        IDLpro.extend(["    q = p"       ,"    z = d"       ,"    bi = zero_array"       ,"    p2 = zero_array"       ,"    d2 = zero_array"       ])
        MATLAB.extend(["    q = p;"      ,"    z = d;"      ,"    bi = zero_array;"      ,"    p2 = zero_array;"      ,"    d2 = zero_array;"      ])
        PYTHON.extend(["    q = p.copy()","    z = d.copy()","    bi = zero_array.copy()","    p2 = zero_array.copy()","    d2 = zero_array.copy()"])

        IDLpro.extend(["    FOR n = m,k DO BEGIN"])
        MATLAB.extend(["    for n = m:k"])
        PYTHON.extend(["    for n in range(m, k_plus1):"])

        IDLpro.extend([  "      mn = "+  "n*(n-1)/2 + m"  ,  "      w  = g[mn]*y + h[mn]*x" ])
        MATLAB.extend(["        mn = "+  "n*(n-1)/2 + m;" ,"        w  = g(mn)*y + h(mn)*x;"])
        PYTHON.extend(["        mn = int( n*(n-1)/2 + m )","        w  = g[mn]*y + h[mn]*x" ])


# Don't need these?
#        MATLAB.extend([  "      %s the next two commented out lines assume scalar inputs"%standards['comment_IDLpro']])
#        MATLAB.extend(["        %s the next two commented out lines assume scalar inputs"%standards['comment_Matlab']])
#        IDLpro.extend([  "      %s IF (abs(p2) LT 1d-38) THEN p2 = 0d %s RJW - Why have these lines?"%(standards['comment_IDLpro'],standards['comment_IDLpro'])]) #note 1d-38 (double), 1e-38 is a float
#        MATLAB.extend(["        %s IF (abs(p2) LT 1e-38) THEN p2 = 0d ;% RJW - Why have these lines?"%(standards['comment_Matlab'],standards['comment_Matlab'])])
#        IDLpro.extend([  "      %s IF (abs(q)  LT 1d-38) THEN q  = 0d %s RJW - Why have these lines?"%(standards['comment_IDLpro'],standards['comment_IDLpro'])])#note 1d_38 (double), 1e-38 is a float
#        MATLAB.extend(["        %s;IF (abs(q)  LT 1e-38) THEN q  = 0d ;% RJW - Why have these lines?"%(standards['comment_Matlab'],standards['comment_Matlab'])])

# Next session is differently formatted for IDL and Matlab, for speed reasons.
        IDLpro.extend([
"      IF scalar_input THEN BEGIN",
"        bbr += b[  n]*w*q",
"        bbt -= a[  n]*w*z",
"        IF bm THEN BEGIN",
"          IF not_bk THEN bi += a[n] * (g[mn]*x-h[mn]*y) * q  %s"%standards['IDLpro_contline'],
"          ELSE           bi += a[n] * (g[mn]*x-h[mn]*y) * z",
"        ENDIF",
"      ENDIF ELSE BEGIN",
"        bbr += b[*,n]*w*q",
"        bbt -= a[*,n]*w*z",
"        IF bm THEN BEGIN",
"          qq = q",
# Fastest why, or not have n_ind, and do N_ELEMENTS(n_ind)
"          ind = WHERE(bk, n_ind, NULL = 1)",
"          IF (n_ind NE 0) THEN qq[ind] = z[ind]",
"          bi += a[*,n] * (g[mn]*x-h[mn]*y) * qq",
"        ENDIF",
"      END"])

        MATLAB.extend([
"        bbr = bbr + b(:,n).*w.*q;",
"        bbt = bbt - a(:,n).*w.*z;",
"        if bm",
"            if scalar_input",
"                if not_bk",
"                    bi = bi + a(  n)  * (g(mn)*x-h(mn)*y)  * q;",
"                else",
"                    bi = bi + a(  n)  * (g(mn)*x-h(mn)*y)  * z;",
"                end",
"            else",
"                qq = q;",
"                ind = find(bk);",
"                if numel(ind)",
"                    qq(ind) = z(ind);",
"                end",
"                bi = bi + a(:,n) .* (g(mn)*x-h(mn)*y) .* qq;",
"            end",
"        end"])

        PYTHON.extend([
"        if scalar_input:",
"            bbr += b[  n]*w*q",
"            bbt -= a[  n]*w*z",
"            if bm:",
"                if not_bk:",
"                    bi += a[n] * (g[mn]*x-h[mn]*y) * q",
"                else:",
"                    bi += a[n] * (g[mn]*x-h[mn]*y) * z",
"        else:",
"            bbr += b[:,n]*w*q",
"            bbt -= a[:,n]*w*z",
"            if bm:",
"                qq = q.copy()",
"                ind = np.where(bk)[0]",
"                if (len(ind) != 0):",
"                    qq[ind] = z[ind]",
"                bi += a[:,n] * (g[mn]*x-h[mn]*y) * qq"])


        IDLpro.extend([  "      xk = rec[mn] %s in IDL 8.4 it's faster to write this to xk, to use below twice.  In Matlab quicker to just use rec(nm)"%standards['comment_IDLpro']])
        PYTHON.extend(["        xk = rec[mn] # faster to write this to xk, to use below twice"])
        IDLpro.extend([  "      dp = cos_theta *z - sin_theta *q - d2*xk"      ])
        MATLAB.extend(["        dp = cos_theta.*z - sin_theta.*q - d2*rec(mn);"])
        PYTHON.extend(["        dp = cos_theta *z - sin_theta *q - d2*xk"      ])


        IDLpro.extend([  "      pm = cos_theta *q                - p2*xk"      ])
        MATLAB.extend(["        pm = cos_theta.*q                - p2*rec(mn);"])
        PYTHON.extend(["        pm = cos_theta *q                - p2*xk"])

        IDLpro.extend([  "      d2 = z"       ,  "      p2 = q"       ,  "      z = dp"       ,  "      q = pm"       ])
        MATLAB.extend(["        d2 = z;"      ,"        p2 = q;"      ,"        z = dp;"      ,"        q = pm;"      ])
        PYTHON.extend(["        d2 = z.copy()","        p2 = q.copy()","        z = dp.copy()","        q = pm.copy()"])

        IDLpro.extend(["    ENDFOR"])
        MATLAB.extend(["    end"])

        IDLpro.extend(["    d = sin_theta *d + cos_theta *p" ,"    p = sin_theta *p" ])
        MATLAB.extend(["    d = sin_theta.*d + cos_theta.*p;","    p = sin_theta.*p;"])
        PYTHON.extend(["    d = sin_theta *d + cos_theta *p" ,"    p = sin_theta *p" ])

        IDLpro.extend([
"    IF bm THEN BEGIN",
"      bi  *= m_minus_1",
"      bbf += bi",
"    ENDIF"])
        MATLAB.extend([
"    if bm",
"        bi  = bi  * m_minus_1;",
"        bbf = bbf + bi;",
"    end"])
        PYTHON.extend([
"    if bm:",
"        bi  *= m_minus_1",
"        bbf += bi"])

        IDLpro.extend(["  ENDFOR"])
        MATLAB.extend([  "end"])

        # Add Blank Line
        (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)

        IDLpro.extend(["  %s br = bbr  %s This doesn't change again"%(standards['comment_IDLpro'],standards['comment_IDLpro'])])
        MATLAB.extend([  "%s br = bbr; %s This doesn't change again"%(standards['comment_Matlab'],standards['comment_Matlab'])])
        PYTHON.extend([  "%s br = bbr  %s This doesn't change again"%(standards['comment_Python'],standards['comment_Python'])])
        IDLpro.extend(["  %s bt = bbt  %s This doesn't change again"%(standards['comment_IDLpro'],standards['comment_IDLpro'])])
        MATLAB.extend([  "%s bt = bbt; %s This doesn't change again"%(standards['comment_Matlab'],standards['comment_Matlab'])])
        PYTHON.extend([  "%s bt = bbt  %s This doesn't change again"%(standards['comment_Python'],standards['comment_Python'])])



        IDLpro.extend(["  IF scalar_input THEN BEGIN"])
        MATLAB.extend(["if scalar_input" ])
        PYTHON.extend(["if scalar_input:"])

        IDLpro.extend(["    IF not_bk THEN"+   " bf = bbf/sin_theta $","    ELSE BEGIN"])
        MATLAB.extend(["    if not_bk" ,"        bf = bbf/sin_theta;" ,"    else"      ])
        PYTHON.extend(["    if not_bk:","        bf = bbf/sin_theta"  ,"    else:"     ])

        IDLpro.extend([  "      IF (cos_theta GE 0d) THEN "+      "bf =  bbf "+              "ELSE "+            "bf = -bbf"])
        MATLAB.extend(["        if (cos_theta >= 0)" ,"            bf =  bbf;"      ,"        else" ,"            bf = -bbf;","        end"])
        PYTHON.extend(["        if (cos_theta >= 0):","            bf =  bbf.copy()","        else:","            bf = np.float(-1)*bbf"])


        IDLpro.extend(["    ENDELSE"])
        MATLAB.extend(["    end"])

        IDLpro.extend(["  ENDIF ELSE BEGIN"])
        MATLAB.extend([  "else" ])
        PYTHON.extend([  "else:"])

        IDLpro.extend([
"    bf = bbf %s set size of array and do the 3rd case"%standards['comment_IDLpro'],
"    ind = WHERE(bk AND (cos_theta LT 0d), NULL = 1)",
"    IF (N_ELEMENTS(ind) NE 0) THEN bf[ind] = -bbf[ind]",
"    ind = WHERE(bk EQ 0, NULL = 1)",
"    IF (N_ELEMENTS(ind) NE 0) THEN bf[ind] =  bbf[ind]/sin_theta[ind]"])

        MATLAB.extend([
"    bf = bbf; %s set size of array and do the 3rd case"%standards['comment_Matlab'],
"    ind = find(bk & (cos_theta <  0));",
"    if numel(ind)",
"        bf(ind) = -bbf(ind);",
"    end",
"    ind = find(~bk); % find(bk == 0)",
"    if numel(ind)",
"        bf(ind) = bbf(ind)./sin_theta(ind);",
"    end"])

        PYTHON.extend([
"    bf = bbf.copy() %s set size of array and do the 3rd case"%standards['comment_Python'],
"    ind = np.where((bk == 1) & (cos_theta < 0))[0]",
"    if (len(ind) != 0):",
"        bf[ind] = -bbf[ind]",
"    ind = np.where(bk == 0)[0]",
"    if (len(ind) != 0):",
"        bf[ind] =  bbf[ind]/sin_theta[ind]"])

        IDLpro.extend(["  ENDELSE"])
        MATLAB.extend([  "end"])

#        IDLpro.extend(["  %s Brtp = [[br],[bt],[bf]]"%( standards['comment_IDLpro'])])
#        MATLAB.extend([  "%s Brtp = [ br , bt , bf ];"%(standards['comment_Matlab'])])

        if (coord == 'rtp'):
            # Add Blank Line
            (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)
            MATLAB.extend(["%%"])

#            IDLpro.extend(["  %s RETURN,Brtp"%standards['comment_IDLpro']])
            IDLpro.extend(["  RETURN,[[bbr],[bbt],[bf]]"])
            MATLAB.extend([  "Brtp = [ bbr , bbt , bf ];"])
            #PYTHON.extend(["return np.transpose(np.array([bbr,bbt,bf]))"])
            # Above output returned a size 3 for scalar, not size 1 by 3, so fix
            PYTHON.extend(["if scalar_input:"])
            PYTHON.extend(["    return             np.array([[bbr,bbt,bf]])"])
            PYTHON.extend(["else:"])
            PYTHON.extend(["    return np.transpose(np.array([bbr,bbt,bf]))"])

        elif (coord == 'xyz'):

            # Add Blank Line
            (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)
            MATLAB.extend(["%%"])

            readme = [
"######################################################################",
"End of RTP code.",
"######################################################################"]
            IDL_indent = '  '
            (IDLpro,MATLAB,PYTHON) = add_commented_line(IDLpro,MATLAB,PYTHON,readme,standards, IDL_indent)

            IDLpro.extend(["  %s Brtp = [[bbr],[bbt],[bf]]"%(  standards['comment_IDLpro'])])
            MATLAB.extend([  "%s Brtp = [ bbr , bbt , bf ];"%( standards['comment_Matlab'])])
            #PYTHON.extend([  "%s Brtp = np.transpose(np.array([bbr , bbt , bf ]))"%(standards['comment_Python'])])
            PYTHON.extend(["%s if scalar_input:"%(standards['comment_Python'])])
            PYTHON.extend(["%s     return             np.array([[bbr,bbt,bf]])"%(standards['comment_Python'])])
            PYTHON.extend(["%s else:"%(standards['comment_Python'])])
            PYTHON.extend(["%s     return np.transpose(np.array([bbr,bbt,bf]))"%(standards['comment_Python'])])


            # Add Blank Line
            (IDLpro,MATLAB,PYTHON) = add_blank_line(IDLpro,MATLAB,PYTHON)
            MATLAB.extend(["%%"])

            IDLpro.extend(["  %s Convert to cartesian coordinates"%standards['comment_IDLpro']])
            MATLAB.extend([  "%s Convert to cartesian coordinates"%standards['comment_Matlab']])
            PYTHON.extend([  "%s Convert to cartesian coordinates"%standards['comment_Python']])
            PYTHON.extend([  "%s Each line is one component, Bx, By then Bz"%standards['comment_Python']])

            IDLpro.extend(["  Bxyz = [ %s"%standards['IDLpro_contline']])
            MATLAB.extend([  "Bxyz = [ %s"%standards['Matlab_contline']])
            #PYTHON.extend([  "Bxyz = np.transpose(np.array([ %s"%standards['Python_contline']])
            # Above output returned a size 3 for scalar, not size 1 by 3, so fix
            PYTHON.extend([  "Bxyz = np.array([ %s"%standards['Python_contline']])

            IDLpro.extend(["    [bbr*sin_theta *cos_phi + bbt *cos_theta *cos_phi - bf *sin_phi], %s %s Bx"%(standards['IDLpro_contline'],standards['comment_IDLpro'])])
            MATLAB.extend(["    bbr.*sin_theta.*cos_phi + bbt.*cos_theta.*cos_phi - bf.*sin_phi   %s %s Bx"%(standards['Matlab_contline'],standards['comment_Matlab'])])
            PYTHON.extend(["    bbr *sin_theta *cos_phi + bbt *cos_theta *cos_phi - bf *sin_phi , %s"%(      standards['Python_contline'])])

            IDLpro.extend(["    [bbr*sin_theta *sin_phi + bbt *cos_theta *sin_phi + bf *cos_phi], %s %s By"%(standards['IDLpro_contline'],standards['comment_IDLpro'])])
            MATLAB.extend(["    bbr.*sin_theta.*sin_phi + bbt.*cos_theta.*sin_phi + bf.*cos_phi   %s %s By"%(standards['Matlab_contline'],standards['comment_Matlab'])])
            PYTHON.extend(["    bbr *sin_theta *sin_phi + bbt *cos_theta *sin_phi + bf *cos_phi , %s"%(      standards['Python_contline'])])

            IDLpro.extend(["    [bbr*cos_theta          - bbt *sin_theta                       ]  %s %s Bz"%(standards['IDLpro_contline'],standards['comment_IDLpro'])])
            MATLAB.extend(["    bbr.*cos_theta          - bbt.*sin_theta                          %s %s Bz"%(standards['Matlab_contline'],standards['comment_Matlab'])])
            PYTHON.extend(["    bbr *cos_theta          - bbt *sin_theta                          %s"%(      standards['Python_contline'])])

            IDLpro.extend(["    ]   %s size n x 3"%standards['comment_IDLpro']])
            MATLAB.extend(["    ];  %s size n x 3"%standards['comment_Matlab']])
            #PYTHON.extend(["    ])) %s size n x 3"%standards['comment_Python']])
            PYTHON.extend(["    ]) %s size 3 x n, or just size 3 if scalar"%standards['comment_Python']]) # since no transpose now

            IDLpro.extend(["  RETURN, Bxyz"])
            #PYTHON.extend([  "return  Bxyz"])
            # Above output returned a size 3 for scalar, not size 1 by 3, so fix
            PYTHON.extend([  "if scalar_input:"])
            PYTHON.extend([  "    return    np.array([Bxyz]) %s size 1 x 3"%standards['comment_Python']])
            PYTHON.extend([  "else:"])
            PYTHON.extend([  "    return np.transpose(Bxyz)  %s size n x 3"%standards['comment_Python']])

            # Matlab doesn't need a return


        # IDL only requires this line to end the function
        IDLpro.extend(["END"])



        if 0: # print to screen
            for i in range(len(IDLpro)):
                print(IDLpro[i])
        if 0: # print to screen
            for i in range(len(MATLAB)):
                print(MATLAB[i])
        if 0: # print to screen
            for i in range(len(PYTHON)):
                print(PYTHON[i])

        if 1: # print to file
            for i in range(len(IDLpro)):
                IDLpro[i] = '%s%s'%(IDLpro[i],line_end) # add line breaks
            file = open(outfile_IDLpro,'w')
            file.writelines(IDLpro)
            file.close()

            for i in range(len(MATLAB)):
                MATLAB[i] = '%s%s'%(MATLAB[i],line_end) # add line breaks
            file = open(outfile_MATLAB,'w')
            file.writelines(MATLAB)
            file.close()

            indent_yet = 0
            for i in range(len(PYTHON)):
                if (indent_yet): # add the initial indent, except first 3 lines
                    PYTHON[i] = '    %s'%PYTHON[i] 
                if (PYTHON[i][0:3] == 'def') and (indent_yet == 0):
                    indent_yet = 1
                PYTHON[i] = '%s%s'%(PYTHON[i],line_end) # add line breaks
            file = open(outfile_PYTHON,'w')
            file.writelines(PYTHON)
            file.close()
