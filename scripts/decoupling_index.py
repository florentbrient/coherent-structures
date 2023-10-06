#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 14:44:58 2023

@author: fbrient

Calculate decoupling index

"""
import sys
import netCDF4 as nc
from collections import OrderedDict
import numpy as np
import Constants as CC
import tools as tl


def format1(value):
    return "%13s" % value
def format2(value):
    return "%10.5f" % value

def _findlcl(datatmp,z):
    # Find lcl (from plot_LES_3D)
    THT         = datatmp['THT'][:]  
    P           = datatmp['PABST'][:]  
    TA          = THT*pow(100000./P,-1*CC.gammafix)
    Qtmp        = datatmp['RVT'][:] 

    Q0          = Qtmp[0,:,:]
    TA0         = TA[0,:,:]
    idxlcl,lcl  = tl.findlcl(Q0,TA0,P,z)
    idxlfc,lfc,idxlnb,lnb = tl.findlfc(idxlcl,Qtmp,TA,z)
    
    del Q0,TA0
    del THT,P,Qtmp,TA
    #ss          = Q.shape
    #ZZ          = np.repeat(np.repeat(z[ :,np.newaxis, np.newaxis],ss[1],axis=1),ss[2],axis=2)
    return idxlcl,idxlfc,idxlnb


def createindex(data,z,index='None'):
    # D1: zb-zLCL (Jones et al 2011 ACP) 
    # D2: (Park et al. 2004, De Roode et al 2016)
    if index=='D1':
        idxlcl,idxlnb,idxlfc=_findlcl(data,z)
#        idxlcl,idxlnb,idxlfc = [int(round(np.nanmean(idxlcl))) for ij in [idxlcl,idxlnb,idxlfc]]
        zlcl,zlnb,zlfc = [np.nanmean(z[idxlcl]) for ij in [idxlcl,idxlnb,idxlfc]]
        
        # Define as the first layer where RCT (ql) > epsilon (1e-6 by default)
        rct  = data['RCT'][:] 
        epsilon = 1e-6
        cloudbase,cloudmiddle,cloudtop,zb,zi = tl.cloudinfo(z,rct,epsilon)
        
        zb0     = zb.flatten()
        alt0    = [z[int(ij)] for ij in zb0 if ~np.isnan(ij)]
        #first choice
        zb_mean = np.mean(alt0) #np.nanmean(zz[zb])
        # second choice (from domain mean): More consistent with article
        zb_mean = z[cloudbase]
        
        print('lcl : ',zlcl)
        print('cloud base 1 (from 3D) : ',np.mean(alt0))
        print('cloud base 2 ',z[cloudbase])
        value =  [zb_mean-zlcl]
    elif 'D2' in index:
        # (Phi_cld - Phi_ML)/(Phi_ziplus - Phi_ML)
        # find inversion
        inv    = 'THLM'
        offset = 0.50 # stronger that useful
        idxzi,toppbl,grad = tl.findpbltop(inv,data,var1D,offset=offset)
        
        # Define as the first layer where RCT (ql) > epsilon (1e-6 by default)
        rct  = data['RCT'][:] 
        epsilon = 1e-6
        cloudbase,cloudmiddle,cloudtop,zb,zi = tl.cloudinfo(z,rct,epsilon)
        
        if 'q' in index:
            var_list=['RNPM',]
        elif 't' in index:
            var_list=['THLM',]
        else:
            print('Problem index: ',index)   
        r = {}
        for var in var_list:
            xx_mean = np.nanmean(data[var],axis=(1,2))
            print(xx_mean)
            xx_ML     = np.nanmean(xx_mean[0:cloudbase])
            xx_cld    = np.nanmean(xx_mean[cloudbase:cloudtop])
            xx_ziplus = xx_mean[idxzi+1]
            print(xx_cld,xx_ML,xx_ziplus)
            r[var] = (xx_cld-xx_ML)/(xx_ziplus-xx_ML)
    
        value = [r[var] for var in var_list]        
    else:
       print('Probleme index')
    return value

simu     = sys.argv[1] #'L25.6'
hour     = sys.argv[2] #'003'
cas      = sys.argv[3] #'FIRE'
name     = sys.argv[4] #'V0301'

# open file
file_dir  = "../data/"+cas+"/"+simu+"/"
file      = "sel_"+simu+".1."+name+"."+hour+".nc4"
file      = file_dir + file

# Name of the output file
dir_out   = "../data/d_mean/"+cas+"/";tl.mkdir(dir_out)
fileout0  = dir_out
#fileout0 += 'stats.'+jointyp.join(typ)+'.'+''.join(selSVT)+'.LL.'+filename_input+'.3'+testch+'.txt'
fileout0 += 'mean_'+cas+'_'+simu+'_'+hour+'_'+'DD.txt'

print(file)

# Open data
DATA    = nc.Dataset(file,'r')

# Dimensions
var1D  = ['vertical_levels','W_E_direction','S_N_direction'] #Z,Y,X
namez  = var1D[0]

data1D,nzyx,sizezyx = [OrderedDict() for ij in range(3)]
for ij in var1D:
  data1D[ij]  = DATA[ij][:]/1000. #km
  nzyx[ij]    = data1D[ij][1]-data1D[ij][0]
  sizezyx[ij] = len(data1D[ij])
  
nxny   = nzyx[var1D[1]]*nzyx[var1D[2]] #km^2
ALT    = data1D[namez]
dz     = [0.5*(ALT[ij+1]-ALT[ij-1]) for ij in range(1,len(ALT)-1)]
dz.insert(0,ALT[1]-ALT[0])
dz.insert(-1,ALT[-1]-ALT[-2])
nxnynz = np.array([nxny*ij for ij in dz]) # volume of each level
nxnynz = np.repeat(np.repeat(nxnynz[:, np.newaxis, np.newaxis]
                  ,sizezyx[var1D[1]], axis=1), sizezyx[var1D[2]], axis=2)

z     = ALT*1000. # m

varindex = 'D2t'
tmp   = createindex(DATA,z,index=varindex)
print(tmp)


# Writing text to save index
fileout = fileout0.replace('DD',varindex)
f = open(fileout, 'w')
f.write("File name : "+file+"\n") 
f.write(format1(varindex)+"\n")
for ij in tmp:
  f.write(format2(ij))
f.close()

