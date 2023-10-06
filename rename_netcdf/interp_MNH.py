#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 15:35:48 2023

@author: fbrient
"""

import sys
import netCDF4 as nc
#from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import RectBivariateSpline

file = sys.argv[1]
print('FILE TO INTERPOLATE : ',file)

fl_dia = nc.Dataset(file, 'r+' ,clobber=True)

#print(fl_dia.variables.items())
vars = fl_dia.variables.items()
#print(vars)
for name, variable in vars:
    if name=='RVT':
        dims = variable.dimensions
    if name=='WT':
        dataWT = variable.datatype
    if name=='UT':
        dataUT = variable.datatype


name_new='UTnew'
x = fl_dia.createVariable(name_new, dataUT, dims)
name_new='VTnew'
x = fl_dia.createVariable(name_new, dataUT, dims)
name_new='WTnew'
x = fl_dia.createVariable(name_new, dataWT, dims)

ni_u=np.array(fl_dia['ni_u'][:])
ni_v=np.array(fl_dia['ni_v'][:])
nj_u=np.array(fl_dia['nj_u'][:])
nj_v=np.array(fl_dia['nj_v'][:])

level_w=fl_dia['level_w'][:]
x=fl_dia['W_E_direction'][:]
y=fl_dia['S_N_direction'][:]
z=fl_dia['vertical_levels'][:]

UT = fl_dia['UT'][:].squeeze()
VT = fl_dia['VT'][:].squeeze()
WT = fl_dia['WT'][:].squeeze()


#zi,yi,xi = np.meshgrid(level_w,y,x)
#zj,yj,xj = np.meshgrid(z,y,x)

# Interpolation that work (WT)
print('********** Start WT interpolation **************')
Source = WT[:]
fit    = interp1d(np.array(level_w), np.array(Source), axis=0)
WT_new = fit(np.array(z))

# Interpolation that work (UT)
# 2D -> 3D
print('********** Start UT/VT interpolation **************')
UT_new = np.zeros((len(z),len(y),len(x)))
VT_new = np.zeros((len(z),len(y),len(x)))
for idxz,zz in enumerate(z):
    print('Interpolation... ',idxz,zz)
    UT2D = UT[idxz,:,:]
    interp_spline = RectBivariateSpline(nj_u, ni_u, UT2D)
    UT_new[idxz,:,:] = interp_spline(y, x)
    
    VT2D = VT[idxz,:,:]
    interp_spline = RectBivariateSpline(nj_u, ni_u, VT2D)
    VT_new[idxz,:,:] = interp_spline(y, x)


fl_dia['UTnew'][:] = UT_new
fl_dia['VTnew'][:] = VT_new
fl_dia['WTnew'][:] = WT_new
fl_dia.close()

#fun  = interp2d(nj_u, ni_u, Source2D, kind='linear')
#znew = fun(y, x)



#ii=156;jj=156
#plt.plot(var[0:5,ii,jj],level_w[0:5],'o');plt.plot(Target[0:5,ii,jj],z[0:5],'o',color='Red')


