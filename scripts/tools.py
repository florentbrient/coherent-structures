#import netCDF4 as nc
from netCDF4 import Dataset
import numpy as np
import os
import Constants as CC
#import scipy as sp
from scipy import interpolate
from scipy import ndimage
from scipy import integrate
import time
#import pylab as plt
#from string import maketrans
from copy import deepcopy
from collections import OrderedDict
import gc

def mkdir(path):
   try:
     os.mkdir(path)
   except:
     pass

def openfilestat(file):
  ijmin = 4
  f = open(file, 'r')
  tab     = [line.rstrip('\n').split() for line in f]
  nametab = tab[2]
  Nmin    = len(tab[ijmin:])
  Nval    = len(tab[ijmin])  
  #print Nmin,Nval
  f.close()
  return tab,nametab,Nmin,Nval,ijmin


def opennc(path,info,name,hour,vtype='V0301',var=None,nc=None,rmsimu=False):
    file0 = path+"{cas}/{simu}/"
    namefile = "sel_{simu}.1.{vtype}.TTTT.nc4"
    if 'short' in info.keys():
       namefile = 'Short_'+namefile
    if nc is not None:
       namefile = namefile.replace('nc4','nc')
    if rmsimu:
       file0 = file0.replace('{simu}/','')
    file0 += namefile 

    cas   = name
    if 'cas' in info.keys():
       cas  = info['cas']
    simu = info['sens']
    file = file0.format(cas=cas,simu=simu,vtype=vtype).replace('TTTT',hour)
    print('1: ',file)
    Data = Dataset(file)
    print('2: ',file)
    if var is not None:
       Data=Data[var]
       try:
        var  = var.replace('PROC1','DATIM')
        Data = Data[var]
       except:
        pass
    return Data,file,cas,simu

def do_unique(tmp):
    tmp[tmp>0]=1
    return tmp

def resiz(tmp): # should be removed by pre-treatment
    return np.squeeze(tmp)

def replacename(nameobj):
  svt0 = ['SVT001','SVT002','SVT003']
  svt1 = ['SVT004','SVT005','SVT006']
  for ss,svt in enumerate(svt0):
    nameobj=nameobj.replace(svt,svt1[ss])
  return nameobj


def removebounds(tmp):
   tmp = resiz(tmp)
   if len(tmp.shape) == 1:
     tmp = tmp[1:-1]
   elif len(tmp.shape) == 2:
     tmp = tmp[1:-1,1:-1]
   elif len(tmp.shape) == 3:
     tmp = tmp[1:-1,1:-1,1:-1]
   else:
     print('Problem removebounds')
   return tmp

# Find nearrest point 
def near(array,value):
  idx=(abs(array-value)).argmin()
  return idx

# Anomalies for 2D and 3D fields
def anomcalc(tmp,ip=0):
    mean  = np.mean(np.mean(tmp,axis=ip+2),axis=ip+1)
    if ip == 0:
      mean    = np.repeat(np.repeat(mean[ :,np.newaxis, np.newaxis],tmp.shape[ip+1],axis=ip+1),tmp.shape[ip+2],axis=ip+2)
    data = tmp-mean
    return data

def equation(a, b):
    """Return a 1D polynomial."""
    return np.polyval(a, b)

# adjust spines in figures
def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 0))  # outward by 10 points
            #print(dir(spine))
            # FB : Need to update with Python3
            #spine.set_smart_bounds(True)
        else:
            spine.set_color('none')  # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])

def tryopen(vv,DATA):
   try:
     tmp = resiz(DATA[vv][:])
   except:
     print('Error in opening ',vv)
     tmp = None
   return tmp

def _nanargmin(arr, axis):
   try:
     return np.nanargmin(arr, axis)
   except ValueError:
     return np.nan

def _nanargmax(arr, axis):
   try:
     return np.nanargmax(arr, axis)
   except ValueError:
     return np.nan

def suffix_name(subcloud,nocloud):
    if subcloud=='1':
        subcloud=True
        subprefix='.scld'
    else:
        subcloud=False
        subprefix=''

    if nocloud=='1':
        nocloud=True
        nocldprefix='.nocld'
    else:
        nocloud=False
        nocldprefix=''
    return subcloud,subprefix,nocloud,nocldprefix
  
def divergence(field, dx, axis=0):
    "return the divergence of a n-D field"
    #data3 = np.sum(data2,axis=0)
    return np.gradient(field, dx, axis=axis)

def qs(temp,pres,order=0) :
  esat  = CC.es0 * np.exp(-(CC.RLVTT*CC.mv/CC.RU)*(1./temp - 1./CC.T0) )
  qstmp = CC.epsilon * esat/ (pres - esat)
  if order >= 1 :
    dqsP = -1*qstmp/pres # change of qsat for a change of P, with fixed T
    dqsT = CC.RLVTT*CC.epsilon*qstmp/ (CC.RD*temp*temp)# change of qsat for a change of T, with fixed P
  else :
    dqsP = 0; dqsT =0
  if order >=2 :
    d2qsP = 2*qstmp/(pres*pres)  # change of dqsat/dP for a change of P, with fixed T
    d2qsT = CC.RLVTT *CC.epsilon*CC.epsilon*qstmp*(2/temp  -1 )/(CC.RD*temp*temp*temp ) # change of dqsat/dT for a change of T, with fixed P
  else :
    d2qsP = 0; d2qsT =0

  return qstmp,dqsP,dqsT,d2qsP,d2qsT

def createrho(T,P):
   ss   = T.shape
   RR   = np.ones(ss)*CC.RD
   rho  = P/(RR*T)
   return rho

def tht2temp(THT,P):
   ss   = P.shape
   p0   = 100000. #np.ones(ss)*100000. #p0
   exner= np.power(P/p0,CC.RD/CC.RCP)
   temp = THT*exner
   return temp

def findTHV(DATA,offset=0.25):
   RVT = np.nanmean(DATA['RVT'],axis=(1,2,))
   THT = np.nanmean(DATA['THT'],axis=(1,2,))
   try:
     RCT = np.nanmean(DATA['RCT'],axis=(1,2,))
   except:
     RCT = np.zeros(len(THT))
   a1  = 0.61
   THV = THT * (np.ones(THT.shape) +a1*RVT - RCT)
   THVint = integrate.cumtrapz(THV)/np.arange(1,len(THV))
   #offset = 0.25
   DT  = THV[:-1]-(THVint+offset)
   idxpbl = np.argmax(DT>0)
   print(THV,THVint,DT)
   #toppbl  = THV[idxpbl]
   return idxpbl

def findTHV3D(ZZ,THV,offset=0.25):
   y_int = integrate.cumtrapz(THV,x=ZZ,axis=0) #initial?
   ss    = ZZ.shape
   #z     = np.arange(1,ss[0])
   #weight= np.repeat(np.repeat(z[ :,np.newaxis, np.newaxis],ss[1],axis=1),ss[2],axis=2)
   #y_int/= weight
   y_int/= ZZ[1:] #weight
   DT    = THV[1:,:,:]-(y_int+offset)
   #print ZZ[:,15,15]
   #print THV[:,15,15]
   #print y_int[:,15,15]
   #print DT[:,15,15]

   #y_int/= ZZ[1:] #weight
   #print y_int[:,15,15]
   idxzi= _nanargmax(DT>0,axis=0)

   #plt.plot(THV[:,15,15],ZZ[:,15,15],'k')
   #plt.plot(y_int[:,15,15],ZZ[1:,15,15],'r')
   #idxpbl15 = idxpbl[15,15]
   #plt.plot(THV[idxpbl15,15,15],ZZ[idxpbl15,15,15],'bo')
   #plt.show()
   return idxzi

def findrh(RVT,T,P):
    return RVT/qs(T,P)[0]

def findmax(data):
    SVT = np.nanmean(data,axis=(1,2,))
    return np.argmax(SVT)
  
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    #print y, box, y_smooth
    return y_smooth

def _openvar(data,var,var1D):
    if var not in data.variables.keys():
        tmp    = createnew(var,data,var1D)
    else:
        tmp    = data[var][:]

    if tmp is not None:
        tmp = np.nanmean(tmp,axis=(1,2,))
    #except:
    #    tmp = None
    return tmp

def findpbltop(typ,DATA,var1D,idx=None,offset=0.25):
    RVT  = _openvar(DATA,'RVT',var1D)
    RNPM = _openvar(DATA,'RNPM',var1D)
    THT  = _openvar(DATA,'THT',var1D)
    THLM = _openvar(DATA,'THLM',var1D)
    P    = _openvar(DATA,'PABST',var1D)
    data = None
    if typ=='gradRH':
     data   = findrh(RVT,tht2temp(THT,P),P)
     #data    = RVT/qs(tht2temp(THT,P),P)[0]
    elif typ=='gradQT':
     data   = RNPM
    elif typ=='gradTHT':
     data   = THT
    elif typ=='gradTHLM':
     data   = THLM
    elif typ == 'PBLhRH':
     data   = findrh(RVT,tht2temp(THT,P),P)
     idx    = np.argmax(data)
    elif typ == 'THV':
     idx    = findTHV(DATA,offset=offset)
    elif typ == 'THLM' or typ=='THT':
     if typ == 'THLM':
        # May changed ! Be caution
        var1D  = ['vertical_levels','W_E_direction','S_N_direction'] #Z,Y,X
        temp = np.nanmean(createnew(typ,DATA,var1D),axis=(1,2,))
     else:
        temp = THT
     THLMint = integrate.cumtrapz(temp)/np.arange(1,len(temp))
     DT      = temp[:-1]-(THLMint+offset)
     idx     = np.argmax(DT>0)

    del RVT,RNPM,THT,THLM,P

    if 'grad' in typ: # Not the best because of the first layer
     Z   = DATA['ZHAT'][:]/1000.
     grad    = abs(np.gradient(data,Z))
     if idx == None:
       idx  = _nanargmax(grad, axis=0) #local gradient
       #idx  = findTHV(DATA) # grad from THV (lowest)
     dz      = 2
     toppbl  = np.mean(grad[idx-dz:idx+dz])#*data[0]
     print(idx,toppbl,Z[idx])
     #plt.plot(grad,Z);plt.show()
     #plt.plot(data,Z);plt.show()
    else:
     toppbl = np.nan #Z[idx]
     grad   = np.nan #Z[idx]
     
    del data

    return idx,toppbl,grad

# Create new variable
def createnew(vv,DATA,var1D):
    # List of variables available for computation
    vc   = {'THV'  :('THT','RVT','RCT'),\
            'THLM' :('THT','PABST','RCT'),\
            'THS1' :('THLM','RNPM'),\
            'THS2' :('THLM','RVT','RNPM'),\
            'TA'   :('THT','PABST'),\
            'DIVUV':('UT','VT',var1D[1],var1D[2]),\
            'REHU' :('RVT','THT','PABST'),\
            'MSE'  :('RVT','THT','PABST',var1D[0]),\
            'WINDSHEAR' :('UT','VT',var1D[0]),\
            'WINDSPEED' :('UT','VT'),\
            'RHO'  :('THT','PABST'),\
            'W2'   :('WT',),\
            'TKEW' :('WT',),\
            'CLDFR':('RCT',),\
            'TKE'  :('UT','VT','WT'),\
#            'LCL'  :('THT','PABST','RVT','ZHAT'),\
           }
    [vc.update({ij:('THT',var1D[0])}) for ij in ['N','NN','PN','CN']]
    [vc.update({ij:('THT','RVT','RCT',var1D[0])}) for ij in ['NNV']]
    [vc.update({ij:('THT','PABST',var1D[0])}) for ij in ['PRW','LWP','Reflectance']]
    [vc.update({ij:('THT','PABST','RVT',var1D[0])}) for ij in ['LCL','LFC','LNB']]

    data = []
    if vv in vc.keys(): 
      data      = [tryopen(ij,DATA) for ij in vc[vv]]
      if vv == 'THV' or vv == 'NNV': # THV = THT * (1 + 0.61 RVT - RCT)
        a1      = 0.61
        if data[1] is None:
          data[1] = np.zeros(data[0].shape)
        if data[2] is None:
          data[2] = np.zeros(data[0].shape)
        tmp     = data[0] * (np.ones(data[0].shape) +a1*data[1] - data[2] )
      if vv == 'THLM':
         tmp = tryopen('THLM',DATA)
         if tmp is None:
             tmp = data[0] -(\
                     (data[0]/tht2temp(data[0],data[1]))\
                    *(CC.RLVTT/CC.RCP)*data[2])
      if vv == 'THS1':
        # THS1 = theta_l * exp( 5.87 * qt ) with qt = rt / (1 + rt)
        # with theta_l --> THLM and QT --> RNPM
        a1      = 5.87
        tmp     = data[0] * np.exp(a1 * data[1]/(np.ones(data[0].shape)+data[1]) )
      if vv == 'THS2':
        # THS2 = theta_l * exp( 5.87 * qt ) * exp( -0.46 * log (rv/r_star)*qt )
        #                * exp( -0.46 * (qt-qv))
        # with qt = rt / (1 + rt), and theta_l --> THLM
        a1    = 5.87; a2 = -0.46
        rstar = np.ones(data[0].shape)*0.0124
        QV    = data[1]/(np.ones(data[0].shape)+data[2])
        QT    = data[2]/(np.ones(data[0].shape)+data[2])
        tmp   = data[0] * np.exp(a1 * QT)\
                 * np.exp( a2 * np.log(data[1]/rstar) * QT )\
                 * np.exp( a2 * (QT-QV) )
      if vv == 'TA': # TA = THT * (P/P0)^(R/CP)
        tmp   = tht2temp(data[0],data[1])
      if vv == 'CLDFR':
        ss    = data[0].shape
        tmp   = data[0].flatten()
        eps   = 1e-6
        tmp   = np.where(tmp < eps, 0., 1.)
        tmp   = np.reshape(tmp,ss)
      if vv == 'DIVUV':
        # DIVUV = DU/DX + DV/DY
       dx   = data[2][2]-data[2][1]; print(dx)
       dy   = data[3][2]-data[3][1]; print(dy)
       tmpU = divergence(data[0], dx, axis=-1) # x-axis
       tmpV = divergence(data[1], dy, axis=-2) # y-axis
       tmp  = tmpU + tmpV
      if vv == 'WINDSHEAR': # WS = DWIND/DZ
       tmp   = np.squeeze(np.sqrt(data[0]*data[0]+data[1]*data[1]))
       ss    = tmp.shape
       zz    = np.repeat(np.repeat(data[2][ :,np.newaxis, np.newaxis],ss[1],axis=1),ss[2],axis=2)
       print(tmp.shape)
       print(np.gradient(zz)[0].shape)
       tmp   = np.gradient(tmp)[0]/np.gradient(zz)[0]; #print tmp #divergence(tmp, zz, axis=0) #
      if vv == 'WINDSPEED':
       tmp   = np.sqrt(data[0]**2.+data[1]**2.)
      if vv == 'RHO': # 
       TA   = tht2temp(data[0],data[1])
       tmp  = createrho(TA,data[1])
      if vv == 'W2': #
       print(vc[vv])
       tmp  = anomcalc(data[0])**2.
      if vv == 'TKEW':
       tmp  = 0.5*anomcalc(data[0])**2.
      if vv == 'NN' or vv == 'NNV' or vv == 'PN' or vv == 'CN':
        # N2 = G/theta * dTH/dZ
       TT  = data[0] ; ZZ = data[-1]
       if vv=='NNV':
         TT = tmp
       ss    = TT.shape
       zz    = np.repeat(np.repeat(ZZ[ :,np.newaxis, np.newaxis],ss[1],axis=1),ss[2],axis=2)
       #print data,data[0],ss,len(data[0]),np.gradient(data[0])[0]
       #print np.gradient(zz)[0]
       dTHdz = np.gradient(TT)[0]/np.gradient(zz)[0]; print(dTHdz)
       tmp   = CC.RG*dTHdz/TT
       if vv == 'PN':
         tmp   = 2*np.pi/np.sqrt(tmp)
       if vv == 'CN':
         tmp   = np.sqrt(tmp)*zz/(2*np.pi)
      if vv == 'LCL' or vv == 'LFC' or vv == 'LNB': #('THT','PABST','RVT','ZHAT')
       Q   = data[2]
       Q0  = Q[0,:,:]
       P   = data[1]
       TA  = tht2temp(data[0],P) #*pow(100000./P,-1*CC.gammafix)
       TA0 = TA[0,:,:]
       z   = data[3]
       idxtmp,tmp = findlcl(Q0,TA0,P,z)
       if vv == 'LFC' or vv == 'LNB':
         idxlfc,lfc,idxlnb,lnb = findlfc(idxtmp,Q,TA,z)
         if vv == 'LFC':
           tmp = lfc
         if vv == 'LNB':
           tmp = lnb
      if vv == 'REHU': #('RVT','THT','PABST')
        TA  = tht2temp(data[1],data[2]) #*pow(100000./data[2],-1*CC.gammafix)
        tmp = data[0]/(qs(TA,data[2])[0])
      if vv == 'MSE': #(('RVT','THT','PABST','ZHAT'))
        TA  = tht2temp(data[1],data[2])
        ss  = TA.shape
        zz  = np.repeat(np.repeat(data[3][ :,np.newaxis, np.newaxis],ss[1],axis=1),ss[2],axis=2)
        tmp = CC.RCP*TA+CC.RG*zz+CC.RLVTT*data[0]
      if vv == 'TKE':
        tmp = 0.5*(anomcalc(data[0])**2.\
            + anomcalc(data[1])**2\
            + anomcalc(data[2])**2)
      if vv == 'PRW' or vv == 'LWP' or vv == 'Reflectance':
        TA   = tht2temp(data[0],data[1])
        rho  = createrho(TA,data[1])
        name= 'RVT'
        if vv == 'LWP' or vv == 'Reflectance':
          name = 'RCT'
        RCT = tryopen(name,DATA)
        if RCT is not None:
          ss  = RCT.shape 
          tmp = np.zeros((ss[1],ss[2]))
          zz  = np.repeat(np.repeat(data[2][ :,np.newaxis, np.newaxis],ss[1],axis=1),ss[2],axis=2)
          for  ij in range(len(data[2])-1):
            tmp[:,:] += rho[ij,:,:]*RCT[ij,:,:]*(zz[ij+1,:,:]-zz[ij,:,:])
        else:
          tmp = None
        if vv == 'Reflectance' and tmp is not None:
          rho_eau  = 1.e+6 #g/m3
          reff     = 5.*1e-6  #10.e-9 #m #
          g        = 0.85
          tmp      = tmp*1000. #kg/m2  --> g/m2
          tmp      = 1.5*tmp/(reff*rho_eau)
          trans    = 1.0/(1.0+0.75*tmp*(1.0-g))  # Quentin Libois
          tmp      = 1.-trans
    else:
       tmp = None
    return tmp

# Defining object
#def def_object(thrs,nbplus=0,AddWT=0): # old version - update routine that use it (Mars 2023)
def def_object(nbplus=0,AddWT=0):
  # Routine that discriminate object in terms of
  # conditional sampling m (m=2 -> 1 s.t.d.)
  # minimum volume Vmin (Vmin = )
  #thrs   = 1
  #thch   = str(thrs).zfill(2)
  #nbmin  = 100 #100 #1000
  objtyp = ['updr','down','down','down']
  objnb  = ['001' ,'001' ,'002' ,'003' ]
  if nbplus == 1:
     objnb = [str(int(ij)+3).zfill(3) for ij in objnb]
     
  WTchar=['' for ij in range(len(objtyp))]
  if AddWT  == 1:
    WTchar = ['_WT','_WT','_WT','_WT']
  
  typs = [field+'_SVT'+objnb[ij]+WTchar[ij] for ij,field in enumerate(objtyp)]
  print(typs)
  return typs, objtyp

# Offset
def findoffset(var):
    offset = {}
    var1000=['LWP','IWP','RC','RT','RVT','RCT','RNPM']
    for ij in var1000:
       offset[ij] = 1000.
    var100 =['lcc','mcc']
    for ij in var100:
       offset[ij] = 100.
    
    rho0=1.14; RLVTT=2.5e6;RCP=1004.
    offset['E0']=rho0*RLVTT
    offset['Q0']=rho0*RCP
    offset['DTHRAD']=86400.

    off    = 1.
    if var in offset.keys():
       off = offset[var]
    return off

def svttyp(case,sens):
   nbplus     = 0
   svt        = {}
   svt['FIRE']={'All'}
   svt['IHOP']={'IHODC','IHOP5'}
   svt['IHOPNW']={'All'}
   svt['AYOTTE']={'All'}
   svt['FIREWIND']={'All'}
   svt['FIRENOWIND']={'All'}
   if case in svt.keys():
     if sens in svt[case] or 'All' in svt[case]:
       nbplus=1
   return nbplus

# Find hatch type for figure
def findhatch(objname,typ='hatch'):
   if objname != 'All': #is not 'All':
     objsplit = objname.split('_')
     objtyp   = objsplit[0]
     svt      = objtyp+objsplit[1]
     if len(objsplit)>2:
       svt+=objsplit[2]
   else:
     svt      = objname

   keys  = OrderedDict()
   keys  = ['updrSVT001WT','downSVT001WT','downSVT002WT','downSVT003WT','downSVT003','All']
   keysd = [string.translate(string.maketrans('123','456')) for string in keys] #({'1': '4', '2': '5', '3': '6'}))


   hatch={} #hatch = ['//','.','\\','+']
   hatch.update(dict.fromkeys([keys[0],keysd[0]],{'hatch':'//' ,'color':'Red'}))
   hatch.update(dict.fromkeys([keys[1],keysd[1]],{'hatch':'.'  ,'color':'Purple'}))
   hatch.update(dict.fromkeys([keys[2],keysd[2]],{'hatch':'+'  ,'color':'Green'}))
   hatch.update(dict.fromkeys([keys[3],keysd[3]],{'hatch':'\\' ,'color':'Blue'}))
   hatch.update(dict.fromkeys([keys[4],keysd[4]],{'hatch':'\\' ,'color':'Blue'}))
   hatch.update(dict.fromkeys([keys[5],keysd[5]],{'hatch':'-'  ,'color':'Black'}))
   #hatch.update([keys[1]]={'hatch':'.' ,'color':'Purple'}
   #hatch.update([keys[2]]={'hatch':'+','color':'Green'}
   #hatch.update([keys[3]]={'hatch':'\\' ,'color':'Blue'}
   #hatch.update([keys[4]]={'hatch':'\\' ,'color':'Blue'}
   #hatch.update(dict.fromkeys(['b', 'e'], 20))
   return hatch[svt][typ]

def cond2plot(case,prefix,boxch=False):
   boxes = None
   xydef = ([0,200],[10,10])
   if boxch:
     offbox = 15
     boxes  = [0,200,0,200] # x1,x2,y1,y2
     boxes  = [ij+offbox for ij in boxes]

   xycas = {}
#   xycas['IHOP']  = { '006':([150,350],[440,440]), # '006':([50,250],[410,410]
#                      '008':([50,250],[410,410])} #[[x1,x2],[y1,y2]]
   xycas['BOMEX'] = { '006':([300,420],[459,459]), 
                      '008':([160,320],[455,455]),  
                      '010':([420,510],[455,455]),  #([270,350],[474,474]),
                    } 
   xycas['FIRE'] = {'012':([420,600],[455,455]), #([270,350],[474,474]),
                    '021':([100,190],[390,390]), 
                    '006':([100,190],[390,390]), 
                    } 
   if case == 'IHOP':
     #xydef=([150,350],[440,440])
     #xydef=([440,440],[150,350])
     #xydef=([140,140],[150,350])
     #xydef=([150,350],[240,240])
     xydef=([100,300],[150,150])
   try:
     xy = np.array(xycas[case][prefix])
   except:
     print('Absence in xy for ',case,prefix,' default used')
     xy = np.array(xydef)
   print(type(np.array(xycas['FIRE']['012'])))
   print(type(xy))
   #stop
   return boxes,xy

# find subdomain
def findsubdom(cas,hour):
   subdom = {}
   subdom = {'IHOP':  {'006':([40,200],[280,440])},\
             'BOMEX': {'006':([40,200],[280,440]),\
#                       '008':([40,200],[280,440])}\
#                       '008':([200,360],[280,440])},\
#                       '008':([240,360],[300,420])},\ default Ru0NW
                       '008':([260,380],[240,360]),\
                       },\

# No clouds
#                        '008':([260,420],[40,200])},\
             'FIRE': {'003':([220,380],[160,320]),\
                      '006':([220,380],[160,320]),\
                      '009':([220,380],[160,320]),\
                      '012':([160,320],[120,280]),\
                      '015':([180,340],[40,200]),\
                      '018':([160,320],[60,220]),\
                      '021':([120,280],[220,380]),\
                      '024':([40,200],[200,360])\
                      },\
             'ARMCU':  {'002':([40,200],[280,440]),\
                        '004':([40,200],[280,440]),\
                        '006':([40,200],[280,440]),\
                        '008':([40,200],[280,440]),\
                        '010':([40,200],[280,440])\
                 }
             }
   yydef = ([40,120],[200,380])
   if cas in subdom.keys():
       xx = subdom[cas]
       if hour in xx.keys():
           yy = np.array(xx[hour])
   else:
       yy = np.array(yydef) # by default
   return yy  

# Infos for figures
def infosfigures(cas, var, mtyp='Mean',relative=False):
   cmaps   = {'Mean':'Greys_r','Anom':'RdBu_r'} # by default
#   cmaps   = {'Mean':'PuBu_r','Anom':'Greys_r'}
   cmaps   = {'Mean':'Greys_r','Anom':'Greys_r'}

   # Switch color label
   switch = []
   # switch  = ['RNPM','RVT']
   if var in switch:
       cmaps['Anom'] = cmaps['Anom'].split('_')[0]
   # Grey color
   greyvar = ['Reflectance']
   if var in greyvar:
     cmaps['Mean'] = 'Greys_r'

   zmin  = 0
   # km
   zmax  = {'IHOP':2, 'FIRE':1.0, 'BOMEX':2, 'ARMCU':2}
   if relative:
       zmax = {'IHOP':2, 'FIRE':1.5, 'BOMEX':3, 'ARMCU':3}

   vrang = {}
   vrang['Mean'] = {'WT':[-5.0,5.0,0.2],
                    'LWP':[0.00,0.14,0.002],
                    'RNPM'  :[0.005,0.01,0.0002],
                    'WINDSHEAR' :[0.,1.,0.02],
                    'REHU':[0.5,0.8,0.01],
                    'Reflectance' :[0.1,1.,0.01],
                    'THLM' : [298,306,2],
                    'RCT' : [0,0.0007,1e-5],
                    'DIVUV': [-0.006,0.009,0.003],
                    }
   vrang['Anom'] = {'WT'   :[-0.8,0.8,0.05],
                    'DIVUV':[-0.05,0.05,0.005],
                    'THV'  :[-1.0,0.6,0.02],
                    'THLM' :[-1.0,1.0,0.02],
                    'RNPM' :[-0.002,0.002,0.0001],
                    'REHU' :[-0.15,0.15,0.005],
                    'PABST':[-2,2,0.1]
                    }

   # Modified range for some case
   modif = {}
   modif['BOMEX'] = {'THV':0.2,'PABST':0.2}
   modif['FIRE']  = {'RNPM':0.2}

   zminmax = None
   if cas in zmax.keys():
     zminmax=[zmin,zmax[cas]]

   levels = None
   if var in vrang[mtyp].keys():
     vmin,vmax,vdiff = vrang[mtyp][var][:]
     mr = 1
     if cas in modif.keys():
         if var in modif[cas].keys():
             mr=modif[cas][var]
     print('Range : ',vmin,vmax,vdiff,var,findoffset(var),mr)
     vmin,vmax,vdiff = [ij*findoffset(var)*mr for ij in (vmin,vmax,vdiff)]
     nb              = abs(vmax-vmin)/vdiff+1
     levels          = [vmin + float(ij)*vdiff for ij in range(int(nb))]

   infofig = {}
   infofig['cmap']    = cmaps #[mtyp]
   infofig['zminmax'] = zminmax
   infofig['levels']  = levels
   return infofig
   

# Axis
def selectdata(zyx,data,slct=None,avg=0):
   datatmp = data 
   if slct is not None:
       
     if slct[0]==0: #xy
       xp,yp = 0,0
       xy    = slct[1]
       if xy[0,0]==xy[0,1]: #same x
         xp = avg
       if xy[1,0]==xy[1,1]: #same y
         yp = avg

       xy    = selectdiag(zyx,xy[:,0],xy[:,1])
       axis  = [zyx[0],np.arange(xy.shape[1])]

       tmp  = np.zeros((len(axis[0]),len(axis[1])))
       for ij in range(xy.shape[1]):
        tmp0      = datatmp[:,xy[1,ij]-yp:xy[1,ij]+yp+1:,xy[0,ij]-xp:xy[0,ij]+xp+1]
        tmp0      = resiz(tmp0)
        if len(tmp0.shape)>1:
          tmp[:,ij] = np.nanmean(tmp0,axis=1)
        else:
          tmp[:,ij] = tmp0
          
     #if slct[0]==2: #useless
     #    xy    = slct[1]
     #    idxx  = np.arange(xy[1,0],xy[1,1])
     #    idxy  = np.arange(xy[0,0],xy[0,1])
     #    axis = [zyx[0],zyx[1][idxy],zyx[2][idxx]]
     #    print('axis 0 ',axis)
     #    idxx, idxy = np.meshgrid(idxx, idxy)
     #    datatmp  = datatmp[:,idxy, idxx]
     # #print('SIZE ',np.shape(datatmp),axis)


     elif slct[0]==1: # or slct[0]==2: #zview
       zp   = avg
       zz   = slct[-1]
       #print zz,slct
       #if slct[0]!=2:
       axis = [zyx[1],zyx[2]]
       #else:
       #axis = [axis[1],axis[2]]
       #  print('axis 1 ',axis)
       #if slct[0]==2:
       #  axis  = [xy[:,0],xy[:,1]]  
       if len(datatmp.shape)==3:
         tmp0 = datatmp[zz-zp:zz+zp+1,:,:]
         tmp  = np.nanmean(tmp0,axis=0)
       else:
         tmp  = datatmp
     else: 
       print('Problem with slct in selectdata')
   #stop
   return tmp,axis

def gridsize(DATA,var1D):
#var1D  = ['vertical_levels','W_E_direction','S_N_direction'] #Z,Y,X
    data1D,nzyx,sizezyx = [OrderedDict() for ij in range(3)]
    for ij in var1D:
        data1D[ij]  = DATA[ij][:]/1000. #km #tl.removebounds(DATA[ij][:])
        nzyx[ij]    = data1D[ij][1]-data1D[ij][0]
        sizezyx[ij] = len(data1D[ij])

    nxny   = nzyx[var1D[1]]*nzyx[var1D[2]] #km^2
    ALT    = data1D[var1D[0]]
    dz     = [0.5*(ALT[ij+1]-ALT[ij-1]) for ij in range(1,len(ALT)-1)]
    dz.insert(0,ALT[1]-ALT[0])
    dz.insert(-1,ALT[-1]-ALT[-2])
    nxnynz = np.array([nxny*ij for ij in dz]) # volume of each level
    nxnynz = np.repeat(np.repeat(nxnynz[:, np.newaxis, np.newaxis]
                             , sizezyx[var1D[1]], axis=1), sizezyx[var1D[2]], axis=2)

    return data1D,ALT,dz,nzyx,nxnynz,sizezyx
     

def selectdiag(zyx,xy1,xy2):
    # diagonal that start on (x1,y1) and ends on (x2,y2)
    # xy1 : [x1,y1]
    # make diagonal
    #print  xy2[0]==xy1[0],xy2[0],xy1[0],xy2[1],xy1[1]
    xx   = np.mgrid[xy1[0]:xy2[0]:len(zyx[-1])*2*1j] # a lot of points
    if xy2[0]==xy1[0]:
      aa = 0
      yy = np.mgrid[xy1[1]:xy2[1]:len(zyx[-2])*2*1j] # a lot of points
    else:
      aa   = float(xy2[1]-xy1[1])/float(xy2[0]-xy1[0])
      bb   = xy1[1]-aa*xy1[0]
      yy   = aa*xx+bb
    #print aa
    xy    = np.zeros((2,len(xx)),dtype='int')#*np.nan
    xorig = np.arange(len(zyx[2]),dtype='int')
    yorig = np.arange(len(zyx[1]),dtype='int')
    ik    = -1
    #print xx,yy
    for ij in range(len(xx)):
      ik += 1
      xy[0,ik] = xorig[near(xorig,xx[ij])] #x
      #print yorig,ij,yy[ij]
      xy[1,ik] = yorig[near(yorig,yy[ij])] #y
      if ik>0:
        if (xy[0,ik]==xy[0,ik-1]) and (xy[1,ik]==xy[1,ik-1]):
          ik = ik-1
    xy = xy[:,:ik+1]
    return xy

def findbasetop(data,epsilon):
    # Find cloud base and cloud top
    # Define as the first layer where RCT (ql) > epsilon
    cloudbase = np.nan; cloudtop = np.nan;
    for ij in range(len(data)):
      if np.isnan(cloudbase)  and data[ij] > epsilon:
         cloudbase = ij
      if ~np.isnan(cloudbase) and np.isnan(cloudtop)  and data[ij] < epsilon:
         cloudtop  = ij
    #print 'base,top : ',cloudbase,cloudtop,data*1000.
    return cloudbase,cloudtop


def findbasetop2(data,epsilon):
    # Find cloud base and cloud top
    # Define as the first layer where RCT (ql) > epsilon
    cloudbase = np.nan; cloudtop = np.nan;
    idx       = np.where(data>epsilon)[0]
    if len(idx):
      #print idx
      cloudbase = idx[0]
      cloudtop  = idx[-1]
    return cloudbase,cloudtop

def cloudinfo(z,rct,epsilon):
    # Find cloud base and cloud top and middle of cloud
    sz  = rct.shape

    # Choice the right computation
    # 1: the first layer > epsilon using the mean profile
    # 2: Every layer > epsilon, and computing the mean of it
    # 3: The minimum between the two of them
    baseandtop = 1

    # First computation : Use the mean 
    meanql             = np.mean(np.mean(rct,axis=2),axis=1)
    cloudbase,cloudtop = findbasetop2(meanql,epsilon)
    print('Compute first :')
    print('cloudbase : ij='+str(cloudbase)+' for z='+str(z[cloudbase]))
    print('cloudtop  : ij='+str(cloudtop) +' for z='+str(z[cloudtop]) )

    # Second computation : Use all data and compute its mean
    cloudbase0 = np.zeros(sz[1:])*np.nan
    cloudtop0  = np.zeros(sz[1:])*np.nan
    idxtab     = np.argwhere(rct[:]>epsilon)
    
    for idx in idxtab:
      cloudtop0[idx[1],idx[2]] = idx[0]
      if np.isnan(cloudbase0[idx[1],idx[2]]):
        cloudbase0[idx[1],idx[2]] = idx[0]

    cloudbase2 = int(np.nanmean(cloudbase0[:]))   # real -> int
    cloudtop2  = int(np.nanmean(cloudtop0[:]))+1  # real -> int
    print ('Compute second :')
    print ('cloudbase : ij='+str(cloudbase2)+' for z='+str(z[cloudbase2]))
    print ('cloudtop  : ij='+str(cloudtop2) +' for z='+str(z[cloudtop2]) )
    
    # Choice the second computation
    if baseandtop == 2:
       cloudbase = cloudbase2
       cloudtop  = cloudtop2
    elif baseandtop == 3:
       cloudbase = min(cloudbase,cloudbase2)
       cloudtop  = min(cloudtop ,cloudtop2) 

    cloudmiddle  = int(round((cloudbase+cloudtop)/2))

    print('Compute final :')
    print('cloudbase : ij='+str(cloudbase)+' for z='+str(z[cloudbase]))
    print('cloudtop  : ij='+str(cloudtop)+' for z='+str(z[cloudtop]))

    return cloudbase,cloudmiddle,cloudtop,cloudbase0,cloudtop0


def findlcl(Q0,TA0,pres,z):
    # Find Lifting Condensation Level
    gradsec = -0.0098 # K/m z in meter
    ss      = pres.shape
    zz      = np.repeat(np.repeat(z[ :,np.newaxis, np.newaxis],ss[1],axis=1),ss[2],axis=2)
    q0      = np.repeat(Q0[np.newaxis,:,:] ,ss[0],axis=0)
    t0      = np.repeat(TA0[np.newaxis,:,:],ss[0],axis=0)
    t1      = t0+gradsec*zz #3D
    dqsat   = (qs(t1,pres)[0] < q0)*zz
    dqsat[dqsat==0]=np.nan

    del q0,t0,t1,zz
    #iw      = 45
    #print t1[:,iw,iw],z[:],qs(t1,pres)[0][:,iw,iw],dqsat[:,iw,iw],q0[0,iw,iw],np.isnan(dqsat).all(),~np.isnan(dqsat).all()
    lcl,idxlcl = [[np.nan] for ij in range(2)]
    if ~np.isnan(dqsat).all():
      lcl      = np.nanmin(dqsat,axis=0)
      idxlcl   = _nanargmin(dqsat,axis=0)
      #print lcl,idxlcl
    del dqsat
    return idxlcl,lcl

def gradhum(RVT,T):
    # 2D
    gradhum = -CC.RG*(1.0+(CC.RLVTT*RVT)/(CC.RD*T)) \
             / (CC.RCP+(pow(CC.RLVTT,2.0)*RVT)/(CC.RV*pow(T,2.0)))
    return gradhum

   
def findlfc(idxlcl,RVT,T,z):
    # Find Level of Free Convection
    gradsec = -0.0098 # K/m z in meter
    ss  = RVT.shape
    zz  = np.repeat(np.repeat(z[ :,np.newaxis, np.newaxis],ss[1],axis=1),ss[2],axis=2)
    dt  = [np.nan]

    t1 = np.zeros(ss)
    if ~np.isnan(idxlcl).all():
      m,n = idxlcl.shape
      I,J = np.ogrid[:m,:n]
      t1  = T[0,:,:]+gradsec*zz[idxlcl, I, J] #3D
      t1  = np.repeat(t1[np.newaxis,:,:],ss[0],axis=0)
      #iw  = 45
      for ij,field in enumerate(z):
        #print iw,ij,t1[ij-1,iw,iw],gradhum(RVT[ij-1,iw,iw],t1[ij-1,iw,iw]),(idxlcl[iw,iw]<ij),(z[ij]-z[ij-1])
        t1[ij,:,:]  = t1[ij-1,:,:] + gradhum(RVT[ij-1,:,:],t1[ij-1,:,:])*(idxlcl<ij)*(z[ij]-z[ij-1])
        #print iw,ij,t1[ij,iw,iw],T[ij,iw,iw]
        #print  gradhum(RVT[ij-1,:,:],t1[ij-1,:,:])*(idxlcl<ij)*(z[ij]-z[ij-1])
      dt         = (t1 > T)*zz #(qs(t1,pres)[0] < q0)*zz
      dt[dt==0.] = np.nan
      
    del RVT, T, t1, zz
    gc.collect()

    lfc,idxlfc,lnb,idxlnb = [np.nan for ij in range(4)]
    #print dt.any()
    #print ~np.isnan(dt).any()
    if ~np.isnan(dt).all():
      #print iw,t1[:,iw,iw],T[:,iw,iw],dt[:,iw,iw]
      lfc        = np.nanmin(dt,axis=0)
      lnb        = np.nanmax(dt,axis=0)

      #dt         = np.array(dt,dtype=object)
      dt1        = deepcopy(dt)
      dt1[np.isnan(dt1)] = np.nanmax(dt)
      idxlfc     = np.nanargmin(dt1,axis=0).astype(float)

      dt2        = deepcopy(dt)
      dt2[np.isnan(dt2)] = np.nanmin(dt)
      idxlnb     = np.nanargmax(dt2,axis=0).astype(float)

      idxlfc[idxlfc<idxlcl]  = np.nan
      idxlnb[((idxlnb<idxlfc) | (idxlnb<idxlcl))] = np.nan

      del dt1,dt2
      
    del dt
    return idxlfc,lfc,idxlnb,lnb

def findinv(DATA,inv,var1D,namez,offset=0.25):
    # Routine to find inversion
    # Goal : add multiple definition
    # namez = vertical_levels
    #DATA    = nc.Dataset(file,'r')
    if inv=='THV' or inv=='THLM' or inv=='THT':
        z           = DATA[namez][:]  # m
        if inv in DATA.variables.keys(): 
            data = DATA[inv] 
        else:
            data = createnew(inv,DATA,var1D)
        ss          = data.shape
        ZZ          = np.repeat(np.repeat(z[ :,np.newaxis, np.newaxis],ss[1],axis=1),ss[2],axis=2)
        idxzi       = findTHV3D(ZZ,data,offset)
        del data,ZZ
    return idxzi

def Interp(lat,lon,data,latI,lonI):
    ip = interpolate.interp2d(lon[:], lat[:], data)
    zi = ip(lonI[:], latI[:])#.conj().transpose()
    zi = zi[::-1,:]
    return zi

def do_delete2(objects,mask,nbmin,rename=True,clouds=None):
    nbmax   = np.max(objects)
    print(nbmax,nbmin)
    time1 = time.time()
    objects = delete_smaller_than(mask,objects,nbmin)
    time2 = time.time()
    print('%s function took %0.3f ms' % ("delete smaller", (time2-time1)*1000.0))
    if clouds is not None:
        print('filter clouds not None')
        objects = delete_clouds(objects,clouds)

    #print np.max(objects),len(np.unique(objects))
    if rename :
        labs = np.unique(objects)
        objects = np.searchsorted(labs, objects)
    nbr = len(np.unique(objects))-1 # except 0
    print('\t', nbmax - nbr, 'objects were too small')
    return objects,nbr

def delete_smaller_than_old(mask,obj,minval):
  #print np.max(mask)
  #print np.max(obj)
  sizes = ndimage.sum(mask,obj,np.unique(obj[obj!=0]))
  del_sizes = sizes < minval
  print(sizes,del_sizes)
  del_cells = np.unique(obj[obj!=0])[del_sizes]
  print(del_cells)
  for cell in del_cells :
    obj[obj==cell] = 0
  return obj

def delete_smaller_than(mask,obj,minval):
  sizes = ndimage.sum(mask,obj,np.unique(obj[obj!=0]))
  del_sizes = sizes < minval
  del_cells = np.unique(obj[obj!=0])[del_sizes]
  # new version
  ss        = obj.shape
  objf      = obj.flatten()
  ind       = np.in1d(objf,del_cells)
  objf[ind] = 0
  obj       = objf.reshape(ss)
  return obj

def delete_clouds(obj,cld,min=1):
    #mask       = do_unique(deepcopy(obj))
    maskclouds = do_unique(deepcopy(cld))
    maskclouds *= obj
    del_cells  = np.unique(maskclouds[maskclouds!=0])
    print('del cells : ', del_cells)
    print('unique : ',np.unique(obj[obj!=0]))
    
    # Remove all object that have clouds
    ss        = obj.shape
    objf      = obj.flatten()
    ind       = np.in1d(objf,del_cells)
    objf[ind] = 0
    obj       = objf.reshape(ss)
    
    return obj

