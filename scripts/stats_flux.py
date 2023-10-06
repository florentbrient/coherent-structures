# -*- coding: utf-8 -*-
"""
Created on Aug 20
"""
import sys
#sys.path.append('/home/brientf/MESONH/scripts/')
from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
import tools as tl
#from itertools import chain
from copy import deepcopy
#import scipy as sp
#from scipy import ndimage as ndi
import time
from collections import OrderedDict
import Constants as CC


def format(value):
    return "%13.3f" % value
def format2(value):
    return "%13s" % value

def plottophat(alt,tot,Fi,tophat,covar,subcomp,var,filesave='output',title='No',legend=1): #tot,Fi,tophat,covar,subcomp
  fig = plt.figure()
  size= [6.0,10.0]; lw=2; fts=15
  ax  = fig.add_subplot(111)
  tl.adjust_spines(ax, ['left', 'bottom'])
  plt.plot(tot,alt,'k',label='Tot',lw=lw)
  plt.plot(Fi,alt,'r',label='Fi',lw=lw)
  plt.plot(tophat,alt,'r',linestyle='--',label='Tophat',lw=lw)
  plt.plot(covar,alt,'r',linestyle=':',label='Covar',lw=lw)
  plt.plot(subcomp,alt,'b',linestyle='-.',label='SubComp',lw=lw)
  plt.plot(tophat+subcomp,alt,'k',linestyle='-.',label='Top-hat + SubComp',lw=lw)
  #ax.set_xlim([0.,100.]) 
  #ax.set_ylim([0.,40.])
  xlabel='Fluxes'
  filesave=filesave.replace('XX',var)
  ax.set_xlabel(xlabel,fontsize=fts)
  ax.set_ylabel('Altitude (km)',fontsize=fts)
  ax.set_title(title)
  plt.legend(loc=legend, prop={'size': fts-5})
  plt.xticks(size=fts)
  plt.yticks(size=fts)
  fig.savefig(filesave + '.pdf')
  plt.close()

# call  : python stats_objects.py downdraft SVT006 L25.6 003
typ      = sys.argv[1].split(',') # downdraft
selSVT   = sys.argv[2].split(',') # SVT006
simu     = sys.argv[3] #'L25.6'
hour     = sys.argv[4] #'003'
cas      = sys.argv[5] #'FIRE'
name     = sys.argv[6] #'V0301'
subcloud = sys.argv[7] #subcloud
nocloud  = sys.argv[8] #cloud
try:
  cond = sys.argv[9]
except:
  cond = ''
  
# Plot tophat
makefigures=True

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
  
print('subcloud ',subcloud)
print('subprefix ',subprefix)
print('nocloud ',nocloud)
print('nocldprefix ',nocldprefix)

# open file
file_dir       = "../data/"+cas+"/"+simu+"/"
filename_input = "sel_"+simu+".1."+name+"."+hour+".nc4"

# out files
dir_out        = "../data/d_stats/"+cas+"/";tl.mkdir(dir_out)
dir_out       += simu+"/";tl.mkdir(dir_out)
dir_decomp     = dir_out+'d_decomp/';tl.mkdir(dir_decomp)

liste_fichiers = file_dir+filename_input
print(liste_fichiers)

# Open data
DATA = Dataset(liste_fichiers)

#################
# Dimensions (similar to plot_LES_3D.py)
var1D  = ['vertical_levels','W_E_direction','S_N_direction'] #Z,Y,X
data1D,nzyx,sizezyx = [OrderedDict() for ij in range(3)]
for ij in var1D:
  data1D[ij]  = DATA[ij][:]/1000. #km #tl.removebounds(DATA[ij][:])
  nzyx[ij]    = data1D[ij][1]-data1D[ij][0]
  sizezyx[ij] = len(data1D[ij])

nxny   = nzyx[var1D[1]]*nzyx[var1D[2]] #km^2

X    = DATA.variables[var1D[1]][:]
Y    = DATA.variables[var1D[2]][:]
ALT    = data1D[var1D[0]]

dz     = [0.5*(ALT[ij+1]-ALT[ij-1]) for ij in range(1,len(ALT)-1)]
dz.insert(0,ALT[1]-ALT[0])
dz.insert(-1,ALT[-1]-ALT[-2])
nxnynz = np.array([nxny*ij for ij in dz]) # volume of each level
nxnynz = np.repeat(np.repeat(nxnynz[:, np.newaxis, np.newaxis]
                             ,len(Y), axis=1), len(X), axis=2)

nx,ny,nz = (X[1]-X[0])/1e3,(Y[1]-Y[0])/1e3,(ALT[1]-ALT[0])/1e3
#################

# Variable to be analyzed/saved at the same time
namevar  = ['THLM','RNPM']
Const    = [CC.RCP,CC.RLVTT] # Constants used for quantified fluxes in W/m2
nbvar    = len(namevar) # number of variables

WT   = DATA.variables['WT'][:]
THLM = DATA.variables['THLM'][:]
RNPM = DATA.variables['RNPM'][:]
P    = DATA.variables['PABST'][:]
THT  = DATA.variables['THT'][:]
#print(WT.shape)
#create rho
TT   = tl.tht2temp(THT,P)
RHO  = tl.createrho(TT,P)

RCTflag=False
try:
  RCT  = DATA.variables['RCT'][:]
  if RCT.max()!=0:
    RCTflag=True
except:
  pass

try:
  CLDFR  = DATA.variables['CLDFR'][:]
  print('CLDFR : ',CLDFR.max())
except:
  CLDFR  = np.nan

# Find maximum with tracer concentration
nameSVT=['SVT006','SVT003']
if subcloud and RCTflag: #RCT.max()!=0:
  nameSVT=['SVT005','SVT002']
try:
  SVTtop = DATA.variables[nameSVT[0]][:]
except:
  SVTtop = DATA.variables[nameSVT[1]][:]
SVTmean  = SVTtop.mean(axis=(1,2))
#print(SVTmean)
idxsvt = SVTmean.argmax()
print('idx, alt(max SVT), max SVT: ',nameSVT,idxsvt,ALT[idxsvt],SVTmean[idxsvt])

# Very important
# Importance of zmax : Maximul altitude for averaging !!

# zmax = Altitude maximum to average fluxes !!
#if cas == 'FIRE':
#  altmax = 0.62 # in km
#elif cas == 'BOMEX' or  cas == 'IHOP':
#  altmax = 2. # in km
#else:
#  altmax = None
#zmax   = tl.near(ALT,altmax)
#print('ALT ',altmax,zmax,ALT)

# Impose idx from highest SVT (cloud top or inversion height)
#zmax  = idx 
# Old version

# Second try: Consistent with plot3D_multicase
# Find PBL top if no clouds
if not RCTflag:
   zmax = idxsvt
# Find PBL top based on theta_l gradient if sub-cloud only
elif subcloud:
   offset = 0.25
   idxzi,toppbl,grad = tl.findpbltop('THLM',DATA,offset=offset)
   zmax   = idxzi
   print("########")
   print('Top PBL : ',idxzi,idxsvt)
   print('height PBL (m) : ',ALT[idxzi],ALT[idxsvt])
   print("########")
# Find PBL top by cloud top
else:
   # Define as the first layer where RCT (ql) > epsilon (1e-6 by default)
   epsilon = 1e-6
   cloudbase,cloudmiddle,cloudtop,zb,zi = tl.cloudinfo(ALT,RCT,epsilon)
   zmax    = cloudtop

# Altitude constrained with zmax
alts = ALT[:zmax+1]
print('*******')
print('Altitude for compuation idx,zmax ',zmax,alts[-1])


# Resize
WT       = WT[:zmax+1,:,:]
THLM     = THLM[:zmax+1,:,:]
RNPM     = RNPM[:zmax+1,:,:] #*1000.
dz       = dz[:zmax+1]
nxnynz   = nxnynz[:zmax+1,:,:]
if nocloud and CLDFR.max()!=0:
   CLDFR = CLDFR[:zmax+1,:,:]

# Calculate anomalies
ATHLM    = tl.anomcalc(THLM)
ARNPM    = tl.anomcalc(RNPM)
# Domain-mean fluxes
WTTHLM   = WT*ATHLM;print(np.mean(WTTHLM))
WTRNPM   = WT*ARNPM;print(np.mean(WTRNPM))
# Domain-mean average
WTM      = np.nanmean(WT,axis=(1,2))
THM      = np.nanmean(THLM,axis=(1,2))
RTM      = np.nanmean(RNPM,axis=(1,2))
# Total volume
vtot     = np.sum(nxnynz)

##### Object selection
# m value (strength of the conditional sampling)
listthrs = [0,1,2,3,4] #[0]
N1       = len(listthrs)
# Vmin (can be defined in units or volume in km^3)
minchar = 'volume' #'volume' #unit
if minchar == 'volume':
    nbmins   = [0.005,0.01,0.02,0.05,0.10,0.25,0.50]
else:
    nbmins   = [0,10,50,100,500,1000,5000,10000,20000]
N2 = len(nbmins) #,50,100,500,1000,5000,10000,20000];
#####

nbr,surf,volume,alpha,altmin,altmax  = [np.zeros((N1,N2)) for ij in range(6)]
rfluxVAR,rfremVAR,fluxVAR2,fremVAR2  = [np.zeros((nbvar,N1,N2)) for ij in range(4)]
fluxVARpos,fluxVARneg                = [np.zeros((nbvar,N1,N2)) for ij in range(2)]
Massflux = np.zeros((N1,N2))

# Test with one value
maketest = True
testch   = ''
if maketest:
  listthrs = [1] #[2]
  #nbmins   = [nbmins[4]] #[0,10,50,100,500,1000,5000,10000,20000] #[10000]
  testch   = '.test'

### Name the output file
# if condition AND/OR
jointyp  = cond
# First letter
typ_first = [x[0] for x in typ]
# info tracer
sel_first = [x[-4:] for x in selSVT]
# Important info from filename
print(filename_input.split('.'))
file_short='_'.join(np.array(filename_input.split('.'))[[0,-2]])
file_short=file_short.replace('sel_', '') 
# Name of the output file
fileout0  = dir_out
#fileout0 += 'stats.'+jointyp.join(typ)+'.'+''.join(selSVT)+'.LL.'+filename_input+'.3'+testch+'.txt'
fileout0 += 'stats.'+''.join(typ_first)+'.'+''.join(sel_first)+\
        '.LL.'+file_short+testch+subprefix+nocldprefix+\
        '.txt'

# Run along the number of m values (conditional sampling)
#id1=-1
for id1,ll in enumerate(listthrs):
  OBJ      = []
  nameobjs = []
  #id1     +=1
  for ik in range(len(typ)):
    nameobjs.append(typ[ik][0:4]+'_'+selSVT[ik]+'_'+str(ll).zfill(2))
    print('nameobjs ',ll,nameobjs[ik])
    OBJ.append(DATA.variables[nameobjs[ik]][:][:zmax+1,:,:])
  sz        = OBJ[0].shape

  fileout = fileout0.replace('LL',str(ll)); 
  f = open(fileout, 'w')
  f.write("File name : "+liste_fichiers+"\n")
  f.write("Object types : "+''.join(selSVT)+"\n")
  f.write(jointyp.join(typ)+" threshold : "+str(ll*0.5)+"\n")
#  names  = ['nbmin',' nb',' surf',' vol',' rvol',' rTHLMflux',' rRNPMflux',' altmin',' altmax\n']
#  unity  = ['-','-','km2','km3','%','%','%','km','km\n']  
  names  = ['nbmin ','nb ','surf ','vol ','rvol ','rTHLMflux ','rRNPMflux','altmin ','altmax', \
            'pTHLMfl','nTHLMfl','pRNPMfl','nRNPMfl','MassFlux','WTmax''\n']
  unity  = ['-','-','km2','km3','%','%','%','km','km','W/m2','W/m2','W/m2','W/m2','kg/m2/s','m/s''\n']
  
  for v in names:
    f.write(format2(v))
  for v in unity:
    f.write(format2(v))
  for id2,bin in enumerate(nbmins):
    print('ib, bin : ',id2,bin)
    mask   = tl.do_unique(deepcopy(OBJ[0]))
    if minchar == 'volume':
        mask = mask*nxnynz # Not good... 2D not 3D
    time1 = time.time()
    OBJ[0],nbr[id1,id2]  = tl.do_delete2(OBJ[0],mask,bin,rename=True)
    NT    = int(np.nanmax(nbr))
    print(nbr[id1,id2],NT)
    time2 = time.time()
    print('%s function 1 took %0.3f s for Vmin=%7.3f' % ("delete", (time2-time1)*1.0, bin))
    # Check volume or number
    #labs = np.unique(OBJ[0])
    #print('labs ',labs)
    #for lab in labs:
    #  print(lab,np.sum(mask[OBJ[0]==lab]))
        
    if nbr[id1,id2] >-1:
     condition         = (OBJ[0]!=0)
     for ik in range(len(typ)-1):
       ikk              = ik+1
       mask             = tl.do_unique(deepcopy(OBJ[ikk]))
       if minchar == 'volume':
         mask = mask*nxnynz # Not good... 2D not 3D
       OBJ[ikk],nbrtmp  = tl.do_delete2(OBJ[ikk],mask,bin,rename=True)
       nbr[id1,id2]    += nbrtmp
       if jointyp  == 'and':
         condition      = (condition & (OBJ[ikk]!=0))
       elif jointyp  == 'or':
         condition      = (condition | (OBJ[ikk]!=0))
         
     if nocloud and CLDFR.max()!=0:
       for inn in range(1,NT+1):
         idx = np.sum(CLDFR[OBJ[0]==inn])
         idxlc = 27 # 27 pixels for a cloud
         print(inn,idx)
         if idx>=idxlc:
           print(CLDFR[OBJ[0]==inn])
           print(inn, np.sum(CLDFR[OBJ[0]==inn]))
           condition[OBJ[0]==inn]=0.
           nbr[id1,id2]=nbr[id1,id2]-1
       print('NT new: ',NT,int(np.max(nbr)))

     surf[id1,id2]     = np.sum(condition)*nxny #km2
     volume[id1,id2]   = np.sum(nxnynz*condition)#*nxnynz #km3
     alpha[id1,id2]    = volume[id1,id2]/vtot
     
     # Initialisation
     #fTHLM,fTHLMfr,fTHLMrem,fTHLMtot,
     tmpflux,tmpvar = [np.zeros((nbvar,sz[1],sz[2])) for ij in range(2)]
     fVAR,fVARfr,fVARrem,fVARtot,VARM= [np.zeros((nbvar,sz[0])) for ij in range(5)]
     
     frac, RHOCD, WTI, THI, RTI = [np.zeros(sz[0]) for ij in range(5)]
     tophatVAR, intraVAR = [np.zeros((nbvar,sz[0])) for ij in range(2)]
     MM = np.zeros(sz[0])
     
     # ratio of fluxes (second try -> good)
     for ik in range(sz[0]):
       cd       = condition[ik,:,:]
       frac[ik] = np.sum(cd)/float(sz[1]*sz[2])
       WT0      = WT[ik,:,:]
       RHO0     = np.copy(RHO[ik,:,:])
       RHOCD[ik]= np.nanmean(RHO0[cd])       
       # Mass flux
       MM[ik]  = frac[ik]*RHOCD[ik]*np.nanmean(WT0[cd])
       # top-hat approx
       WTI[ik]        = np.nanmean(WT0[cd])
       
       # tmpvar and tmpflux
       # Var 0 = THLM
       tmpflux[0,:,:] = WTTHLM[ik,:,:]
       tmpvar[0,:,:]  = THLM[ik,:,:]
       VARM[0,:]      = THM
       # Var 1 = RNPM
       tmpflux[1,:,:] = WTRNPM[ik,:,:]
       tmpvar[1,:,:]  = RNPM[ik,:,:]
       VARM[1,:]      = RTM
       
       #tmp0     = ; tmp1 = WTRNPM[ik,:,:]
       #TH0      = THLM[ik,:,:];    RT0 = RNPM[ik,:,:]
       
       for iv in range(nbvar):
         tmpflux0 = tmpflux[iv,:,:]
         tmpvar0  = tmpvar[iv,:,:]
         
         # CS flux contribution
         fVARfr[iv,ik]  = frac[ik]*np.nanmean(tmpflux0[cd])
         fVAR[iv,ik]    = np.nanmean(tmpflux0[cd])
         # remaining part
         fVARrem[iv,ik]   = (1.0-frac[ik])*np.nanmean(tmpflux0[~cd])
         fVARtot[iv,ik]   = np.mean(tmpflux0)
         plt.plot()
       
         # internal variability
         VARi             = np.nanmean(tmpvar0[cd])
         tophatVAR[iv,ik] = (WTI[ik]-WTM[ik])*(VARi-VARM[iv,ik])
         intraVAR[iv,ik]  = np.nanmean((tmpvar[iv,cd]-VARi)*(WT0[cd]-WTI[ik]))
         
       #intraTHLM[ik]  = np.nanmean((WT0[cd]-WTI[ik])*(TH0[cd]-THI[ik]))
       # CS flux contribution
       #fRNPMfr[ik]  = frac[ik]*np.nanmean(tmp1[cd])
       #fRNPM[ik]    = np.nanmean(tmp1[cd])
       # remaining part
       #fRNPMrem[ik] = (1.0-frac[ik])*np.nanmean(tmp1[~cd])
       #fRNPMtot[ik] = np.mean(tmp1)
       # top-hat approx
       #RTI[ik]        = np.nanmean(RT0[cd])
       #tophatRNPM[ik] = (WTI[ik]-WTM[ik])*(RTI[ik]-RTM[ik])
       # internal variability
       #intraRNPM[ik]  = np.nanmean((RT0[cd]-RTI[ik])*(WT0[cd]-WTI[ik]))

        # Flux totaux, Fi, tophat, intra-inter, sub compensatoire
     #filesave='decomposition_flux_'+var+'_'+typ
     filesave  = dir_decomp+'decomposition_flux_'+'XX'+'_'+jointyp.join(typ)+'_'+''.join(selSVT)
     filesave += '_'+simu+'_'+hour+'_'+str(ll)+'_'+str(bin)
     if name != 'V0301':
        filesave += '_'+name
     title     = 'Flux '+cas+' '+simu+' '+hour \
               +' '+jointyp.join(typ)+' (m='+str(float(ll)/2.)+',Vmin='+str(bin)+')'

     # Areal fraction of the object (ij)
     alphai = frac
     # maximum vertical velocity
     WTmax     = np.nan
     if ~np.isnan(WTI).all():
       idx       = np.nanargmax(abs(WTI))
       WTmax     = WTI[idx]
     
     # Mass flux in kg/m2/s
     Massflux[id1,id2]  = np.nansum(MM*dz)/np.nansum(dz)

     
     for iv in range(nbvar):
       tot=fVARtot[iv,:];Fi=alphai*fVAR[iv,:]; tophat=alphai*tophatVAR[iv,:]
       covar=alphai*intraVAR[iv,:]; subcomp=(alphai*alphai/(1.0-alphai))*tophatVAR[iv,:]
       if makefigures:
         plottophat(alts,tot,Fi,tophat,covar,subcomp,namevar[iv],filesave,title,legend=3)

#     tot=fRNPMtot; Fi=alphai*fRNPM; tophat=alphai*tophatRNPM; 
#     covar=alphai*intraRNPM; subcomp=(alphai*alphai/(1.0-alphai))*tophatRNPM
#     plottophat(alts,tot,Fi,tophat,covar,subcomp,'RNPM',filesave,title,legend=1) #4

       # Vertical averaging
       # flux in W/m2
       fluxVARpos[iv,id1,id2]  = Const[iv]*np.nansum(fVARfr[iv,:]*RHOCD*dz*(fVARfr[iv,:]>0))/np.nansum(dz)
       fluxVARneg[iv,id1,id2]  = Const[iv]*np.nansum(fVARfr[iv,:]*RHOCD*dz*(fVARfr[iv,:]<=0))/np.nansum(dz)

       rfluxVAR[iv,id1,id2]    = 100.*np.nansum(abs(fVARfr[iv,:]*dz))/np.nansum(abs(fVARtot[iv,:]*dz))
       #fluxRNPM[id1,id2] = 100.*np.nansum(abs(fRNPMfr*dz))/np.nansum(abs(fRNPMtot*dz))
       print('*****')
       print(iv,id1,id2,np.nansum(abs(fVARfr[iv,:])),np.nansum(abs(fVARtot[iv,:])),\
             np.nansum(abs(fVARfr[iv,:]*dz))/np.nansum(abs(fVARtot[iv,:]*dz)))
       print(iv,id1,id2,np.nansum((fVARfr[iv,:])),np.nansum((fVARtot[iv,:])),\
             np.nansum((fVARfr[iv,:]*dz))/np.nansum((fVARtot[iv,:]*dz)))
       #print(id1,id2,np.nansum(abs(fRNPMfr)),np.nansum(abs(fRNPMtot)),np.nansum(abs(fRNPMfr*dz))/np.nansum(abs(fRNPMtot*dz)))

       rfremVAR[iv,id1,id2]    = 100.*np.nansum(abs(fVARrem[iv,:]*dz))/np.nansum(abs(fVARtot[iv,:]*dz))
       #fremRNPM[id1,id2] = 100.*np.nansum(abs(fRNPMrem*dz))/np.nansum(abs(fRNPMtot*dz))

     # size of the objects
     if len(typ)==1:
       listObjs          = np.unique(OBJ[0][OBJ[0]!=0]);#print listObjss
       min0,max0         = [],[]
       for il in listObjs: #range(max(listObjs)):
          cd = np.where(OBJ[0].reshape((sz[0],sz[1]*sz[2]))==il)
          # Error? ? si nz change with height???
          #min0tmp = np.min(cd[0])*nz*1e3
          # Fixed?
          min0tmp = ALT[round(np.min(cd[0]))]
          min0.append(min0tmp);
          max0tmp = ALT[round(np.max(cd[0]))] 
          max0.append(max0tmp)
          #print il,cd,altmin,altmax
       altmin[id1,id2]   = np.mean(min0)
       altmax[id1,id2]   = np.mean(max0)
    else: # all objects
       altmin[id1,id2]   = np.nan
       altmax[id1,id2]   = np.nan


#    tab               = [bin,nbr[id1,id2],surf[id1,id2],volume[id1,id2],100.*alpha[id1,id2],fluxTHLM[id1,id2],fluxRNPM[id1,id2],altmin[id1,id2],altmax[id1,id2]]
    tab               = [bin,nbr[id1,id2],surf[id1,id2],volume[id1,id2],100.*alpha[id1,id2] \
                        ,rfluxVAR[0,id1,id2],rfluxVAR[1,id1,id2],altmin[id1,id2],altmax[id1,id2] \
                        ,fluxVARpos[0,id1,id2],fluxVARneg[0,id1,id2] \
                        ,fluxVARpos[1,id1,id2],fluxVARneg[1,id1,id2] \
                        ,Massflux[id1,id2],WTmax
                        ]
    print(tab)
    for v in tab:
      f.write(str(format(v)))
    f.write("\n")
  f.close()




