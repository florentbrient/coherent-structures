#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 14:19:50 2022

@author: fbrient

# Goal:
# Plot temporal evolution of fluxes

"""
import sys
import numpy as np
from matplotlib import pyplot as plt
import tools as tl
import makefigures as mf

# Read txt files
def openfilestat(file):
  ijmin   = 5
  f       = open(file, 'r')
  tab     = [line.rstrip('\n').split() for line in f]
  nametab = tab[3]
  Nmin    = len(tab[ijmin:])
  Nval    = len(tab[ijmin])  
  print(Nmin,Nval)
  f.close()
  return tab,nametab,Nmin,Nval,ijmin


# Read txt files
def openfilemean(file):
  ijmin   = 2 # line where values are
  f       = open(file, 'r')
  tab     = [line.rstrip('\n').split() for line in f]
  print(tab)
  value   = tab[ijmin];print(len(value))
  value   = float(value[0])
  f.close()
  return value

# 
def _plotscatter(ylab,extray,nameobjs,hours,filesave,\
                 tabij=[0,3],fts=15,size=[8.0,8.0],\
                 xlab='',regress=True):
    #tabij= [0,3] # updraft and downdraft
        
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    for ij_tab in tabij:
        tmp = abs(extray[ij_tab].squeeze())
        color = tl.findhatch(nameobjs[ij_tab],typ='color')
        
        ax.plot(value_mean,tmp,'.',color=color)
        for ih,hh in enumerate(hours):
            ihour = str(int(hh))
            ax.text(value_mean[ih],tmp[ih],ihour,color=color,fontsize=12)
        if regress:
            coef = np.polyfit(value_mean,tmp,1)
            poly1d_fn = np.poly1d(coef)
            corr = np.corrcoef(value_mean,tmp)[0,1]
            xx = np.linspace(value_mean.min(),value_mean.max())
            plt.plot(xx, poly1d_fn(xx),'--',\
                     label=f'$r = {corr:.2f}$',\
                     color=color)
            
    plt.legend(fontsize="15")
    ax.locator_params(axis='y', nbins=6)
    ax.axvline(x=0, color='k', linewidth=1)
    tl.adjust_spines(ax, ['left', 'bottom'])
    ax.get_yaxis().set_tick_params(direction='out')
    ax.get_xaxis().set_tick_params(direction='out')
    ax.set_xlabel(xlab,fontsize=fts)
    ax.set_ylabel(ylab,fontsize=fts)
    plt.xticks(size=fts)
    plt.yticks(size=fts)
    fig.set_size_inches(size[0], size[1])
    fig.savefig(filesave + '.png')
    fig.savefig(filesave + '.pdf')   
    plt.close()
    return None


# Same opening as stats_flux.py
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
  
#if cas=='IHOP':
#   typ     = [typ[i] for i, e in enumerate(selSVT) if '002' not in e]
#   selSVT  = [ij for ij in selSVT if '002' not in ij]

# Plot scatter
plotscat_cond  = True
varindex       = 'D2t'
scat_cond_Nmin = [0.02]

print('iCI: ',selSVT,typ)
subcloud,subprefix,nocloud,nocldprefix = tl.suffix_name(subcloud,nocloud)

#hours     = ['002','004','006','008']
hours     = [ij for ij in hour.split(',')]
print('hours ',hours)
typs      = ['updraft','downdraft','downdraft','downdraft']
typs.append('or'.join(typs))

listthrs  = [0,1,2,3,4] ; 

NH        = len(hours)
NL        = len(listthrs)

# input directory
dir_in    = "../data/d_stats/"+cas+"/"+simu+"/"

### Name the input file
# if condition AND/OR
jointyp  = cond
# First letter
typ_first = [x[0] for x in typ] # updraft etc...
# info tracer
sel_first = [x[-4:] for x in selSVT] # SVT001 etc...
print(typ_first,selSVT)
if len(typ_first)>2:
  typ_first.append(''.join(typ_first))
  sel_first.append(''.join(sel_first))
  print(selSVT,sel_first)
  
NT        = len(typ_first)

# Important info from filename

filein0   = dir_in
filein0  += 'stats.'+'TT'+'.'+'SS'+\
        '.LL.'+simu+'_HH'+subprefix+nocldprefix+\
        '.txt'

# out files
dir_out     = "../data/figures/"+cas+"/";tl.mkdir(dir_out)
dir_out    += simu+"/";tl.mkdir(dir_out)
dir_out    += 'd_time/';tl.mkdir(dir_out)


# Output filename
namefile0=dir_out+'SSSS_VVVV_MM_NNNN'+subprefix+nocldprefix 
# time_Var_Nmin_subcloud_nocloud

#####################
# Goal : reprodude nameobjs 
# Is it relevant to do it here?
# Need to be improved
# Only for used for having typs2
#
# Objects characteristics (From plot_LES_3D.py)
#thrs   = 2
#thch   = str(thrs).zfill(2)
# Minimal volume for object detection vmin
# Brient et. al 19 (GRL): 0.25 km^3
# By default : 0.02 km^3
  
# Select by? 
#minchar = 'volume' #unit
#if minchar == 'volume':
#    vmin   = 0.02
#    suffixmin = '_vol'+str(vmin) 
#elif minchar == 'unit':
#    nbmin  = 100 #100 #1000
#    suffixmin = '_nb'+str(nbmin)
nbplus = tl.svttyp(cas,simu) #1
# Object based on tracer concentration AND vertical velocity?
AddWT  = 1
typs2, objtyp = tl.def_object(nbplus=nbplus,AddWT=AddWT)
#typs2 = [ij+'_'+thch  for ij in typs2] 
#print(typs2,objtyp)
typs2.append('All')
print('TT ',typs2, typs)

############################

#print(typs,typs2)updr_SVT004_WT
data     = {}
nameobjs = []
nameidx  = []
if cas=='IHOP':
  nameidx = [ij for ij in typs2 if '002' in ij]

for it,typ in enumerate(typs2):
  print('color : ',typ)
  #print(tl.findhatch(typs2[it],typ='color'))
  
  filein1a = filein0.replace('TT',typ_first[it])
  filein1b = filein1a.replace('SS',sel_first[it])   
  for ih,hh in enumerate(hours):
    filein1c = filein1b.replace('HH',hh); 
    for il,ll in enumerate(listthrs):
      filein = filein1c.replace('LL',str(ll));
      print(filein)
      tab,nametab,Nmin,Nval,ijmin  = openfilestat(filein)
            
      if il==0 and ih==0:
          data[typ] = np.zeros((NH,NL,Nmin,Nval))
      #if it==0 and il==0 and ih==0:
      #  data=np.zeros((NT,NH,NL,Nmin,Nval)) #[typs],hours,m,Vmin,Variable
        
      for ij in range(Nmin):
        for ik in range(Nval):
          if typ not in nameidx:
            data[typ][ih,il,ij,ik]=float(tab[ijmin+ij][ik])
      
  # replace typs2
  if typ not in nameidx:
    nameobjs  += [typ]
    
if cas=='IHOP':
    for ij in nameobjs:
       data[ij.replace('003','002')] = data.pop(ij)
    nameobjs = [ij.replace('003','002') for ij in nameobjs]

print('HERE : ',nameobjs,data.keys())

# Info figure
xxname=['time']
unitx =['hours']
ymin = dict(); ymax = dict()
ymin['rvol']=0.
#ymin['rTHLMflux']=0.;ymin['rRNPMflux']=0.
#ymax['rTHLMflux']=109.;ymax['rRNPMflux']=109.
levely=dict()
levely['rvol']=range(0,40,10)
levely['rTHLMflux']=range(0,120,10)
levely['rRNPMflux']=range(0,120,10)
levelx=dict()
xx   = 2
levelx['time']=range(0,int(hours[-1])+2*xx,xx)
#if case == 'IHOP':
#  levelx['time']=range(0,12+2,2)
if cas == 'BOMEX':
  levelx['time']=range(0,14+2,2)
if cas == 'ARMCU':
  levelx['time']=range(0,12+2,2)

DH    = float(hours[1])-float(hours[0]) # bof !
datax = np.arange(1,NH+1)*DH 
datax = np.repeat(datax,1,axis=0).reshape((1,NH))


#Nmax   = len(typs2)-1 #typs.index('downdraftorupdraft')
#idx     = Nmax
Nmax   = nameobjs[:-1]
idx    = nameobjs[-1]
print('Nmax ',Nmax)

# Which m you want to plot (m=1, index 2)
listplot = 2
if len(listthrs)==1:
  listplot=listthrs[0]
#idxlist  = np.argmin(listplot==listthrs); 
idxlist = listthrs.index(listplot)
print('idxlist : ',idxlist)

# Open domain mean characteristics of the simulation
dir_xx   = "../data/d_mean/"+cas+"/"
file_xx  = dir_xx+'mean_'+cas+'_'+simu+'_HH_DD.txt'
file_xx   = file_xx.replace('DD',varindex)
value_mean = np.zeros(len(hours))
for ih,hh in enumerate(hours):
    filein   = file_xx.replace('HH',hh)
    value_mean[ih]    = openfilemean(filein)
print('value_mean : ',value_mean)
    


namefile0 = namefile0.replace('MM',str(listplot))
for xx,xlab in enumerate(xxname):
  print('xlab',xx,xlab)
  namefile1 = namefile0.replace('SSSS',xlab)
  Nmintab   = data[idx][0,idxlist,:,0] # List of Vmin
  #print('N: ',Nmintab) 
  for ij,Nminx in enumerate(Nmintab):  # Vmin (0.005 for example)
     print('Nmin ',ij,Nminx)
     namefile2 = namefile1.replace('NNNN',str(Nminx))
     ysave1, ysave2 = None, None
     for ik,ylab in enumerate(nametab): #[1:]): #range(1,Nval):
        print('nametab ',ik,ylab) # ylab is the name of the value (surf, vol...)
        namefile3 = namefile2.replace('VVVV',ylab)
        
        datay = data[idx][:,idxlist,ij,ik]
        datay = datay.reshape((1,NH))
        #print('y ',datay)
        
        miny  = np.min(data[idx][:,idxlist,:,ik],axis=1) 
        maxy  = np.max(data[idx][:,idxlist,:,ik],axis=1)
        mins  = np.reshape((miny,maxy),(2,len(miny)))
        
        if len(typs)>1:
          extrax,extray = [],[] 
          #mins1         = []
          for it in Nmax: #range(Nmax):
            extrax.append(datax)
            tmp = data[it][:,idxlist,ij,ik]
            tmp = np.repeat(tmp,1,axis=0).reshape((1,NH))
            extray.append(tmp)
            #print('y ',it,extray)
            #print(typs[it],tmp)
            #miny1  = np.min(data[it,:,idxlist,:,ik],axis=1) 
            #maxy1  = np.max(data[it,:,idxlist,:,ik],axis=1)
            #mins1.append(np.reshape((miny1,maxy1),(2,len(miny1))))

        #print(datay)
        #plt.plot(datax,datay);plt.show()
        
        if plotscat_cond and (Nminx in scat_cond_Nmin):
            filesave = dir_xx+'scatter_'+varindex+'_'+ylab+'_'+str(listplot)+'_'+str(Nminx)
            print('filesave scatter: ',filesave)
            _plotscatter(ylab,extray,nameobjs,hours,filesave,\
                     xlab=varindex)
            if ylab == 'altmax':
                ysave1=extray[0] # updraft
            if ylab == 'altmin':
                ysave2=extray[3] #downdraft
            if ysave1 is not None and ysave2 is not None:
                filesave = dir_xx+'scatter_'+varindex+'_altmax_altmin_'+str(listplot)+'_'+str(Nminx)
                _plotscatter('',[ysave1,ysave2],
                     [nameobjs[0],nameobjs[3]],hours,filesave,\
                     xlab=varindex,tabij=[0,1])
                
        
        if ylab == 'altmax':
            datay = datay*np.nan         

        
        filesave = namefile3
        print('!! ',nameobjs)
        mf.plotstats(nameobjs,datax,datay,filesave=filesave,xlab=xlab,ylab=ylab,\
                     unitx=unitx[xx],xaxis=levelx,yaxis=levely,\
                     mins=mins,extrax=extrax,extray=extray)

        


