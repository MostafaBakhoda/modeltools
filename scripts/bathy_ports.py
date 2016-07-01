#!/usr/bin/env python
import modeltools.hycom
import argparse
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import modeltools.forcing.bathy
#import modeltools.hycom.io
import abfile
import modeltools.cice.io
import numpy
from mpl_toolkits.basemap import Basemap
import netCDF4
import logging
import scipy.ndimage.measurements

# Set up logger
_loglevel=logging.DEBUG
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate=False


def write_port_location(fid,kdport,ifport,ilport,jfport,jlport) :
   #c ---     'kdport' = port orientation (1=N, 2=S, 3=E, 4=W)
   #c ---     'ifport' = first i-index
   #c ---     'ilport' = last  i-index (=ifport for N or S orientation)
   #c ---     'jfport' = first j-index
   #c ---     'jlport' = last  j-index (=jfport for E or W orientation)
   #c ---     'lnport' = port length (calculated, not input)
   fid.write("%6d  'kdport' = port orientation (1=N, 2=S, 3=E, 4=W)\n"%kdport)
   fid.write("%6d  'ifport' = first i-index\n"%ifport)
   fid.write("%6d  'ilport' = last  i-index (=ifport for N or S orientation)\n"%ilport)
   fid.write("%6d  'jfport' = first j-index\n"%jfport)
   fid.write("%6d  'jlport' = last  j-index (=jfport for E or W orientation)\n"%jlport)
   return


def main(infile):

   bathy_threshold=0. # TODO

   # Read plon,plat
   gfile=abfile.ABFileGrid("regional.grid","r")
   plon=gfile.read_field("plon")
   plat=gfile.read_field("plat")
   gfile.close()

   # Read input bathymetri
   bfile=abfile.ABFileBathy(infile,"r",idm=gfile.idm,jdm=gfile.jdm,mask=True)
   in_depth_m=bfile.read_field("depth")
   bfile.close()
   in_depth=numpy.ma.filled(in_depth_m,bathy_threshold)

   L=numpy.zeros(in_depth.shape)

   ifports=[]
   ilports=[]
   jfports=[]
   jlports=[]
   kdports=[]

   process_south=True
   process_north=True
   process_west=True
   process_east=True


   # ports on  "southern edge"
   if process_south :
      jfport=1 # starts from 0
      jlport=1 # starts from 0
      S_label, S_nf = scipy.ndimage.measurements.label(~in_depth_m.mask[jfport,:])
      L[jfport,:] = S_label[:]
      kdport=2
      tmp=[str(elem) for elem in S_label]
      tmp="".join(tmp)
      logger.info("Southern edge segments:%s"%tmp)
      for i in range(S_nf ) :
         I,=numpy.where(S_label==i+1)
         ifport=I[0]
         ilport=I[-1]
         #print I
         #print i, ifport,ilport
         logger.info("Segment %5d on southern edge: First and last i (from 0): %5d %5d"%( i+1, ifport+1,ilport+1))
         ifports.append(ifport)
         ilports.append(ilport)
         jfports.append(jfport)
         jlports.append(jlport)
         kdports.append(kdport)

   # ports on  "northern edge"
   if process_north :
      jfport=in_depth.shape[0]-2 # starts from 0
      jlport =jfport             # starts from 0
      N_label, N_nf = scipy.ndimage.measurements.label(~in_depth_m.mask[jfport,:])
      L[jfport,N_label>0] = N_label[N_label>0] + L.max()
      kdport=1
      tmp=[str(elem) for elem in N_label]
      tmp="".join(tmp)
      logger.info("Northern edge segments:%s"%tmp)
      for i in range(N_nf ) :
         I,=numpy.where(N_label==i+1)
         ifport=I[0]
         ilport=I[-1]
         #print I
         #print i, ifport,ilport
         logger.info("Segment %5d on Northern edge: First and last i (from 0): %5d %5d"%( i+1, ifport+1,ilport+1))
         ifports.append(ifport)
         ilports.append(ilport)
         jfports.append(jfport)
         jlports.append(jlport)
         kdports.append(kdport)

   # ports on  "eastern edge"
   if process_east :
      ifport=in_depth_m.shape[1]-2 # starts from 0
      ilport=1                     # starts from 0
      E_label, E_nf = scipy.ndimage.measurements.label(~in_depth_m.mask[:,ifport])
      L[E_label>0,ifport] = E_label[E_label>0] + L.max()
      kdport=3


      tmp=[str(elem) for elem in E_label]
      tmp="".join(tmp)
      logger.info("Easterm  edge segments:%s"%tmp)
      for i in range(E_nf ) :
         I,=numpy.where(E_label==i+1)
         jfport=I[0]
         jlport=I[-1]
         #print I
         #print i, jfport,jlport
         logger.info("Segment %5d on Eastern  edge: First and last j (from 0): %5d %5d"%( i+1, jfport,jlport))
         ifports.append(ifport)
         ilports.append(ilport)
         jfports.append(jfport)
         jlports.append(jlport)
         kdports.append(kdport)

   # ports on  "western edge"
   if process_west : 
      ifport=1
      ilport=1
      W_label, W_nf = scipy.ndimage.measurements.label(~in_depth_m.mask[:,ifport])
      L[W_label>0,ifport] = W_label[W_label>0] + L.max()
      kdport=4
      tmp=[str(elem) for elem in W_label]
      tmp="".join(tmp)
      logger.info("Western  edge segments:%s"%tmp)
      for i in range(W_nf ) :
         I,=numpy.where(W_label==i+1)
         jfport=I[0]
         jlport=I[-1]
         #print I
         #print i, jfport,jlport
         logger.info("Segment %5d on Western  edge: First and last j (from 0): %5d %5d"%( i+1, jfport,jlport))
         ifports.append(ifport)
         ilports.append(ilport)
         jfports.append(jfport)
         jlports.append(jlport)
         kdports.append(kdport)


   # Open port output file
   logger.info("Writing to ports.input.tmp (NB: Transforming from C to Fortran indexing)")
   fid=open("ports.input.tmp","w")
   fid.write("%6d  'nports' = Number of ports \n"%len(kdports))
   for i in range(len(kdports)) :
      write_port_location(fid,kdports[i],ifports[i]+1,ilports[i]+1,jfports[i]+1,jlports[i]+1)
   fid.close()




   
   logger.info("Writing bathymetry plot to file ports_all.png")
   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   cmap=matplotlib.pyplot.get_cmap("Greys_r")
   ax.add_patch(matplotlib.patches.Rectangle((1,1),in_depth.shape[0],in_depth.shape[0],color=".5",alpha=.5))
   P=ax.pcolormesh(in_depth_m,cmap=cmap)
   #figure.colorbar(P,norm=matplotlib.colors.LogNorm(vmin=in_depth_m.min(), vmax=in_depth_m.max()),cmap=cmap)
   I,J=numpy.where(L>0)
   S=ax.scatter(J,I,50,L[I,J],edgecolor='none')
   CB=ax.figure.colorbar(S)
   ax.set_xlim(0,in_depth.shape[1])
   ax.set_ylim(0,in_depth.shape[0])
   CB.ax.set_title("Port number")
   logger.info("Writing to ports_all.png")
   figure.canvas.print_figure("ports_all.png")



   # Port diagnostics plot
   for i in range(len(kdports)) :
      figure.clf()
      ax=figure.add_subplot(111)
      cmap=matplotlib.pyplot.get_cmap("Greys_r")
      P=ax.pcolormesh(in_depth_m,cmap=cmap,edgecolor=".4",alpha=.5,linewidth=.05)
      iwidth = ilports[i] - ifports[i]+1
      jwidth = jlports[i] - jfports[i]+1
      #print ifports[i],jfports[i],iwidth,jwidth
      d=1
      if kdports[i] == 1 :
         xy = (ifports[i],jfports[i])
         jwidth=d
         c="r"
      elif kdports[i] == 2 :
         xy = (ifports[i],jfports[i])
         jwidth=d
         c="g"
      elif kdports[i] == 3 :
         xy = (ifports[i],jfports[i])
         iwidth=d
         c="b"
      elif kdports[i] == 4 :
         xy = (ifports[i],jfports[i])
         iwidth=d
         c="m"
      ax.add_patch(matplotlib.patches.Rectangle(xy,iwidth,jwidth, color=c,alpha=.5))
      ax.grid()

      ax.set_xlim(xy[0]-20,xy[0]+iwidth+20)
      ax.set_ylim(xy[1]-20,xy[1]+jwidth+20)
      ax.set_title("Port number %d" %i)
      fname="port_%03d.png"%(i+1)
      logger.info("Writing Diagnostics to %s"%fname)
      figure.canvas.print_figure(fname,bbox_inches='tight',dpi=200)




if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='Find ports')
   parser.add_argument('infile', type=str)
   args = parser.parse_args()

   main(args.infile)
