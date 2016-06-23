#!/usr/bin/env python
##!/usr/bin/python -E
import modeltools.hycom
import logging
import argparse
import datetime
import numpy
import os

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


def main(blkdat_file):

   bp=modeltools.hycom.BlkdatParser(blkdat_file)

   kdm=bp["kdm"]
   nhybrd=bp["nhybrd"]
   nsigma=bp["nsigma"]
   logger.info("%d layers, %d hybrid layers, %d sigma layers"%(kdm,nhybrd,nsigma))

   # Check for "dp0k" 
   if bp["dp0k"] :
      dp0k=bp["dp0k"]
   else  :
      # Create dp0k (deep z-level) from parameters
      dp00=bp["dp00"]
      dp00x=bp["dp00x"]
      dp00f=bp["dp00f"]
      dp0k=[]
      for k in range(kdm) :
         dp0k.append(dp00*dp00f**k)
      dp0k=[min(elem,dp00x) for elem in dp0k]
      dp0k=numpy.array(dp0k)

   # Check for "ds0k" 
   if bp["ds0k"] :
      ds0k=bp["ds0k"]
   else  :
      # Create ds0k (shallow z-level)  from parameters
      ds00=bp["ds00"]
      ds00x=bp["ds00x"]
      ds00f=bp["ds00f"]
      ds0k=[]
      for k in range(kdm) :
         ds0k.append(ds00*ds00f**k)
      ds0k=[min(elem,ds00x) for elem in ds0k]
      ds0k=numpy.array(ds0k)
    
   for i in dp0k : print "%12.3f "%i,

   # Interface from dp
   intf0k=numpy.zeros((kdm+1))
   intf0s=numpy.zeros((kdm+1))
   for i in range(nsigma) :
      intf0s[i+1] = intf0s[i] + ds0k[i]
   for i in range(nhybrd) :
      intf0k[i+1] = intf0k[i] + dp0k[i]
   print "intf0k:",intf0k
   print "intf0s:",intf0s

   # x and depth arrays
   x      = numpy.linspace(0., 30.*1000,120)
   nx=x.size
   bottom=numpy.zeros((nx))
   bottom[:nx/2] = numpy.linspace(5., 100.,nx/2) + 80.*numpy.sin(x[0:nx/2] * 2*numpy.pi / 20000.)
   bottom[nx/2:] = numpy.linspace(bottom[nx/2-1], 2000.,nx/2)
   bottom[nx/2:] = bottom[nx/2:] + 250.*numpy.sin((x[nx/2:]-x[nx/2-1]) * 2*numpy.pi / 20000.)
   #bottom = numpy.linspace(5., 500.,x.shape[0])

   #Initialize interfaces
   intf=numpy.zeros((bottom.shape[0],intf0k.shape[0]))

   # Depth of deepest sigma interface, compared to bottom
   ideep      = intf0k[nsigma]
   ishallow   = intf0s[nsigma]
   f_ishallow = ishallow/bottom
   f_ideep    = ideep/bottom

   # Fixed shallow z_level (nsigma shallow z levels extend beyond ocean floor)
   Imask=f_ishallow>1.
   I=numpy.where(Imask)
   print "Imask",Imask
   intf[I[0],:nsigma+1] = intf0s[:nsigma+1]

   # Fixed deep z_level (nsigma deep z levels extend beyond ocean floor)
   Jmask=f_ideep>1.
   J=numpy.where(~Jmask)
   print "Jmask",Jmask
   #intf[J[0],:nsigma+1] = intf0k[:nsigma+1]
   intf[J[0],:nhybrd+1] = intf0k[:nhybrd+1]

   # Sigma coordinates where sigma-th deep z levels is below ocean floor, and sigma-th shallow z level is above ocean floor
   Kmask=numpy.logical_and(Jmask,numpy.logical_not(Imask))
   print "Kmask",Kmask
   K=numpy.where(Kmask)
   intf[K[0],:nsigma+1] = intf0k[:nsigma+1]
   tmp        = numpy.transpose(intf[K[0],:nsigma+1])*bottom[K[0]]/ideep
   intf[K[0],:nsigma+1] = tmp.transpose()

   #print intf[:,-1] 
   #print intf[:,nhybrd-1]


   intf = numpy.transpose(numpy.minimum(numpy.transpose(intf),bottom))
   #print x.shape
   #print intf.shape
   import matplotlib.pyplot as plt
   ax=plt.gca()
   ax.hold(True)
   #print x
   #print intf[Imask,-1],
   #print intf[Imask,kdm-nsigma+1]

   ax.fill_between(x,intf[:,nsigma]*-1.,0.,color="r",interpolate=False,where=Imask,label="Shallow z")
   ax.fill_between(x,intf[:,nhybrd]*-1.,0.,color="b",interpolate=False,where=~Jmask,label="Deep z")
   ax.fill_between(x,intf[:,nsigma]*-1.,0.,color="g",interpolate=False,where=Kmask,label="Sigma")

   for k in range(intf.shape[1]) :
      if (k+1)%5 == 0 :
         plt.plot(x,-intf[:,k],color="k",linestyle="--",label=str(k+1))
      else:
         plt.plot(x,-intf[:,k],color=".5")

   ax.plot(x,-bottom,lw=4,color="k")
   ax.legend(fontsize=6)
   plt.gcf().savefig("vcoord.png")
      
   ax.set_ylim(-200,0)
   plt.gcf().savefig("vcoord200.png")

   






   dp00i=bp["dp00i"]

   isotop=bp["isotop"]




if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='')
   parser.add_argument('file' , help="blkdat file")
   args = parser.parse_args()
   main(args.file)



