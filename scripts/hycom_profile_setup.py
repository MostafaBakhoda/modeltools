#!/usr/bin/env python
##!/usr/bin/python -E
import modeltools.hycom
import argparse
import datetime
import numpy
import os


def main(blkdat_file):

   bp=modeltools.hycom.BlkdatParser(blkdat_file)

   kdm=bp["kdm"]
   nhybrd=bp["nhybrd"]
   nsigma=bp["nsigma"]

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
      intf0k[i+1] = intf0k[i] + dp0k[i]
      intf0s[i+1] = intf0s[i] + ds0k[i]
   #
   for i in range(nsigma,kdm) :
      intf0k[i+1] = intf0k[i] + dp0k[nsigma-1]
      intf0s[i+1] = intf0s[i] + ds0k[nsigma-1]
   #print intf0k
   #print intf0s

   #ideep = numpy.sum(intf0k[0:nsigma+1]) # Starts
   #ishallow = numpy.sum(intf0s[0:nsigma+1]) # Ends
   ideep    = intf0k[nsigma+1]
   ishallow = intf0s[nsigma+1]

   maxdeep = max(numpy.max(ishallow),numpy.max(ideep))
   #print "maxdeep=",maxdeep
#   print ishallow

   x      = numpy.linspace(0., 30.*1000,120)
   nx=x.size
   bottom=numpy.zeros((nx))
   bottom[:nx/2] = numpy.linspace(5., 100.,nx/2) + 80.*numpy.sin(x[0:nx/2] * 2*numpy.pi / 20000.)

   bottom[nx/2:] = numpy.linspace(bottom[nx/2-1], 600.,nx/2)
   bottom[nx/2:] = bottom[nx/2:] + 250.*numpy.sin((x[nx/2:]-x[nx/2-1]) * 2*numpy.pi / 20000.)
   #bottom = numpy.hstack((bottom.transpose(),bottom2.transpose())).transpose()
   #print bottom.shape
   ##raise NameError,"test"

   #bottom = numpy.linspace(5., 500.,x.shape[0])
   intf=numpy.zeros((bottom.shape[0],intf0k.shape[0]))

   f_ishallow = ishallow/bottom
   f_ideep = ideep/bottom

   # Fixed shallow z_level
   Imask=f_ishallow>=1.
   I=numpy.where(Imask)
   intf[I[0],:] = intf0s

   # Fixed deep z_level
   Jmask=f_ideep<=1.
   J=numpy.where(Jmask)
   intf[J[0],:] = intf0k

   #In between
   Kmask=numpy.logical_and(f_ishallow<1.,f_ideep>1.)
   K=numpy.where(Kmask)
   intf[K[0],:] = intf0k
   tmp        = numpy.transpose(intf[K[0],:])*bottom[K[0]]/ideep
   intf[K[0],:] = tmp.transpose()

   intf = numpy.transpose(numpy.minimum(numpy.transpose(intf),bottom))
   #print x.shape
   #print intf.shape
   import matplotlib.pyplot as plt
   ax=plt.gca()
   ax.hold(True)
   #print x
   #print intf[Imask,-1],
   #print intf[Imask,kdm-nsigma+1]

   #ax.fill_between(x,intf[:,-1]*-1.,0.,color=".1",interpolate=False,where=Imask)
   ax.fill_between(x,intf[:,nsigma+1]*-1.,0.,color="r",interpolate=False,where=Imask,label="Shallow z")
   ax.fill_between(x,intf[:,nsigma+1]*-1.,0.,color="g",interpolate=False,where=Jmask,label="Sigma")
   ax.fill_between(x,intf[:,nsigma+1]*-1.,0.,color="b",interpolate=False,where=Kmask,label="Deep z")
   #ax.fill_between(x,intf[:,-1]*-1.,0.,color=".8",interpolate=False,where=Jmask)

   for k in range(1,intf.shape[1]) :
      #plt.plot(x,-intf[:,k],label=str(k))
      if k%5 == 0 :
         plt.plot(x,-intf[:,k],color=".5",linestyle="--",label=str(k+1))
      else:
         plt.plot(x,-intf[:,k],color=".5")

   ax.plot(x,-bottom,lw=4,color="k")
   ax.legend(fontsize=6)
   ax.set_ylim(-50,0)
   plt.gcf().savefig("tst.png")
      

   






   dp00i=bp["dp00i"]

   isotop=bp["isotop"]




if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='')
   parser.add_argument('file' , help="blkdat file")
   args = parser.parse_args()
   main(args.file)



