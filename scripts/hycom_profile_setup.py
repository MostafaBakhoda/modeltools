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

   dp00=bp["dp00"]
   dp00x=bp["dp00x"]
   dp00f=bp["dp00f"]

   ds00=bp["ds00"]
   ds00x=bp["ds00x"]
   ds00f=bp["ds00f"]

   # Create dp0k (shallow z-level) and ds0k (deep z-level)
   dp0k=[]
   ds0k=[]
   for k in range(1,kdm+1) :
      dp0k.append(dp00*dp00f**(k-1))
      ds0k.append(ds00*ds00f**(k-1))
   dp0k=[min(elem,dp00x) for elem in dp0k]
   ds0k=[min(elem,ds00x) for elem in ds0k]
   dp0k=numpy.array(dp0k)
   ds0k=numpy.array(ds0k)

   intf0k=numpy.zeros((kdm+1))
   intf0s=numpy.zeros((kdm+1))
   for i in range(1,kdm) :
      print i
      intf0k[i] = intf0k[i-1] + dp0k[i]
      intf0s[i] = intf0s[i-1] + ds0k[i]
   print intf0k
   print intf0s

   ideep = numpy.sum(intf0k[0:nsigma]) # Starts
   ishallow = numpy.sum(intf0s[0:nsigma]) # Ends

   x      = numpy.linspace(0., 30.*1000,60)
   bottom = numpy.linspace(5., max(ishallow,ideep)+500.,x.shape[0]) + 100.*numpy.sin(x * 2*numpy.pi / 20000.)
   #dp=numpy.zeros((intf0k.shape[0],bottom.shape[0]))
   intf=numpy.zeros((bottom.shape[0],intf0k.shape[0]))
   print intf.shape

   print ideep
   print ishallow
   print bottom
   f_ishallow = ishallow/bottom
   f_ideep = ideep/bottom

   # Fixed shallow z_level
   print f_ishallow.shape
   I=numpy.where(f_ishallow>=1.)
   print intf[I[0],:].shape
   print I[0]
   intf[I[0],:] = intf0s


   # Fixed deep z_level
   I=numpy.where(f_ideep<=1.)
   print I[0]
   intf[I[0],:] = intf0k

   #In between
   I=numpy.where(numpy.logical_and(f_ishallow<1.,f_ideep>1.))
   print I[0]
   intf[I[0],:] = intf0k
   tmp        = numpy.transpose(intf[I[0],:])*bottom[I[0]]/ideep
   intf[I[0],:] = tmp.transpose()

   intf = numpy.transpose(numpy.minimum(numpy.transpose(intf),bottom))



   import matplotlib.pyplot as plt
   plt.plot(-intf)
   ax=plt.gca()
   ax.plot(-bottom,lw=4)
   print ax
   #plt.ylabel('some numbers')
   plt.gcf().savefig("tst.png")
      

   






   dp00i=bp["dp00i"]

   isotop=bp["isotop"]




if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='')
   parser.add_argument('file' , help="blkdat file")
   args = parser.parse_args()
   main(args.file)



