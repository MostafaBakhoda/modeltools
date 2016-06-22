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


def main(infile,blo,bla) :

   bathy_threshold=0. # TODO

   # Read plon,plat
   gfile=abfile.ABFileGrid("regional.grid","r")
   plon=gfile.read_field("plon")
   plat=gfile.read_field("plat")
   gfile.close()

   # Read input bathymetri
   bfile=abfile.ABFileBathy(infile,"r",idm=gfile.idm,jdm=gfile.jdm)
   in_depth=bfile.read_field("depth",None)
   bfile.close()



   # Modify basin 
   depth=-1.*numpy.copy(in_depth)
   in_depth_m=numpy.ma.masked_where(depth>=bathy_threshold,in_depth)
   it=1
   while it==1 or numpy.count_nonzero(depth-depth_old) > 0 :
      depth_old = numpy.copy(depth)
      logger.info("Basin modifications ... pass %d"%(it))
      depth=modeltools.tools.remove_isolated_basins(plon,plat,depth,blo,bla,threshold=bathy_threshold)
      depth=modeltools.tools.remove_islets(depth,threshold=bathy_threshold)
      depth=modeltools.tools.remove_one_neighbour_cells(depth,threshold=bathy_threshold)
      logger.info("Modified %d points "%numpy.count_nonzero(depth-depth_old) )
      it+=1
   w5=numpy.copy(depth)


   w5[:,0]=bathy_threshold
   w5[:,-1]=bathy_threshold
   w5[0,:]=bathy_threshold
   w5[-1,:]=bathy_threshold

   # Mask data where depth below threshold
   w5=numpy.ma.masked_where(w5>=bathy_threshold,w5)

   # Print to HYCOM and CICE bathymetry files
   abfile.write_bathymetry("NEW",0,-w5,-bathy_threshold)
   kmt=numpy.where(~w5.mask,1.,0.)
   logger.info("Writing cice mask to file cice_kmt.nc")
   modeltools.cice.io.write_netcdf_kmt(kmt,"cice_kmt.nc")

   # Create netcdf file with all  stages for analysis
   logger.info("Writing bathymetry to file bathy_consistency.nc")
   ncid = netCDF4.Dataset("bathy_consistency.nc","w")
   ncid.createDimension("idm",w5.shape[1])
   ncid.createDimension("jdm",w5.shape[0])
   ncid.createVariable("lon","f8",("jdm","idm"))
   ncid.createVariable("lat","f8",("jdm","idm"))
   ncid.createVariable("old","f8",("jdm","idm"))
   ncid.createVariable("old_masked","f8",("jdm","idm"))
   ncid.createVariable("new","f8",("jdm","idm"))
   ncid.createVariable("new_masked","f8",("jdm","idm"))
   ncid.createVariable("modified","i4",("jdm","idm"))
   ncid.variables["lon"][:]=plon
   ncid.variables["lat"][:]=plat
   ncid.variables["old"][:]=in_depth
   ncid.variables["old_masked"][:]=in_depth_m
   ncid.variables["new"][:]=-depth
   ncid.variables["new_masked"][:]=-w5
   modmask=numpy.abs(in_depth-(-depth))>.1
   ncid.variables["modified"][:]=modmask.astype("i4")
   ncid.close()
   
   logger.info("Writing bathymetry plot to file newbathy.png")
   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   P=ax.pcolormesh(w5)
   figure.colorbar(P,norm=matplotlib.colors.LogNorm(vmin=w5.min(), vmax=w5.max()))
   I,J=numpy.where(numpy.abs(in_depth-(-depth))>.1)
   ax.scatter(J,I,20,"r")
   figure.canvas.print_figure("newbathy.png")





if __name__ == "__main__" :
   class PointParseAction(argparse.Action) :
     def __call__(self, parser, args, values, option_string=None):
       tmp = values[0].split(",")
       tmp = [float(elem) for elem in tmp[0:2]]
       tmp1= getattr(args, self.dest)
       tmp1.append(tmp)
       setattr(args, self.dest, tmp1)

   parser = argparse.ArgumentParser(description='Ensure consistenct of HYCOM bathy files')
   parser.add_argument('--basin_point', nargs="*", action=PointParseAction,default=[])
   parser.add_argument('infile', type=str)
   args = parser.parse_args()

   if args.basin_point :
      blo = [elem[0] for elem in args.basin_point]
      bla = [elem[1] for elem in args.basin_point]
   else :
      blo=[]
      bla=[]
   main(args.infile,blo,bla)
