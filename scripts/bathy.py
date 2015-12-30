#!/usr/bin/env python
import modeltools.hycom
import argparse
import datetime
import modeltools.grid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import modeltools.forcing.bathy
import modeltools.hycom.io
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

if __name__ == "__main__" :
   class PointParseAction(argparse.Action) :
     def __call__(self, parser, args, values, option_string=None):
       tmp = values[0].split(",")
       tmp = [float(elem) for elem in tmp[0:2]]
       tmp1= getattr(args, self.dest)
       tmp1.append(tmp)
       setattr(args, self.dest, tmp1)

   parser = argparse.ArgumentParser(description='Prepare HYCOM bathy files')
   parser.add_argument('proj4_string', help="proj4 string ")
   parser.add_argument('--basin_point', nargs="*", action=PointParseAction,default=[])
   parser.add_argument('--shapiro_passes', type=int, default=1, help="Number of shapiro passes to apply")
   parser.add_argument('ll_lon', type=float, help='lower left corner longitude')
   parser.add_argument('ll_lat', type=float, help='lower left corner latitude')
   parser.add_argument('dx',     type=int, help='Grid spacing 1st index [m]')
   parser.add_argument('dy',     type=int, help='Grid spacing 2nd index [m]')
   parser.add_argument('nx',     type=int, help='Grid dimension 1st index []')
   parser.add_argument('ny',     type=int, help='Grid dimension 2nd index []')
   args = parser.parse_args()

   if args.basin_point :
      blo = [elem[0] for elem in args.basin_point]
      bla = [elem[1] for elem in args.basin_point]
   else :
      blo=[]
      bla=[]
   #print blo,bla
   #raise NameError,"test"

   bathy_threshold=-5.

   # Create grids and write to file
   #proj4_string="+proj=stere  +lon_0=-45 +lat_0=90 +lat_ts=80 +ellps=sphere"
   #grid1=modeltools.grid.Proj4Grid("+proj=stere  +lon_0=-45 +lat_0=90 +lat_ts=80 +ellps=sphere",
   #                              -89.5,45.5,20000,20000,400,300)
   grid1=modeltools.grid.Proj4Grid(args.proj4_string,args.ll_lon,args.ll_lat,args.dx,args.dy,
         args.nx,args.ny)
   modeltools.hycom.io.write_regional_grid(grid1)
   modeltools.cice.io.write_netcdf_grid(grid1,"cice_grid.nc")

   # Interpolation of bathymetry
   #gebco = modeltools.bathy.GEBCO2014("/Users/knutal/Bathymetry/GEBCO/GEBCO_2014_2D.nc")
   gebco = modeltools.forcing.bathy.GEBCO2014(filename="/work/shared/nersc/msc/ModelInput/bathymetry/GEBCO_2014/GEBCO_2014_2D_median8km.nc") 
   lon,lat=grid1.pgrid()
   w2=gebco.regrid(lon,lat,width=grid1.dx)
   #w2=numpy.ma.masked_where(w2>=bathy_threshold,w2)

   # Run shapiro filter on interpolated data to remove 2 DeltaX noise
   w3=numpy.copy(w2)
   for i in range(args.shapiro_passes):
      logger.info("Shapiro filter ... pass %d"%(i+1))
      w3=modeltools.tools.shapiro_filter(w3,threshold=bathy_threshold)
   #print w3[0:10,200]

   # Modify basin 
   w4=numpy.copy(w3)
   it=1
   while it==1 or numpy.count_nonzero(w4-w4old) > 0 :
      w4old = numpy.copy(w4)
      logger.info("Basin modifications ... pass %d"%(it))
      w4=modeltools.tools.remove_isolated_basins(lon,lat,w4,blo,bla,threshold=bathy_threshold)
      w4=modeltools.tools.remove_islets(w4,threshold=bathy_threshold)
      w4=modeltools.tools.remove_one_neighbour_cells(w4,threshold=bathy_threshold)
      logger.info("Modified %d points "%numpy.count_nonzero(w4-w4old) )
      it+=1
   w5=numpy.copy(w4)

   # Mask data where depth below threshold
   w5=numpy.ma.masked_where(w5>=bathy_threshold,w5)
   #print w5[0:10,200]

   # Print to HYCOM and CICE bathymetry files
   # TODO: Find nice generic name for hycom
   modeltools.hycom.io.write_bathymetry("TPTa0.20",1,-w5,-bathy_threshold)
   kmt=numpy.where(~w5.mask,1.,0.)
   modeltools.cice.io.write_netcdf_kmt(kmt,"cice_kmt.nc")

   # Create netcdf file with all  stages for analysis
   logger.info("Writing bathymetry to file gridgen_stages.nc")
   ncid = netCDF4.Dataset("gridgen_stages.nc","w")
   ncid.createDimension("idm",w5.shape[1])
   ncid.createDimension("jdm",w5.shape[0])
   ncid.createVariable("lon","f8",("jdm","idm"))
   ncid.createVariable("lat","f8",("jdm","idm"))
   ncid.createVariable("h_1","f8",("jdm","idm"))
   ncid.createVariable("h_2","f8",("jdm","idm"))
   ncid.createVariable("h_3","f8",("jdm","idm"))
   ncid.variables["lon"][:]=lon
   ncid.variables["lat"][:]=lat
   ncid.variables["h_1"][:]=w2
   ncid.variables["h_2"][:]=w3
   ncid.variables["h_3"][:]=w5
   ncid.close()
   
   # Show grid on a map
   logger.info("grid shown in grid.png")
   grid1.plotgrid(2.).canvas.print_figure("grid.png")

   # Show some grid statistics 
   scpx=grid1.scpx()
   scpy=grid1.scpy()

   sx = (w2[1:,:-1]-w2[:-1,:-1])/scpy[1:,:-1]
   sy = (w2[:-1,1:]-w2[:-1,:-1])/scpx[:-1,1:]
   grad = sx + sy
   slopefac=numpy.sqrt(sx**2+sy**2)
   logger.info("Maximum slope factor after interpolation =%.4f"%slopefac.max())

   sx = (w5[1:,:-1]-w5[:-1,:-1])/scpy[1:,:-1]
   sy = (w5[:-1,1:]-w5[:-1,:-1])/scpx[:-1,1:]
   grad = sx + sy
   slopefac=numpy.sqrt(sx**2+sy**2)
   logger.info("Maximum slope factor after smoothing    =%.4f"%slopefac.max())

   
   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   P=ax.pcolor(slopefac)
   figure.colorbar(P,norm=matplotlib.colors.LogNorm(vmin=w5.min(), vmax=w5.max()))
   #ax.contour(w5)#,[-10.,-100.,-500.,-1000.])
   #ax.set_title("Slope fac in color, depth contours in black")
   logger.info("Slope factor in slopefac.png")
   figure.canvas.print_figure("slopefac.png")




