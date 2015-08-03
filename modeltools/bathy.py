import scipy.io.netcdf
import numpy
class GEBCO2014(object) :

   def __init__(self,filename="/work/shared/nersc/msc/ModelInput/bathymetry/GEBCO_2014/GEBCO_2014_2D.nc") :
      self._filename = filename
      self._nc = scipy.io.netcdf.netcdf_file(self._filename)

      self._lon = self._nc.variables["lon"][:]
      self._lat = self._nc.variables["lat"][:]
      self._elevation = self._nc.variables["elevation"][:]


   def regrid(self,grid) :
      pass



