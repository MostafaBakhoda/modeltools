import scipy.io.netcdf
import numpy
class GEBCO2014(object) :

   def __init__(self,filename="/work/shared/nersc/msc/ModelInput/bathymetry/GEBCO_2014/GEBCO_2014_2D.nc") :
      self._filename = filename
      self._nc = scipy.io.netcdf.netcdf_file(self._filename)

      self._lon = self._nc.variables["lon"][:]
      self._lat = self._nc.variables["lat"][:]
      self._elevation = self._nc.variables["elevation"][:]


      self._dlon= ( self._lon[-1] - self._lon[0] ) / (self._lon.size-1)
      self._dlat= ( self._lat[-1] - self._lat[0] ) / (self._lat.size-1)

      print self._dlon, self._lon[1] - self._lon[0]
      self._dlat= ( self._lat[-1] - self._lat[0] ) / (self._lat.size-1)
      print self._dlat, self._lat[1] - self._lat[0]


   def regrid(self,lon,lat) :

      # As of now, only nearest 

      i = numpy.floor(numpy.mod(lon-self._lon[0],360.) / self._dlon).astype("int")
      j = numpy.floor((lat-self._lat[0]) / self._dlat).astype("int")
      print i.min(),i.max(),self._lon.shape
      print j.min(),j.max(),self._lat.shape
      print i.shape

      w = self._elevation[j,i]
      return w



