#import scipy.io.netcdf
import netCDF4 
import numpy
import logging
import re  
import xml.etree.ElementTree
import datetime
import cfunits
import netcdftime
import scipy

class FieldReaderError(Exception):
    """Base class for exceptions in this module."""
    pass

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)




class FieldReader(object) :
   def __init__(self,filenametemplate) :
      self._filenametemplate = filenametemplate
      self._filename         = None


   def file_is_open(self,newfilename) :
      if self._filename is None :
         return False
      elif newfilename == self._filename :
         return True
      else :
         return False

   def find_timestep(self,dt,varname) :
      tmp = self._coordmap[varname]["time"]-dt
      tmp = [ i.days*86400 + i.seconds  for i in tmp]
      I = numpy.where(numpy.array(tmp)==0)[0]
      return I


   def get_grid(self,varname,dt) :
      self.open_if_needed(dt)
      if self._coordrank[varname]["lon"] > self._coordrank[varname]["lat"] :
         return numpy.meshgrid(self._coordvar["lon"],self._coordvar["lat"])
      else:
         return numpy.meshgrid(self._coordvar["lat"],self._coordvar["lon"])


   def get_coords(self,varname,dt) :
      self.open_if_needed(dt)
      return self._coordvar["lon"],self._coordvar["lat"]


   def get_proj4grid(self,varname,dt) :
      raise NotImplementedError,"get_proj4grid not implemented"



   @classmethod
   def get_field_reader(cls,filenametemplate,format) :
      if format == "netcdf" :
         return NetcdfFieldReader(filenametemplate)
      else :
         raise FieldReaderError,"Only netcdf supported at the moment"


class NetcdfFieldReader(FieldReader) :
   def __init__(self,filenametemplate) :
      super(NetcdfFieldReader,self).__init__(filenametemplate)

   def open(self) :
      #self._nc = scipy.io.netcdf.netcdf_file(self._filename,"r")
      #print "open started"
      self._nc = netCDF4.Dataset(self._filename,"r")

      # Read and map coordinate variables of all input vars
      self._coordvar={}
      self._coordmap={}
      self._coordrank={}
      for varname,var in self._nc.variables.items() :

         #print varname
         self._coordmap[varname] = {}
         self._coordrank[varname] = {}

         for inumber,i in enumerate(var.dimensions) :
            coordvar = self._nc.variables[i]
            unit = cfunits.Units(coordvar.units)

            if i not in self._coordvar.keys() :
               # Convert to datetime. use netcdftime as handling is better. 
               coordvals = numpy.array(self._nc.variables[i][:])
               if unit.isreftime :
                  tmp=netcdftime.utime(coordvar.units,calendar=coordvar.calendar)
                  self._coordvar["time"]=tmp.num2date(coordvals)
               elif unit.islongitude :
                  self._coordvar["lon"] = coordvar[:]
               elif unit.islatitude :
                  self._coordvar["lat"] = coordvar[:]
               else :
                  raise FieldReaderError,"Dont know how to handle coordinate variable %s"%i

            if unit.isreftime :
               self._coordmap[varname]["time"] = self._coordvar["time"]
               self._coordrank[varname]["time"] = inumber
            elif unit.islongitude :
               self._coordmap[varname]["lon"] = self._coordvar["lon"]
               self._coordrank[varname]["lon"] = inumber
            elif unit.islatitude :
               self._coordmap[varname]["lat"] = self._coordvar["lat"]
               self._coordrank[varname]["lat"] = inumber
            else :
               raise FieldReaderError,"Dont know how to handle coordinate variable %s"%i

         #print self._coordmap[varname].keys()
      #print "open finished"

   def close(self) :
      self._nc.close()
      self._filename=None

   def open_if_needed(self,dt) :
      # Open file if necessary
      newfilename=dt.strftime(self._filenametemplate)
      if not self.file_is_open(newfilename) : 
         if self._filename is not None : self._nc.close()
         self._filename = newfilename
         self.open()


   def get_timestep(self,varname,dt) :
      # Open file if necessary
      self.open_if_needed(dt)

      #TODO: Get rank of time variable in field. For now its assumed to be first (which is the normal)

      #Find timestep
      I = self.find_timestep(dt,varname)
      return self._nc.variables[varname][I,:] # Will actually read 3D, 4D etc




# TODO:
# To implement other readers:
# Subclass FieldReader
# on open (or init) : Define self._coordvar (gets coordinate variables in file)
#                     Time coordinate must be defined as datetime
# on open (or init) : Define self._coordmap (links variables to coord variables)
# Implement init, open, get_timestep and close






      
