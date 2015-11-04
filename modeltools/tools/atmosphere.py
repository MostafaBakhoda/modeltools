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

class AtmosphericForcingError(Exception):
    """Base class for exceptions in this module."""
    pass

class FieldReaderError(Exception):
    """Base class for exceptions in this module."""
    pass

class FieldInterpolatorError(object) :
    """Base class for exceptions in this module."""
    pass

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


default_threshold=-5.

class FieldInterpolator(object) :
   def __init__(self,x,y,targetx,targety,x_is_longitude=True) :
      self._x=x
      self._y=y


      dx=self._x[1:]-self._x[:-1]
      absdiffx= numpy.abs((numpy.max(dx)-numpy.min(dx))/numpy.abs(numpy.max(dx)))
      dy=self._y[1:]-self._y[:-1]
      absdiffy= numpy.abs((numpy.max(dy)-numpy.min(dy))/numpy.abs(numpy.max(dy)))

      if len(self._x.shape) <> 1 or len(self._x.shape) <> 1 :
         raise FieldInterpolatorError,"Field coordinates (x) must be 1D and regular"
      elif absdiffx>1e-5 :
         raise FieldInterpolatorError,"Field coordinates (x) must be 1D and regular"

      if len(self._y.shape) <> 1 or len(self._y.shape) <> 1 :
         raise FieldInterpolatorError,"Field coordinates (y) must be 1D and regular"
      elif absdiffy>1e-5 :
         raise FieldInterpolatorError,"Field coordinates (y) must be 1D and regular"


      # Simplify internal calculations by flipping coordinate directions to be increasing
      self._dx=numpy.abs(dx[0])
      self._nx=self._x.size
      if dx[0] < 0 :
         self._x0=x[-1]
         self._flipx=True
         self._x=-self._x
      else :
         self._x0=x[0]
         self._flipx=False
      self._dy=numpy.abs(dy[0])
      self._ny=self._y.size
      if dy[0] < 0 :
         self._y0=y[-1]
         self._flipy=True
         self._y=-self._y
      else :
         self._y0=y[0]
         self._flipy=False


      if x_is_longitude and abs(self._dx*self._nx - 360.0) < 1e-4:
         self._wrapmode=("raise","wrap")
      self._targetx=targetx
      self._targety=targety





class FieldInterpolatorBilinear(FieldInterpolator)  :
   def __init__(self,x,y,fld,targetx,targety,x_is_longitude=True) :
      super(FieldInterpolatorBilinear,self).__init__(x,y,targetx,targety,x_is_longitude=x_is_longitude)
      self._calc_bilinear_weights()


   def _calc_pivot_points(self) :

      float_ip = (self._targetx-self._x0)/self._dx
      float_jp = (self._targety-self._y0)/self._dy

      # Pivot points
      self._ip=numpy.floor(float_ip).astype("int")
      self._jp=numpy.floor(float_jp).astype("int")

      # Remainder used in weights
      self._irem = numpy.remainder(float_ip,1.)
      self._jrem = numpy.remainder(float_jp,1.)

      # Pivot points displaced
      self._ip1=self._ip+1
      self._jp1=self._jp+1

      # Indexes into flattened array
      self._Ill = numpy.ravel_multi_index([self._jp ,self._ip  ],(self._ny,self._nx),mode=self._wrapmode)
      self._Ilr = numpy.ravel_multi_index([self._jp ,self._ip1 ],(self._ny,self._nx),mode=self._wrapmode)
      self._Iur = numpy.ravel_multi_index([self._jp1,self._ip1 ],(self._ny,self._nx),mode=self._wrapmode)
      self._Iul = numpy.ravel_multi_index([self._jp1,self._ip  ],(self._ny,self._nx),mode=self._wrapmode)


   def _calc_bilinear_weights(self) :
      self._calc_pivot_points()
      # Weights of points used in bilinear calculation
      self._wll = (1.-self._irem)*(1.-self._jrem) # Lower left
      self._wlr = (   self._irem)*(1.-self._jrem) # Lower right
      self._wur = (   self._irem)*(   self._jrem) # Upper right
      self._wul = (1.-self._irem)*(   self._jrem) # Upper left


   def interpolate(self,infld) :
      tmp = numpy.array(infld)

      # Flip array if needed for the interpolation (see init)
      if self._flipx :
         tmp = numpy.fliplr(tmp)
      if self._flipy :
         tmp = numpy.flipud(tmp)

      tmp=tmp.flatten()

      newfld = tmp[self._Ill]*self._wll 
      newfld+= tmp[self._Ilr]*self._wlr 
      newfld+= tmp[self._Iur]*self._wur
      newfld+= tmp[self._Iul]*self._wul

      newfld.shape=self._targetx.shape
      return newfld




class FieldInterpolatorRectBivariateSpline(FieldInterpolator)  :
   def __init__(self,x,y,field,targetx,targety,x_is_longitude=True) :
      super(FieldInterpolatorRectBivariateSpline,self).__init__(x,y,targetx,targety,x_is_longitude=x_is_longitude)
      tmp=numpy.array(field)
      if self._flipx :
         tmp=numpy.fliplr(tmp)
      if self._flipy :
         tmp=numpy.flipud(tmp)
      self._intobj = scipy.interpolate.RectBivariateSpline(self._x,self._y,tmp.transpose())

   def interpolate(self,fld) :
      return self._intobj(self._targetx,self._targety,grid=False)





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
         raise AtmosphericForcingError,"Only netcdf supported at the moment"


class NetcdfFieldReader(FieldReader) :
   def __init__(self,filenametemplate) :
      super(NetcdfFieldReader,self).__init__(filenametemplate)

   def open(self) :
      #self._nc = scipy.io.netcdf.netcdf_file(self._filename,"r")
      print "open started"
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

         print self._coordmap[varname].keys()
      print "open finished"

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


   def get_timestep(self,dt,varname) :
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






      


      

class AtmosphericForcing(object) :
   # These are the fileds this routine knows about, and can use to calculate fields
   _known_fields = [
         "10u",
         "10v",
         "2t",
         "2d",
         "msl",
         "ci",
         "tcc",
         "tp",
         "ro",
         "ssrd"]

   def __init__(self,configfile,forcing_dataset) :
      self._configfile = configfile
      self._forcing_dataset=forcing_dataset
      self._tree=xml.etree.ElementTree.ElementTree(file=configfile)

      elements = self._tree.findall('forcing_datasets/forcing_dataset[@name="%s"]'%forcing_dataset)
      if len(elements) > 1 : 
         msg = "Could not find unique dataset %s"%forcing_dataset
         raise AtmosphericForcingError,msg
      if len(elements) == 0 : 
         msg = "Could not find dataset %s"%forcing_dataset
         raise AtmosphericForcingError,msg

      # Parse the  Self._element attributes 
      self._rootPath = None
      self._element=elements[0]
      self._rootPath = self._element.attrib["rootPath"]
      if "rootPath" in self._element.attrib.keys(): self._rootPath = self._element.attrib["rootPath"]

      # Get format - only netcdf currently supported 
      if "format" in self._element.attrib.keys(): 
         self._format=self._element.attrib["format"]
         if self._format <> "netcdf" :
            raise AtmosphericForcingError,"Only netcdf supported at the moment"
      else :
         raise AtmosphericForcingError,"Format must be specified"

      # Parse the available fields
      elements = self._element.findall('field')
      self._fields   = {} 
      self._varnames = {}
      for i in elements :
         name             = i.attrib["name"]    # Variable names known to this module
         filenametemplate = i.attrib["file"]    # File known to this module
         varname          = i.attrib["varname"] # Variable name in file
         if self._rootPath is not None :
            filenametemplate = filenametemplate.replace("[rootPath]",self._rootPath)
         if name not in self._known_fields :
            msg = "Unknown field with name %s"%name
            raise AtmosphericForcingError,msg
         self._fields    [name] = FieldReader.get_field_reader(filenametemplate,self._format) 
         self._varnames[name] = varname

   def get_timestep(self,dt,varnames=None) :
      flddict={}
      for k,v in self._fields.items() :
         if varnames is None or k in varnames :
            flddict[k]=numpy.squeeze(v.get_timestep(dt,self._varnames[k]))
      return flddict


   def get_grid(self,dt,varnames=None) :
      flddict={}
      for k,v in self._fields.items() :
         if varnames is None or k in varnames :
            flddict[k]=v.get_grid(self._varnames[k],dt)
      return flddict


   def get_coords(self,dt,varnames=None) :
      flddict={}
      for k,v in self._fields.items() :
         if varnames is None or k in varnames :
            flddict[k]=v.get_coords(self._varnames[k],dt)
      return flddict


#   def get_proj4grid(self,dt,varnames=None) :
#      flddict={}
#      for k,v in self._fields.items() :
#         if varnames is None or k in varnames :
#            flddict[k]=v.get_proj4grid(self._varnames[k],dt)
#      return flddict


def qsw0(qswtime,daysinyear,cc,plat,plon) :
   #
   # --- -------------------------------------------------------------------
   # --- compute 24 hrs mean solar irrradiance at the marine surface layer
   # --- (unit: w/m^2)
   # --- -------------------------------------------------------------------
   #
   # --- Average number of days in year over a 400-year cycle (Gregorian Calendar)
   daysinyear400=365.2425
   #c --- set various quantities
   pi2=8.*numpy.arctan(1.)          #        2 times pi
   deg=360./pi2             #        convert from radians to degrees
   rad=pi2/360.             #        convert from degrees to radians
   eepsil=1.e-9             #        small number
   ifrac=24                 #        split each 12 hrs day into ifrac parts
   fraci=1./ifrac           #        1 over ifrac
   absh2o=0.09              # ---    absorption of water and ozone
   s0=1365.                 # w/m^2  solar constant
   radian=rad
#c
#c --- -------------------------------------------------------------------
#c --- compute 24 hrs mean solar radiation at the marine surface layer 
#c --- -------------------------------------------------------------------
#C --- KAL: TODO - adhere to hycom time setup
   day=numpy.mod(qswtime,daysinyear)    #0 < day < 364
   day=numpy.floor(day)
#c
   dangle=pi2*day/float(daysinyear)   #day-number-angle, in radians 
   if day<0. or day>daysinyear+1 :
      print 'qsw0: Error in day for day angle'
      print 'Day angle is ',day,daysinyear,qswtime
      raise NameError,"test"
      
# --- compute astronomic quantities -- 
   decli=.006918+.070257*numpy.sin(dangle)   -.399912*numpy.cos(dangle)      \
                +.000907*numpy.sin(2.*dangle)-.006758*numpy.cos(2.*dangle)   \
                +.001480*numpy.sin(3.*dangle)-.002697*numpy.cos(3.*dangle)

   sundv=1.00011+.001280*numpy.sin(dangle)   +.034221*numpy.cos(dangle)      \
                +.000077*numpy.sin(2.*dangle)+.000719*numpy.cos(2.*dangle)

# --- compute astronomic quantities

   sin2=numpy.sin(plat*radian)*numpy.sin(decli)
   cos2=numpy.cos(plat*radian)*numpy.cos(decli)
#
# --- split each day into ifrac parts, and compute the solar radiance for 
# --- each part. by assuming symmetry of the irradiance about noon, it
# --- is sufficient to compute the irradiance for the first 12 hrs of
# --- the (24 hrs) day (mean for the first 12 hrs equals then the mean
# --- for the last 12 hrs)
#
# --- TODO - This routine can also return daily varying solar heat flux
   scosz=0.
   stot=0.
   for npart in range(1,25) :
      bioday=day+(npart-.5)*fraci*.5
      biohr=bioday*86400.                #hour of day in seconds
      biohr=numpy.mod(biohr+43200.,86400.)    #hour of day;  biohr=0  at noon
      hangle=pi2*biohr/86400.            #hour angle, in radians
#
      cosz=numpy.maximum(0.,sin2+cos2*numpy.cos(hangle)) #cosine of the zenith angle
      scosz=scosz+cosz                     #  ..accumulated..
      srad =s0*sundv*cosz                  #extraterrestrial radiation
#
#         sdir=srad*0.7**(1./(cosz+eepsil))    #direct radiation component
#         sdir=srad * exp(-0.356674943938732447/(cosz+eepsil))         
# ---    KAL prevent underflow - .7^100 = 3x10^-16 
      sdir=srad*0.7**(numpy.minimum(100.,1./(cosz+eepsil)))    #direct radiation component
#
      sdif=((1.-absh2o)*srad-sdir)*.5               #diffusive radiation component
      altdeg=numpy.maximum(0.,numpy.arcsin(numpy.minimum(1.0,sin2+cos2)))*deg #solar noon altitude in degrees
      cfac=(1.-0.62*cc+0.0019*altdeg)               #cloudiness correction 
      ssurf=(sdir+sdif)*cfac
      stot=stot+ssurf

#     enddo
   scosz=scosz*fraci               #24-hrs mean of  cosz
   radfl0=stot*fraci               #24-hrs mean shortw rad in w/m^2
#
# --  Original formula was wrong ...
#     !cawdir(i,j)=1.-numpy.maximum(0.15,0.05/(scosz+0.15)) 
#     !cawdir(i,j)=1.-numpy.maximum(0.03,0.05/(scosz+0.15))  !Correction   - Mats
   cawdir=1.-numpy.minimum(0.15,0.05/(scosz+0.15))   #Correction 2 - KAL
#     enddo
#     enddo
#     enddo
#$OMP END PARALLEL DO
#
#     end subroutine qsw0

   return radfl0,cawdir





