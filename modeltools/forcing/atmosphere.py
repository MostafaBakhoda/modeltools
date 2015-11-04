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
import modeltools.tools.indata

_all_known_names = [
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

# Units used by internal calculations
_assumed_units = {
      "10u":"m s**-1",
      "10v":"m s**-1",
      "2t":"K",
      "2d":"K",
      "msl":"Pa",
      "ci":"1",
      "tcc":"1",
      "tp":"m s**-1",
      "ro":"1",
      "ssrd":"W m**2 s**-1"
      }

class AtmosphericForcingError(Exception):
    """Base class for exceptions in this module."""
    pass

class ForcingField(object) :
   def __init__(self,filenametemplate,varname,unit,format) :
      self._filenametemplate= filenametemplate
      self._varname         = varname
      self._unit            = unit
      self._cfunit          = cfunits.Units(units=self._unit)
      self._format          = format
      self._fieldreader = modeltools.tools.indata.FieldReader.get_field_reader(filenametemplate,format) 

   @property
   def filenametemplate(self) : 
      return self._filenametemplate

   @property
   def varname(self) : 
      return self._varname

   @property
   def unit(self) : 
      return self._unit

   @property
   def format(self) : 
      return self._format

   def get_timestep(self,dt,unit=None) : 

      # Do unit conversion to get correct output unit
      if unit is not None :
         mycfunit = cfunits.Units(unit)
      tmp = numpy.squeeze(self._fieldreader.get_timestep(self._varname,dt))
      if not self._cfunit.equals(mycfunit) :
         print "Unit conversion:",self.varname,"unit=",self._cfunit, "targetunit=", mycfunit
         print "Unit conversion:max=",tmp.max()
         tmp=cfunits.Units.conform(tmp,self._cfunit,mycfunit)
         print "Unit conversion:max after=",tmp.max()
      return tmp

   def get_coords(self,dt) : 
      return self._fieldreader.get_coords(self._varname,dt)

   def get_grid(self,dt) : 
      return self._fieldreader.get_grid(self._varname,dt)


class AtmosphericForcing(object) :
   # These are the fileds this routine knows about, and can use to calculate new fields

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
      self._timestep = self._element.attrib["timestep"]
      if "rootPath" in self._element.attrib.keys(): self._rootPath = self._element.attrib["rootPath"]
      if self._timestep[-1] == "h" :
         self._timestep = datetime.timedelta(hours=int(self._timestep[:-1]))
      else :
         raise AtmosphericForcingError,"time step must be specified in hours (hours + letter 'h')"

      # Get format - only netcdf currently supported 
      if "format" in self._element.attrib.keys(): 
         self._format=self._element.attrib["format"]
         if self._format <> "netcdf" :
            raise AtmosphericForcingError,"Only netcdf supported at the moment"
      else :
         raise AtmosphericForcingError,"Format must be specified"

      # Parse the available fields and create forcingfield class
      elements = self._element.findall('field')
      self._fields   = {} 
      for i in elements :
         name             = i.attrib["known_name"]    # Variable names known to this module
         filenametemplate = i.attrib["file"]    # File known to this module
         varname          = i.attrib["varname"] # Variable name in file
         unit             = i.attrib["unit"]    # 
         if self._rootPath is not None :
            filenametemplate = filenametemplate.replace("[rootPath]",self._rootPath)
         if name not in _all_known_names :
            msg = "Unknown field with name %s"%name
            raise AtmosphericForcingError,msg

         self._fields[name] = ForcingField(filenametemplate,varname,unit,self._format) 

   def get_timestep(self,dt,varnames=None) :
      flddict={}
      for k,v in self._fields.items() :
         if varnames is None or k in varnames :
            flddict[k]=v.get_timestep(dt,unit=_assumed_units[k])
            #print k,v.varname,v.unit,numpy.max(flddict[k]),_assumed_units[k]
      self._fielddata=flddict
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
            flddict[k]=v.get_coords(dt)
      return flddict


#   def get_proj4grid(self,dt,varnames=None) :
#      flddict={}
#      for k,v in self._fields.items() :
#         if varnames is None or k in varnames :
#            flddict[k]=v.get_proj4grid(self._varnames[k],dt)
#      return flddict

   @property
   def timestep(self) :
      return self._timestep

   @property
   def varnames(self) :
      return [ elem.varname for elem in self._fields.values ]

   @property
   def known_names(self) :
      return self._fields.keys()


   def calculate_windstress(self) :
      if "10u" in self.known_names and "10v" in self.known_names :
         self._fielddata["taux"], self._fielddata["tauy"] = calculate_windstress(self["10u"],self["10v"])
      else :
         raise AtmosphericForcingError,"Can not calculate wind stress without 10 meter winds"


   def calculate_windspeed(self) :
      if "10u" in self.known_names and "10v" in self.known_names :
         self._fielddata["wspd"] = numpy.sqrt(self["10u"]**2+self["10v"]**2)
      else :
         raise AtmosphericForcingError,"Can not calculate wind speed without 10 meter winds"

   def calculate_ustar(self) :
      if "taux" in self.known_names and "taux" in self.known_names :
         self._fielddata["ustar"] = numpy.sqrt((self["taux"]**2+self["tauy"]**2)*1e-3)
      else :
         raise AtmosphericForcingError,"Can not calculate wind stress without 10 meter winds"

   def calculate_vapmix(self) :
      if "2t" in self.known_names and "msl" in self.known_names and "2d" in self.known_names:
         e = satvap(self["2d"])
         self._fielddata["vapmix"] = vapmix(e,self["msl"])
      else :
         raise AtmosphericForcingError,"Can not calculate wind stress without 10 meter winds"
     
     
     
     
     
     
     
     

   def __getitem__(self,i) :
      return self._fielddata[i]





def windstress(uwind,vwind) :
   ws=numpy.sqrt(uwind**2+vwind**2)
   if karalight :
      wndfac=numpy.max(2.5,numpy.min(32.5,ws))
      cd_new = 1.0E-3*(.692 + .0710*wndfac - .000700*wndfac**2)
   else :
      wndfac=(1.+sign(1.,ws-11.))*.5
      cd_new=(0.49+0.065*ws)*1.0e-3*wndfac+cd*(1.-wndfac)
   wfact=ws*airdns*cd_new
   taux = uwind*wfact
   tauy = vwind*wfact
   ustar = numpy.sqrt((taux**2+tauy**2)*1e-3)
   return taux, tauy



def vapmix(e,p) :
   # Input is
   # e = vapour pressure = saturation vapor pressure at dewpoint temperature 
   # p = air pressure
   vapmix = 0.622 * e / (p-e)
   return vapmix


def satvap(t) :
   # This function calculates the saturation vapour pressure
   # [Pa] from the temperature [deg K].
   # Modified: Anita Jacob, June '97
   #
   # Input: t: temperature [deg K]
   # Output: satvap: saturation vapour pressure at temp. t
   #
   # es(T) = C1 * exp(C3*(T - T0)/(T - C4)) from ECMWF manual
   #data c1/610.78/,t00/273.16/
   c1=610.78
   t00=273.16

   #if (t < t00) then
   #   c3 = 21.875
   #   c4 = 7.66
   #else
   #   c3 = 17.269
   #   c4 = 35.86
   #endif
   c3 = numpy.where(t < t00,21.875,7.66)
   c4 = numpy.where(t < t00,17.269,35.86)
   aa = c3 * (t - t00)
   bb = t - c4
   cc=aa/bb

   #if (cc < -20.0) then
   #   satvap=0.0
   #else
   #   satvap = c1 * exp(aa/bb)
   satvap=numpy.where(cc<-20.0,0.0,c1 * numpy.exp(aa/bb))
   return satvap

def  relhumid(sva,svd,msl) :
   # This routine calculates the relative humidity by the 
   # dew point temperature and the mean sea level pressure.
   # Modified: Anita Jacob, June '97
   # Input:
   #    sva: saturatn vapour press at air temp [K]
   #    svd: saturatn vapour press at dew pt temp [K]
   #    msl: pressure at mean sea level [Pa]
   # Output: 
   #   relhumid: Relative Humidity

   # We use the Tetens formula:
   # es(T) = C1 * exp(C3*(T - T0)/(T - C4)) from ECMWF manual
   #              es(Tdew)        p - es(Tair)
   # RH = 100 *  -----------   *  ------------
   #             p - es(tdew)       es(Tair)
   aaa=msl - svd
   aaa = svd/aaa
   bbb = (msl - sva)/sva
   relhumid = 100. * aaa * bbb
   return relhumid



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





