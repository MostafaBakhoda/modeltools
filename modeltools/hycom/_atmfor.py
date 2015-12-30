import datetime
import cfunits
import numpy
#import modeltools.grid
import modeltools.hycom
import modeltools.tools
import modeltools.forcing.atmosphere
import modeltools.tools
#from mpl_toolkits.basemap import Basemap, shiftgrid
import logging


_loglevel=logging.DEBUG
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate=False # Dont propagate to parent in hierarchy (determined by "." in __name__)


# Plot scalar field on model grid
def plot_fig(fld,dt,varname,filename) :
   from matplotlib.backends.backend_agg import FigureCanvasAgg
   from matplotlib.figure import Figure
   fig = Figure(figsize=(5,4), dpi=100)
   ax = fig.add_subplot(111)
   canvas = FigureCanvasAgg(fig)
   P=ax.pcolor(fld)
   fig.colorbar(P)
   ax.set_title("%s at %s"%(varname,str(dt)))
   canvas.print_figure(filename)

def atmfor(start,end,af,grid_file="regional.grid",blkdat_file="blkdat.input",plot_diag=False) :

   # Open hycom grid file, read longitude and latitude@
   # TODO: HYCOM-specific
   za = modeltools.hycom.io.ABFileRegionalGrid(grid_file,"r")
   mlon = za.readfield("plon")
   mlat = za.readfield("plat")
   Nx=mlon.shape[1]
   Ny=mlon.shape[0]
   za.close()

   # parse blkdat to get yearflag
   # TODO: HYCOM-specific
   blkd = modeltools.hycom.BlkdatParser(blkdat_file)
   yrflag = blkd["yrflag"]
   wndflg = blkd["wndflg"]
   lwflag  = blkd["lwflag"]

   # Main loop 
   ffiles={}
   dt = start
   while dt <= end :
       
       logger.info("Reading at %s"%str(dt))
       #print af.known_names

       # Read variables
       af.get_timestep(dt)

       print "msl shape",af["msl"].data.shape
       print "2d.shape ",af["2d"].data.shape

       # Estimate dependent variable on native grid
       # radflx is downwelling longwave radiation
       # TODO: HYCOM-specific
       if wndflg in [1,2,3] :
           af.calculate_windstress()
           af.calculate_windspeed()
           af.calculate_ustar()
       
       if "vapmix" not in af.known_names_explicit : af.calculate_vapmix()
       if "ssrd"   not in af.known_names_explicit : af.calculate_ssrd()
       if lwflag == -1 :
           if "strd"   not in af.known_names_explicit : af.calculate_strd()
       else :
           raise ValueError,"TODO: lwflag<>-1 not supported"

       # Open output files. Dict uses "known name" when mapping to file object
       # TODO: HYCOM-specific
       print modeltools.hycom.atmosphere_variable_names.keys()
       if dt == start :
           for k,vname in modeltools.hycom.atmosphere_variable_names.items() :
               if k in af.known_names :
                   ffiles[k]=modeltools.hycom.io.ABFileForcing("forcing.%s"%vname,"w",idm=Nx, jdm=Ny,
                                                               cline1="ERA Interim",cline2=vname)
                   
       # Interpolation of all fields and unit conversion
       newfld={}
       for kn in [elem for elem in af.known_names if elem in ffiles.keys()] :

           # Coordinates
           lo,la=af[kn].coords

           # Read and convert field to units used by HYCOM
          # TODO: HYCOM-specific
           fld=af[kn].data_to_unit(modeltools.hycom.atmosphere_variable_cfunit[kn])

           #TODO: : Possible to choos interpolation here
           fi=modeltools.tools.FieldInterpolatorBilinear(lo,la,fld,mlon,mlat)
           newfld[kn]=fi.interpolate(fld)

       # Do rotation of u and v components if this the first component of a vector field
       for kn in af.known_names :
           if kn in modeltools.hycom.atmosphere_variable_vectors.keys() :
               knu,knv = modeltools.hycom.atmosphere_variable_vectors[kn]
               ur,vr=modeltools.tools.rotate_vector(newfld[knu],newfld[knv],mlon,mlat)
               newfld[knu]=ur
               newfld[knv]=vr

       # Loop over open files 
       # TODO: HYCOM specific
       for kn in ffiles.keys() :
               
           # Name used by hycom
           vname=modeltools.hycom.atmosphere_variable_names[kn]

           # Write to hycom file
           newdt=af[kn].time
           ord_day,hour,isec=modeltools.hycom.datetime_to_ordinal(newdt,yrflag)
           dtime=modeltools.hycom.dayfor(newdt.year,ord_day,hour,yrflag)
           ffiles[kn].writefield(newfld[kn],newfld[kn],vname,dtime,af.timestep_in_days)

           # Write diagnostics, if requested
           if plot_diag :
              tmp="forcing.%s."%vname 
              tmp=tmp+"%Y%m%d%H.png"
              tmp=newdt.strftime(tmp)
              print "plotting %s"%tmp
              plot_fig(newfld[kn],newdt,vname,tmp)

       # Increase time
       dt = dt + af.timestep
       

               
   for kn in ffiles.keys() :
       ffiles[kn].close()
   af=[]



