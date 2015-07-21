""" Module for generating grids, mainly for use by models.  """

import pyproj
import numpy
import logging
import re 
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


class Grid(object) :
   pass

class ConformalGrid(Grid) :
   pass

class Proj4Grid(Grid) :

   def __init__(self,proj4string,ll_lon,ll_lat,dx,dy,Nx,Ny) :

      


      self._proj4string=proj4string
      self._proj=pyproj.Proj(proj4string)
      self._initgrid(ll_lon,ll_lat,dx,dy,Nx,Ny)


   def _initgrid(self,ll_lon,ll_lat,dx,dy,Nx,Ny) :
      self._dx=dx
      self._dy=dy
      self._Nx=Nx
      self._Ny=Ny
      self._ll_lon=ll_lon
      self._ll_lat=ll_lat

      # Calculate 0,0 (LL=Lower Left) in projection coordinates. 
      self._ll_x,self._ll_y = self._proj(ll_lon,ll_lat)

      # Create grid. This is the P-grid. Note increase of stencil - used to calculate grid sizes
      self._x=self._ll_x + numpy.linspace(-dx,dx*(Nx),Nx+2)
      self._y=self._ll_y + numpy.linspace(-dy,dy*(Ny),Ny+2)
      self._X,self._Y = numpy.meshgrid(self._x,self._y)
      #print self._X.shape,self._x.shape,self._y.shape

      tmp= self._proj(self._X,self._Y,inverse=True)
      logger.debug("Initialized P-grid using projection %s"%self._proj4string)
      logger.debug("Lower left corner lon/lat of grid: (%.3g,%.3g)" % (ll_lon,ll_lat))
      logger.debug("Grid spacing in projection coords: (%.3g,%.3g)" % (self._dx,self._dy))
      logger.debug("Number of grid Nodes in x/y      : (%5d,%5d)"   % (self._Nx,self._Ny))

      logger.debug("Min   x projection coordinate = %.3g"%self._x.min())
      logger.debug("Max   x projection coordinate = %.3g"%self._x.max())
      logger.debug("Min   y projection coordinate = %.3g"%self._y.min())
      logger.debug("Max   y projection coordinate = %.3g"%self._y.max())

      logger.debug("Min lon = %.3g"%tmp[0].min())
      logger.debug("Max lon = %.3g"%tmp[0].max())
      logger.debug("Min lat = %.3g"%tmp[1].min())
      logger.debug("Max lat = %.3g"%tmp[1].max())

      # Should insist that the ellipse is spherical
      #if bla bla bla :
      #   msg = "The geoid used in the projection is not spherical. Aborting"
      #   logger.error(msg)
      #   raise ValueError,msg

   

   def pgrid(self,extended=False) : return self._grid(0.,0.,extended)

   def ugrid(self,extended=False) : return self._grid(-0.5*self._dx,0.,extended)

   def vgrid(self,extended=False) : return self._grid(0.,-0.5*self._dy,extended)

   def qgrid(self,extended=False) : return self._grid(-0.5*self._dx,-0.5*self._dy,extended)
      
      
   def _grid(self,deltax,deltay,extended=False) : 
      if extended :
         return self._proj(self._X-deltax,self._Y-deltay,inverse=True)
      else :
         return self._proj(self._X[1:-1,1:-1]-deltax,self._Y[1:-1,1:-1]-deltay,inverse=True)

   
   @property
   def width(self) : return self._Nx*self._dx

   @property
   def height(self) : return self._Ny*self._dy

   @property
   def Nx(self) : return self._Nx

   @property
   def Ny(self) : return self._Ny

   def plotgrid(self,fac=1.) :
      return plotgrid(*self.pgrid(),width=self.width*fac,height=self.height*fac)

# --- ------------------------------------------------------------------
# --- Calc scuy and scvx from qlat, qlon
# ---               *--*--*
# ---               |  |  |
# --- stencil:      u--p--*
# ---               |  |  |
# --- We have q   ->q--v--*
# --- ------------------------------------------------------------------

   def scuy(self) :
      qlon,qlat = self.qgrid(extended=True)
      return actual_grid_spacing(qlon[1:-1,2:],qlat[1:-1,2:], qlon[1:-1,1:-1], qlat[1:-1,1:-1])

   def scvx(self) :
      qlon,qlat = self.qgrid(extended=True)
      return actual_grid_spacing(qlon[2:,1:-1],qlat[2:,1:-1], qlon[1:-1,1:-1],qlat[1:-1,1:-1])

# --- ------------------------------------------------------------------
# --- Calc scvy and scux from plat, plon
# ---               *--*--*
# ---               |  |  |
# --- stencil:      u--p--*
# ---               |  |  |
# --- We have p   ->q--v--*
# --- ------------------------------------------------------------------

   def scvy(self) :
      plon,plat = self.pgrid(extended=True)
      return actual_grid_spacing(plon[1:-1,1:-1],plat[1:-1,1:-1], plon[1:-1,:-2], plat[1:-1,:-2])

   def scux(self) :
      plon,plat = self.pgrid(extended=True)
      return actual_grid_spacing(plon[1:-1,1:-1],plat[1:-1,1:-1], plon[:-2,1:-1],plat[:-2,1:-1])

# --- ------------------------------------------------------------------
# --- Calc scpx and scqy from ulat, ulon
# ---               *--*--*
# ---               |  |  |
# --- stencil:      u--p--*
# ---               |  |  |
# --- We have u   ->q--v--*
# --- ------------------------------------------------------------------

   def scpx(self) :
      ulon,ulat = self.ugrid(extended=True)
      return actual_grid_spacing(ulon[1:-1,1:-1],ulat[1:-1,1:-1], ulon[2:,1:-1], ulat[2:,1:-1])

   def scqy(self) :
      ulon,ulat = self.ugrid(extended=True)
      return actual_grid_spacing(ulon[1:-1,1:-1],ulat[1:-1,1:-1], ulon[1:-1,:-2],ulat[1:-1,:-2])

# --- ------------------------------------------------------------------
# --- Calc scpy and scqx from vlat,vlon
# ---               *--*--*
# ---               |  |  |
# --- stencil:      u--p--*
# ---               |  |  |
# --- We have v   ->q--v--*
# --- ------------------------------------------------------------------

   def scpy(self) :
      vlon,vlat = self.vgrid(extended=True)
      return actual_grid_spacing(vlon[1:-1,1:-1],vlat[1:-1,1:-1], vlon[1:-1,2:], vlat[1:-1,2:])

   def scqx(self) :
      vlon,vlat = self.vgrid(extended=True)
      return actual_grid_spacing(vlon[1:-1,1:-1],vlat[1:-1,1:-1], vlon[1:-1,:-2],vlat[1:-1,:-2])

   def p_azimuth(self) :
      plon,plat = self.pgrid(extended=True)
      return fwd_azimuth( plon[1:-1,1:-1],plat[1:-1,1:-1], 
                          plon[2:,1:-1]   ,plat[2:,1:-1])

      
   def corio(self) :
      qlon,qlat = self.qgrid()
      return numpy.radians(qlat) * 4. * numpy.pi / 86164.0 # Sidereal day
      return 


   def aspect_ratio(self) :
      scpx=self.scpx()
      scpy=self.scpy()
      asp = numpy.where(scpy==0.,99.0,scpx/scpy)
      return asp


   def write_my_projection_info(self) :
      self.write_projection_info(self._proj4string,self._ll_lon,self._ll_lat,self._dx,self._dy,
         self._Nx,self._Ny)

   @classmethod
   def write_projection_info(cls,p4s,ll_lon,ll_lat,dx,dy,nx,ny) :
      fid=open("proj.info","w")
      fid.write("File Version         : %.1f\n"%1.0)
      fid.write("Proj4 String         : %s\n"%p4s)
      fid.write("Lower Left Longitude : %16.10f\n"%ll_lon)
      fid.write("Lower Left Latitude  : %16.10f\n"%ll_lat)
      fid.write("Projection Delta X   : %16.10g\n"%dx)
      fid.write("Projection Delta Y   : %16.10g\n"%dy)
      fid.write("Projection Nodes X   : %6i\n"%nx)
      fid.write("Projection Nodes Y   : %6i\n"%ny)
      fid.close()
      
   @classmethod
   def grid_from_file(cls,filename="proj.info") :
      p4s,ll_lon,ll_lat,dx,dy,nx,ny = cls.read_projection_info(filename)
      return cls(p4s,ll_lon,ll_lat,dx,dy,nx,ny)

   
   @classmethod
   def read_projection_info(cls,filename="proj.info") :
      fid=open(filename,"r")

      res = True

      line=fid.readline().strip()
      m = re.search("File Version[ ]*:(.*)",line)
      if m :
         #print m.group(1)
         v = float(m.group(1))
      else :
         msg = "Failed to read version number from projection file %s"%filename
         logger.error(msg)
         raise ValueError,msg


      line=fid.readline()
      m = re.search("Proj4 String[ ]*:(.*)",line)
      if m :
         p4s = m.group(1)
      else :
         msg = "Failed to read proj 4 string from projection file %s"%filename
         logger.error(msg)
         raise ValueError,msg

      line=fid.readline()
      m = re.search("Lower Left Longitude[ ]*:(.*)",line)
      if m :
         ll_lon = float(m.group(1))
      else :
         msg = "Failed to read lower left longitude from projection file %s"%filename
         logger.error(msg)
         raise ValueError,msg

      line=fid.readline()
      m = re.search("Lower Left Latitude[ ]*:(.*)",line)
      if m :
         ll_lat = float(m.group(1))
      else :
         msg = "Failed to read lower left latitude number from projection file %s"%filename
         logger.error(msg)
         raise ValueError,msg

      line=fid.readline()
      m = re.search("Projection Delta X[ ]*:(.*)",line)
      if m :
         dx = float(m.group(1))
      else :
         msg = "Failed to read delta x from projection file %s"%filename
         logger.error(msg)
         raise ValueError,msg

      line=fid.readline()
      m = re.search("Projection Delta Y[ ]*:(.*)",line)
      if m :
         dy = float(m.group(1))
      else :
         msg = "Failed to read delta y from projection file %s"%filename
         logger.error(msg)
         raise ValueError,msg

      line=fid.readline()
      m = re.search("Projection Nodes X[ ]*:(.*)",line)
      if m :
         nx = int(m.group(1))
      else :
         msg = "Failed to read number of x nodes from projection file %s"%filename
         logger.error(msg)
         raise ValueError,msg

      line=fid.readline()
      m = re.search("Projection Nodes Y[ ]*:(.*)",line)
      if m :
         ny = int(m.group(1))
      else :
         msg = "Failed to read number of y nodes from projection file %s"%filename
         logger.error(msg)
         raise ValueError,msg
      fid.close()

      return  p4s,ll_lon,ll_lat,dx,dy,nx,ny



   def raster(self) :
      pass



def actual_grid_spacing(lon1,lat1,lon2,lat2) : 
   geod=pyproj.Geod(ellps="sphere")
   az1,az2,dist = geod.inv(lon1,lat1,lon2,lat2)
   return dist

def fwd_azimuth(lon1,lat1,lon2,lat2) : 
   geod=pyproj.Geod(ellps="sphere")
   az1,az2,dist = geod.inv(lon1,lat1,lon2,lat2)
   return az1

   


def plotgrid(lon,lat,width=3000000,height=3000000) :
   import matplotlib
   from matplotlib.figure import Figure
   from matplotlib.backends.backend_agg import FigureCanvasAgg
   from mpl_toolkits.basemap import Basemap

   #figure = Figure()
   #ax     = figure.add_subplot(111)
   #canvas = FigureCanvasAgg(figure)

   figure = matplotlib.pyplot.figure()
   ax=figure.add_subplot(111)

   #Pick center longitude
   ix,iy=[elem/2 for elem in lon.shape]
   clon=lon[ix,iy]
   clat=lat[ix,iy]

   # Probably a way of estimating the width here...
   m = Basemap(projection='stere',lon_0=clon,lat_0=clat,resolution='l',width=width,height=height,ax=ax)
   x,y = m(lon,lat)

   # Pick a suitable set of grid lines
   stepx,stepy=[elem/10 for elem in lon.shape]

   x=x[::stepx,::stepy]
   y=y[::stepx,::stepy]


   m.drawcoastlines()
   m.drawmapboundary() # draw a line around the map region
   m.drawparallels(numpy.arange(-90.,120.,30.),labels=[1,0,0,0]) # draw parallels
   m.drawmeridians(numpy.arange(0.,420.,60.),labels=[0,0,0,1]) # draw meridians
   v=numpy.zeros(x.shape)
   col=".8"
   cmap=matplotlib.colors.ListedColormap([col,col])
   #m.pcolormesh(x,y,v,ax=ax,edgecolor="k",cmap=cmap)
   m.pcolormesh(x,y,v,ax=ax,edgecolor="none",cmap=cmap)
   for j in range(y.shape[1]) :
      m.plot(x[:,j],y[:,j],color="b",lw=2)
   for i in range(y.shape[0]) :
      m.plot(x[i,:],y[i,:],color="r",lw=2)
   ax.set_title("Every %d in x(blue) and every %d in y(red) shown"%(stepy,stepx))

   return figure



