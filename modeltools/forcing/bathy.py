import scipy.io.netcdf
import numpy
import logging
import re 



# Set up logger
_loglevel=logging.DEBUG
_logger = logging.getLogger(__name__)
_logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
_logger.addHandler(ch)


default_threshold=-5.


class GEBCO2014(object) :

   def __init__(self,filename="/work/shared/nersc/msc/ModelInput/bathymetry/GEBCO_2014/GEBCO_2014_2D.nc") :
      self._filename = filename
      self._nc = scipy.io.netcdf.netcdf_file(self._filename)

      self._lon = self._nc.variables["lon"][:]
      self._lat = self._nc.variables["lat"][:]
      self._elevation = self._nc.variables["elevation"][:]


      self._dlon= ( self._lon[-1] - self._lon[0] ) / (self._lon.size-1)
      self._dlat= ( self._lat[-1] - self._lat[0] ) / (self._lat.size-1)


   def regrid(self,lon,lat,width=None,myfilter="median"):

      # Pivot points
      wi= numpy.mod(lon-self._lon[0],360.) / self._dlon
      wj= (lat-self._lat[0]) / self._dlat

      i = numpy.floor(wi).astype("int")
      j = numpy.floor(wj).astype("int")
      ip1 = numpy.mod(i+1,self._lon.size)
      jp1 = numpy.minimum(j+1,self._lat.size-1)

      #print i
      #print wi
      wi = wi - i
      wj = wj - j



      if width is not None :

         numpts = int(numpy.round(width / (self._dlat*111*1000 )))
         if numpts < 3 :
            logger.debug("Too few points to do filtering. Using raw bathymetry")
            w = (1.-wi)*(1.-wj)*self._elevation[j  ,i  ] + \
                (   wi)*(1.-wj)*self._elevation[j  ,ip1] + \
                (   wi)*(   wj)*self._elevation[jp1,ip1] + \
                (1.-wi)*(   wj)*self._elevation[jp1,i  ]

         else :
            numpts = numpts/2 # -> half width
            logger.debug("Filtering width half-width = %d grid cells"%numpts)
            wtmp = numpy.zeros( ( (2*numpts)**2,lon.shape[0],lon.shape[1] ) )
            cnt=0
            for k in range(-numpts,numpts) :
               for l in range(-numpts,numpts) :
                  tmpi = i + k
                  tmpi = numpy.mod(tmpi,self._lon.shape[0])
                  tmpj = j + l
                  I = numpy.where(tmpj < 0)
                  if I[0].size>0 : 
                     tmpj[I] = numpy.abs(tmpj[I])
                     tmpi[I] = numpy.mod(tmpi[I]+ self._lon.size/2,self._lon.shape[0])
                  I = numpy.where(tmpj >= self._lat.size)
                  if I[0].size>0 : 
                     tmpj[I] = self._lat.size - 1 -  (tmpj[I] - (self._lat.size - 1))
                     tmpi[I] = numpy.mod(tmpi[I]+self._lon.size/2,self._lon.shape[0])
                  #print cnt,wtmp.shape
                  #print k,l,numpts,cnt,wtmp.shape
                  wtmp[cnt,:] = self._elevation[tmpj,tmpi]
                  cnt+=1
            if myfilter == "median" :
               w = numpy.median(wtmp,axis=0)
            elif myfilter == "average" :
               w = numpy.mean(wtmp,axis=0)
            else :
               w = numpy.mean(wtmp,axis=0)


      else :
         logger.debug("Using raw bathymetry")
         w = (1.-wi)*(1.-wj)*self._elevation[j  ,i  ] +\
             (   wi)*(1.-wj)*self._elevation[j  ,ip1] +\
             (   wi)*(   wj)*self._elevation[jp1,ip1] +\
             (1.-wi)*(   wj)*self._elevation[jp1,i  ]


      return w


def remove_one_neighbour_cells(inv,threshold=default_threshold) :
   v = numpy.copy(inv) 
   v[-1,:] = threshold
   v[1,:] = threshold
   v[:,-1] = threshold
   v[:,1] = threshold
   I=[[-1]]
   while len(I[0])>0 :
      v0   = v[1:-1,1:-1] < threshold
      vim1 = numpy.where(v[:-2,1:-1]< threshold,1,0)
      vip1 = numpy.where(v[2:,1:-1] < threshold,1,0)
      vjm1 = numpy.where(v[1:-1,:-2]< threshold,1,0)
      vjp1 = numpy.where(v[1:-1,2:] < threshold,1,0)
      tmp = vim1 + vip1 + vjm1 + vjp1
      I = numpy.where(numpy.logical_and(tmp<=1, v0  ))
      print "Found %d one neighbour cells"%I[0].size
      if I[0].size > 0 :
         #print v[I]
         #print I[0].size,tmp.min(),tmp.max()
         #print v0[I][0], vim1[I][0], vip1[I][0], vjm1[I][0], vjp1[I][0]
         v[1:-1,1:-1][I] = threshold
   return v

def remove_islets(inv,threshold=default_threshold) :
   v = numpy.copy(inv) 
   v[-1,:] = threshold
   v[1,:] = threshold
   v[:,-1] = threshold
   v[:,1] = threshold
   I=[[-1]]
   while len(I[0])>0 :
      v0   = v[1:-1,1:-1] >= threshold

      fim1 = numpy.where(v[:-2,1:-1]< threshold,1,0)
      fip1 = numpy.where(v[2:,1:-1] < threshold,1,0)
      fjm1 = numpy.where(v[1:-1,:-2]< threshold,1,0)
      fjp1 = numpy.where(v[1:-1,2:] < threshold,1,0)


      tmp = fim1 + fip1 + fjm1 + fjp1
      I = numpy.where(numpy.logical_and(tmp>=3, v0  ))
      print "Found %d islets"%I[0].size


      if I[0].size > 0 :
         vim1 = v[:-2,1:-1][I] * fim1[I]
         vip1 = v[2:,1:-1] [I] * fip1[I]
         vjm1 = v[1:-1,:-2][I] * fjm1[I]
         vjp1 = v[1:-1,2:] [I] * fjp1[I]
         v[1:-1,1:-1][I] = (vim1 + vip1 + vjm1 + vjp1)/tmp[I]
   return v


def remove_isolated_basins(lon,lat,inv,lon0,lat0,threshold=default_threshold) :
    import scipy.ndimage.measurements

    mask = inv < threshold
    labarray,num_features=scipy.ndimage.measurements.label(mask)
    feat_count = {}
    for i in range(num_features):
       feat_count[i+1] = numpy.count_nonzero(labarray==i+1)

    # Order by feature count
    tmp = sorted(feat_count.items(), key=lambda x:x[1],reverse=True)
    #for i in range(num_features) :
    #   print "Feature %03d: %d cells" %tuple(tmp[i])
    main_feature = tmp[0][0] 
    print "Main feature in terms of cells is feature %d"%main_feature

    # Find points nearest to lon0, lat0
    outv=numpy.ones(inv.shape)*threshold
    if len(lon0)> 0 :
       for lo,la in zip(lon0,lat0) :
         dist = numpy.sqrt((lo-lon)**2+(lat-la)**2) # close enough for this purpose
         I = numpy.argmin(dist)
         if I and inv.flatten()[I] < threshold :
            i,j=numpy.unravel_index(I,inv.shape)
            feature=labarray[i,j]
            print "Position (%7.3f,%7.3f) : Feature %d is used"%(lo,la,feature)
            outv = numpy.where(labarray==feature,inv,outv)
    else :
       print "No control points given, using main feature"
       outv = numpy.where(labarray==main_feature,inv,outv)
    return outv


def shapiro_filter(w,threshold=default_threshold,S=0.25) :
   myw = numpy.copy(w)
   # TODO: Account for periodic grids

   # x - direction
   wm1=numpy.copy(myw[1:-1,1:-3])
   w0 =numpy.copy(myw[1:-1,2:-2])
   wp1=numpy.copy(myw[1:-1,3:-1])
   fwm1=numpy.where(wm1<threshold,1,0)
   fw  =numpy.where(w0 <threshold,1,0)
   fwp1=numpy.where(wp1<threshold,1,0)
   I = numpy.where(fwm1+fw+fwp1==3)
   if I[0].size > 0 :
      w0[I] = w0[I] + 0.5*S * (wm1[I]+wp1[I]-2*w0[I])
   myw[1:-1,2:-2] = w0


   # x - direction
   wm1=numpy.copy(myw[1:-3,1:-1])
   w0 =numpy.copy(myw[2:-2,1:-1])
   wp1=numpy.copy(myw[3:-1,1:-1])
   fwm1=numpy.where(wm1<threshold,1,0)
   fw  =numpy.where(w0 <threshold,1,0)
   fwp1=numpy.where(wp1<threshold,1,0)
   I = numpy.where(fwm1+fw+fwp1==3)
   if I[0].size > 0 :
      w0[I] = w0[I] + 0.5*S * (wm1[I]+wp1[I]-2*w0[I])
   myw[2:-2,1:-1] = w0

   return myw




   



    
   
