import logging
import numpy

_default_threshold=-5.
_default_S=0.25

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


def shapiro_filter(w,threshold=_default_threshold,S=_default_S) :
   logger.info("Applying shapiro filter (1 pass, S=%.3f)"%S)
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




def remove_one_neighbour_cells(inv,threshold=_default_threshold) :
   v = numpy.copy(inv) 
   v[-1,:] = threshold
   v[0,:] = threshold
   v[:,-1] = threshold
   v[:,0] = threshold
   I=[[-1]]
   while len(I[0])>0 :
      v0   = v[1:-1,1:-1] < threshold
      vim1 = numpy.where(v[:-2,1:-1]< threshold,1,0)
      vip1 = numpy.where(v[2:,1:-1] < threshold,1,0)
      vjm1 = numpy.where(v[1:-1,:-2]< threshold,1,0)
      vjp1 = numpy.where(v[1:-1,2:] < threshold,1,0)
      tmp = vim1 + vip1 + vjm1 + vjp1
      I = numpy.where(numpy.logical_and(tmp<=1, v0  ))
      logger.info("Found %d one neighbour cells"%I[0].size)
      if I[0].size > 0 :
         #print v[I]
         #print I[0].size,tmp.min(),tmp.max()
         #print v0[I][0], vim1[I][0], vip1[I][0], vjm1[I][0], vjp1[I][0]
         v[1:-1,1:-1][I] = threshold
   return v

def remove_islets(inv,threshold=_default_threshold) :
   v = numpy.copy(inv) 
   v[-1,:] = threshold
   v[0,:] = threshold
   v[:,-1] = threshold
   v[:,0] = threshold
   I=[[-1]]
   while len(I[0])>0 :
      v0   = v[1:-1,1:-1] >= threshold

      fim1 = numpy.where(v[:-2,1:-1]< threshold,1,0)
      fip1 = numpy.where(v[2:,1:-1] < threshold,1,0)
      fjm1 = numpy.where(v[1:-1,:-2]< threshold,1,0)
      fjp1 = numpy.where(v[1:-1,2:] < threshold,1,0)


      tmp = fim1 + fip1 + fjm1 + fjp1
      I = numpy.where(numpy.logical_and(tmp>=3, v0  ))
      logger.info("Found %d islets"%I[0].size)


      if I[0].size > 0 :
         vim1 = v[:-2,1:-1][I] * fim1[I]
         vip1 = v[2:,1:-1] [I] * fip1[I]
         vjm1 = v[1:-1,:-2][I] * fjm1[I]
         vjp1 = v[1:-1,2:] [I] * fjp1[I]
         v[1:-1,1:-1][I] = (vim1 + vip1 + vjm1 + vjp1)/tmp[I]
   return v


def remove_isolated_basins(lon,lat,inv,lon0,lat0,threshold=_default_threshold) :
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
    logger.info("Main feature in terms of cells is feature %d"%main_feature)

    # Find points nearest to lon0, lat0
    outv=numpy.ones(inv.shape)*threshold
    if len(lon0)> 0 :
       for lo,la in zip(lon0,lat0) :
         dist = numpy.sqrt((lo-lon)**2+(lat-la)**2) # close enough for this purpose
         I = numpy.argmin(dist)
         if I and inv.flatten()[I] < threshold :
            i,j=numpy.unravel_index(I,inv.shape)
            feature=labarray[i,j]
            logger.info( "Position (%7.3f,%7.3f) : Feature %d is used"%(lo,la,feature))
            outv = numpy.where(labarray==feature,inv,outv)
    else :
       logger.info( "No control points given, using main feature")
       outv = numpy.where(labarray==main_feature,inv,outv)
    return outv



   



    
   
