""" MNodule for doing IO on files used by hycom """
import numpy 
import struct
import sys
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



class ZAFileError(Exception) :
   pass

class ZAFile(object) :
   def __init__(self,idm,jdm,filename,action,mask=False,real4=True,endian="native") :
      self._idm = idm
      self._jdm = jdm 
      self._filename = filename 
      self._action = action 
      self._mask = mask 
      self._real4= real4
      self._endian= endian
      if self._action.lower() not in ["r","w"]  :
         raise ZAFileError("action argument must be either r(ead) or w(rite)")
      if self._endian.lower() not in ["little","big","native"]  :
         raise ZAFileError("action argument must be either native, little ort big")

      
      if self._endian.lower() == "native" :
         self._swap_endian = False
      else :
         self._swap_endian = str(self._endian) <> str(sys.byteorder)

      logging.debug("Endianness set to %s",self._endian)
      logging.debug("Byteswapping is   %s"%str(self._swap_endian))


      self._zaiost()
      self._zaiopf()



   def _zaiost(self) :
      # Size of output 2D Array
      self._n2drec= ((self._idm*self._jdm+4095)/4096)*4096

      # Init sequenctial record counter
      self._iarec = 0
      self._spval = 2**100.


   def _zaiopf(self) :
      # Open .a and .b file
      self._filea = open(self._filename+".a",self._action+"b")
      self._fileb = open(self._filename+".b",self._action)



   def zaiowr_a(self,h,mask) :
      w=numpy.ones(self._n2drec)*self._spval
      w[0:self._idm*self._jdm] = h.flatten('F') # Fortran order

      # Calc min and mask
      if self._mask :
         tmp = mask>0
         I=numpy.where(tmp)
         hmax=h(I).max()
         hmin=h(I).min()
         J=numpy.where(~tmp)
         w[J] = self._spval

      else :
         hmax=h.max()
         hmin=h.min()

      if self._real4 :
         w = w.astype(numpy.float32)
         struct_fmt="f"
      else :
         w = w.astype(numpy.float64)
         struct_fmt="d"



      # Write to file, handle endianness ?
      if self._swap_endian :
         w = w.byteswap()
      #print w[0]

      # This is okay for sequential access. For accessing specific records, use fseek
      #self._filea.write(w)
      #
      #binpack=struct.pack("=%d%s"%(w.size,struct_fmt),*w[:])
      #w.write(binpack)
      #
      w.tofile(self._filea)

      return hmin,hmax


   def close(self) :
      self._filea.close()
      self._fileb.close()

   @property
   def n2drec(self):
      return self._n2drec



class ZAFileRegionalGrid(ZAFile) :
   def __init__(self,idm,jdm,filename,action,mask=False,real4=True,endian="native",mapflg=-1) :
      super(ZAFileRegionalGrid,self).__init__(idm,jdm,filename,action,mask=mask,real4=real4,endian=endian)

      # Regional file .b header
      self._fileb.write("%5d    %s\n"%(self._idm,"'idm   '"))
      self._fileb.write("%5d    %s\n"%(self._idm,"'jdm   '"))
      self._fileb.write("%5d    %s\n"%(mapflg,"'mapflg'"))



   def zaiowr(self,field,mask,fieldname) :
      hmin,hmax = self.zaiowr_a(field,mask)
      self.zaiowr_b(fieldname,hmin,hmax)


   def zaiowr_b(self,fieldname,hmin,hmax) :
      self._fileb.write("%4s:  min,max =%16.5f%16.5f\n"%(fieldname,hmin,hmax))
      



def write_regional_grid(grid) :
   plon,plat=grid.pgrid()
   ulon,ulat=grid.ugrid()
   vlon,vlat=grid.vgrid()
   qlon,qlat=grid.qgrid()


   regf = ZAFileRegionalGrid(grid.Nx,grid.Ny,"regional.grid","w",endian="native")
   regf.zaiowr(plon,plon,"plon")
   regf.zaiowr(plat,plat,"plat")
   regf.zaiowr(qlon,qlon,"qlon")
   regf.zaiowr(qlat,qlat,"qlat")
   regf.zaiowr(ulon,ulon,"ulon")
   regf.zaiowr(ulat,ulat,"ulat")
   regf.zaiowr(vlon,vlon,"vlon")
   regf.zaiowr(vlat,vlat,"vlat")

   regf.zaiowr(grid.p_azimuth(),vlat,"pang")

   regf.zaiowr(grid.scpx(),plon,"scpx")
   regf.zaiowr(grid.scpy(),plon,"scpy")
   regf.zaiowr(grid.scqx(),plon,"scqx")
   regf.zaiowr(grid.scqy(),plon,"scqy")
   regf.zaiowr(grid.scux(),plon,"scux")
   regf.zaiowr(grid.scuy(),plon,"scuy")
   regf.zaiowr(grid.scvx(),plon,"scvx")
   regf.zaiowr(grid.scvy(),plon,"scvy")

   regf.zaiowr(grid.corio(),plon,"cori")
   regf.zaiowr(grid.aspect_ratio(),plon,"pasp")

   #print regf.n2drec * 19 * 4

   regf.close()



def write_newpos(grid,endian="native") :
   if endian.lower() == "native" :
      swap_endian = False
   else :
      swap_endian = str(endian) <> str(sys.byteorder)
   logging.debug("Endianness set to %s",endian)
   logging.debug("Byteswapping is   %s"%str(swap_endian))

   # Number of records
   plon,plat=grid.pgrid()
   nrec = numpy.array([plon.size+plat.size])
   nrec = nrec.astype("i4")

   # Write to file, handle endianness ?
   plon = plon.astype(numpy.float64)
   plat = plat.astype(numpy.float64)
   if swap_endian :
      plon = plon.byteswap()
      plat = plat.byteswap()
      nrec = nrec.byteswap()

   print nrec,nrec.byteswap()
   print map(hex,nrec)
   print map(hex,nrec.byteswap())
   print map(hex,nrec.byteswap().byteswap())


   fid=open("newpos.uf","wb")

   # Fortran sequential unformatted access requires number of elements as first 
   nrec.tofile(fid)
   plat.tofile(fid)
   plon.tofile(fid)
   nrec.tofile(fid)
   #fid.write(sdruct.pack("i4",nrec))
   fid.close()


   
