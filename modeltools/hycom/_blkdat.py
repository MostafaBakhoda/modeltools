""" MNodule for parsing hycom blkdat """
import numpy 
import sys
import logging
import re

class BlkdatError(Exception):
    """Base class for exceptions in this module."""
    pass

# Set up logger
_loglevel=logging.DEBUG
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch) 



class BlkdatParser(object) :
   """ Class for doing binary input/output on hycom .a files """
   _integer_fields = [ "iversn","iexpt","idm","jdm","itest","jtest","kdm",
      "nhybrd","nsigma","kapref","thflag","vsigma","iniflg","jerlv0","yrflag",
      "sshflg","cplifq","incflg","incstp","incupf","hybmap","hybflg","advflg",
      "advtyp","momtyp","slip","ishelf","ntracr","trcflg","tsofrq","mlflag","pensol",
      "dypflg","bblkpp","shinst","dbdiff","nonloc","bodiw","difout","difsmo","hblflg",
      "niter","langmr","clmflg","fltflg","nfladv","wndflg","ustflg","flxflg","empflg",
      "dswflg","albflg","sssflg","sstflg","lwflag","icmflg",",slprs","stroff","flxoff",
      "flxsmo","relax","trcrlx","priver","epmass"]
   _float_fields = []
      
   def __init__(self,filename) :

      self._sigma    = []
      self._datadict = {}
      self._descdict = {}
      fid=open(filename,"r")

      # 4 first line is header
      self.header=[]
      self.header.append(fid.readline())
      self.header.append(fid.readline())
      self.header.append(fid.readline())
      self.header.append(fid.readline())


      for line in fid :

         # general pattern to search for
         #m = re.match("^(.*)'([a-z0-9]{6})'[ ]*=[ ]*(.*$)",line.strip())
         m = re.match("^[\s]*(.*)[\s]*'([a-z0-9 _]{6})'[\s]*=[\s]*(.*$)",line.strip())
         if not m :
            print m,line
            raise BlkdatError,"Error when parsing file. Line does not match required pattern"

         value=m.group(1).strip()
         key  =m.group(2).strip()
         desc =m.group(3).strip()

         #print key
         if key == "sigma" :
            self._sigma.append(float(value))
         else  :
            self._datadict[key] = value

         if key in  self._integer_fields :
            self._datadict[key] = int(self._datadict[key])
         elif key in  self._float_fields :
            self._datadict[key] = float(self._datadict[key])

   def __getitem__(self,i) :
      if i=="sigma" :
         return self._sigma
      elif i in self._integer_fields :
         return int(self._datadict[i])
      else :
         return self._datadict[i]











