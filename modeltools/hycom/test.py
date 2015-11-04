import unittest
import timetools
import random

class TimeTest(unittest.TestCase):

    def test_time_conversion(self):
       """ Tests that forday and dayfor are consistent """
       numtest=2000000
       for k in range(numtest) :
          yrflag = random.randrange(0,4)
          selector = random.randrange(0,3)
          dtime = random.randrange(1,200000)
          
          # Tests wnen close to integer, and normal
          if selector==0 :
             dtime = dtime + random.random()/100.
          elif selector==1 :
             dtime = dtime + 1.-random.random()/100.
          elif selector==2 :
             dtime = dtime + random.random()   

          if k%100000 == 0 : 
              print "Doing test %10d of %10d. yrflag=%d, dtime=%20.16g"%(k,numtest,yrflag,dtime)

          try :
              #print "dtime,yrflag=",dtime,yrflag
              iyear,iday,ihour = timetools.forday(dtime,yrflag)
              #print "iyear,iday,ihour=",iyear,iday,ihour
              dtimeout=timetools.dayfor(iyear,iday,ihour,yrflag)
              #print "dtime, reversed:",dtime
          except:
              unittest.fail("Caught error. dtime = %20.10g,yrflag=%d"%(dtime,yrflag))
          
          #print yrflag,dtime,dtimeout-dtime,3600./86400.
          
          #dtimeout is "floored" to nearest hour. So error should never be greater than this:
          if (dtime-dtimeout)>3600./86400:
              unittest.fail("Error: inn forday or dayfor, dtime=%20.10g, yrflag%d="%(dtime,yrflag))

             
if __name__ == "__main__" :
   unittest.main()
