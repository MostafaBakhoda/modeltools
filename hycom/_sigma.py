class Sigma(object) :

   DTHIRD=1.0/3.0


   def __init__(self,flag) :
# --- coefficients for sigma-0 (based on Brydon & Sun fit)
      self._sigma=flag
      if flag == 0 :
         self.C1=-1.36471E-01
         self.C2= 4.68181E-02
         self.C3= 8.07004E-01
         self.C4=-7.45353E-03
         self.C5=-2.94418E-03
         self.C6= 3.43570E-05
         self.C7= 3.48658E-05
      elif flag == 2 :
         self.C1=-1.36471E-01
         self.C2= 4.68181E-02
         self.C3= 8.07004E-01
         self.C4=-7.45353E-03
         self.C5=-2.94418E-03,
         self.C6= 3.43570E-05
         self.C7= 3.48658E-05
      else :
         raise ValueError,"flag<>0 not implemented"

   def A0(self,S) :
      return (self.C1+self.C3*S)/self.C6

   def A1(self,S) :
      return (self.C2+self.C5*S)/self.C6

   def A2(self,S) :
      return (self.C4+self.C7*S)/self.C6

   def CUBQ(self,S) :
      return self.DTHIRD*A1(S)-(self.DTHIRD*A2(S))**2

   def CUBR(self,R,S) :
      return self.DTHIRD*(0.50*A1(S)*A2(S)-1.50*(A0(S)-R/self.C6)) -(self.DTHIRD*A2(S))**3

   def CUBAN(self,R,S) :
      return self.DTHIRD*ATAN2(SQRT(MAX(DZERO,-(self.CUBQ(S)**3+CUBR(R,S)**2))),CUBR(R,S))

   def CUBRL(self,R,S) :
      return SQRT(-self.CUBQ(S))*COS(self.CUBAN(R,S))

   def CUBIM(self,R,S) :
      return SQRT(-self.CUBQ(S))*SIN(self.CUBAN(R,S))

# --- temp (deg c) as a function of sigma and salinity (mil)
   def TOFSIG(self,R,S) :
      return -self.CUBRL(R,S)+SQRT(3.)*self.CUBIM(R,S)-self.DTHIRD*self.A2(S)

# --- salinity (mil) as a function of sigma and temperature (deg c)
   def SOFSIG(self,R,T) :
      return (R-self.C1-T*(self.C2+T*(self.C4+self.C6*T)))/(self.C3+T*(self.C5+self.C7*T))

# --- sigma-theta as a function of temp (deg c) and salinity (mil)
# --- (friedrich-levitus 3rd degree polynomial fit)
   def SIG(self,T,S) :
      return (self.C1+self.C3*S+T*(self.C2+self.C5*S+T*(self.C4+self.C7*S+self.C6*T)))


# --- auxiliary statements for finding root of 3rd degree polynomial
#     A0(S)=(C1+C3*S)/C6
#     A1(S)=(C2+C5*S)/C6
#     A2(S)=(C4+C7*S)/C6
#     CUBQ(S)=DTHIRD*A1(S)-(DTHIRD*A2(S))**2
#     CUBR(R,S)=DTHIRD*(0.5D0*A1(S)*A2(S)-1.5D0*(A0(S)-R/C6))
#    .   -(DTHIRD*A2(S))**3
# --- if q**3+r**2>0, water is too dense to yield real root at given
# --- salinitiy. setting q**3+r**2=0 in that case is equivalent to
# --- lowering sigma until a double real root is obtained.
#     CUBAN(R,S)=DTHIRD*ATAN2(SQRT(MAX(DZERO,
#    .   -(CUBQ(S)**3+CUBR(R,S)**2))),CUBR(R,S))
#     CUBRL(R,S)=SQRT(-CUBQ(S))*COS(CUBAN(R,S))
#     CUBIM(R,S)=SQRT(-CUBQ(S))*SIN(CUBAN(R,S))
#
# --- temp (deg c) as a function of sigma and salinity (mil)
#     TOFSIG(R,S)=-CUBRL(R,S)+SQRT(3.)*CUBIM(R,S)-DTHIRD*A2(S)
#
# --- salinity (mil) as a function of sigma and temperature (deg c)
#     SOFSIG(R,T)=(R-C1-T*(C2+T*(C4+C6*T)))/(C3+T*(C5+C7*T))
#
# --- sigma-theta as a function of temp (deg c) and salinity (mil)
# --- (friedrich-levitus 3rd degree polynomial fit)
#     SIG(T,S)=(C1+C3*S+T*(C2+C5*S+T*(C4+C7*S+C6*T)))

   @property
   def sigma(self) : return self._sigma
