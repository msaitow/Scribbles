// Copyright (C) 2005-2014 by Takeshi Yanai (yanait@gmail.com)
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA

#include <tensor/tensor.h>
#include <sci/lct/lct_entry.h>
#include <sci/lct/lct_pno.h>
#include <sci/lct/lct_fmatrix.h>
#include <sci/lct/lct_linalg.h>
namespace orz { namespace lct {
    void FMatrixOrthogonalizer::FormLHS_Vm1_Vm2_CEd(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &c2, const orz::DTensor &d3, DensityPack &dPack,
                                                    TensorContainer &ovl, TensorContainer &PNO4, const bool do4EXTinPNO,
                                                    const double QE0, const double Ecas, const double h0_energy, TensorContainer &DebugInts, const orz::DTensor &casfock, const orz::DTensor &FV,
                                                    TensorContainer &T2, TensorContainer &R2){

#include "lct_preparations_orthdens.cpp"      
      
      // Back-transform PNO ampitude into canonical basis
      orz::DTensor T_m2_(_Nvirt,_Nvirt,Nrhom2);
      for(auto &&P : irange(0L, Nrhom2)){
        const long I_P = TActList.GetIndex(P);
        if (I_P < 0) continue;
        orz::DTensor U2 = T2.GetTensor(P);
        TActList.GetPair(P).BackTransform2EXT(U2);
        for(auto &&a : irange(0L, _Nvirt))
          for(auto &&b : irange(0L, _Nvirt))
            T_m2_(a,b,P) = U2(a,b);
      }

#include "lct_integrals4EXT.cpp"

 //Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn
  orz::DTensor tRFiC(_Nvirt,_Nvirt,Nrhom1,_Ndocc);
  {
    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Ndocc,_Nact,_Nact,_Nact);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(p*_Nact*_Nact+a1*_Nact+a0)*tmp1.shape(1)+(a2)] = d2.cptr()[p*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a2
          }//a0
        }//a1
      }//p
      orz::DTensor tmp2(_Nact,_Ndocc);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          tmp2.cptr()[(a2)*tmp2.shape(1)+(i)] = CFip.cptr()[i*CFip.shape(1)+a2];
        }//i
      }//a2
      // @W0(i,p,a1,a0) <<= +1 @D2(p,a1,a2,a0) @Fcore(i,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&i : irange(0L, _Ndocc)){
              W0.cptr()[i*_Nact*_Nact*_Nact+p*_Nact*_Nact+a1*_Nact+a0] += tmp3.cptr()[p*_Nact*_Nact*_Ndocc+a1*_Nact*_Ndocc+a0*_Ndocc+i];
            }//i
          }//a0
        }//a1
      }//p
    }
    orz::DTensor W1(_Ndocc,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Ndocc*_Nact,_Nact*_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&p : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nact+p)*tmp1.shape(1)+(a1*_Nact+a0)] = W0.cptr()[i*W0.shape(1)*W0.shape(2)*W0.shape(3)+p*W0.shape(2)*W0.shape(3)+a1*W0.shape(3)+a0];
            }//a0
          }//a1
        }//p
      }//i
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a1*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//_P_
        }//a0
      }//a1
      // @W1(i,p,_P_) <<= +1 @W0(i,p,a1,a0) @X(-2)(_P_,a1,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&p : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            W1.cptr()[i*_Nact*Nrhom2+p*Nrhom2+_P_] += tmp3.cptr()[i*_Nact*Nrhom2+p*Nrhom2+_P_];
          }//_P_
        }//p
      }//i
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(a*_Nvirt+b)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[a*T_m2_.shape(1)*T_m2_.shape(2)+b*T_m2_.shape(2)+_P_];
          }//_P_
        }//b
      }//a
      orz::DTensor tmp2(Nrhom2,_Ndocc*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Nact+p)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)+p*W1.shape(2)+_P_];
          }//p
        }//i
      }//_P_
      // @W2(i,p,a,b) <<= +1 @T(-2)(a,b,_P_) @W1(i,p,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&p : irange(0L, _Nact)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Nvirt*_Ndocc*_Nact+b*_Ndocc*_Nact+i*_Nact+p];
            }//p
          }//i
        }//b
      }//a
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+a*_Nvirt+b)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+b];
            }//p
          }//b
        }//a
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,a,b) @X(-1)(+)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Ndocc,_Nact,_Nact,_Nact);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(p*_Nact*_Nact+a1*_Nact+a0)*tmp1.shape(1)+(a2)] = d2.cptr()[p*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a2
          }//a0
        }//a1
      }//p
      orz::DTensor tmp2(_Nact,_Ndocc);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          tmp2.cptr()[(a2)*tmp2.shape(1)+(i)] = CFip.cptr()[i*CFip.shape(1)+a2];
        }//i
      }//a2
      // @W0(i,p,a1,a0) <<= +1 @D2(p,a1,a2,a0) @Fcore(i,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&i : irange(0L, _Ndocc)){
              W0.cptr()[i*_Nact*_Nact*_Nact+p*_Nact*_Nact+a1*_Nact+a0] += tmp3.cptr()[p*_Nact*_Nact*_Ndocc+a1*_Nact*_Ndocc+a0*_Ndocc+i];
            }//i
          }//a0
        }//a1
      }//p
    }
    orz::DTensor W1(_Ndocc,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Ndocc*_Nact,_Nact*_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&p : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nact+p)*tmp1.shape(1)+(a1*_Nact+a0)] = W0.cptr()[i*W0.shape(1)*W0.shape(2)*W0.shape(3)+p*W0.shape(2)*W0.shape(3)+a1*W0.shape(3)+a0];
            }//a0
          }//a1
        }//p
      }//i
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a1*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//_P_
        }//a0
      }//a1
      // @W1(i,p,_P_) <<= +1 @W0(i,p,a1,a0) @X(-2)(_P_,a1,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&p : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            W1.cptr()[i*_Nact*Nrhom2+p*Nrhom2+_P_] += tmp3.cptr()[i*_Nact*Nrhom2+p*Nrhom2+_P_];
          }//_P_
        }//p
      }//i
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(b*_Nvirt+a)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[b*T_m2_.shape(1)*T_m2_.shape(2)+a*T_m2_.shape(2)+_P_];
          }//_P_
        }//a
      }//b
      orz::DTensor tmp2(Nrhom2,_Ndocc*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Nact+p)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)+p*W1.shape(2)+_P_];
          }//p
        }//i
      }//_P_
      // @W2(i,p,b,a) <<= +1 @T(-2)(b,a,_P_) @W1(i,p,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&p : irange(0L, _Nact)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[b*_Nvirt*_Ndocc*_Nact+a*_Ndocc*_Nact+i*_Nact+p];
            }//p
          }//i
        }//a
      }//b
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+a];
            }//p
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,b,a) @X(-1)(-)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//b
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Ndocc,_Nact,_Nact,_Nact);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(p*_Nact*_Nact+a1*_Nact+a0)*tmp1.shape(1)+(a2)] = d2.cptr()[p*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a2
          }//a0
        }//a1
      }//p
      orz::DTensor tmp2(_Nact,_Ndocc);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          tmp2.cptr()[(a2)*tmp2.shape(1)+(i)] = CFip.cptr()[i*CFip.shape(1)+a2];
        }//i
      }//a2
      // @W0(i,p,a1,a0) <<= +1 @D2(p,a1,a2,a0) @Fcore(i,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&i : irange(0L, _Ndocc)){
              W0.cptr()[i*_Nact*_Nact*_Nact+p*_Nact*_Nact+a1*_Nact+a0] += tmp3.cptr()[p*_Nact*_Nact*_Ndocc+a1*_Nact*_Ndocc+a0*_Ndocc+i];
            }//i
          }//a0
        }//a1
      }//p
    }
    orz::DTensor W1(_Ndocc,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Ndocc*_Nact,_Nact*_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&p : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nact+p)*tmp1.shape(1)+(a1*_Nact+a0)] = W0.cptr()[i*W0.shape(1)*W0.shape(2)*W0.shape(3)+p*W0.shape(2)*W0.shape(3)+a1*W0.shape(3)+a0];
            }//a0
          }//a1
        }//p
      }//i
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+_P_];
          }//_P_
        }//a0
      }//a1
      // @W1(i,p,_P_) <<= +1 @W0(i,p,a1,a0) @X(-2)(_P_,a0,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&p : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            W1.cptr()[i*_Nact*Nrhom2+p*Nrhom2+_P_] += tmp3.cptr()[i*_Nact*Nrhom2+p*Nrhom2+_P_];
          }//_P_
        }//p
      }//i
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(b*_Nvirt+a)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[b*T_m2_.shape(1)*T_m2_.shape(2)+a*T_m2_.shape(2)+_P_];
          }//_P_
        }//a
      }//b
      orz::DTensor tmp2(Nrhom2,_Ndocc*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Nact+p)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)+p*W1.shape(2)+_P_];
          }//p
        }//i
      }//_P_
      // @W2(i,p,b,a) <<= +1 @T(-2)(b,a,_P_) @W1(i,p,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&p : irange(0L, _Nact)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[b*_Nvirt*_Ndocc*_Nact+a*_Ndocc*_Nact+i*_Nact+p];
            }//p
          }//i
        }//a
      }//b
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+a];
            }//p
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,b,a) @X(-1)(+)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//b
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Ndocc,_Nact,_Nact,_Nact);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(p*_Nact*_Nact+a1*_Nact+a0)*tmp1.shape(1)+(a2)] = d2.cptr()[p*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a2
          }//a0
        }//a1
      }//p
      orz::DTensor tmp2(_Nact,_Ndocc);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          tmp2.cptr()[(a2)*tmp2.shape(1)+(i)] = CFip.cptr()[i*CFip.shape(1)+a2];
        }//i
      }//a2
      // @W0(i,p,a1,a0) <<= +1 @D2(p,a1,a2,a0) @Fcore(i,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&i : irange(0L, _Ndocc)){
              W0.cptr()[i*_Nact*_Nact*_Nact+p*_Nact*_Nact+a1*_Nact+a0] += tmp3.cptr()[p*_Nact*_Nact*_Ndocc+a1*_Nact*_Ndocc+a0*_Ndocc+i];
            }//i
          }//a0
        }//a1
      }//p
    }
    orz::DTensor W1(_Ndocc,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Ndocc*_Nact,_Nact*_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&p : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nact+p)*tmp1.shape(1)+(a1*_Nact+a0)] = W0.cptr()[i*W0.shape(1)*W0.shape(2)*W0.shape(3)+p*W0.shape(2)*W0.shape(3)+a1*W0.shape(3)+a0];
            }//a0
          }//a1
        }//p
      }//i
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+_P_];
          }//_P_
        }//a0
      }//a1
      // @W1(i,p,_P_) <<= +1 @W0(i,p,a1,a0) @X(-2)(_P_,a0,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&p : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            W1.cptr()[i*_Nact*Nrhom2+p*Nrhom2+_P_] += tmp3.cptr()[i*_Nact*Nrhom2+p*Nrhom2+_P_];
          }//_P_
        }//p
      }//i
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(a*_Nvirt+b)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[a*T_m2_.shape(1)*T_m2_.shape(2)+b*T_m2_.shape(2)+_P_];
          }//_P_
        }//b
      }//a
      orz::DTensor tmp2(Nrhom2,_Ndocc*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Nact+p)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)+p*W1.shape(2)+_P_];
          }//p
        }//i
      }//_P_
      // @W2(i,p,a,b) <<= +1 @T(-2)(a,b,_P_) @W1(i,p,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&p : irange(0L, _Nact)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Nvirt*_Ndocc*_Nact+b*_Ndocc*_Nact+i*_Nact+p];
            }//p
          }//i
        }//b
      }//a
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+a*_Nvirt+b)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+b];
            }//p
          }//b
        }//a
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,a,b) @X(-1)(-)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Ndocc,_Nact,_Nact,_Nact);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact*_Nact*_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a4 : irange(0L, _Nact)){
              for(auto &&a3 : irange(0L, _Nact)){
                for(auto &&a0 : irange(0L, _Nact)){
                  tmp1.cptr()[(p*_Nact*_Nact+a2*_Nact+a1)*tmp1.shape(1)+(a4*_Nact*_Nact+a3*_Nact+a0)] = d3.cptr()[p*d3.shape(1)*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a2*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a4*d3.shape(3)*d3.shape(4)*d3.shape(5)+a1*d3.shape(4)*d3.shape(5)+a3*d3.shape(5)+a0];
                }//a0
              }//a3
            }//a4
          }//a1
        }//a2
      }//p
      orz::DTensor tmp2(_Nact*_Nact*_Nact,_Ndocc);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&i : irange(0L, _Ndocc)){
              tmp2.cptr()[(a4*_Nact*_Nact+a3*_Nact+a0)*tmp2.shape(1)+(i)] = V2_FiC_aaac.cptr()[a0*V2_FiC_aaac.shape(1)*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a3*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a4*V2_FiC_aaac.shape(3)+i];
            }//i
          }//a0
        }//a3
      }//a4
      // @W0(i,p,a2,a1) <<= +1 @D3(p,a2,a4,a1,a3,a0) @V2_FiC_aaac(a0,a3,a4,i)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&i : irange(0L, _Ndocc)){
              W0.cptr()[i*_Nact*_Nact*_Nact+p*_Nact*_Nact+a2*_Nact+a1] += tmp3.cptr()[p*_Nact*_Nact*_Ndocc+a2*_Nact*_Ndocc+a1*_Ndocc+i];
            }//i
          }//a1
        }//a2
      }//p
    }
    orz::DTensor W1(_Ndocc,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Ndocc*_Nact,_Nact*_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&p : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nact+p)*tmp1.shape(1)+(a2*_Nact+a1)] = W0.cptr()[i*W0.shape(1)*W0.shape(2)*W0.shape(3)+p*W0.shape(2)*W0.shape(3)+a2*W0.shape(3)+a1];
            }//a1
          }//a2
        }//p
      }//i
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a1)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a2*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+_P_];
          }//_P_
        }//a1
      }//a2
      // @W1(i,p,_P_) <<= +1 @W0(i,p,a2,a1) @X(-2)(_P_,a2,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&p : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            W1.cptr()[i*_Nact*Nrhom2+p*Nrhom2+_P_] += tmp3.cptr()[i*_Nact*Nrhom2+p*Nrhom2+_P_];
          }//_P_
        }//p
      }//i
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(a*_Nvirt+b)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[a*T_m2_.shape(1)*T_m2_.shape(2)+b*T_m2_.shape(2)+_P_];
          }//_P_
        }//b
      }//a
      orz::DTensor tmp2(Nrhom2,_Ndocc*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Nact+p)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)+p*W1.shape(2)+_P_];
          }//p
        }//i
      }//_P_
      // @W2(i,p,a,b) <<= +1 @T(-2)(a,b,_P_) @W1(i,p,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&p : irange(0L, _Nact)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Nvirt*_Ndocc*_Nact+b*_Ndocc*_Nact+i*_Nact+p];
            }//p
          }//i
        }//b
      }//a
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+a*_Nvirt+b)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+b];
            }//p
          }//b
        }//a
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,a,b) @X(-1)(+)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Ndocc,_Nact,_Nact,_Nact);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact*_Nact*_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a4 : irange(0L, _Nact)){
              for(auto &&a3 : irange(0L, _Nact)){
                for(auto &&a0 : irange(0L, _Nact)){
                  tmp1.cptr()[(p*_Nact*_Nact+a2*_Nact+a1)*tmp1.shape(1)+(a4*_Nact*_Nact+a3*_Nact+a0)] = d3.cptr()[p*d3.shape(1)*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a2*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a4*d3.shape(3)*d3.shape(4)*d3.shape(5)+a1*d3.shape(4)*d3.shape(5)+a3*d3.shape(5)+a0];
                }//a0
              }//a3
            }//a4
          }//a1
        }//a2
      }//p
      orz::DTensor tmp2(_Nact*_Nact*_Nact,_Ndocc);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&i : irange(0L, _Ndocc)){
              tmp2.cptr()[(a4*_Nact*_Nact+a3*_Nact+a0)*tmp2.shape(1)+(i)] = V2_FiC_aaac.cptr()[a0*V2_FiC_aaac.shape(1)*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a3*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a4*V2_FiC_aaac.shape(3)+i];
            }//i
          }//a0
        }//a3
      }//a4
      // @W0(i,p,a2,a1) <<= +1 @D3(p,a2,a4,a1,a3,a0) @V2_FiC_aaac(a0,a3,a4,i)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&i : irange(0L, _Ndocc)){
              W0.cptr()[i*_Nact*_Nact*_Nact+p*_Nact*_Nact+a2*_Nact+a1] += tmp3.cptr()[p*_Nact*_Nact*_Ndocc+a2*_Nact*_Ndocc+a1*_Ndocc+i];
            }//i
          }//a1
        }//a2
      }//p
    }
    orz::DTensor W1(_Ndocc,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Ndocc*_Nact,_Nact*_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&p : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nact+p)*tmp1.shape(1)+(a2*_Nact+a1)] = W0.cptr()[i*W0.shape(1)*W0.shape(2)*W0.shape(3)+p*W0.shape(2)*W0.shape(3)+a2*W0.shape(3)+a1];
            }//a1
          }//a2
        }//p
      }//i
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a1)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a2*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+_P_];
          }//_P_
        }//a1
      }//a2
      // @W1(i,p,_P_) <<= +1 @W0(i,p,a2,a1) @X(-2)(_P_,a2,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&p : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            W1.cptr()[i*_Nact*Nrhom2+p*Nrhom2+_P_] += tmp3.cptr()[i*_Nact*Nrhom2+p*Nrhom2+_P_];
          }//_P_
        }//p
      }//i
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(b*_Nvirt+a)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[b*T_m2_.shape(1)*T_m2_.shape(2)+a*T_m2_.shape(2)+_P_];
          }//_P_
        }//a
      }//b
      orz::DTensor tmp2(Nrhom2,_Ndocc*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Nact+p)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)+p*W1.shape(2)+_P_];
          }//p
        }//i
      }//_P_
      // @W2(i,p,b,a) <<= +1 @T(-2)(b,a,_P_) @W1(i,p,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&p : irange(0L, _Nact)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[b*_Nvirt*_Ndocc*_Nact+a*_Ndocc*_Nact+i*_Nact+p];
            }//p
          }//i
        }//a
      }//b
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+a];
            }//p
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,b,a) @X(-1)(-)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//b
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Ndocc,_Nact,_Nact,_Nact);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact*_Nact*_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a4 : irange(0L, _Nact)){
              for(auto &&a3 : irange(0L, _Nact)){
                for(auto &&a0 : irange(0L, _Nact)){
                  tmp1.cptr()[(p*_Nact*_Nact+a2*_Nact+a1)*tmp1.shape(1)+(a4*_Nact*_Nact+a3*_Nact+a0)] = d3.cptr()[p*d3.shape(1)*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a2*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a4*d3.shape(3)*d3.shape(4)*d3.shape(5)+a1*d3.shape(4)*d3.shape(5)+a3*d3.shape(5)+a0];
                }//a0
              }//a3
            }//a4
          }//a1
        }//a2
      }//p
      orz::DTensor tmp2(_Nact*_Nact*_Nact,_Ndocc);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&i : irange(0L, _Ndocc)){
              tmp2.cptr()[(a4*_Nact*_Nact+a3*_Nact+a0)*tmp2.shape(1)+(i)] = V2_FiC_aaac.cptr()[a0*V2_FiC_aaac.shape(1)*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a3*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a4*V2_FiC_aaac.shape(3)+i];
            }//i
          }//a0
        }//a3
      }//a4
      // @W0(i,p,a2,a1) <<= +1 @D3(p,a2,a4,a1,a3,a0) @V2_FiC_aaac(a0,a3,a4,i)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&i : irange(0L, _Ndocc)){
              W0.cptr()[i*_Nact*_Nact*_Nact+p*_Nact*_Nact+a2*_Nact+a1] += tmp3.cptr()[p*_Nact*_Nact*_Ndocc+a2*_Nact*_Ndocc+a1*_Ndocc+i];
            }//i
          }//a1
        }//a2
      }//p
    }
    orz::DTensor W1(_Ndocc,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Ndocc*_Nact,_Nact*_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&p : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nact+p)*tmp1.shape(1)+(a2*_Nact+a1)] = W0.cptr()[i*W0.shape(1)*W0.shape(2)*W0.shape(3)+p*W0.shape(2)*W0.shape(3)+a2*W0.shape(3)+a1];
            }//a1
          }//a2
        }//p
      }//i
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a1)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a1*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+_P_];
          }//_P_
        }//a1
      }//a2
      // @W1(i,p,_P_) <<= +1 @W0(i,p,a2,a1) @X(-2)(_P_,a1,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&p : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            W1.cptr()[i*_Nact*Nrhom2+p*Nrhom2+_P_] += tmp3.cptr()[i*_Nact*Nrhom2+p*Nrhom2+_P_];
          }//_P_
        }//p
      }//i
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(b*_Nvirt+a)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[b*T_m2_.shape(1)*T_m2_.shape(2)+a*T_m2_.shape(2)+_P_];
          }//_P_
        }//a
      }//b
      orz::DTensor tmp2(Nrhom2,_Ndocc*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Nact+p)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)+p*W1.shape(2)+_P_];
          }//p
        }//i
      }//_P_
      // @W2(i,p,b,a) <<= +1 @T(-2)(b,a,_P_) @W1(i,p,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&p : irange(0L, _Nact)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[b*_Nvirt*_Ndocc*_Nact+a*_Ndocc*_Nact+i*_Nact+p];
            }//p
          }//i
        }//a
      }//b
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+a];
            }//p
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,b,a) @X(-1)(+)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//b
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Ndocc,_Nact,_Nact,_Nact);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact*_Nact*_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a4 : irange(0L, _Nact)){
              for(auto &&a3 : irange(0L, _Nact)){
                for(auto &&a0 : irange(0L, _Nact)){
                  tmp1.cptr()[(p*_Nact*_Nact+a2*_Nact+a1)*tmp1.shape(1)+(a4*_Nact*_Nact+a3*_Nact+a0)] = d3.cptr()[p*d3.shape(1)*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a2*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a4*d3.shape(3)*d3.shape(4)*d3.shape(5)+a1*d3.shape(4)*d3.shape(5)+a3*d3.shape(5)+a0];
                }//a0
              }//a3
            }//a4
          }//a1
        }//a2
      }//p
      orz::DTensor tmp2(_Nact*_Nact*_Nact,_Ndocc);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&i : irange(0L, _Ndocc)){
              tmp2.cptr()[(a4*_Nact*_Nact+a3*_Nact+a0)*tmp2.shape(1)+(i)] = V2_FiC_aaac.cptr()[a0*V2_FiC_aaac.shape(1)*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a3*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a4*V2_FiC_aaac.shape(3)+i];
            }//i
          }//a0
        }//a3
      }//a4
      // @W0(i,p,a2,a1) <<= +1 @D3(p,a2,a4,a1,a3,a0) @V2_FiC_aaac(a0,a3,a4,i)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&i : irange(0L, _Ndocc)){
              W0.cptr()[i*_Nact*_Nact*_Nact+p*_Nact*_Nact+a2*_Nact+a1] += tmp3.cptr()[p*_Nact*_Nact*_Ndocc+a2*_Nact*_Ndocc+a1*_Ndocc+i];
            }//i
          }//a1
        }//a2
      }//p
    }
    orz::DTensor W1(_Ndocc,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Ndocc*_Nact,_Nact*_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&p : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nact+p)*tmp1.shape(1)+(a2*_Nact+a1)] = W0.cptr()[i*W0.shape(1)*W0.shape(2)*W0.shape(3)+p*W0.shape(2)*W0.shape(3)+a2*W0.shape(3)+a1];
            }//a1
          }//a2
        }//p
      }//i
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a1)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a1*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+_P_];
          }//_P_
        }//a1
      }//a2
      // @W1(i,p,_P_) <<= +1 @W0(i,p,a2,a1) @X(-2)(_P_,a1,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&p : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            W1.cptr()[i*_Nact*Nrhom2+p*Nrhom2+_P_] += tmp3.cptr()[i*_Nact*Nrhom2+p*Nrhom2+_P_];
          }//_P_
        }//p
      }//i
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(a*_Nvirt+b)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[a*T_m2_.shape(1)*T_m2_.shape(2)+b*T_m2_.shape(2)+_P_];
          }//_P_
        }//b
      }//a
      orz::DTensor tmp2(Nrhom2,_Ndocc*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Nact+p)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)+p*W1.shape(2)+_P_];
          }//p
        }//i
      }//_P_
      // @W2(i,p,a,b) <<= +1 @T(-2)(a,b,_P_) @W1(i,p,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&p : irange(0L, _Nact)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Nvirt*_Ndocc*_Nact+b*_Ndocc*_Nact+i*_Nact+p];
            }//p
          }//i
        }//b
      }//a
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+a*_Nvirt+b)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+b];
            }//p
          }//b
        }//a
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,a,b) @X(-1)(-)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(v0*_Nvirt+b)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[v0*T_m2_.shape(1)*T_m2_.shape(2)+b*T_m2_.shape(2)+_P_];
          }//_P_
        }//b
      }//v0
      orz::DTensor tmp2(Nrhom2,_Nact*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(a0*_Nact+a1)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+_P_];
          }//a1
        }//a0
      }//_P_
      // @W0(a0,a1,v0,b) <<= +1 @T(-2)(v0,b,_P_) @X(-2)(_P_,a0,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              W0.cptr()[a0*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+v0*_Nvirt+b] += tmp3.cptr()[v0*_Nvirt*_Nact*_Nact+b*_Nact*_Nact+a0*_Nact+a1];
            }//a1
          }//a0
        }//b
      }//v0
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(p*_Nact+a2)*tmp1.shape(1)+(a1*_Nact+a0)] = d2.cptr()[p*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a2
      }//p
      orz::DTensor tmp2(_Nact*_Nact,_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(v0*_Nvirt+b)] = W0.cptr()[a0*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+v0*W0.shape(3)+b];
            }//b
          }//v0
        }//a0
      }//a1
      // @W1(p,a2,v0,b) <<= +1 @D2(p,a1,a2,a0) @W0(a0,a1,v0,b)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              W1.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+v0*_Nvirt+b] += tmp3.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+v0*_Nvirt+b];
            }//b
          }//v0
        }//a2
      }//p
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a2)] = V2_FiC_vavc.cptr()[v0*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a2*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a*V2_FiC_vavc.shape(3)+i];
            }//a2
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*_Nact,_Nact*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a2)*tmp2.shape(1)+(p*_Nvirt+b)] = W1.cptr()[p*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+v0*W1.shape(3)+b];
            }//b
          }//p
        }//a2
      }//v0
      // @W2(i,p,a,b) <<= +1 @V2_FiC_vavc(v0,a2,a,i) @W1(p,a2,v0,b)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&b : irange(0L, _Nvirt)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Ndocc*_Nact*_Nvirt+i*_Nact*_Nvirt+p*_Nvirt+b];
            }//b
          }//p
        }//i
      }//a
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+a*_Nvirt+b)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+b];
            }//p
          }//b
        }//a
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,a,b) @X(-1)(+)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(v0*_Nvirt+a)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[v0*T_m2_.shape(1)*T_m2_.shape(2)+a*T_m2_.shape(2)+_P_];
          }//_P_
        }//a
      }//v0
      orz::DTensor tmp2(Nrhom2,_Nact*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(a0*_Nact+a1)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+_P_];
          }//a1
        }//a0
      }//_P_
      // @W0(a0,a1,v0,a) <<= +1 @T(-2)(v0,a,_P_) @X(-2)(_P_,a0,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              W0.cptr()[a0*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+v0*_Nvirt+a] += tmp3.cptr()[v0*_Nvirt*_Nact*_Nact+a*_Nact*_Nact+a0*_Nact+a1];
            }//a1
          }//a0
        }//a
      }//v0
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(p*_Nact+a2)*tmp1.shape(1)+(a1*_Nact+a0)] = d2.cptr()[p*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a2
      }//p
      orz::DTensor tmp2(_Nact*_Nact,_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(v0*_Nvirt+a)] = W0.cptr()[a0*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+v0*W0.shape(3)+a];
            }//a
          }//v0
        }//a0
      }//a1
      // @W1(p,a2,v0,a) <<= +1 @D2(p,a1,a2,a0) @W0(a0,a1,v0,a)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              W1.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+v0*_Nvirt+a] += tmp3.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+v0*_Nvirt+a];
            }//a
          }//v0
        }//a2
      }//p
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a2)] = V2_FiC_vavc.cptr()[v0*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a2*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+b*V2_FiC_vavc.shape(3)+i];
            }//a2
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*_Nact,_Nact*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a2)*tmp2.shape(1)+(p*_Nvirt+a)] = W1.cptr()[p*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+v0*W1.shape(3)+a];
            }//a
          }//p
        }//a2
      }//v0
      // @W2(i,p,b,a) <<= +1 @V2_FiC_vavc(v0,a2,b,i) @W1(p,a2,v0,a)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&a : irange(0L, _Nvirt)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[b*_Ndocc*_Nact*_Nvirt+i*_Nact*_Nvirt+p*_Nvirt+a];
            }//a
          }//p
        }//i
      }//b
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+a];
            }//p
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,b,a) @X(-1)(-)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//b
      }//i
    }
  }//End Contra
  {
    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(v0*_Nvirt+a)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[v0*T_m2_.shape(1)*T_m2_.shape(2)+a*T_m2_.shape(2)+_P_];
          }//_P_
        }//a
      }//v0
      orz::DTensor tmp2(Nrhom2,_Nact*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(a0*_Nact+a1)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+_P_];
          }//a1
        }//a0
      }//_P_
      // @W0(a0,a1,v0,a) <<= +1 @T(-2)(v0,a,_P_) @X(-2)(_P_,a0,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              W0.cptr()[a0*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+v0*_Nvirt+a] += tmp3.cptr()[v0*_Nvirt*_Nact*_Nact+a*_Nact*_Nact+a0*_Nact+a1];
            }//a1
          }//a0
        }//a
      }//v0
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(p*_Nact+a2)*tmp1.shape(1)+(a1*_Nact+a0)] = d2.cptr()[p*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a2
      }//p
      orz::DTensor tmp2(_Nact*_Nact,_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(v0*_Nvirt+a)] = W0.cptr()[a0*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+v0*W0.shape(3)+a];
            }//a
          }//v0
        }//a0
      }//a1
      // @W1(p,a2,v0,a) <<= +1 @D2(p,a1,a2,a0) @W0(a0,a1,v0,a)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              W1.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+v0*_Nvirt+a] += tmp3.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+v0*_Nvirt+a];
            }//a
          }//v0
        }//a2
      }//p
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a2)] = V2_FiC_vavc.cptr()[v0*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a2*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+b*V2_FiC_vavc.shape(3)+i];
            }//a2
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*_Nact,_Nact*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a2)*tmp2.shape(1)+(p*_Nvirt+a)] = W1.cptr()[p*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+v0*W1.shape(3)+a];
            }//a
          }//p
        }//a2
      }//v0
      // @W2(i,p,b,a) <<= +1 @V2_FiC_vavc(v0,a2,b,i) @W1(p,a2,v0,a)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&a : irange(0L, _Nvirt)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[b*_Ndocc*_Nact*_Nvirt+i*_Nact*_Nvirt+p*_Nvirt+a];
            }//a
          }//p
        }//i
      }//b
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+a];
            }//p
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,b,a) @X(-1)(+)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//b
      }//i
    }
  }//End Contra
  {
    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(v0*_Nvirt+b)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[v0*T_m2_.shape(1)*T_m2_.shape(2)+b*T_m2_.shape(2)+_P_];
          }//_P_
        }//b
      }//v0
      orz::DTensor tmp2(Nrhom2,_Nact*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(a0*_Nact+a1)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+_P_];
          }//a1
        }//a0
      }//_P_
      // @W0(a0,a1,v0,b) <<= +1 @T(-2)(v0,b,_P_) @X(-2)(_P_,a0,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              W0.cptr()[a0*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+v0*_Nvirt+b] += tmp3.cptr()[v0*_Nvirt*_Nact*_Nact+b*_Nact*_Nact+a0*_Nact+a1];
            }//a1
          }//a0
        }//b
      }//v0
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(p*_Nact+a2)*tmp1.shape(1)+(a1*_Nact+a0)] = d2.cptr()[p*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a2
      }//p
      orz::DTensor tmp2(_Nact*_Nact,_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(v0*_Nvirt+b)] = W0.cptr()[a0*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+v0*W0.shape(3)+b];
            }//b
          }//v0
        }//a0
      }//a1
      // @W1(p,a2,v0,b) <<= +1 @D2(p,a1,a2,a0) @W0(a0,a1,v0,b)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              W1.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+v0*_Nvirt+b] += tmp3.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+v0*_Nvirt+b];
            }//b
          }//v0
        }//a2
      }//p
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a2)] = V2_FiC_vavc.cptr()[v0*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a2*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a*V2_FiC_vavc.shape(3)+i];
            }//a2
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*_Nact,_Nact*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a2)*tmp2.shape(1)+(p*_Nvirt+b)] = W1.cptr()[p*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+v0*W1.shape(3)+b];
            }//b
          }//p
        }//a2
      }//v0
      // @W2(i,p,a,b) <<= +1 @V2_FiC_vavc(v0,a2,a,i) @W1(p,a2,v0,b)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&b : irange(0L, _Nvirt)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Ndocc*_Nact*_Nvirt+i*_Nact*_Nvirt+p*_Nvirt+b];
            }//b
          }//p
        }//i
      }//a
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+a*_Nvirt+b)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+b];
            }//p
          }//b
        }//a
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,a,b) @X(-1)(-)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(v0*_Nvirt+a)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[v0*T_m2_.shape(1)*T_m2_.shape(2)+a*T_m2_.shape(2)+_P_];
          }//_P_
        }//a
      }//v0
      orz::DTensor tmp2(Nrhom2,_Nact*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(a0*_Nact+a1)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+_P_];
          }//a1
        }//a0
      }//_P_
      // @W0(a0,a1,v0,a) <<= +1 @T(-2)(v0,a,_P_) @X(-2)(_P_,a0,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              W0.cptr()[a0*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+v0*_Nvirt+a] += tmp3.cptr()[v0*_Nvirt*_Nact*_Nact+a*_Nact*_Nact+a0*_Nact+a1];
            }//a1
          }//a0
        }//a
      }//v0
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(p*_Nact+a2)*tmp1.shape(1)+(a1*_Nact+a0)] = d2.cptr()[p*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a2
      }//p
      orz::DTensor tmp2(_Nact*_Nact,_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(v0*_Nvirt+a)] = W0.cptr()[a0*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+v0*W0.shape(3)+a];
            }//a
          }//v0
        }//a0
      }//a1
      // @W1(p,a2,v0,a) <<= +1 @D2(p,a1,a2,a0) @W0(a0,a1,v0,a)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              W1.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+v0*_Nvirt+a] += tmp3.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+v0*_Nvirt+a];
            }//a
          }//v0
        }//a2
      }//p
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a2)] = V2_FiC_vvac.cptr()[b*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a2*V2_FiC_vvac.shape(3)+i];
            }//a2
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*_Nact,_Nact*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a2)*tmp2.shape(1)+(p*_Nvirt+a)] = W1.cptr()[p*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+v0*W1.shape(3)+a];
            }//a
          }//p
        }//a2
      }//v0
      // @W2(i,p,b,a) <<= +1 @V2_FiC_vvac(b,v0,a2,i) @W1(p,a2,v0,a)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&a : irange(0L, _Nvirt)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[b*_Ndocc*_Nact*_Nvirt+i*_Nact*_Nvirt+p*_Nvirt+a];
            }//a
          }//p
        }//i
      }//b
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+a];
            }//p
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,b,a) @X(-1)(+)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//b
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(v0*_Nvirt+b)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[v0*T_m2_.shape(1)*T_m2_.shape(2)+b*T_m2_.shape(2)+_P_];
          }//_P_
        }//b
      }//v0
      orz::DTensor tmp2(Nrhom2,_Nact*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(a0*_Nact+a1)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+_P_];
          }//a1
        }//a0
      }//_P_
      // @W0(a0,a1,v0,b) <<= +1 @T(-2)(v0,b,_P_) @X(-2)(_P_,a0,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              W0.cptr()[a0*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+v0*_Nvirt+b] += tmp3.cptr()[v0*_Nvirt*_Nact*_Nact+b*_Nact*_Nact+a0*_Nact+a1];
            }//a1
          }//a0
        }//b
      }//v0
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(p*_Nact+a2)*tmp1.shape(1)+(a1*_Nact+a0)] = d2.cptr()[p*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a2
      }//p
      orz::DTensor tmp2(_Nact*_Nact,_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(v0*_Nvirt+b)] = W0.cptr()[a0*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+v0*W0.shape(3)+b];
            }//b
          }//v0
        }//a0
      }//a1
      // @W1(p,a2,v0,b) <<= +1 @D2(p,a1,a2,a0) @W0(a0,a1,v0,b)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              W1.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+v0*_Nvirt+b] += tmp3.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+v0*_Nvirt+b];
            }//b
          }//v0
        }//a2
      }//p
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a2)] = V2_FiC_vvac.cptr()[a*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a2*V2_FiC_vvac.shape(3)+i];
            }//a2
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*_Nact,_Nact*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a2)*tmp2.shape(1)+(p*_Nvirt+b)] = W1.cptr()[p*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+v0*W1.shape(3)+b];
            }//b
          }//p
        }//a2
      }//v0
      // @W2(i,p,a,b) <<= +1 @V2_FiC_vvac(a,v0,a2,i) @W1(p,a2,v0,b)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&b : irange(0L, _Nvirt)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Ndocc*_Nact*_Nvirt+i*_Nact*_Nvirt+p*_Nvirt+b];
            }//b
          }//p
        }//i
      }//a
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+a*_Nvirt+b)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+b];
            }//p
          }//b
        }//a
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,a,b) @X(-1)(-)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(v0*_Nvirt+b)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[v0*T_m2_.shape(1)*T_m2_.shape(2)+b*T_m2_.shape(2)+_P_];
          }//_P_
        }//b
      }//v0
      orz::DTensor tmp2(Nrhom2,_Nact*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1*_Nact+a0)] = Cm2.cptr()[a1*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//a0
        }//a1
      }//_P_
      // @W0(a1,a0,v0,b) <<= +1 @T(-2)(v0,b,_P_) @X(-2)(_P_,a1,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              W0.cptr()[a1*_Nact*_Nvirt*_Nvirt+a0*_Nvirt*_Nvirt+v0*_Nvirt+b] += tmp3.cptr()[v0*_Nvirt*_Nact*_Nact+b*_Nact*_Nact+a1*_Nact+a0];
            }//a0
          }//a1
        }//b
      }//v0
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(p*_Nact+a2)*tmp1.shape(1)+(a1*_Nact+a0)] = d2.cptr()[p*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a2
      }//p
      orz::DTensor tmp2(_Nact*_Nact,_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(v0*_Nvirt+b)] = W0.cptr()[a1*W0.shape(1)*W0.shape(2)*W0.shape(3)+a0*W0.shape(2)*W0.shape(3)+v0*W0.shape(3)+b];
            }//b
          }//v0
        }//a0
      }//a1
      // @W1(p,a2,v0,b) <<= +1 @D2(p,a1,a2,a0) @W0(a1,a0,v0,b)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              W1.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+v0*_Nvirt+b] += tmp3.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+v0*_Nvirt+b];
            }//b
          }//v0
        }//a2
      }//p
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a2)] = V2_FiC_vvac.cptr()[a*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a2*V2_FiC_vvac.shape(3)+i];
            }//a2
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*_Nact,_Nact*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a2)*tmp2.shape(1)+(p*_Nvirt+b)] = W1.cptr()[p*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+v0*W1.shape(3)+b];
            }//b
          }//p
        }//a2
      }//v0
      // @W2(i,p,a,b) <<= +1 @V2_FiC_vvac(a,v0,a2,i) @W1(p,a2,v0,b)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&b : irange(0L, _Nvirt)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Ndocc*_Nact*_Nvirt+i*_Nact*_Nvirt+p*_Nvirt+b];
            }//b
          }//p
        }//i
      }//a
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+a*_Nvirt+b)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+b];
            }//p
          }//b
        }//a
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,a,b) @X(-1)(+)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(v0*_Nvirt+a)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[v0*T_m2_.shape(1)*T_m2_.shape(2)+a*T_m2_.shape(2)+_P_];
          }//_P_
        }//a
      }//v0
      orz::DTensor tmp2(Nrhom2,_Nact*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1*_Nact+a0)] = Cm2.cptr()[a1*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//a0
        }//a1
      }//_P_
      // @W0(a1,a0,v0,a) <<= +1 @T(-2)(v0,a,_P_) @X(-2)(_P_,a1,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              W0.cptr()[a1*_Nact*_Nvirt*_Nvirt+a0*_Nvirt*_Nvirt+v0*_Nvirt+a] += tmp3.cptr()[v0*_Nvirt*_Nact*_Nact+a*_Nact*_Nact+a1*_Nact+a0];
            }//a0
          }//a1
        }//a
      }//v0
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(p*_Nact+a2)*tmp1.shape(1)+(a1*_Nact+a0)] = d2.cptr()[p*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a2
      }//p
      orz::DTensor tmp2(_Nact*_Nact,_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(v0*_Nvirt+a)] = W0.cptr()[a1*W0.shape(1)*W0.shape(2)*W0.shape(3)+a0*W0.shape(2)*W0.shape(3)+v0*W0.shape(3)+a];
            }//a
          }//v0
        }//a0
      }//a1
      // @W1(p,a2,v0,a) <<= +1 @D2(p,a1,a2,a0) @W0(a1,a0,v0,a)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              W1.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+v0*_Nvirt+a] += tmp3.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+v0*_Nvirt+a];
            }//a
          }//v0
        }//a2
      }//p
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a2)] = V2_FiC_vvac.cptr()[b*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a2*V2_FiC_vvac.shape(3)+i];
            }//a2
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*_Nact,_Nact*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a2)*tmp2.shape(1)+(p*_Nvirt+a)] = W1.cptr()[p*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+v0*W1.shape(3)+a];
            }//a
          }//p
        }//a2
      }//v0
      // @W2(i,p,b,a) <<= +1 @V2_FiC_vvac(b,v0,a2,i) @W1(p,a2,v0,a)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&a : irange(0L, _Nvirt)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[b*_Ndocc*_Nact*_Nvirt+i*_Nact*_Nvirt+p*_Nvirt+a];
            }//a
          }//p
        }//i
      }//b
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+a];
            }//p
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,b,a) @X(-1)(-)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//b
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(b*_Nvirt+v0)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[b*T_m2_.shape(1)*T_m2_.shape(2)+v0*T_m2_.shape(2)+_P_];
          }//_P_
        }//v0
      }//b
      orz::DTensor tmp2(Nrhom2,_Nact*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1*_Nact+a0)] = Cm2.cptr()[a1*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//a0
        }//a1
      }//_P_
      // @W0(a1,a0,b,v0) <<= +1 @T(-2)(b,v0,_P_) @X(-2)(_P_,a1,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              W0.cptr()[a1*_Nact*_Nvirt*_Nvirt+a0*_Nvirt*_Nvirt+b*_Nvirt+v0] += tmp3.cptr()[b*_Nvirt*_Nact*_Nact+v0*_Nact*_Nact+a1*_Nact+a0];
            }//a0
          }//a1
        }//v0
      }//b
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(p*_Nact+a2)*tmp1.shape(1)+(a1*_Nact+a0)] = d2.cptr()[p*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a2
      }//p
      orz::DTensor tmp2(_Nact*_Nact,_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(b*_Nvirt+v0)] = W0.cptr()[a1*W0.shape(1)*W0.shape(2)*W0.shape(3)+a0*W0.shape(2)*W0.shape(3)+b*W0.shape(3)+v0];
            }//v0
          }//b
        }//a0
      }//a1
      // @W1(p,a2,b,v0) <<= +1 @D2(p,a1,a2,a0) @W0(a1,a0,b,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              W1.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+b*_Nvirt+v0] += tmp3.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+b*_Nvirt+v0];
            }//v0
          }//b
        }//a2
      }//p
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a2)] = V2_FiC_vavc.cptr()[v0*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a2*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a*V2_FiC_vavc.shape(3)+i];
            }//a2
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*_Nact,_Nact*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a2)*tmp2.shape(1)+(p*_Nvirt+b)] = W1.cptr()[p*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+b*W1.shape(3)+v0];
            }//b
          }//p
        }//a2
      }//v0
      // @W2(i,p,a,b) <<= +1 @V2_FiC_vavc(v0,a2,a,i) @W1(p,a2,b,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&b : irange(0L, _Nvirt)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Ndocc*_Nact*_Nvirt+i*_Nact*_Nvirt+p*_Nvirt+b];
            }//b
          }//p
        }//i
      }//a
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+a*_Nvirt+b)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+b];
            }//p
          }//b
        }//a
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,a,b) @X(-1)(+)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(a*_Nvirt+v0)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[a*T_m2_.shape(1)*T_m2_.shape(2)+v0*T_m2_.shape(2)+_P_];
          }//_P_
        }//v0
      }//a
      orz::DTensor tmp2(Nrhom2,_Nact*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1*_Nact+a0)] = Cm2.cptr()[a1*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//a0
        }//a1
      }//_P_
      // @W0(a1,a0,a,v0) <<= +1 @T(-2)(a,v0,_P_) @X(-2)(_P_,a1,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              W0.cptr()[a1*_Nact*_Nvirt*_Nvirt+a0*_Nvirt*_Nvirt+a*_Nvirt+v0] += tmp3.cptr()[a*_Nvirt*_Nact*_Nact+v0*_Nact*_Nact+a1*_Nact+a0];
            }//a0
          }//a1
        }//v0
      }//a
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(p*_Nact+a2)*tmp1.shape(1)+(a1*_Nact+a0)] = d2.cptr()[p*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a2
      }//p
      orz::DTensor tmp2(_Nact*_Nact,_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(a*_Nvirt+v0)] = W0.cptr()[a1*W0.shape(1)*W0.shape(2)*W0.shape(3)+a0*W0.shape(2)*W0.shape(3)+a*W0.shape(3)+v0];
            }//v0
          }//a
        }//a0
      }//a1
      // @W1(p,a2,a,v0) <<= +1 @D2(p,a1,a2,a0) @W0(a1,a0,a,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              W1.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+a*_Nvirt+v0] += tmp3.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+a*_Nvirt+v0];
            }//v0
          }//a
        }//a2
      }//p
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a2)] = V2_FiC_vavc.cptr()[v0*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a2*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+b*V2_FiC_vavc.shape(3)+i];
            }//a2
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*_Nact,_Nact*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a2)*tmp2.shape(1)+(p*_Nvirt+a)] = W1.cptr()[p*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+a*W1.shape(3)+v0];
            }//a
          }//p
        }//a2
      }//v0
      // @W2(i,p,b,a) <<= +1 @V2_FiC_vavc(v0,a2,b,i) @W1(p,a2,a,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&a : irange(0L, _Nvirt)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[b*_Ndocc*_Nact*_Nvirt+i*_Nact*_Nvirt+p*_Nvirt+a];
            }//a
          }//p
        }//i
      }//b
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+a];
            }//p
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,b,a) @X(-1)(-)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//b
      }//i
    }
  }//End Contra
  {
    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(a*_Nvirt+v0)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[a*T_m2_.shape(1)*T_m2_.shape(2)+v0*T_m2_.shape(2)+_P_];
          }//_P_
        }//v0
      }//a
      orz::DTensor tmp2(Nrhom2,_Nact*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1*_Nact+a0)] = Cm2.cptr()[a1*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//a0
        }//a1
      }//_P_
      // @W0(a1,a0,a,v0) <<= +1 @T(-2)(a,v0,_P_) @X(-2)(_P_,a1,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              W0.cptr()[a1*_Nact*_Nvirt*_Nvirt+a0*_Nvirt*_Nvirt+a*_Nvirt+v0] += tmp3.cptr()[a*_Nvirt*_Nact*_Nact+v0*_Nact*_Nact+a1*_Nact+a0];
            }//a0
          }//a1
        }//v0
      }//a
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(p*_Nact+a2)*tmp1.shape(1)+(a1*_Nact+a0)] = d2.cptr()[p*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a2
      }//p
      orz::DTensor tmp2(_Nact*_Nact,_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(a*_Nvirt+v0)] = W0.cptr()[a1*W0.shape(1)*W0.shape(2)*W0.shape(3)+a0*W0.shape(2)*W0.shape(3)+a*W0.shape(3)+v0];
            }//v0
          }//a
        }//a0
      }//a1
      // @W1(p,a2,a,v0) <<= +1 @D2(p,a1,a2,a0) @W0(a1,a0,a,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              W1.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+a*_Nvirt+v0] += tmp3.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+a*_Nvirt+v0];
            }//v0
          }//a
        }//a2
      }//p
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a2)] = V2_FiC_vavc.cptr()[v0*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a2*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+b*V2_FiC_vavc.shape(3)+i];
            }//a2
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*_Nact,_Nact*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a2)*tmp2.shape(1)+(p*_Nvirt+a)] = W1.cptr()[p*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+a*W1.shape(3)+v0];
            }//a
          }//p
        }//a2
      }//v0
      // @W2(i,p,b,a) <<= +1 @V2_FiC_vavc(v0,a2,b,i) @W1(p,a2,a,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&a : irange(0L, _Nvirt)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[b*_Ndocc*_Nact*_Nvirt+i*_Nact*_Nvirt+p*_Nvirt+a];
            }//a
          }//p
        }//i
      }//b
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+a];
            }//p
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,b,a) @X(-1)(+)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//b
      }//i
    }
  }//End Contra
  {
    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(b*_Nvirt+v0)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[b*T_m2_.shape(1)*T_m2_.shape(2)+v0*T_m2_.shape(2)+_P_];
          }//_P_
        }//v0
      }//b
      orz::DTensor tmp2(Nrhom2,_Nact*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1*_Nact+a0)] = Cm2.cptr()[a1*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//a0
        }//a1
      }//_P_
      // @W0(a1,a0,b,v0) <<= +1 @T(-2)(b,v0,_P_) @X(-2)(_P_,a1,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              W0.cptr()[a1*_Nact*_Nvirt*_Nvirt+a0*_Nvirt*_Nvirt+b*_Nvirt+v0] += tmp3.cptr()[b*_Nvirt*_Nact*_Nact+v0*_Nact*_Nact+a1*_Nact+a0];
            }//a0
          }//a1
        }//v0
      }//b
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(p*_Nact+a2)*tmp1.shape(1)+(a1*_Nact+a0)] = d2.cptr()[p*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a2
      }//p
      orz::DTensor tmp2(_Nact*_Nact,_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(b*_Nvirt+v0)] = W0.cptr()[a1*W0.shape(1)*W0.shape(2)*W0.shape(3)+a0*W0.shape(2)*W0.shape(3)+b*W0.shape(3)+v0];
            }//v0
          }//b
        }//a0
      }//a1
      // @W1(p,a2,b,v0) <<= +1 @D2(p,a1,a2,a0) @W0(a1,a0,b,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              W1.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+b*_Nvirt+v0] += tmp3.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+b*_Nvirt+v0];
            }//v0
          }//b
        }//a2
      }//p
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a2)] = V2_FiC_vavc.cptr()[v0*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a2*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a*V2_FiC_vavc.shape(3)+i];
            }//a2
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*_Nact,_Nact*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a2)*tmp2.shape(1)+(p*_Nvirt+b)] = W1.cptr()[p*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+b*W1.shape(3)+v0];
            }//b
          }//p
        }//a2
      }//v0
      // @W2(i,p,a,b) <<= +1 @V2_FiC_vavc(v0,a2,a,i) @W1(p,a2,b,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&b : irange(0L, _Nvirt)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Ndocc*_Nact*_Nvirt+i*_Nact*_Nvirt+p*_Nvirt+b];
            }//b
          }//p
        }//i
      }//a
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+a*_Nvirt+b)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+b];
            }//p
          }//b
        }//a
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,a,b) @X(-1)(-)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(a*_Nvirt+v0)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[a*T_m2_.shape(1)*T_m2_.shape(2)+v0*T_m2_.shape(2)+_P_];
          }//_P_
        }//v0
      }//a
      orz::DTensor tmp2(Nrhom2,_Nact*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1*_Nact+a0)] = Cm2.cptr()[a1*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//a0
        }//a1
      }//_P_
      // @W0(a1,a0,a,v0) <<= +1 @T(-2)(a,v0,_P_) @X(-2)(_P_,a1,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              W0.cptr()[a1*_Nact*_Nvirt*_Nvirt+a0*_Nvirt*_Nvirt+a*_Nvirt+v0] += tmp3.cptr()[a*_Nvirt*_Nact*_Nact+v0*_Nact*_Nact+a1*_Nact+a0];
            }//a0
          }//a1
        }//v0
      }//a
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(p*_Nact+a2)*tmp1.shape(1)+(a1*_Nact+a0)] = d2.cptr()[p*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a2
      }//p
      orz::DTensor tmp2(_Nact*_Nact,_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(a*_Nvirt+v0)] = W0.cptr()[a1*W0.shape(1)*W0.shape(2)*W0.shape(3)+a0*W0.shape(2)*W0.shape(3)+a*W0.shape(3)+v0];
            }//v0
          }//a
        }//a0
      }//a1
      // @W1(p,a2,a,v0) <<= +1 @D2(p,a1,a2,a0) @W0(a1,a0,a,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              W1.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+a*_Nvirt+v0] += tmp3.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+a*_Nvirt+v0];
            }//v0
          }//a
        }//a2
      }//p
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a2)] = V2_FiC_vvac.cptr()[b*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a2*V2_FiC_vvac.shape(3)+i];
            }//a2
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*_Nact,_Nact*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a2)*tmp2.shape(1)+(p*_Nvirt+a)] = W1.cptr()[p*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+a*W1.shape(3)+v0];
            }//a
          }//p
        }//a2
      }//v0
      // @W2(i,p,b,a) <<= +1 @V2_FiC_vvac(b,v0,a2,i) @W1(p,a2,a,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&a : irange(0L, _Nvirt)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[b*_Ndocc*_Nact*_Nvirt+i*_Nact*_Nvirt+p*_Nvirt+a];
            }//a
          }//p
        }//i
      }//b
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+a];
            }//p
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,b,a) @X(-1)(+)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//b
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(b*_Nvirt+v0)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[b*T_m2_.shape(1)*T_m2_.shape(2)+v0*T_m2_.shape(2)+_P_];
          }//_P_
        }//v0
      }//b
      orz::DTensor tmp2(Nrhom2,_Nact*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1*_Nact+a0)] = Cm2.cptr()[a1*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//a0
        }//a1
      }//_P_
      // @W0(a1,a0,b,v0) <<= +1 @T(-2)(b,v0,_P_) @X(-2)(_P_,a1,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              W0.cptr()[a1*_Nact*_Nvirt*_Nvirt+a0*_Nvirt*_Nvirt+b*_Nvirt+v0] += tmp3.cptr()[b*_Nvirt*_Nact*_Nact+v0*_Nact*_Nact+a1*_Nact+a0];
            }//a0
          }//a1
        }//v0
      }//b
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(p*_Nact+a2)*tmp1.shape(1)+(a1*_Nact+a0)] = d2.cptr()[p*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a2
      }//p
      orz::DTensor tmp2(_Nact*_Nact,_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(b*_Nvirt+v0)] = W0.cptr()[a1*W0.shape(1)*W0.shape(2)*W0.shape(3)+a0*W0.shape(2)*W0.shape(3)+b*W0.shape(3)+v0];
            }//v0
          }//b
        }//a0
      }//a1
      // @W1(p,a2,b,v0) <<= +1 @D2(p,a1,a2,a0) @W0(a1,a0,b,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              W1.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+b*_Nvirt+v0] += tmp3.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+b*_Nvirt+v0];
            }//v0
          }//b
        }//a2
      }//p
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a2)] = V2_FiC_vvac.cptr()[a*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a2*V2_FiC_vvac.shape(3)+i];
            }//a2
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*_Nact,_Nact*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a2)*tmp2.shape(1)+(p*_Nvirt+b)] = W1.cptr()[p*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+b*W1.shape(3)+v0];
            }//b
          }//p
        }//a2
      }//v0
      // @W2(i,p,a,b) <<= +1 @V2_FiC_vvac(a,v0,a2,i) @W1(p,a2,b,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&b : irange(0L, _Nvirt)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Ndocc*_Nact*_Nvirt+i*_Nact*_Nvirt+p*_Nvirt+b];
            }//b
          }//p
        }//i
      }//a
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+a*_Nvirt+b)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+b];
            }//p
          }//b
        }//a
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,a,b) @X(-1)(-)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(b*_Nvirt+v0)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[b*T_m2_.shape(1)*T_m2_.shape(2)+v0*T_m2_.shape(2)+_P_];
          }//_P_
        }//v0
      }//b
      orz::DTensor tmp2(Nrhom2,_Nact*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(a0*_Nact+a1)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+_P_];
          }//a1
        }//a0
      }//_P_
      // @W0(a0,a1,b,v0) <<= +1 @T(-2)(b,v0,_P_) @X(-2)(_P_,a0,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              W0.cptr()[a0*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+b*_Nvirt+v0] += tmp3.cptr()[b*_Nvirt*_Nact*_Nact+v0*_Nact*_Nact+a0*_Nact+a1];
            }//a1
          }//a0
        }//v0
      }//b
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(p*_Nact+a2)*tmp1.shape(1)+(a1*_Nact+a0)] = d2.cptr()[p*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a2
      }//p
      orz::DTensor tmp2(_Nact*_Nact,_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(b*_Nvirt+v0)] = W0.cptr()[a0*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+b*W0.shape(3)+v0];
            }//v0
          }//b
        }//a0
      }//a1
      // @W1(p,a2,b,v0) <<= +1 @D2(p,a1,a2,a0) @W0(a0,a1,b,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              W1.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+b*_Nvirt+v0] += tmp3.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+b*_Nvirt+v0];
            }//v0
          }//b
        }//a2
      }//p
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a2)] = V2_FiC_vvac.cptr()[a*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a2*V2_FiC_vvac.shape(3)+i];
            }//a2
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*_Nact,_Nact*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a2)*tmp2.shape(1)+(p*_Nvirt+b)] = W1.cptr()[p*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+b*W1.shape(3)+v0];
            }//b
          }//p
        }//a2
      }//v0
      // @W2(i,p,a,b) <<= +1 @V2_FiC_vvac(a,v0,a2,i) @W1(p,a2,b,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&b : irange(0L, _Nvirt)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Ndocc*_Nact*_Nvirt+i*_Nact*_Nvirt+p*_Nvirt+b];
            }//b
          }//p
        }//i
      }//a
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+a*_Nvirt+b)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+b];
            }//p
          }//b
        }//a
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,a,b) @X(-1)(+)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(a*_Nvirt+v0)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[a*T_m2_.shape(1)*T_m2_.shape(2)+v0*T_m2_.shape(2)+_P_];
          }//_P_
        }//v0
      }//a
      orz::DTensor tmp2(Nrhom2,_Nact*_Nact);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(a0*_Nact+a1)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+_P_];
          }//a1
        }//a0
      }//_P_
      // @W0(a0,a1,a,v0) <<= +1 @T(-2)(a,v0,_P_) @X(-2)(_P_,a0,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              W0.cptr()[a0*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+a*_Nvirt+v0] += tmp3.cptr()[a*_Nvirt*_Nact*_Nact+v0*_Nact*_Nact+a0*_Nact+a1];
            }//a1
          }//a0
        }//v0
      }//a
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(p*_Nact+a2)*tmp1.shape(1)+(a1*_Nact+a0)] = d2.cptr()[p*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a2
      }//p
      orz::DTensor tmp2(_Nact*_Nact,_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(a*_Nvirt+v0)] = W0.cptr()[a0*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+a*W0.shape(3)+v0];
            }//v0
          }//a
        }//a0
      }//a1
      // @W1(p,a2,a,v0) <<= +1 @D2(p,a1,a2,a0) @W0(a0,a1,a,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              W1.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+a*_Nvirt+v0] += tmp3.cptr()[p*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+a*_Nvirt+v0];
            }//v0
          }//a
        }//a2
      }//p
    }
    orz::DTensor W2(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a2)] = V2_FiC_vvac.cptr()[b*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a2*V2_FiC_vvac.shape(3)+i];
            }//a2
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*_Nact,_Nact*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a2)*tmp2.shape(1)+(p*_Nvirt+a)] = W1.cptr()[p*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+a*W1.shape(3)+v0];
            }//a
          }//p
        }//a2
      }//v0
      // @W2(i,p,b,a) <<= +1 @V2_FiC_vvac(b,v0,a2,i) @W1(p,a2,a,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&p : irange(0L, _Nact)){
            for(auto &&a : irange(0L, _Nvirt)){
              W2.cptr()[i*_Nact*_Nvirt*_Nvirt+p*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[b*_Ndocc*_Nact*_Nvirt+i*_Nact*_Nvirt+p*_Nvirt+a];
            }//a
          }//p
        }//i
      }//b
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&p : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(p)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+p*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+a];
            }//p
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&p : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(p)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+p];
        }//P
      }//p
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W2(i,p,b,a) @X(-1)(-)(P,p)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//b
      }//i
    }
  }//End Contra
  //genCode.insert(std::make_pair("tRFiC", boost::shared_ptr<orz::DTensor>( new orz::DTensor(tRFiC) )));
 //Code-End-Code-End-Code-End-Code-End-Code-End-Code-End-Code-End
      
      // Transform tRFiC into PNO basis
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          const long I_Pi = RActList.GetIndex(P,i);
          if (I_Pi < 0) continue;
          orz::DTensor RP(_Nvirt, _Nvirt);
          for(auto &&a : irange(0L, _Nvirt))
            for(auto &&b : irange(0L, _Nvirt))
              RP(a,b) = tRFiC(a,b,P,i);
          RActList.GetPair(P,i).Transform2EXT(RP);
          orz::DTensor Rorig = R2.GetTensor(P,i);
          RP += Rorig;          
          R2.PutTensor(P,i, RP);
        }
      }
                  
      return;
    }//End func
  } }
