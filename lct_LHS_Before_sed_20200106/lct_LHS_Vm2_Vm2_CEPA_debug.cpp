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
#include <chrono>

using namespace std::chrono;

namespace orz { namespace lct {
    void FMatrixOrthogonalizer::FormLHS_Vm2_Vm2_CEd(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &c2, const orz::DTensor &d3, DensityPack &dPack,
                                                    TensorContainer &ovl, TensorContainer &PNO4, const bool do4EXTinPNO,
                                                    const double QE0, const double Ecas, const double h0_energy, TensorContainer &DebugInts, const orz::DTensor &casfock, const orz::DTensor &FV,
                                                    TensorContainer &T2, TensorContainer &R2){

      orz::ProgressTimer pt("* V(-2)/V(-2)");
      
#include "lct_preparations_orthdens.cpp"      
      R2.initContainer(Nrhom2);                   
      for(auto &&P : irange(0L, Nrhom2)){
        const long I_P = TActList.GetIndex(P);
        if (I_P < 0) continue;        
        const long Nvirt = RActList.GetPair(0,P).GetNPNO();
        orz::DTensor R_P(Nvirt, Nvirt);                    
        R2.PutTensor(P, R_P);                              
      }                                                    
      
      // Back-transform PNO ampitude into canonical basis
      orz::DTensor T_m2_(_Nvirt,_Nvirt,Nrhom2);
      for(auto &&P : irange(0L, Nrhom2)){
        const long I_P = RActList.GetIndex(P);
        if (I_P < 0) continue;
        orz::DTensor U2 = T2.GetTensor(P);
        TActList.GetPair(P).BackTransform2EXT(U2);
        for(auto &&a : irange(0L, _Nvirt))
          for(auto &&b : irange(0L, _Nvirt))
            T_m2_(a,b,P) = U2(a,b);
      }
        
#include "lct_integrals4EXT.cpp"

#define _PROFILE_LHS
      
 //Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn
  orz::DTensor tRFiC(_Nvirt,_Nvirt,Nrhom2);
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= +1/4 @W1(P,_P_) @W2(b,a,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a2*_Nact+a0)*tmp1.shape(1)+(a3*_Nact+a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a3
        }//a0
      }//a2
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a3*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-2)(P,a3,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W0.cptr()[a2*_Nact*Nrhom2+a0*Nrhom2+P] += tmp3.cptr()[a2*_Nact*Nrhom2+a0*Nrhom2+P];
          }//P
        }//a0
      }//a2
    }
    orz::DTensor W1(Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(Nrhom2,_Nact*_Nact);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a2*W0.shape(1)*W0.shape(2)+a0*W0.shape(2)+P];
          }//a0
        }//a2
      }//P
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+_P_];
          }//_P_
        }//a0
      }//a2
      // @W1(P,_P_) <<= +1 @W0(a2,a0,P) @X(-2)(_P_,a0,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          W1.cptr()[P*Nrhom2+_P_] += tmp3.cptr()[P*Nrhom2+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt,_Nvirt);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          tmp1.cptr()[(b)*tmp1.shape(1)+(v0)] = CFab.cptr()[v0*CFab.shape(1)+b];
        }//v0
      }//b
      orz::DTensor tmp2(_Nvirt,_Nvirt*Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(v0)*tmp2.shape(1)+(a*Nrhom2+_P_)] = T_m2_.cptr()[v0*T_m2_.shape(1)*T_m2_.shape(2)+a*T_m2_.shape(2)+_P_];
          }//_P_
        }//a
      }//v0
      // @W2(b,a,_P_) <<= +1 @Fcore(v0,b) @T(-2)(v0,a,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            W2.cptr()[b*_Nvirt*Nrhom2+a*Nrhom2+_P_] += tmp3.cptr()[b*_Nvirt*Nrhom2+a*Nrhom2+_P_];
          }//_P_
        }//a
      }//b
    }
    {
      orz::DTensor tmp1(Nrhom2,Nrhom2);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(_P_)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//_P_
      }//P
      orz::DTensor tmp2(Nrhom2,_Nvirt*_Nvirt);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(b*_Nvirt+a)] = W2.cptr()[b*W2.shape(1)*W2.shape(2)+a*W2.shape(2)+_P_];
          }//a
        }//b
      }//_P_
      // @tRFiC(a,b,P) <<= +1 _X_ @W1(P,_P_) @W2(b,a,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[P*_Nvirt*_Nvirt+b*_Nvirt+a];
          }//a
        }//b
      }//P
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= +1/4 @W1(P,_P_) @W2(a,b,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a2*_Nact+a0)*tmp1.shape(1)+(a3*_Nact+a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a3
        }//a0
      }//a2
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a3*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-2)(P,a3,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W0.cptr()[a2*_Nact*Nrhom2+a0*Nrhom2+P] += tmp3.cptr()[a2*_Nact*Nrhom2+a0*Nrhom2+P];
          }//P
        }//a0
      }//a2
    }
    orz::DTensor W1(Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(Nrhom2,_Nact*_Nact);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a2*W0.shape(1)*W0.shape(2)+a0*W0.shape(2)+P];
          }//a0
        }//a2
      }//P
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a2*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//_P_
        }//a0
      }//a2
      // @W1(P,_P_) <<= +1 @W0(a2,a0,P) @X(-2)(_P_,a2,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          W1.cptr()[P*Nrhom2+_P_] += tmp3.cptr()[P*Nrhom2+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt,_Nvirt);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          tmp1.cptr()[(a)*tmp1.shape(1)+(v0)] = CFab.cptr()[v0*CFab.shape(1)+a];
        }//v0
      }//a
      orz::DTensor tmp2(_Nvirt,_Nvirt*Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(v0)*tmp2.shape(1)+(b*Nrhom2+_P_)] = T_m2_.cptr()[v0*T_m2_.shape(1)*T_m2_.shape(2)+b*T_m2_.shape(2)+_P_];
          }//_P_
        }//b
      }//v0
      // @W2(a,b,_P_) <<= +1 @Fcore(v0,a) @T(-2)(v0,b,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            W2.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+_P_] += tmp3.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+_P_];
          }//_P_
        }//b
      }//a
    }
    {
      orz::DTensor tmp1(Nrhom2,Nrhom2);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(_P_)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//_P_
      }//P
      orz::DTensor tmp2(Nrhom2,_Nvirt*_Nvirt);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(a*_Nvirt+b)] = W2.cptr()[a*W2.shape(1)*W2.shape(2)+b*W2.shape(2)+_P_];
          }//b
        }//a
      }//_P_
      // @tRFiC(a,b,P) <<= +1 _X_ @W1(P,_P_) @W2(a,b,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[P*_Nvirt*_Nvirt+a*_Nvirt+b];
          }//b
        }//a
      }//P
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= +1/4 @W1(P,_P_) @W2(b,a,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a2*_Nact+a0)*tmp1.shape(1)+(a3*_Nact+a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a3
        }//a0
      }//a2
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a3*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-2)(P,a3,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W0.cptr()[a2*_Nact*Nrhom2+a0*Nrhom2+P] += tmp3.cptr()[a2*_Nact*Nrhom2+a0*Nrhom2+P];
          }//P
        }//a0
      }//a2
    }
    orz::DTensor W1(Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(Nrhom2,_Nact*_Nact);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a2*W0.shape(1)*W0.shape(2)+a0*W0.shape(2)+P];
          }//a0
        }//a2
      }//P
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a2*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//_P_
        }//a0
      }//a2
      // @W1(P,_P_) <<= +1 @W0(a2,a0,P) @X(-2)(_P_,a2,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          W1.cptr()[P*Nrhom2+_P_] += tmp3.cptr()[P*Nrhom2+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt,_Nvirt);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          tmp1.cptr()[(b)*tmp1.shape(1)+(v0)] = CFab.cptr()[v0*CFab.shape(1)+b];
        }//v0
      }//b
      orz::DTensor tmp2(_Nvirt,_Nvirt*Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(v0)*tmp2.shape(1)+(a*Nrhom2+_P_)] = T_m2_.cptr()[a*T_m2_.shape(1)*T_m2_.shape(2)+v0*T_m2_.shape(2)+_P_];
          }//_P_
        }//a
      }//v0
      // @W2(b,a,_P_) <<= +1 @Fcore(v0,b) @T(-2)(a,v0,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            W2.cptr()[b*_Nvirt*Nrhom2+a*Nrhom2+_P_] += tmp3.cptr()[b*_Nvirt*Nrhom2+a*Nrhom2+_P_];
          }//_P_
        }//a
      }//b
    }
    {
      orz::DTensor tmp1(Nrhom2,Nrhom2);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(_P_)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//_P_
      }//P
      orz::DTensor tmp2(Nrhom2,_Nvirt*_Nvirt);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(b*_Nvirt+a)] = W2.cptr()[b*W2.shape(1)*W2.shape(2)+a*W2.shape(2)+_P_];
          }//a
        }//b
      }//_P_
      // @tRFiC(a,b,P) <<= +1 _X_ @W1(P,_P_) @W2(b,a,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[P*_Nvirt*_Nvirt+b*_Nvirt+a];
          }//a
        }//b
      }//P
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= +1/4 @W1(P,_P_) @W2(a,b,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a2*_Nact+a0)*tmp1.shape(1)+(a3*_Nact+a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a3
        }//a0
      }//a2
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a3*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-2)(P,a3,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W0.cptr()[a2*_Nact*Nrhom2+a0*Nrhom2+P] += tmp3.cptr()[a2*_Nact*Nrhom2+a0*Nrhom2+P];
          }//P
        }//a0
      }//a2
    }
    orz::DTensor W1(Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(Nrhom2,_Nact*_Nact);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a2*W0.shape(1)*W0.shape(2)+a0*W0.shape(2)+P];
          }//a0
        }//a2
      }//P
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+_P_];
          }//_P_
        }//a0
      }//a2
      // @W1(P,_P_) <<= +1 @W0(a2,a0,P) @X(-2)(_P_,a0,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          W1.cptr()[P*Nrhom2+_P_] += tmp3.cptr()[P*Nrhom2+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt,_Nvirt);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          tmp1.cptr()[(a)*tmp1.shape(1)+(v0)] = CFab.cptr()[v0*CFab.shape(1)+a];
        }//v0
      }//a
      orz::DTensor tmp2(_Nvirt,_Nvirt*Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(v0)*tmp2.shape(1)+(b*Nrhom2+_P_)] = T_m2_.cptr()[b*T_m2_.shape(1)*T_m2_.shape(2)+v0*T_m2_.shape(2)+_P_];
          }//_P_
        }//b
      }//v0
      // @W2(a,b,_P_) <<= +1 @Fcore(v0,a) @T(-2)(b,v0,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            W2.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+_P_] += tmp3.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+_P_];
          }//_P_
        }//b
      }//a
    }
    {
      orz::DTensor tmp1(Nrhom2,Nrhom2);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(_P_)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//_P_
      }//P
      orz::DTensor tmp2(Nrhom2,_Nvirt*_Nvirt);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(a*_Nvirt+b)] = W2.cptr()[a*W2.shape(1)*W2.shape(2)+b*W2.shape(2)+_P_];
          }//b
        }//a
      }//_P_
      // @tRFiC(a,b,P) <<= +1 _X_ @W1(P,_P_) @W2(a,b,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[P*_Nvirt*_Nvirt+a*_Nvirt+b];
          }//b
        }//a
      }//P
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/4 @T(-2)(b,a,_P_) @W2(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a4*_Nact+a2)*tmp1.shape(1)+(a3*_Nact+a1)] = d2.cptr()[a4*d2.shape(1)*d2.shape(2)*d2.shape(3)+a3*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a1];
            }//a1
          }//a3
        }//a2
      }//a4
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a1*Cm2.shape(1)*Cm2.shape(2)+a3*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a4,a2,P) <<= +1 @D2(a4,a3,a2,a1) @X(-2)(P,a1,a3)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W0.cptr()[a4*_Nact*Nrhom2+a2*Nrhom2+P] += tmp3.cptr()[a4*_Nact*Nrhom2+a2*Nrhom2+P];
          }//P
        }//a2
      }//a4
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a4)] = CFpq.cptr()[a4*CFpq.shape(1)+a0];
        }//a4
      }//a0
      orz::DTensor tmp2(_Nact,_Nact*Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a4)*tmp2.shape(1)+(a2*Nrhom2+P)] = W0.cptr()[a4*W0.shape(1)*W0.shape(2)+a2*W0.shape(2)+P];
          }//P
        }//a2
      }//a4
      // @W1(a0,a2,P) <<= +1 @Fcore(a4,a0) @W0(a4,a2,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a0*_Nact*Nrhom2+a2*Nrhom2+P] += tmp3.cptr()[a0*_Nact*Nrhom2+a2*Nrhom2+P];
          }//P
        }//a2
      }//a0
    }
    orz::DTensor W2(Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(Nrhom2,_Nact*_Nact);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a0*_Nact+a2)] = W1.cptr()[a0*W1.shape(1)*W1.shape(2)+a2*W1.shape(2)+P];
          }//a2
        }//a0
      }//P
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a0*_Nact+a2)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+_P_];
          }//_P_
        }//a2
      }//a0
      // @W2(P,_P_) <<= +1 @W1(a0,a2,P) @X(-2)(_P_,a0,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          W2.cptr()[P*Nrhom2+_P_] += tmp3.cptr()[P*Nrhom2+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(b*_Nvirt+a)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[b*T_m2_.shape(1)*T_m2_.shape(2)+a*T_m2_.shape(2)+_P_];
          }//_P_
        }//a
      }//b
      orz::DTensor tmp2(Nrhom2,Nrhom2);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&P : irange(0L, Nrhom2)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W2.cptr()[P*W2.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-2)(b,a,_P_) @W2(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[b*_Nvirt*Nrhom2+a*Nrhom2+P];
          }//P
        }//a
      }//b
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/4 @T(-2)(a,b,_P_) @W2(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a4*_Nact+a2)*tmp1.shape(1)+(a3*_Nact+a1)] = d2.cptr()[a4*d2.shape(1)*d2.shape(2)*d2.shape(3)+a3*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a1];
            }//a1
          }//a3
        }//a2
      }//a4
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a3*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a4,a2,P) <<= +1 @D2(a4,a3,a2,a1) @X(-2)(P,a3,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W0.cptr()[a4*_Nact*Nrhom2+a2*Nrhom2+P] += tmp3.cptr()[a4*_Nact*Nrhom2+a2*Nrhom2+P];
          }//P
        }//a2
      }//a4
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a4)] = CFpq.cptr()[a4*CFpq.shape(1)+a0];
        }//a4
      }//a0
      orz::DTensor tmp2(_Nact,_Nact*Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a4)*tmp2.shape(1)+(a2*Nrhom2+P)] = W0.cptr()[a4*W0.shape(1)*W0.shape(2)+a2*W0.shape(2)+P];
          }//P
        }//a2
      }//a4
      // @W1(a0,a2,P) <<= +1 @Fcore(a4,a0) @W0(a4,a2,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a0*_Nact*Nrhom2+a2*Nrhom2+P] += tmp3.cptr()[a0*_Nact*Nrhom2+a2*Nrhom2+P];
          }//P
        }//a2
      }//a0
    }
    orz::DTensor W2(Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(Nrhom2,_Nact*_Nact);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a0*_Nact+a2)] = W1.cptr()[a0*W1.shape(1)*W1.shape(2)+a2*W1.shape(2)+P];
          }//a2
        }//a0
      }//P
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a0*_Nact+a2)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+_P_];
          }//_P_
        }//a2
      }//a0
      // @W2(P,_P_) <<= +1 @W1(a0,a2,P) @X(-2)(_P_,a0,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          W2.cptr()[P*Nrhom2+_P_] += tmp3.cptr()[P*Nrhom2+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(a*_Nvirt+b)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[a*T_m2_.shape(1)*T_m2_.shape(2)+b*T_m2_.shape(2)+_P_];
          }//_P_
        }//b
      }//a
      orz::DTensor tmp2(Nrhom2,Nrhom2);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&P : irange(0L, Nrhom2)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W2.cptr()[P*W2.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-2)(a,b,_P_) @W2(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P];
          }//P
        }//b
      }//a
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/4 @T(-2)(a,b,_P_) @W2(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a4*_Nact+a2)*tmp1.shape(1)+(a3*_Nact+a1)] = d2.cptr()[a4*d2.shape(1)*d2.shape(2)*d2.shape(3)+a3*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a1];
            }//a1
          }//a3
        }//a2
      }//a4
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a1*Cm2.shape(1)*Cm2.shape(2)+a3*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a4,a2,P) <<= +1 @D2(a4,a3,a2,a1) @X(-2)(P,a1,a3)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W0.cptr()[a4*_Nact*Nrhom2+a2*Nrhom2+P] += tmp3.cptr()[a4*_Nact*Nrhom2+a2*Nrhom2+P];
          }//P
        }//a2
      }//a4
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a4)] = CFpq.cptr()[a4*CFpq.shape(1)+a0];
        }//a4
      }//a0
      orz::DTensor tmp2(_Nact,_Nact*Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a4)*tmp2.shape(1)+(a2*Nrhom2+P)] = W0.cptr()[a4*W0.shape(1)*W0.shape(2)+a2*W0.shape(2)+P];
          }//P
        }//a2
      }//a4
      // @W1(a0,a2,P) <<= +1 @Fcore(a4,a0) @W0(a4,a2,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a0*_Nact*Nrhom2+a2*Nrhom2+P] += tmp3.cptr()[a0*_Nact*Nrhom2+a2*Nrhom2+P];
          }//P
        }//a2
      }//a0
    }
    orz::DTensor W2(Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(Nrhom2,_Nact*_Nact);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a0*_Nact+a2)] = W1.cptr()[a0*W1.shape(1)*W1.shape(2)+a2*W1.shape(2)+P];
          }//a2
        }//a0
      }//P
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a0*_Nact+a2)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a2*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//_P_
        }//a2
      }//a0
      // @W2(P,_P_) <<= +1 @W1(a0,a2,P) @X(-2)(_P_,a2,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          W2.cptr()[P*Nrhom2+_P_] += tmp3.cptr()[P*Nrhom2+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(a*_Nvirt+b)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[a*T_m2_.shape(1)*T_m2_.shape(2)+b*T_m2_.shape(2)+_P_];
          }//_P_
        }//b
      }//a
      orz::DTensor tmp2(Nrhom2,Nrhom2);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&P : irange(0L, Nrhom2)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W2.cptr()[P*W2.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-2)(a,b,_P_) @W2(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P];
          }//P
        }//b
      }//a
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/4 @T(-2)(b,a,_P_) @W2(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a4*_Nact+a2)*tmp1.shape(1)+(a3*_Nact+a1)] = d2.cptr()[a4*d2.shape(1)*d2.shape(2)*d2.shape(3)+a3*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a1];
            }//a1
          }//a3
        }//a2
      }//a4
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a3*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a4,a2,P) <<= +1 @D2(a4,a3,a2,a1) @X(-2)(P,a3,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W0.cptr()[a4*_Nact*Nrhom2+a2*Nrhom2+P] += tmp3.cptr()[a4*_Nact*Nrhom2+a2*Nrhom2+P];
          }//P
        }//a2
      }//a4
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a4)] = CFpq.cptr()[a4*CFpq.shape(1)+a0];
        }//a4
      }//a0
      orz::DTensor tmp2(_Nact,_Nact*Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a4)*tmp2.shape(1)+(a2*Nrhom2+P)] = W0.cptr()[a4*W0.shape(1)*W0.shape(2)+a2*W0.shape(2)+P];
          }//P
        }//a2
      }//a4
      // @W1(a0,a2,P) <<= +1 @Fcore(a4,a0) @W0(a4,a2,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a0*_Nact*Nrhom2+a2*Nrhom2+P] += tmp3.cptr()[a0*_Nact*Nrhom2+a2*Nrhom2+P];
          }//P
        }//a2
      }//a0
    }
    orz::DTensor W2(Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(Nrhom2,_Nact*_Nact);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a0*_Nact+a2)] = W1.cptr()[a0*W1.shape(1)*W1.shape(2)+a2*W1.shape(2)+P];
          }//a2
        }//a0
      }//P
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a0*_Nact+a2)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a2*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//_P_
        }//a2
      }//a0
      // @W2(P,_P_) <<= +1 @W1(a0,a2,P) @X(-2)(_P_,a2,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          W2.cptr()[P*Nrhom2+_P_] += tmp3.cptr()[P*Nrhom2+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(b*_Nvirt+a)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[b*T_m2_.shape(1)*T_m2_.shape(2)+a*T_m2_.shape(2)+_P_];
          }//_P_
        }//a
      }//b
      orz::DTensor tmp2(Nrhom2,Nrhom2);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&P : irange(0L, Nrhom2)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W2.cptr()[P*W2.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-2)(b,a,_P_) @W2(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[b*_Nvirt*Nrhom2+a*Nrhom2+P];
          }//P
        }//a
      }//b
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= +1/4 @T(-2)(v0,b,_P_) @W2(a,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*_Nact,_Nact*_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a3 : irange(0L, _Nact)){
                for(auto &&a1 : irange(0L, _Nact)){
                  tmp1.cptr()[(a5*_Nact*_Nact*_Nact+a4*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a3*_Nact+a1)] = d3.cptr()[a5*d3.shape(1)*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a4*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a3*d3.shape(3)*d3.shape(4)*d3.shape(5)+a2*d3.shape(4)*d3.shape(5)+a1*d3.shape(5)+a0];
                }//a1
              }//a3
            }//a0
          }//a2
        }//a4
      }//a5
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a3*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a5,a4,a2,a0,P) <<= +1 @D3(a5,a4,a3,a2,a1,a0) @X(-2)(P,a3,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                W0.cptr()[a5*_Nact*_Nact*_Nact*Nrhom2+a4*_Nact*_Nact*Nrhom2+a2*_Nact*Nrhom2+a0*Nrhom2+P] += tmp3.cptr()[a5*_Nact*_Nact*_Nact*Nrhom2+a4*_Nact*_Nact*Nrhom2+a2*_Nact*Nrhom2+a0*Nrhom2+P];
              }//P
            }//a0
          }//a2
        }//a4
      }//a5
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Nvirt+v0)*tmp1.shape(1)+(a2*_Nact+a5)] = V2_FiC_vava.cptr()[a*V2_FiC_vava.shape(1)*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+a2*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+v0*V2_FiC_vava.shape(3)+a5];
            }//a5
          }//a2
        }//v0
      }//a
      orz::DTensor tmp2(_Nact*_Nact,_Nact*_Nact*Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                tmp2.cptr()[(a2*_Nact+a5)*tmp2.shape(1)+(a4*_Nact*Nrhom2+a0*Nrhom2+P)] = W0.cptr()[a5*W0.shape(1)*W0.shape(2)*W0.shape(3)*W0.shape(4)+a4*W0.shape(2)*W0.shape(3)*W0.shape(4)+a2*W0.shape(3)*W0.shape(4)+a0*W0.shape(4)+P];
              }//P
            }//a0
          }//a4
        }//a5
      }//a2
      // @W1(a4,a0,a,v0,P) <<= +1 @V2_FiC_vava(a,a2,v0,a5) @W0(a5,a4,a2,a0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                W1.cptr()[a4*_Nact*_Nvirt*_Nvirt*Nrhom2+a0*_Nvirt*_Nvirt*Nrhom2+a*_Nvirt*Nrhom2+v0*Nrhom2+P] += tmp3.cptr()[a*_Nvirt*_Nact*_Nact*Nrhom2+v0*_Nact*_Nact*Nrhom2+a4*_Nact*Nrhom2+a0*Nrhom2+P];
              }//P
            }//a0
          }//a4
        }//v0
      }//a
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*Nrhom2,_Nact*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            for(auto &&a4 : irange(0L, _Nact)){
              for(auto &&a0 : irange(0L, _Nact)){
                tmp1.cptr()[(a*_Nvirt*Nrhom2+v0*Nrhom2+P)*tmp1.shape(1)+(a4*_Nact+a0)] = W1.cptr()[a4*W1.shape(1)*W1.shape(2)*W1.shape(3)*W1.shape(4)+a0*W1.shape(2)*W1.shape(3)*W1.shape(4)+a*W1.shape(3)*W1.shape(4)+v0*W1.shape(4)+P];
              }//a0
            }//a4
          }//P
        }//v0
      }//a
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a4*_Nact+a0)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a4*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//_P_
        }//a0
      }//a4
      // @W2(a,v0,P,_P_) <<= +1 @W1(a4,a0,a,v0,P) @X(-2)(_P_,a4,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            for(auto &&_P_ : irange(0L, Nrhom2)){
              W2.cptr()[a*_Nvirt*Nrhom2*Nrhom2+v0*Nrhom2*Nrhom2+P*Nrhom2+_P_] += tmp3.cptr()[a*_Nvirt*Nrhom2*Nrhom2+v0*Nrhom2*Nrhom2+P*Nrhom2+_P_];
            }//_P_
          }//P
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt,_Nvirt*Nrhom2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(b)*tmp1.shape(1)+(v0*Nrhom2+_P_)] = T_m2_.cptr()[v0*T_m2_.shape(1)*T_m2_.shape(2)+b*T_m2_.shape(2)+_P_];
          }//_P_
        }//v0
      }//b
      orz::DTensor tmp2(_Nvirt*Nrhom2,_Nvirt*Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom2)){
              tmp2.cptr()[(v0*Nrhom2+_P_)*tmp2.shape(1)+(a*Nrhom2+P)] = W2.cptr()[a*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//a
        }//_P_
      }//v0
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-2)(v0,b,_P_) @W2(a,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[b*_Nvirt*Nrhom2+a*Nrhom2+P];
          }//P
        }//a
      }//b
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= +1/4 @T(-2)(v0,a,_P_) @W2(b,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*_Nact,_Nact*_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a3 : irange(0L, _Nact)){
                for(auto &&a1 : irange(0L, _Nact)){
                  tmp1.cptr()[(a5*_Nact*_Nact*_Nact+a4*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a3*_Nact+a1)] = d3.cptr()[a5*d3.shape(1)*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a4*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a3*d3.shape(3)*d3.shape(4)*d3.shape(5)+a2*d3.shape(4)*d3.shape(5)+a1*d3.shape(5)+a0];
                }//a1
              }//a3
            }//a0
          }//a2
        }//a4
      }//a5
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a1*Cm2.shape(1)*Cm2.shape(2)+a3*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a5,a4,a2,a0,P) <<= +1 @D3(a5,a4,a3,a2,a1,a0) @X(-2)(P,a1,a3)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                W0.cptr()[a5*_Nact*_Nact*_Nact*Nrhom2+a4*_Nact*_Nact*Nrhom2+a2*_Nact*Nrhom2+a0*Nrhom2+P] += tmp3.cptr()[a5*_Nact*_Nact*_Nact*Nrhom2+a4*_Nact*_Nact*Nrhom2+a2*_Nact*Nrhom2+a0*Nrhom2+P];
              }//P
            }//a0
          }//a2
        }//a4
      }//a5
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Nvirt+v0)*tmp1.shape(1)+(a2*_Nact+a5)] = V2_FiC_vava.cptr()[b*V2_FiC_vava.shape(1)*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+a2*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+v0*V2_FiC_vava.shape(3)+a5];
            }//a5
          }//a2
        }//v0
      }//b
      orz::DTensor tmp2(_Nact*_Nact,_Nact*_Nact*Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                tmp2.cptr()[(a2*_Nact+a5)*tmp2.shape(1)+(a4*_Nact*Nrhom2+a0*Nrhom2+P)] = W0.cptr()[a5*W0.shape(1)*W0.shape(2)*W0.shape(3)*W0.shape(4)+a4*W0.shape(2)*W0.shape(3)*W0.shape(4)+a2*W0.shape(3)*W0.shape(4)+a0*W0.shape(4)+P];
              }//P
            }//a0
          }//a4
        }//a5
      }//a2
      // @W1(a4,a0,b,v0,P) <<= +1 @V2_FiC_vava(b,a2,v0,a5) @W0(a5,a4,a2,a0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                W1.cptr()[a4*_Nact*_Nvirt*_Nvirt*Nrhom2+a0*_Nvirt*_Nvirt*Nrhom2+b*_Nvirt*Nrhom2+v0*Nrhom2+P] += tmp3.cptr()[b*_Nvirt*_Nact*_Nact*Nrhom2+v0*_Nact*_Nact*Nrhom2+a4*_Nact*Nrhom2+a0*Nrhom2+P];
              }//P
            }//a0
          }//a4
        }//v0
      }//b
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*Nrhom2,_Nact*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            for(auto &&a4 : irange(0L, _Nact)){
              for(auto &&a0 : irange(0L, _Nact)){
                tmp1.cptr()[(b*_Nvirt*Nrhom2+v0*Nrhom2+P)*tmp1.shape(1)+(a4*_Nact+a0)] = W1.cptr()[a4*W1.shape(1)*W1.shape(2)*W1.shape(3)*W1.shape(4)+a0*W1.shape(2)*W1.shape(3)*W1.shape(4)+b*W1.shape(3)*W1.shape(4)+v0*W1.shape(4)+P];
              }//a0
            }//a4
          }//P
        }//v0
      }//b
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a4*_Nact+a0)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a4*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//_P_
        }//a0
      }//a4
      // @W2(b,v0,P,_P_) <<= +1 @W1(a4,a0,b,v0,P) @X(-2)(_P_,a4,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            for(auto &&_P_ : irange(0L, Nrhom2)){
              W2.cptr()[b*_Nvirt*Nrhom2*Nrhom2+v0*Nrhom2*Nrhom2+P*Nrhom2+_P_] += tmp3.cptr()[b*_Nvirt*Nrhom2*Nrhom2+v0*Nrhom2*Nrhom2+P*Nrhom2+_P_];
            }//_P_
          }//P
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt,_Nvirt*Nrhom2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(a)*tmp1.shape(1)+(v0*Nrhom2+_P_)] = T_m2_.cptr()[v0*T_m2_.shape(1)*T_m2_.shape(2)+a*T_m2_.shape(2)+_P_];
          }//_P_
        }//v0
      }//a
      orz::DTensor tmp2(_Nvirt*Nrhom2,_Nvirt*Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom2)){
              tmp2.cptr()[(v0*Nrhom2+_P_)*tmp2.shape(1)+(b*Nrhom2+P)] = W2.cptr()[b*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//b
        }//_P_
      }//v0
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-2)(v0,a,_P_) @W2(b,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P];
          }//P
        }//b
      }//a
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= +1/4 @T(-2)(v0,a,_P_) @W2(b,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*_Nact,_Nact*_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a3 : irange(0L, _Nact)){
                for(auto &&a1 : irange(0L, _Nact)){
                  tmp1.cptr()[(a5*_Nact*_Nact*_Nact+a4*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a3*_Nact+a1)] = d3.cptr()[a5*d3.shape(1)*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a4*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a3*d3.shape(3)*d3.shape(4)*d3.shape(5)+a2*d3.shape(4)*d3.shape(5)+a1*d3.shape(5)+a0];
                }//a1
              }//a3
            }//a0
          }//a2
        }//a4
      }//a5
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a3*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a5,a4,a2,a0,P) <<= +1 @D3(a5,a4,a3,a2,a1,a0) @X(-2)(P,a3,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                W0.cptr()[a5*_Nact*_Nact*_Nact*Nrhom2+a4*_Nact*_Nact*Nrhom2+a2*_Nact*Nrhom2+a0*Nrhom2+P] += tmp3.cptr()[a5*_Nact*_Nact*_Nact*Nrhom2+a4*_Nact*_Nact*Nrhom2+a2*_Nact*Nrhom2+a0*Nrhom2+P];
              }//P
            }//a0
          }//a2
        }//a4
      }//a5
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Nvirt+v0)*tmp1.shape(1)+(a4*_Nact+a5)] = V2_FiC_vvaa.cptr()[b*V2_FiC_vvaa.shape(1)*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+v0*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+a4*V2_FiC_vvaa.shape(3)+a5];
            }//a5
          }//a4
        }//v0
      }//b
      orz::DTensor tmp2(_Nact*_Nact,_Nact*_Nact*Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                tmp2.cptr()[(a4*_Nact+a5)*tmp2.shape(1)+(a2*_Nact*Nrhom2+a0*Nrhom2+P)] = W0.cptr()[a5*W0.shape(1)*W0.shape(2)*W0.shape(3)*W0.shape(4)+a4*W0.shape(2)*W0.shape(3)*W0.shape(4)+a2*W0.shape(3)*W0.shape(4)+a0*W0.shape(4)+P];
              }//P
            }//a0
          }//a2
        }//a5
      }//a4
      // @W1(a2,a0,b,v0,P) <<= +1 @V2_FiC_vvaa(b,v0,a4,a5) @W0(a5,a4,a2,a0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                W1.cptr()[a2*_Nact*_Nvirt*_Nvirt*Nrhom2+a0*_Nvirt*_Nvirt*Nrhom2+b*_Nvirt*Nrhom2+v0*Nrhom2+P] += tmp3.cptr()[b*_Nvirt*_Nact*_Nact*Nrhom2+v0*_Nact*_Nact*Nrhom2+a2*_Nact*Nrhom2+a0*Nrhom2+P];
              }//P
            }//a0
          }//a2
        }//v0
      }//b
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*Nrhom2,_Nact*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&a0 : irange(0L, _Nact)){
                tmp1.cptr()[(b*_Nvirt*Nrhom2+v0*Nrhom2+P)*tmp1.shape(1)+(a2*_Nact+a0)] = W1.cptr()[a2*W1.shape(1)*W1.shape(2)*W1.shape(3)*W1.shape(4)+a0*W1.shape(2)*W1.shape(3)*W1.shape(4)+b*W1.shape(3)*W1.shape(4)+v0*W1.shape(4)+P];
              }//a0
            }//a2
          }//P
        }//v0
      }//b
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+_P_];
          }//_P_
        }//a0
      }//a2
      // @W2(b,v0,P,_P_) <<= +1 @W1(a2,a0,b,v0,P) @X(-2)(_P_,a0,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            for(auto &&_P_ : irange(0L, Nrhom2)){
              W2.cptr()[b*_Nvirt*Nrhom2*Nrhom2+v0*Nrhom2*Nrhom2+P*Nrhom2+_P_] += tmp3.cptr()[b*_Nvirt*Nrhom2*Nrhom2+v0*Nrhom2*Nrhom2+P*Nrhom2+_P_];
            }//_P_
          }//P
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt,_Nvirt*Nrhom2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(a)*tmp1.shape(1)+(v0*Nrhom2+_P_)] = T_m2_.cptr()[v0*T_m2_.shape(1)*T_m2_.shape(2)+a*T_m2_.shape(2)+_P_];
          }//_P_
        }//v0
      }//a
      orz::DTensor tmp2(_Nvirt*Nrhom2,_Nvirt*Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom2)){
              tmp2.cptr()[(v0*Nrhom2+_P_)*tmp2.shape(1)+(b*Nrhom2+P)] = W2.cptr()[b*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//b
        }//_P_
      }//v0
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-2)(v0,a,_P_) @W2(b,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P];
          }//P
        }//b
      }//a
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= +1/4 @T(-2)(v0,b,_P_) @W2(a,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*_Nact,_Nact*_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a3 : irange(0L, _Nact)){
                for(auto &&a1 : irange(0L, _Nact)){
                  tmp1.cptr()[(a5*_Nact*_Nact*_Nact+a4*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a3*_Nact+a1)] = d3.cptr()[a5*d3.shape(1)*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a4*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a3*d3.shape(3)*d3.shape(4)*d3.shape(5)+a2*d3.shape(4)*d3.shape(5)+a1*d3.shape(5)+a0];
                }//a1
              }//a3
            }//a0
          }//a2
        }//a4
      }//a5
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a3*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a5,a4,a2,a0,P) <<= +1 @D3(a5,a4,a3,a2,a1,a0) @X(-2)(P,a3,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                W0.cptr()[a5*_Nact*_Nact*_Nact*Nrhom2+a4*_Nact*_Nact*Nrhom2+a2*_Nact*Nrhom2+a0*Nrhom2+P] += tmp3.cptr()[a5*_Nact*_Nact*_Nact*Nrhom2+a4*_Nact*_Nact*Nrhom2+a2*_Nact*Nrhom2+a0*Nrhom2+P];
              }//P
            }//a0
          }//a2
        }//a4
      }//a5
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Nvirt+v0)*tmp1.shape(1)+(a4*_Nact+a5)] = V2_FiC_vvaa.cptr()[a*V2_FiC_vvaa.shape(1)*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+v0*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+a4*V2_FiC_vvaa.shape(3)+a5];
            }//a5
          }//a4
        }//v0
      }//a
      orz::DTensor tmp2(_Nact*_Nact,_Nact*_Nact*Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                tmp2.cptr()[(a4*_Nact+a5)*tmp2.shape(1)+(a2*_Nact*Nrhom2+a0*Nrhom2+P)] = W0.cptr()[a5*W0.shape(1)*W0.shape(2)*W0.shape(3)*W0.shape(4)+a4*W0.shape(2)*W0.shape(3)*W0.shape(4)+a2*W0.shape(3)*W0.shape(4)+a0*W0.shape(4)+P];
              }//P
            }//a0
          }//a2
        }//a5
      }//a4
      // @W1(a2,a0,a,v0,P) <<= +1 @V2_FiC_vvaa(a,v0,a4,a5) @W0(a5,a4,a2,a0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                W1.cptr()[a2*_Nact*_Nvirt*_Nvirt*Nrhom2+a0*_Nvirt*_Nvirt*Nrhom2+a*_Nvirt*Nrhom2+v0*Nrhom2+P] += tmp3.cptr()[a*_Nvirt*_Nact*_Nact*Nrhom2+v0*_Nact*_Nact*Nrhom2+a2*_Nact*Nrhom2+a0*Nrhom2+P];
              }//P
            }//a0
          }//a2
        }//v0
      }//a
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*Nrhom2,_Nact*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&a0 : irange(0L, _Nact)){
                tmp1.cptr()[(a*_Nvirt*Nrhom2+v0*Nrhom2+P)*tmp1.shape(1)+(a2*_Nact+a0)] = W1.cptr()[a2*W1.shape(1)*W1.shape(2)*W1.shape(3)*W1.shape(4)+a0*W1.shape(2)*W1.shape(3)*W1.shape(4)+a*W1.shape(3)*W1.shape(4)+v0*W1.shape(4)+P];
              }//a0
            }//a2
          }//P
        }//v0
      }//a
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a2*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//_P_
        }//a0
      }//a2
      // @W2(a,v0,P,_P_) <<= +1 @W1(a2,a0,a,v0,P) @X(-2)(_P_,a2,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            for(auto &&_P_ : irange(0L, Nrhom2)){
              W2.cptr()[a*_Nvirt*Nrhom2*Nrhom2+v0*Nrhom2*Nrhom2+P*Nrhom2+_P_] += tmp3.cptr()[a*_Nvirt*Nrhom2*Nrhom2+v0*Nrhom2*Nrhom2+P*Nrhom2+_P_];
            }//_P_
          }//P
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt,_Nvirt*Nrhom2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(b)*tmp1.shape(1)+(v0*Nrhom2+_P_)] = T_m2_.cptr()[v0*T_m2_.shape(1)*T_m2_.shape(2)+b*T_m2_.shape(2)+_P_];
          }//_P_
        }//v0
      }//b
      orz::DTensor tmp2(_Nvirt*Nrhom2,_Nvirt*Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom2)){
              tmp2.cptr()[(v0*Nrhom2+_P_)*tmp2.shape(1)+(a*Nrhom2+P)] = W2.cptr()[a*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//a
        }//_P_
      }//v0
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-2)(v0,b,_P_) @W2(a,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[b*_Nvirt*Nrhom2+a*Nrhom2+P];
          }//P
        }//a
      }//b
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= +1/4 @T(-2)(b,v0,_P_) @W2(a,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*_Nact,_Nact*_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a3 : irange(0L, _Nact)){
                for(auto &&a1 : irange(0L, _Nact)){
                  tmp1.cptr()[(a5*_Nact*_Nact*_Nact+a4*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a3*_Nact+a1)] = d3.cptr()[a5*d3.shape(1)*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a4*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a3*d3.shape(3)*d3.shape(4)*d3.shape(5)+a2*d3.shape(4)*d3.shape(5)+a1*d3.shape(5)+a0];
                }//a1
              }//a3
            }//a0
          }//a2
        }//a4
      }//a5
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a3*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a5,a4,a2,a0,P) <<= +1 @D3(a5,a4,a3,a2,a1,a0) @X(-2)(P,a3,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                W0.cptr()[a5*_Nact*_Nact*_Nact*Nrhom2+a4*_Nact*_Nact*Nrhom2+a2*_Nact*Nrhom2+a0*Nrhom2+P] += tmp3.cptr()[a5*_Nact*_Nact*_Nact*Nrhom2+a4*_Nact*_Nact*Nrhom2+a2*_Nact*Nrhom2+a0*Nrhom2+P];
              }//P
            }//a0
          }//a2
        }//a4
      }//a5
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Nvirt+v0)*tmp1.shape(1)+(a2*_Nact+a5)] = V2_FiC_vava.cptr()[a*V2_FiC_vava.shape(1)*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+a2*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+v0*V2_FiC_vava.shape(3)+a5];
            }//a5
          }//a2
        }//v0
      }//a
      orz::DTensor tmp2(_Nact*_Nact,_Nact*_Nact*Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                tmp2.cptr()[(a2*_Nact+a5)*tmp2.shape(1)+(a4*_Nact*Nrhom2+a0*Nrhom2+P)] = W0.cptr()[a5*W0.shape(1)*W0.shape(2)*W0.shape(3)*W0.shape(4)+a4*W0.shape(2)*W0.shape(3)*W0.shape(4)+a2*W0.shape(3)*W0.shape(4)+a0*W0.shape(4)+P];
              }//P
            }//a0
          }//a4
        }//a5
      }//a2
      // @W1(a4,a0,a,v0,P) <<= +1 @V2_FiC_vava(a,a2,v0,a5) @W0(a5,a4,a2,a0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                W1.cptr()[a4*_Nact*_Nvirt*_Nvirt*Nrhom2+a0*_Nvirt*_Nvirt*Nrhom2+a*_Nvirt*Nrhom2+v0*Nrhom2+P] += tmp3.cptr()[a*_Nvirt*_Nact*_Nact*Nrhom2+v0*_Nact*_Nact*Nrhom2+a4*_Nact*Nrhom2+a0*Nrhom2+P];
              }//P
            }//a0
          }//a4
        }//v0
      }//a
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*Nrhom2,_Nact*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            for(auto &&a4 : irange(0L, _Nact)){
              for(auto &&a0 : irange(0L, _Nact)){
                tmp1.cptr()[(a*_Nvirt*Nrhom2+v0*Nrhom2+P)*tmp1.shape(1)+(a4*_Nact+a0)] = W1.cptr()[a4*W1.shape(1)*W1.shape(2)*W1.shape(3)*W1.shape(4)+a0*W1.shape(2)*W1.shape(3)*W1.shape(4)+a*W1.shape(3)*W1.shape(4)+v0*W1.shape(4)+P];
              }//a0
            }//a4
          }//P
        }//v0
      }//a
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a4*_Nact+a0)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a4*Cm2.shape(2)+_P_];
          }//_P_
        }//a0
      }//a4
      // @W2(a,v0,P,_P_) <<= +1 @W1(a4,a0,a,v0,P) @X(-2)(_P_,a0,a4)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            for(auto &&_P_ : irange(0L, Nrhom2)){
              W2.cptr()[a*_Nvirt*Nrhom2*Nrhom2+v0*Nrhom2*Nrhom2+P*Nrhom2+_P_] += tmp3.cptr()[a*_Nvirt*Nrhom2*Nrhom2+v0*Nrhom2*Nrhom2+P*Nrhom2+_P_];
            }//_P_
          }//P
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt,_Nvirt*Nrhom2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(b)*tmp1.shape(1)+(v0*Nrhom2+_P_)] = T_m2_.cptr()[b*T_m2_.shape(1)*T_m2_.shape(2)+v0*T_m2_.shape(2)+_P_];
          }//_P_
        }//v0
      }//b
      orz::DTensor tmp2(_Nvirt*Nrhom2,_Nvirt*Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom2)){
              tmp2.cptr()[(v0*Nrhom2+_P_)*tmp2.shape(1)+(a*Nrhom2+P)] = W2.cptr()[a*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//a
        }//_P_
      }//v0
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-2)(b,v0,_P_) @W2(a,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[b*_Nvirt*Nrhom2+a*Nrhom2+P];
          }//P
        }//a
      }//b
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= +1/4 @T(-2)(a,v0,_P_) @W2(b,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*_Nact,_Nact*_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a3 : irange(0L, _Nact)){
                for(auto &&a1 : irange(0L, _Nact)){
                  tmp1.cptr()[(a5*_Nact*_Nact*_Nact+a4*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a3*_Nact+a1)] = d3.cptr()[a5*d3.shape(1)*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a4*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a3*d3.shape(3)*d3.shape(4)*d3.shape(5)+a2*d3.shape(4)*d3.shape(5)+a1*d3.shape(5)+a0];
                }//a1
              }//a3
            }//a0
          }//a2
        }//a4
      }//a5
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a1*Cm2.shape(1)*Cm2.shape(2)+a3*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a5,a4,a2,a0,P) <<= +1 @D3(a5,a4,a3,a2,a1,a0) @X(-2)(P,a1,a3)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                W0.cptr()[a5*_Nact*_Nact*_Nact*Nrhom2+a4*_Nact*_Nact*Nrhom2+a2*_Nact*Nrhom2+a0*Nrhom2+P] += tmp3.cptr()[a5*_Nact*_Nact*_Nact*Nrhom2+a4*_Nact*_Nact*Nrhom2+a2*_Nact*Nrhom2+a0*Nrhom2+P];
              }//P
            }//a0
          }//a2
        }//a4
      }//a5
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Nvirt+v0)*tmp1.shape(1)+(a2*_Nact+a5)] = V2_FiC_vava.cptr()[b*V2_FiC_vava.shape(1)*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+a2*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+v0*V2_FiC_vava.shape(3)+a5];
            }//a5
          }//a2
        }//v0
      }//b
      orz::DTensor tmp2(_Nact*_Nact,_Nact*_Nact*Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                tmp2.cptr()[(a2*_Nact+a5)*tmp2.shape(1)+(a4*_Nact*Nrhom2+a0*Nrhom2+P)] = W0.cptr()[a5*W0.shape(1)*W0.shape(2)*W0.shape(3)*W0.shape(4)+a4*W0.shape(2)*W0.shape(3)*W0.shape(4)+a2*W0.shape(3)*W0.shape(4)+a0*W0.shape(4)+P];
              }//P
            }//a0
          }//a4
        }//a5
      }//a2
      // @W1(a4,a0,b,v0,P) <<= +1 @V2_FiC_vava(b,a2,v0,a5) @W0(a5,a4,a2,a0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                W1.cptr()[a4*_Nact*_Nvirt*_Nvirt*Nrhom2+a0*_Nvirt*_Nvirt*Nrhom2+b*_Nvirt*Nrhom2+v0*Nrhom2+P] += tmp3.cptr()[b*_Nvirt*_Nact*_Nact*Nrhom2+v0*_Nact*_Nact*Nrhom2+a4*_Nact*Nrhom2+a0*Nrhom2+P];
              }//P
            }//a0
          }//a4
        }//v0
      }//b
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*Nrhom2,_Nact*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            for(auto &&a4 : irange(0L, _Nact)){
              for(auto &&a0 : irange(0L, _Nact)){
                tmp1.cptr()[(b*_Nvirt*Nrhom2+v0*Nrhom2+P)*tmp1.shape(1)+(a4*_Nact+a0)] = W1.cptr()[a4*W1.shape(1)*W1.shape(2)*W1.shape(3)*W1.shape(4)+a0*W1.shape(2)*W1.shape(3)*W1.shape(4)+b*W1.shape(3)*W1.shape(4)+v0*W1.shape(4)+P];
              }//a0
            }//a4
          }//P
        }//v0
      }//b
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a4*_Nact+a0)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a4*Cm2.shape(2)+_P_];
          }//_P_
        }//a0
      }//a4
      // @W2(b,v0,P,_P_) <<= +1 @W1(a4,a0,b,v0,P) @X(-2)(_P_,a0,a4)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            for(auto &&_P_ : irange(0L, Nrhom2)){
              W2.cptr()[b*_Nvirt*Nrhom2*Nrhom2+v0*Nrhom2*Nrhom2+P*Nrhom2+_P_] += tmp3.cptr()[b*_Nvirt*Nrhom2*Nrhom2+v0*Nrhom2*Nrhom2+P*Nrhom2+_P_];
            }//_P_
          }//P
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt,_Nvirt*Nrhom2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(a)*tmp1.shape(1)+(v0*Nrhom2+_P_)] = T_m2_.cptr()[a*T_m2_.shape(1)*T_m2_.shape(2)+v0*T_m2_.shape(2)+_P_];
          }//_P_
        }//v0
      }//a
      orz::DTensor tmp2(_Nvirt*Nrhom2,_Nvirt*Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom2)){
              tmp2.cptr()[(v0*Nrhom2+_P_)*tmp2.shape(1)+(b*Nrhom2+P)] = W2.cptr()[b*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//b
        }//_P_
      }//v0
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-2)(a,v0,_P_) @W2(b,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P];
          }//P
        }//b
      }//a
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= +1/4 @T(-2)(a,v0,_P_) @W2(b,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*_Nact,_Nact*_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a3 : irange(0L, _Nact)){
                for(auto &&a1 : irange(0L, _Nact)){
                  tmp1.cptr()[(a5*_Nact*_Nact*_Nact+a4*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a3*_Nact+a1)] = d3.cptr()[a5*d3.shape(1)*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a4*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a3*d3.shape(3)*d3.shape(4)*d3.shape(5)+a2*d3.shape(4)*d3.shape(5)+a1*d3.shape(5)+a0];
                }//a1
              }//a3
            }//a0
          }//a2
        }//a4
      }//a5
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a3*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a5,a4,a2,a0,P) <<= +1 @D3(a5,a4,a3,a2,a1,a0) @X(-2)(P,a3,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                W0.cptr()[a5*_Nact*_Nact*_Nact*Nrhom2+a4*_Nact*_Nact*Nrhom2+a2*_Nact*Nrhom2+a0*Nrhom2+P] += tmp3.cptr()[a5*_Nact*_Nact*_Nact*Nrhom2+a4*_Nact*_Nact*Nrhom2+a2*_Nact*Nrhom2+a0*Nrhom2+P];
              }//P
            }//a0
          }//a2
        }//a4
      }//a5
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Nvirt+v0)*tmp1.shape(1)+(a4*_Nact+a5)] = V2_FiC_vvaa.cptr()[b*V2_FiC_vvaa.shape(1)*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+v0*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+a4*V2_FiC_vvaa.shape(3)+a5];
            }//a5
          }//a4
        }//v0
      }//b
      orz::DTensor tmp2(_Nact*_Nact,_Nact*_Nact*Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                tmp2.cptr()[(a4*_Nact+a5)*tmp2.shape(1)+(a2*_Nact*Nrhom2+a0*Nrhom2+P)] = W0.cptr()[a5*W0.shape(1)*W0.shape(2)*W0.shape(3)*W0.shape(4)+a4*W0.shape(2)*W0.shape(3)*W0.shape(4)+a2*W0.shape(3)*W0.shape(4)+a0*W0.shape(4)+P];
              }//P
            }//a0
          }//a2
        }//a5
      }//a4
      // @W1(a2,a0,b,v0,P) <<= +1 @V2_FiC_vvaa(b,v0,a4,a5) @W0(a5,a4,a2,a0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                W1.cptr()[a2*_Nact*_Nvirt*_Nvirt*Nrhom2+a0*_Nvirt*_Nvirt*Nrhom2+b*_Nvirt*Nrhom2+v0*Nrhom2+P] += tmp3.cptr()[b*_Nvirt*_Nact*_Nact*Nrhom2+v0*_Nact*_Nact*Nrhom2+a2*_Nact*Nrhom2+a0*Nrhom2+P];
              }//P
            }//a0
          }//a2
        }//v0
      }//b
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*Nrhom2,_Nact*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&a0 : irange(0L, _Nact)){
                tmp1.cptr()[(b*_Nvirt*Nrhom2+v0*Nrhom2+P)*tmp1.shape(1)+(a2*_Nact+a0)] = W1.cptr()[a2*W1.shape(1)*W1.shape(2)*W1.shape(3)*W1.shape(4)+a0*W1.shape(2)*W1.shape(3)*W1.shape(4)+b*W1.shape(3)*W1.shape(4)+v0*W1.shape(4)+P];
              }//a0
            }//a2
          }//P
        }//v0
      }//b
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a2*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//_P_
        }//a0
      }//a2
      // @W2(b,v0,P,_P_) <<= +1 @W1(a2,a0,b,v0,P) @X(-2)(_P_,a2,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            for(auto &&_P_ : irange(0L, Nrhom2)){
              W2.cptr()[b*_Nvirt*Nrhom2*Nrhom2+v0*Nrhom2*Nrhom2+P*Nrhom2+_P_] += tmp3.cptr()[b*_Nvirt*Nrhom2*Nrhom2+v0*Nrhom2*Nrhom2+P*Nrhom2+_P_];
            }//_P_
          }//P
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt,_Nvirt*Nrhom2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(a)*tmp1.shape(1)+(v0*Nrhom2+_P_)] = T_m2_.cptr()[a*T_m2_.shape(1)*T_m2_.shape(2)+v0*T_m2_.shape(2)+_P_];
          }//_P_
        }//v0
      }//a
      orz::DTensor tmp2(_Nvirt*Nrhom2,_Nvirt*Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom2)){
              tmp2.cptr()[(v0*Nrhom2+_P_)*tmp2.shape(1)+(b*Nrhom2+P)] = W2.cptr()[b*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//b
        }//_P_
      }//v0
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-2)(a,v0,_P_) @W2(b,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P];
          }//P
        }//b
      }//a
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= +1/4 @T(-2)(b,v0,_P_) @W2(a,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*_Nact,_Nact*_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a3 : irange(0L, _Nact)){
                for(auto &&a1 : irange(0L, _Nact)){
                  tmp1.cptr()[(a5*_Nact*_Nact*_Nact+a4*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a3*_Nact+a1)] = d3.cptr()[a5*d3.shape(1)*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a4*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a3*d3.shape(3)*d3.shape(4)*d3.shape(5)+a2*d3.shape(4)*d3.shape(5)+a1*d3.shape(5)+a0];
                }//a1
              }//a3
            }//a0
          }//a2
        }//a4
      }//a5
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a3*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a5,a4,a2,a0,P) <<= +1 @D3(a5,a4,a3,a2,a1,a0) @X(-2)(P,a3,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                W0.cptr()[a5*_Nact*_Nact*_Nact*Nrhom2+a4*_Nact*_Nact*Nrhom2+a2*_Nact*Nrhom2+a0*Nrhom2+P] += tmp3.cptr()[a5*_Nact*_Nact*_Nact*Nrhom2+a4*_Nact*_Nact*Nrhom2+a2*_Nact*Nrhom2+a0*Nrhom2+P];
              }//P
            }//a0
          }//a2
        }//a4
      }//a5
    }
    orz::DTensor W1(_Nact,_Nact,_Nvirt,_Nvirt,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Nvirt+v0)*tmp1.shape(1)+(a4*_Nact+a5)] = V2_FiC_vvaa.cptr()[a*V2_FiC_vvaa.shape(1)*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+v0*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+a4*V2_FiC_vvaa.shape(3)+a5];
            }//a5
          }//a4
        }//v0
      }//a
      orz::DTensor tmp2(_Nact*_Nact,_Nact*_Nact*Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                tmp2.cptr()[(a4*_Nact+a5)*tmp2.shape(1)+(a2*_Nact*Nrhom2+a0*Nrhom2+P)] = W0.cptr()[a5*W0.shape(1)*W0.shape(2)*W0.shape(3)*W0.shape(4)+a4*W0.shape(2)*W0.shape(3)*W0.shape(4)+a2*W0.shape(3)*W0.shape(4)+a0*W0.shape(4)+P];
              }//P
            }//a0
          }//a2
        }//a5
      }//a4
      // @W1(a2,a0,a,v0,P) <<= +1 @V2_FiC_vvaa(a,v0,a4,a5) @W0(a5,a4,a2,a0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                W1.cptr()[a2*_Nact*_Nvirt*_Nvirt*Nrhom2+a0*_Nvirt*_Nvirt*Nrhom2+a*_Nvirt*Nrhom2+v0*Nrhom2+P] += tmp3.cptr()[a*_Nvirt*_Nact*_Nact*Nrhom2+v0*_Nact*_Nact*Nrhom2+a2*_Nact*Nrhom2+a0*Nrhom2+P];
              }//P
            }//a0
          }//a2
        }//v0
      }//a
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*Nrhom2,_Nact*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&a0 : irange(0L, _Nact)){
                tmp1.cptr()[(a*_Nvirt*Nrhom2+v0*Nrhom2+P)*tmp1.shape(1)+(a2*_Nact+a0)] = W1.cptr()[a2*W1.shape(1)*W1.shape(2)*W1.shape(3)*W1.shape(4)+a0*W1.shape(2)*W1.shape(3)*W1.shape(4)+a*W1.shape(3)*W1.shape(4)+v0*W1.shape(4)+P];
              }//a0
            }//a2
          }//P
        }//v0
      }//a
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+_P_];
          }//_P_
        }//a0
      }//a2
      // @W2(a,v0,P,_P_) <<= +1 @W1(a2,a0,a,v0,P) @X(-2)(_P_,a0,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            for(auto &&_P_ : irange(0L, Nrhom2)){
              W2.cptr()[a*_Nvirt*Nrhom2*Nrhom2+v0*Nrhom2*Nrhom2+P*Nrhom2+_P_] += tmp3.cptr()[a*_Nvirt*Nrhom2*Nrhom2+v0*Nrhom2*Nrhom2+P*Nrhom2+_P_];
            }//_P_
          }//P
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt,_Nvirt*Nrhom2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(b)*tmp1.shape(1)+(v0*Nrhom2+_P_)] = T_m2_.cptr()[b*T_m2_.shape(1)*T_m2_.shape(2)+v0*T_m2_.shape(2)+_P_];
          }//_P_
        }//v0
      }//b
      orz::DTensor tmp2(_Nvirt*Nrhom2,_Nvirt*Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom2)){
              tmp2.cptr()[(v0*Nrhom2+_P_)*tmp2.shape(1)+(a*Nrhom2+P)] = W2.cptr()[a*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//a
        }//_P_
      }//v0
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-2)(b,v0,_P_) @W2(a,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[b*_Nvirt*Nrhom2+a*Nrhom2+P];
          }//P
        }//a
      }//b
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= +1/4 @V2_FiC_vvvv(b,v0,a,v1) @W2(v1,v0,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a2*_Nact+a0)*tmp1.shape(1)+(a3*_Nact+a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a3
        }//a0
      }//a2
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a3*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-2)(P,a3,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W0.cptr()[a2*_Nact*Nrhom2+a0*Nrhom2+P] += tmp3.cptr()[a2*_Nact*Nrhom2+a0*Nrhom2+P];
          }//P
        }//a0
      }//a2
    }
    orz::DTensor W1(Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(Nrhom2,_Nact*_Nact);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a2*W0.shape(1)*W0.shape(2)+a0*W0.shape(2)+P];
          }//a0
        }//a2
      }//P
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a2*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//_P_
        }//a0
      }//a2
      // @W1(P,_P_) <<= +1 @W0(a2,a0,P) @X(-2)(_P_,a2,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          W1.cptr()[P*Nrhom2+_P_] += tmp3.cptr()[P*Nrhom2+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&v1 : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(v1*_Nvirt+v0)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[v1*T_m2_.shape(1)*T_m2_.shape(2)+v0*T_m2_.shape(2)+_P_];
          }//_P_
        }//v0
      }//v1
      orz::DTensor tmp2(Nrhom2,Nrhom2);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&P : irange(0L, Nrhom2)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @W2(v1,v0,P) <<= +1 @T(-2)(v1,v0,_P_) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v1 : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            W2.cptr()[v1*_Nvirt*Nrhom2+v0*Nrhom2+P] += tmp3.cptr()[v1*_Nvirt*Nrhom2+v0*Nrhom2+P];
          }//P
        }//v0
      }//v1
    }
    if(!do4EXTinPNO)    
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nvirt*_Nvirt);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&v1 : irange(0L, _Nvirt)){
              tmp1.cptr()[(b*_Nvirt+a)*tmp1.shape(1)+(v0*_Nvirt+v1)] = V2_FiC_vvvv.cptr()[b*V2_FiC_vvvv.shape(1)*V2_FiC_vvvv.shape(2)*V2_FiC_vvvv.shape(3)+v0*V2_FiC_vvvv.shape(2)*V2_FiC_vvvv.shape(3)+a*V2_FiC_vvvv.shape(3)+v1];
            }//v1
          }//v0
        }//a
      }//b
      orz::DTensor tmp2(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&v1 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(v0*_Nvirt+v1)*tmp2.shape(1)+(P)] = W2.cptr()[v1*W2.shape(1)*W2.shape(2)+v0*W2.shape(2)+P];
          }//P
        }//v1
      }//v0
      // @tRFiC(a,b,P) <<= +1 _X_ @V2_FiC_vvvv(b,v0,a,v1) @W2(v1,v0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[b*_Nvirt*Nrhom2+a*Nrhom2+P];
          }//P
        }//a
      }//b
    }
    else
    {
      //orz::ProgressTimer pt("  ==> 4ext contraction (PNO)");
      for(auto &&P : irange(0L, Nrhom2)){
        const long I_P = TActList.GetIndex(P);
        if (I_P < 0) continue;
        PairData P_P = TActList.GetPair(P);
        const long Nvir = TActList.GetNPNO(P);
        orz::DTensor I4;
        PNO4.GetTensor(P,I4);
        // Sort I4        
        orz::DTensor tmp1(Nvir*Nvir,Nvir*Nvir);
        const long Nvir2 = Nvir*(Nvir+1)/2;        
        for(auto &&a : irange(0L, Nvir)){
          for(auto &&b : irange(0L, Nvir)){
            const long idx_ab = (b<a ? a*(a+1)/2+b : b*(b+1)/2+a);  
            for(auto &&c : irange(0L, Nvir)){
              for(auto &&d : irange(0L, Nvir)){
                const long idx_cd = (d<c ? c*(c+1)/2+d : d*(d+1)/2+c);
                const long idx_abcd = (idx_cd<idx_ab ? idx_ab*(idx_ab+1)/2+idx_cd : idx_cd*(idx_cd+1)/2+idx_ab);  
                tmp1.cptr()[b*Nvir*Nvir*Nvir+d*Nvir*Nvir+a*Nvir+c] = I4.cptr()[idx_abcd];
              }
            }
          }
        }
        orz::DTensor tmp2(_Nvirt,_Nvirt);
        for(auto &&v1 : irange(0L, _Nvirt))        
          for(auto &&v0 : irange(0L, _Nvirt))
            tmp2.cptr()[v0*_Nvirt+v1] = W2.cptr()[v1*W2.shape(1)*W2.shape(2)+v0*W2.shape(2)+P];
        P_P.Transform2EXT(tmp2);
        tmp2.reshape_inplace(1,Nvir*Nvir);                
        orz::DTensor tmp3(tmp2*tmp1);
        tmp3 *= _X_;
        orz::DTensor S2 = R2.GetTensor(P);
        for(auto &&a : irange(0L, Nvir))        
          for(auto &&b : irange(0L, Nvir))
            S2.cptr()[a*S2.shape(1)+b] += tmp3.cptr()[b*Nvir+a];
        R2.PutTensor(P,S2);
      }//P  
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= +1/4 @V2_FiC_vvvv(a,v0,b,v1) @W2(v1,v0,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a2*_Nact+a0)*tmp1.shape(1)+(a3*_Nact+a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a3
        }//a0
      }//a2
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a3*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-2)(P,a3,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W0.cptr()[a2*_Nact*Nrhom2+a0*Nrhom2+P] += tmp3.cptr()[a2*_Nact*Nrhom2+a0*Nrhom2+P];
          }//P
        }//a0
      }//a2
    }
    orz::DTensor W1(Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(Nrhom2,_Nact*_Nact);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a2*W0.shape(1)*W0.shape(2)+a0*W0.shape(2)+P];
          }//a0
        }//a2
      }//P
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+_P_];
          }//_P_
        }//a0
      }//a2
      // @W1(P,_P_) <<= +1 @W0(a2,a0,P) @X(-2)(_P_,a0,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          W1.cptr()[P*Nrhom2+_P_] += tmp3.cptr()[P*Nrhom2+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom2);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&v1 : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(v1*_Nvirt+v0)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[v1*T_m2_.shape(1)*T_m2_.shape(2)+v0*T_m2_.shape(2)+_P_];
          }//_P_
        }//v0
      }//v1
      orz::DTensor tmp2(Nrhom2,Nrhom2);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&P : irange(0L, Nrhom2)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @W2(v1,v0,P) <<= +1 @T(-2)(v1,v0,_P_) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v1 : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            W2.cptr()[v1*_Nvirt*Nrhom2+v0*Nrhom2+P] += tmp3.cptr()[v1*_Nvirt*Nrhom2+v0*Nrhom2+P];
          }//P
        }//v0
      }//v1
    }
    if(!do4EXTinPNO)        
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nvirt*_Nvirt);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&v1 : irange(0L, _Nvirt)){
              tmp1.cptr()[(a*_Nvirt+b)*tmp1.shape(1)+(v0*_Nvirt+v1)] = V2_FiC_vvvv.cptr()[a*V2_FiC_vvvv.shape(1)*V2_FiC_vvvv.shape(2)*V2_FiC_vvvv.shape(3)+v0*V2_FiC_vvvv.shape(2)*V2_FiC_vvvv.shape(3)+b*V2_FiC_vvvv.shape(3)+v1];
            }//v1
          }//v0
        }//b
      }//a
      orz::DTensor tmp2(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&v1 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(v0*_Nvirt+v1)*tmp2.shape(1)+(P)] = W2.cptr()[v1*W2.shape(1)*W2.shape(2)+v0*W2.shape(2)+P];
          }//P
        }//v1
      }//v0
      // @tRFiC(a,b,P) <<= +1 _X_ @V2_FiC_vvvv(a,v0,b,v1) @W2(v1,v0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P];
          }//P
        }//b
      }//a
    }
    else
    {
      //orz::ProgressTimer pt("  ==> 4ext contraction (PNO)");
      for(auto &&P : irange(0L, Nrhom2)){
        const long I_P = TActList.GetIndex(P);
        if (I_P < 0) continue;
        PairData P_P = TActList.GetPair(P);
        const long Nvir = TActList.GetNPNO(P);
        orz::DTensor I4;
        PNO4.GetTensor(P,I4);
        // Sort I4        
        orz::DTensor tmp1(Nvir*Nvir,Nvir*Nvir);
        const long Nvir2 = Nvir*(Nvir+1)/2;        
        for(auto &&a : irange(0L, Nvir)){
          for(auto &&b : irange(0L, Nvir)){
            const long idx_ab = (b<a ? a*(a+1)/2+b : b*(b+1)/2+a);  
            for(auto &&c : irange(0L, Nvir)){
              for(auto &&d : irange(0L, Nvir)){
                const long idx_cd = (d<c ? c*(c+1)/2+d : d*(d+1)/2+c);
                const long idx_abcd = (idx_cd<idx_ab ? idx_ab*(idx_ab+1)/2+idx_cd : idx_cd*(idx_cd+1)/2+idx_ab);  
                tmp1.cptr()[b*Nvir*Nvir*Nvir+d*Nvir*Nvir+a*Nvir+c] = I4.cptr()[idx_abcd];
              }
            }
          }
        }
        orz::DTensor tmp2(_Nvirt,_Nvirt);
        for(auto &&v1 : irange(0L, _Nvirt))        
          for(auto &&v0 : irange(0L, _Nvirt))
            tmp2.cptr()[v0*_Nvirt+v1] = W2.cptr()[v1*W2.shape(1)*W2.shape(2)+v0*W2.shape(2)+P];
        P_P.Transform2EXT(tmp2);
        tmp2.reshape_inplace(1,Nvir*Nvir);
        orz::DTensor tmp3(tmp2*tmp1);
        tmp3 *= _X_;
        orz::DTensor S2 = R2.GetTensor(P);
        for(auto &&a : irange(0L, Nvir))        
          for(auto &&b : irange(0L, Nvir))
            S2.cptr()[a*S2.shape(1)+b] += tmp3.cptr()[a*Nvir+b];
        R2.PutTensor(P,S2);
      }//P  
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/4 @T(-2)(b,a,_P_) @W2(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*_Nact,_Nact*_Nact);
      for(auto &&a6 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&a3 : irange(0L, _Nact)){
                for(auto &&a1 : irange(0L, _Nact)){
                  tmp1.cptr()[(a6*_Nact*_Nact*_Nact+a5*_Nact*_Nact+a4*_Nact+a2)*tmp1.shape(1)+(a3*_Nact+a1)] = d3.cptr()[a6*d3.shape(1)*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a5*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a4*d3.shape(3)*d3.shape(4)*d3.shape(5)+a3*d3.shape(4)*d3.shape(5)+a2*d3.shape(5)+a1];
                }//a1
              }//a3
            }//a2
          }//a4
        }//a5
      }//a6
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a1*Cm2.shape(1)*Cm2.shape(2)+a3*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a6,a5,a4,a2,P) <<= +1 @D3(a6,a5,a4,a3,a2,a1) @X(-2)(P,a1,a3)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a6 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                W0.cptr()[a6*_Nact*_Nact*_Nact*Nrhom2+a5*_Nact*_Nact*Nrhom2+a4*_Nact*Nrhom2+a2*Nrhom2+P] += tmp3.cptr()[a6*_Nact*_Nact*_Nact*Nrhom2+a5*_Nact*_Nact*Nrhom2+a4*_Nact*Nrhom2+a2*Nrhom2+P];
              }//P
            }//a2
          }//a4
        }//a5
      }//a6
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact,_Nact*_Nact*_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a5 : irange(0L, _Nact)){
            for(auto &&a6 : irange(0L, _Nact)){
              tmp1.cptr()[(a0)*tmp1.shape(1)+(a4*_Nact*_Nact+a5*_Nact+a6)] = V2_FiC_aaaa.cptr()[a0*V2_FiC_aaaa.shape(1)*V2_FiC_aaaa.shape(2)*V2_FiC_aaaa.shape(3)+a4*V2_FiC_aaaa.shape(2)*V2_FiC_aaaa.shape(3)+a5*V2_FiC_aaaa.shape(3)+a6];
            }//a6
          }//a5
        }//a4
      }//a0
      orz::DTensor tmp2(_Nact*_Nact*_Nact,_Nact*Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&a6 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                tmp2.cptr()[(a4*_Nact*_Nact+a5*_Nact+a6)*tmp2.shape(1)+(a2*Nrhom2+P)] = W0.cptr()[a6*W0.shape(1)*W0.shape(2)*W0.shape(3)*W0.shape(4)+a5*W0.shape(2)*W0.shape(3)*W0.shape(4)+a4*W0.shape(3)*W0.shape(4)+a2*W0.shape(4)+P];
              }//P
            }//a2
          }//a6
        }//a5
      }//a4
      // @W1(a0,a2,P) <<= +1 @V2_FiC_aaaa(a0,a4,a5,a6) @W0(a6,a5,a4,a2,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a0*_Nact*Nrhom2+a2*Nrhom2+P] += tmp3.cptr()[a0*_Nact*Nrhom2+a2*Nrhom2+P];
          }//P
        }//a2
      }//a0
    }
    orz::DTensor W2(Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(Nrhom2,_Nact*_Nact);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a0*_Nact+a2)] = W1.cptr()[a0*W1.shape(1)*W1.shape(2)+a2*W1.shape(2)+P];
          }//a2
        }//a0
      }//P
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a0*_Nact+a2)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+_P_];
          }//_P_
        }//a2
      }//a0
      // @W2(P,_P_) <<= +1 @W1(a0,a2,P) @X(-2)(_P_,a0,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          W2.cptr()[P*Nrhom2+_P_] += tmp3.cptr()[P*Nrhom2+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(b*_Nvirt+a)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[b*T_m2_.shape(1)*T_m2_.shape(2)+a*T_m2_.shape(2)+_P_];
          }//_P_
        }//a
      }//b
      orz::DTensor tmp2(Nrhom2,Nrhom2);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&P : irange(0L, Nrhom2)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W2.cptr()[P*W2.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-2)(b,a,_P_) @W2(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[b*_Nvirt*Nrhom2+a*Nrhom2+P];
          }//P
        }//a
      }//b
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/4 @T(-2)(a,b,_P_) @W2(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*_Nact,_Nact*_Nact);
      for(auto &&a6 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&a3 : irange(0L, _Nact)){
                for(auto &&a1 : irange(0L, _Nact)){
                  tmp1.cptr()[(a6*_Nact*_Nact*_Nact+a5*_Nact*_Nact+a4*_Nact+a2)*tmp1.shape(1)+(a3*_Nact+a1)] = d3.cptr()[a6*d3.shape(1)*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a5*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a4*d3.shape(3)*d3.shape(4)*d3.shape(5)+a3*d3.shape(4)*d3.shape(5)+a2*d3.shape(5)+a1];
                }//a1
              }//a3
            }//a2
          }//a4
        }//a5
      }//a6
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a3*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a6,a5,a4,a2,P) <<= +1 @D3(a6,a5,a4,a3,a2,a1) @X(-2)(P,a3,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a6 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                W0.cptr()[a6*_Nact*_Nact*_Nact*Nrhom2+a5*_Nact*_Nact*Nrhom2+a4*_Nact*Nrhom2+a2*Nrhom2+P] += tmp3.cptr()[a6*_Nact*_Nact*_Nact*Nrhom2+a5*_Nact*_Nact*Nrhom2+a4*_Nact*Nrhom2+a2*Nrhom2+P];
              }//P
            }//a2
          }//a4
        }//a5
      }//a6
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact,_Nact*_Nact*_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a5 : irange(0L, _Nact)){
            for(auto &&a6 : irange(0L, _Nact)){
              tmp1.cptr()[(a0)*tmp1.shape(1)+(a4*_Nact*_Nact+a5*_Nact+a6)] = V2_FiC_aaaa.cptr()[a0*V2_FiC_aaaa.shape(1)*V2_FiC_aaaa.shape(2)*V2_FiC_aaaa.shape(3)+a4*V2_FiC_aaaa.shape(2)*V2_FiC_aaaa.shape(3)+a5*V2_FiC_aaaa.shape(3)+a6];
            }//a6
          }//a5
        }//a4
      }//a0
      orz::DTensor tmp2(_Nact*_Nact*_Nact,_Nact*Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&a6 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                tmp2.cptr()[(a4*_Nact*_Nact+a5*_Nact+a6)*tmp2.shape(1)+(a2*Nrhom2+P)] = W0.cptr()[a6*W0.shape(1)*W0.shape(2)*W0.shape(3)*W0.shape(4)+a5*W0.shape(2)*W0.shape(3)*W0.shape(4)+a4*W0.shape(3)*W0.shape(4)+a2*W0.shape(4)+P];
              }//P
            }//a2
          }//a6
        }//a5
      }//a4
      // @W1(a0,a2,P) <<= +1 @V2_FiC_aaaa(a0,a4,a5,a6) @W0(a6,a5,a4,a2,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a0*_Nact*Nrhom2+a2*Nrhom2+P] += tmp3.cptr()[a0*_Nact*Nrhom2+a2*Nrhom2+P];
          }//P
        }//a2
      }//a0
    }
    orz::DTensor W2(Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(Nrhom2,_Nact*_Nact);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a0*_Nact+a2)] = W1.cptr()[a0*W1.shape(1)*W1.shape(2)+a2*W1.shape(2)+P];
          }//a2
        }//a0
      }//P
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a0*_Nact+a2)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+_P_];
          }//_P_
        }//a2
      }//a0
      // @W2(P,_P_) <<= +1 @W1(a0,a2,P) @X(-2)(_P_,a0,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          W2.cptr()[P*Nrhom2+_P_] += tmp3.cptr()[P*Nrhom2+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(a*_Nvirt+b)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[a*T_m2_.shape(1)*T_m2_.shape(2)+b*T_m2_.shape(2)+_P_];
          }//_P_
        }//b
      }//a
      orz::DTensor tmp2(Nrhom2,Nrhom2);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&P : irange(0L, Nrhom2)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W2.cptr()[P*W2.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-2)(a,b,_P_) @W2(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P];
          }//P
        }//b
      }//a
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/4 @T(-2)(a,b,_P_) @W2(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*_Nact,_Nact*_Nact);
      for(auto &&a6 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&a3 : irange(0L, _Nact)){
                for(auto &&a1 : irange(0L, _Nact)){
                  tmp1.cptr()[(a6*_Nact*_Nact*_Nact+a5*_Nact*_Nact+a4*_Nact+a2)*tmp1.shape(1)+(a3*_Nact+a1)] = d3.cptr()[a6*d3.shape(1)*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a5*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a4*d3.shape(3)*d3.shape(4)*d3.shape(5)+a3*d3.shape(4)*d3.shape(5)+a2*d3.shape(5)+a1];
                }//a1
              }//a3
            }//a2
          }//a4
        }//a5
      }//a6
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a1*Cm2.shape(1)*Cm2.shape(2)+a3*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a6,a5,a4,a2,P) <<= +1 @D3(a6,a5,a4,a3,a2,a1) @X(-2)(P,a1,a3)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a6 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                W0.cptr()[a6*_Nact*_Nact*_Nact*Nrhom2+a5*_Nact*_Nact*Nrhom2+a4*_Nact*Nrhom2+a2*Nrhom2+P] += tmp3.cptr()[a6*_Nact*_Nact*_Nact*Nrhom2+a5*_Nact*_Nact*Nrhom2+a4*_Nact*Nrhom2+a2*Nrhom2+P];
              }//P
            }//a2
          }//a4
        }//a5
      }//a6
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact,_Nact*_Nact*_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a5 : irange(0L, _Nact)){
            for(auto &&a6 : irange(0L, _Nact)){
              tmp1.cptr()[(a0)*tmp1.shape(1)+(a4*_Nact*_Nact+a5*_Nact+a6)] = V2_FiC_aaaa.cptr()[a0*V2_FiC_aaaa.shape(1)*V2_FiC_aaaa.shape(2)*V2_FiC_aaaa.shape(3)+a4*V2_FiC_aaaa.shape(2)*V2_FiC_aaaa.shape(3)+a5*V2_FiC_aaaa.shape(3)+a6];
            }//a6
          }//a5
        }//a4
      }//a0
      orz::DTensor tmp2(_Nact*_Nact*_Nact,_Nact*Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&a6 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                tmp2.cptr()[(a4*_Nact*_Nact+a5*_Nact+a6)*tmp2.shape(1)+(a2*Nrhom2+P)] = W0.cptr()[a6*W0.shape(1)*W0.shape(2)*W0.shape(3)*W0.shape(4)+a5*W0.shape(2)*W0.shape(3)*W0.shape(4)+a4*W0.shape(3)*W0.shape(4)+a2*W0.shape(4)+P];
              }//P
            }//a2
          }//a6
        }//a5
      }//a4
      // @W1(a0,a2,P) <<= +1 @V2_FiC_aaaa(a0,a4,a5,a6) @W0(a6,a5,a4,a2,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a0*_Nact*Nrhom2+a2*Nrhom2+P] += tmp3.cptr()[a0*_Nact*Nrhom2+a2*Nrhom2+P];
          }//P
        }//a2
      }//a0
    }
    orz::DTensor W2(Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(Nrhom2,_Nact*_Nact);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a0*_Nact+a2)] = W1.cptr()[a0*W1.shape(1)*W1.shape(2)+a2*W1.shape(2)+P];
          }//a2
        }//a0
      }//P
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a0*_Nact+a2)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a2*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//_P_
        }//a2
      }//a0
      // @W2(P,_P_) <<= +1 @W1(a0,a2,P) @X(-2)(_P_,a2,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          W2.cptr()[P*Nrhom2+_P_] += tmp3.cptr()[P*Nrhom2+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(a*_Nvirt+b)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[a*T_m2_.shape(1)*T_m2_.shape(2)+b*T_m2_.shape(2)+_P_];
          }//_P_
        }//b
      }//a
      orz::DTensor tmp2(Nrhom2,Nrhom2);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&P : irange(0L, Nrhom2)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W2.cptr()[P*W2.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-2)(a,b,_P_) @W2(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P];
          }//P
        }//b
      }//a
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/4 @T(-2)(b,a,_P_) @W2(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*_Nact,_Nact*_Nact);
      for(auto &&a6 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&a3 : irange(0L, _Nact)){
                for(auto &&a1 : irange(0L, _Nact)){
                  tmp1.cptr()[(a6*_Nact*_Nact*_Nact+a5*_Nact*_Nact+a4*_Nact+a2)*tmp1.shape(1)+(a3*_Nact+a1)] = d3.cptr()[a6*d3.shape(1)*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a5*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a4*d3.shape(3)*d3.shape(4)*d3.shape(5)+a3*d3.shape(4)*d3.shape(5)+a2*d3.shape(5)+a1];
                }//a1
              }//a3
            }//a2
          }//a4
        }//a5
      }//a6
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(P)] = Cm2.cptr()[a3*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+P];
          }//P
        }//a1
      }//a3
      // @W0(a6,a5,a4,a2,P) <<= +1 @D3(a6,a5,a4,a3,a2,a1) @X(-2)(P,a3,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a6 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                W0.cptr()[a6*_Nact*_Nact*_Nact*Nrhom2+a5*_Nact*_Nact*Nrhom2+a4*_Nact*Nrhom2+a2*Nrhom2+P] += tmp3.cptr()[a6*_Nact*_Nact*_Nact*Nrhom2+a5*_Nact*_Nact*Nrhom2+a4*_Nact*Nrhom2+a2*Nrhom2+P];
              }//P
            }//a2
          }//a4
        }//a5
      }//a6
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact,_Nact*_Nact*_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a5 : irange(0L, _Nact)){
            for(auto &&a6 : irange(0L, _Nact)){
              tmp1.cptr()[(a0)*tmp1.shape(1)+(a4*_Nact*_Nact+a5*_Nact+a6)] = V2_FiC_aaaa.cptr()[a0*V2_FiC_aaaa.shape(1)*V2_FiC_aaaa.shape(2)*V2_FiC_aaaa.shape(3)+a4*V2_FiC_aaaa.shape(2)*V2_FiC_aaaa.shape(3)+a5*V2_FiC_aaaa.shape(3)+a6];
            }//a6
          }//a5
        }//a4
      }//a0
      orz::DTensor tmp2(_Nact*_Nact*_Nact,_Nact*Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&a6 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom2)){
                tmp2.cptr()[(a4*_Nact*_Nact+a5*_Nact+a6)*tmp2.shape(1)+(a2*Nrhom2+P)] = W0.cptr()[a6*W0.shape(1)*W0.shape(2)*W0.shape(3)*W0.shape(4)+a5*W0.shape(2)*W0.shape(3)*W0.shape(4)+a4*W0.shape(3)*W0.shape(4)+a2*W0.shape(4)+P];
              }//P
            }//a2
          }//a6
        }//a5
      }//a4
      // @W1(a0,a2,P) <<= +1 @V2_FiC_aaaa(a0,a4,a5,a6) @W0(a6,a5,a4,a2,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a0*_Nact*Nrhom2+a2*Nrhom2+P] += tmp3.cptr()[a0*_Nact*Nrhom2+a2*Nrhom2+P];
          }//P
        }//a2
      }//a0
    }
    orz::DTensor W2(Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(Nrhom2,_Nact*_Nact);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a0*_Nact+a2)] = W1.cptr()[a0*W1.shape(1)*W1.shape(2)+a2*W1.shape(2)+P];
          }//a2
        }//a0
      }//P
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a0*_Nact+a2)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a2*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//_P_
        }//a2
      }//a0
      // @W2(P,_P_) <<= +1 @W1(a0,a2,P) @X(-2)(_P_,a2,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          W2.cptr()[P*Nrhom2+_P_] += tmp3.cptr()[P*Nrhom2+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(b*_Nvirt+a)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[b*T_m2_.shape(1)*T_m2_.shape(2)+a*T_m2_.shape(2)+_P_];
          }//_P_
        }//a
      }//b
      orz::DTensor tmp2(Nrhom2,Nrhom2);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&P : irange(0L, Nrhom2)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W2.cptr()[P*W2.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-2)(b,a,_P_) @W2(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[b*_Nvirt*Nrhom2+a*Nrhom2+P];
          }//P
        }//a
      }//b
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/4 @T(-2)(a,b,_P_) @W2(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a5*_Nact+a3)*tmp1.shape(1)+(a4*_Nact+a2)] = d2.cptr()[a5*d2.shape(1)*d2.shape(2)*d2.shape(3)+a4*d2.shape(2)*d2.shape(3)+a3*d2.shape(3)+a2];
            }//a2
          }//a4
        }//a3
      }//a5
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a4*_Nact+a2)*tmp2.shape(1)+(P)] = Cm2.cptr()[a4*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+P];
          }//P
        }//a2
      }//a4
      // @W0(a5,a3,P) <<= +1 @D2(a5,a4,a3,a2) @X(-2)(P,a4,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W0.cptr()[a5*_Nact*Nrhom2+a3*Nrhom2+P] += tmp3.cptr()[a5*_Nact*Nrhom2+a3*Nrhom2+P];
          }//P
        }//a3
      }//a5
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              tmp1.cptr()[(a0*_Nact+a1)*tmp1.shape(1)+(a3*_Nact+a5)] = V2_FiC_aaaa.cptr()[a0*V2_FiC_aaaa.shape(1)*V2_FiC_aaaa.shape(2)*V2_FiC_aaaa.shape(3)+a3*V2_FiC_aaaa.shape(2)*V2_FiC_aaaa.shape(3)+a1*V2_FiC_aaaa.shape(3)+a5];
            }//a5
          }//a3
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a5)*tmp2.shape(1)+(P)] = W0.cptr()[a5*W0.shape(1)*W0.shape(2)+a3*W0.shape(2)+P];
          }//P
        }//a5
      }//a3
      // @W1(a0,a1,P) <<= +1 @V2_FiC_aaaa(a0,a3,a1,a5) @W0(a5,a3,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a0*_Nact*Nrhom2+a1*Nrhom2+P] += tmp3.cptr()[a0*_Nact*Nrhom2+a1*Nrhom2+P];
          }//P
        }//a1
      }//a0
    }
    orz::DTensor W2(Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(Nrhom2,_Nact*_Nact);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a0*_Nact+a1)] = W1.cptr()[a0*W1.shape(1)*W1.shape(2)+a1*W1.shape(2)+P];
          }//a1
        }//a0
      }//P
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a0*_Nact+a1)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a1*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+_P_];
          }//_P_
        }//a1
      }//a0
      // @W2(P,_P_) <<= +1 @W1(a0,a1,P) @X(-2)(_P_,a1,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          W2.cptr()[P*Nrhom2+_P_] += tmp3.cptr()[P*Nrhom2+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(a*_Nvirt+b)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[a*T_m2_.shape(1)*T_m2_.shape(2)+b*T_m2_.shape(2)+_P_];
          }//_P_
        }//b
      }//a
      orz::DTensor tmp2(Nrhom2,Nrhom2);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&P : irange(0L, Nrhom2)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W2.cptr()[P*W2.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-2)(a,b,_P_) @W2(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P];
          }//P
        }//b
      }//a
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/4 @T(-2)(b,a,_P_) @W2(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.2500000000000000;
    orz::DTensor W0(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a5*_Nact+a3)*tmp1.shape(1)+(a4*_Nact+a2)] = d2.cptr()[a5*d2.shape(1)*d2.shape(2)*d2.shape(3)+a4*d2.shape(2)*d2.shape(3)+a3*d2.shape(3)+a2];
            }//a2
          }//a4
        }//a3
      }//a5
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a4*_Nact+a2)*tmp2.shape(1)+(P)] = Cm2.cptr()[a4*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+P];
          }//P
        }//a2
      }//a4
      // @W0(a5,a3,P) <<= +1 @D2(a5,a4,a3,a2) @X(-2)(P,a4,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W0.cptr()[a5*_Nact*Nrhom2+a3*Nrhom2+P] += tmp3.cptr()[a5*_Nact*Nrhom2+a3*Nrhom2+P];
          }//P
        }//a3
      }//a5
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              tmp1.cptr()[(a0*_Nact+a1)*tmp1.shape(1)+(a3*_Nact+a5)] = V2_FiC_aaaa.cptr()[a0*V2_FiC_aaaa.shape(1)*V2_FiC_aaaa.shape(2)*V2_FiC_aaaa.shape(3)+a3*V2_FiC_aaaa.shape(2)*V2_FiC_aaaa.shape(3)+a1*V2_FiC_aaaa.shape(3)+a5];
            }//a5
          }//a3
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a5 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a5)*tmp2.shape(1)+(P)] = W0.cptr()[a5*W0.shape(1)*W0.shape(2)+a3*W0.shape(2)+P];
          }//P
        }//a5
      }//a3
      // @W1(a0,a1,P) <<= +1 @V2_FiC_aaaa(a0,a3,a1,a5) @W0(a5,a3,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a0*_Nact*Nrhom2+a1*Nrhom2+P] += tmp3.cptr()[a0*_Nact*Nrhom2+a1*Nrhom2+P];
          }//P
        }//a1
      }//a0
    }
    orz::DTensor W2(Nrhom2,Nrhom2);
    {
      orz::DTensor tmp1(Nrhom2,_Nact*_Nact);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a0*_Nact+a1)] = W1.cptr()[a0*W1.shape(1)*W1.shape(2)+a1*W1.shape(2)+P];
          }//a1
        }//a0
      }//P
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a0*_Nact+a1)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a1*Cm2.shape(2)+_P_];
          }//_P_
        }//a1
      }//a0
      // @W2(P,_P_) <<= +1 @W1(a0,a1,P) @X(-2)(_P_,a0,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&_P_ : irange(0L, Nrhom2)){
          W2.cptr()[P*Nrhom2+_P_] += tmp3.cptr()[P*Nrhom2+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(b*_Nvirt+a)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[b*T_m2_.shape(1)*T_m2_.shape(2)+a*T_m2_.shape(2)+_P_];
          }//_P_
        }//a
      }//b
      orz::DTensor tmp2(Nrhom2,Nrhom2);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&P : irange(0L, Nrhom2)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W2.cptr()[P*W2.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-2)(b,a,_P_) @W2(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[b*_Nvirt*Nrhom2+a*Nrhom2+P];
          }//P
        }//a
      }//b
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  //genCode.insert(std::make_pair("tRFiC", boost::shared_ptr<orz::DTensor>( new orz::DTensor(tRFiC) )));
 //Code-End-Code-End-Code-End-Code-End-Code-End-Code-End-Code-End
            
      // Transform tRFiC into PNO basis
      for(auto &&P : irange(0L, Nrhom2)){
        const long I_P = RActList.GetIndex(P);
        if (I_P < 0) continue;
        orz::DTensor RP(_Nvirt, _Nvirt);
        for(auto &&a : irange(0L, _Nvirt))
          for(auto &&b : irange(0L, _Nvirt))
            RP(a,b) = tRFiC(a,b,P);
        orz::DTensor Rorig = R2.GetTensor(P);
        RActList.GetPair(P).Transform2EXT(RP);
        RP += Rorig;
        R2.PutTensor(P, RP);
      }
            
      return;
    }//End func
  } }
