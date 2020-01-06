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
    void FMatrixOrthogonalizer::FormLHS_Vm2_Vm1_CEd(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &c2, const orz::DTensor &d3, DensityPack &dPack,
                                                    TensorContainer &ovl, TensorContainer &PNO4, const bool do4EXTinPNO,
                                                    const double QE0, const double Ecas, const double h0_energy, TensorContainer &DebugInts, const orz::DTensor &casfock, const orz::DTensor &FV,
                                                    TensorContainer &T2, TensorContainer &R2){

#include "lct_preparations_orthdens.cpp"      

      // Back-transform PNO ampitude into canonical basis
      orz::DTensor T_m1_ (_Nvirt,_Nvirt,Nrhom1,_Ndocc);
      orz::DTensor Tc_m1_(_Nvirt,_Nvirt,Nrhom1,_Ndocc);      
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          const long I_Pi = TActList.GetIndex(P,i);
          if (I_Pi < 0) continue;
          orz::DTensor U2 = T2.GetTensor(P,i);
          orz::DTensor U2tilde(U2.copy());
          U2tilde.gaxpy(+2.0, U2.swapdim(1,0), -1.0);          
          TActList.GetPair(P,i).BackTransform2EXT(U2     );
          TActList.GetPair(P,i).BackTransform2EXT(U2tilde);          
          for(auto &&a : irange(0L, _Nvirt))
            for(auto &&b : irange(0L, _Nvirt)){
              T_m1_ (a,b,P,i) = U2     (a,b);
              Tc_m1_(a,b,P,i) = U2tilde(a,b);
            }//ab
        }//i
      }//P

#include "lct_integrals4EXT.cpp"

#define _PROFILE_LHS

 //Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn
  orz::DTensor tRFiC(_Nvirt,_Nvirt,Nrhom2);
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= +1/2 @W1(a3,_P_,P) @W2(a3,a,b,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a0
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a3,a2,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(-)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_];
            }//_P_
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*Nrhom1,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*Nrhom1+_P_)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a0
          }//a2
        }//_P_
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a2
      // @W1(a3,_P_,P) <<= +1 @W0(a3,a2,a0,_P_) @X(-2)(P,a0,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P];
          }//P
        }//_P_
      }//a3
    }
    orz::DTensor W2(_Nact,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*Nrhom1,_Nvirt*_Ndocc);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(a*Nrhom1+_P_)*tmp1.shape(1)+(v0*_Ndocc+c0)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//c0
          }//v0
        }//_P_
      }//a
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(b*_Nact+a3)] = V2_FiC_vavc.cptr()[b*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a3*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+v0*V2_FiC_vavc.shape(3)+c0];
            }//a3
          }//b
        }//c0
      }//v0
      // @W2(a3,a,b,_P_) <<= +1 @Tc(-1)(v0,a,_P_,c0) @V2_FiC_vavc(b,a3,v0,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              W2.cptr()[a3*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+_P_] += tmp3.cptr()[a*Nrhom1*_Nvirt*_Nact+_P_*_Nvirt*_Nact+b*_Nact+a3];
            }//a3
          }//b
        }//_P_
      }//a
    }
    {
      orz::DTensor tmp1(Nrhom2,_Nact*Nrhom1);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a3*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)+_P_*W1.shape(2)+P];
          }//_P_
        }//a3
      }//P
      orz::DTensor tmp2(_Nact*Nrhom1,_Nvirt*_Nvirt);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(a3*Nrhom1+_P_)*tmp2.shape(1)+(a*_Nvirt+b)] = W2.cptr()[a3*W2.shape(1)*W2.shape(2)*W2.shape(3)+a*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+_P_];
            }//b
          }//a
        }//_P_
      }//a3
      // @tRFiC(a,b,P) <<= +1 _X_ @W1(a3,_P_,P) @W2(a3,a,b,_P_)
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= +1/2 @W1(a3,_P_,P) @W2(a3,a,b,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a0
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a3,a2,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_];
            }//_P_
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*Nrhom1,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*Nrhom1+_P_)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a0
          }//a2
        }//_P_
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a2
      // @W1(a3,_P_,P) <<= +1 @W0(a3,a2,a0,_P_) @X(-2)(P,a0,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P];
          }//P
        }//_P_
      }//a3
    }
    orz::DTensor W2(_Nact,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*Nrhom1,_Nvirt*_Ndocc);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(a*Nrhom1+_P_)*tmp1.shape(1)+(v0*_Ndocc+c0)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//c0
          }//v0
        }//_P_
      }//a
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(b*_Nact+a3)] = V2_FiC_vavc.cptr()[b*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a3*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+v0*V2_FiC_vavc.shape(3)+c0];
            }//a3
          }//b
        }//c0
      }//v0
      // @W2(a3,a,b,_P_) <<= +1 @Tc(-1)(a,v0,_P_,c0) @V2_FiC_vavc(b,a3,v0,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              W2.cptr()[a3*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+_P_] += tmp3.cptr()[a*Nrhom1*_Nvirt*_Nact+_P_*_Nvirt*_Nact+b*_Nact+a3];
            }//a3
          }//b
        }//_P_
      }//a
    }
    {
      orz::DTensor tmp1(Nrhom2,_Nact*Nrhom1);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a3*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)+_P_*W1.shape(2)+P];
          }//_P_
        }//a3
      }//P
      orz::DTensor tmp2(_Nact*Nrhom1,_Nvirt*_Nvirt);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(a3*Nrhom1+_P_)*tmp2.shape(1)+(a*_Nvirt+b)] = W2.cptr()[a3*W2.shape(1)*W2.shape(2)*W2.shape(3)+a*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+_P_];
            }//b
          }//a
        }//_P_
      }//a3
      // @tRFiC(a,b,P) <<= +1 _X_ @W1(a3,_P_,P) @W2(a3,a,b,_P_)
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= +1/2 @W1(a3,_P_,P) @W2(a3,b,a,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a0
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a3,a2,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(-)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_];
            }//_P_
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*Nrhom1,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*Nrhom1+_P_)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a0
          }//a2
        }//_P_
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a2*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a2
      // @W1(a3,_P_,P) <<= +1 @W0(a3,a2,a0,_P_) @X(-2)(P,a2,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P];
          }//P
        }//_P_
      }//a3
    }
    orz::DTensor W2(_Nact,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*Nrhom1,_Nvirt*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*Nrhom1+_P_)*tmp1.shape(1)+(v0*_Ndocc+c0)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//c0
          }//v0
        }//_P_
      }//b
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(a*_Nact+a3)] = V2_FiC_vavc.cptr()[a*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a3*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+v0*V2_FiC_vavc.shape(3)+c0];
            }//a3
          }//a
        }//c0
      }//v0
      // @W2(a3,b,a,_P_) <<= +1 @Tc(-1)(v0,b,_P_,c0) @V2_FiC_vavc(a,a3,v0,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              W2.cptr()[a3*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+_P_] += tmp3.cptr()[b*Nrhom1*_Nvirt*_Nact+_P_*_Nvirt*_Nact+a*_Nact+a3];
            }//a3
          }//a
        }//_P_
      }//b
    }
    {
      orz::DTensor tmp1(Nrhom2,_Nact*Nrhom1);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a3*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)+_P_*W1.shape(2)+P];
          }//_P_
        }//a3
      }//P
      orz::DTensor tmp2(_Nact*Nrhom1,_Nvirt*_Nvirt);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(a3*Nrhom1+_P_)*tmp2.shape(1)+(b*_Nvirt+a)] = W2.cptr()[a3*W2.shape(1)*W2.shape(2)*W2.shape(3)+b*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+_P_];
            }//a
          }//b
        }//_P_
      }//a3
      // @tRFiC(a,b,P) <<= +1 _X_ @W1(a3,_P_,P) @W2(a3,b,a,_P_)
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= +1/2 @W1(a3,_P_,P) @W2(a3,b,a,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a0
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a3,a2,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_];
            }//_P_
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*Nrhom1,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*Nrhom1+_P_)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a0
          }//a2
        }//_P_
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a2*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a2
      // @W1(a3,_P_,P) <<= +1 @W0(a3,a2,a0,_P_) @X(-2)(P,a2,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P];
          }//P
        }//_P_
      }//a3
    }
    orz::DTensor W2(_Nact,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*Nrhom1,_Nvirt*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*Nrhom1+_P_)*tmp1.shape(1)+(v0*_Ndocc+c0)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//c0
          }//v0
        }//_P_
      }//b
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(a*_Nact+a3)] = V2_FiC_vavc.cptr()[a*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a3*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+v0*V2_FiC_vavc.shape(3)+c0];
            }//a3
          }//a
        }//c0
      }//v0
      // @W2(a3,b,a,_P_) <<= +1 @Tc(-1)(b,v0,_P_,c0) @V2_FiC_vavc(a,a3,v0,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              W2.cptr()[a3*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+_P_] += tmp3.cptr()[b*Nrhom1*_Nvirt*_Nact+_P_*_Nvirt*_Nact+a*_Nact+a3];
            }//a3
          }//a
        }//_P_
      }//b
    }
    {
      orz::DTensor tmp1(Nrhom2,_Nact*Nrhom1);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a3*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)+_P_*W1.shape(2)+P];
          }//_P_
        }//a3
      }//P
      orz::DTensor tmp2(_Nact*Nrhom1,_Nvirt*_Nvirt);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(a3*Nrhom1+_P_)*tmp2.shape(1)+(b*_Nvirt+a)] = W2.cptr()[a3*W2.shape(1)*W2.shape(2)*W2.shape(3)+b*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+_P_];
            }//a
          }//b
        }//_P_
      }//a3
      // @tRFiC(a,b,P) <<= +1 _X_ @W1(a3,_P_,P) @W2(a3,b,a,_P_)
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/2 @W1(a3,_P_,P) @W2(a3,a,b,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a0
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a3,a2,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_];
            }//_P_
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*Nrhom1,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*Nrhom1+_P_)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a0
          }//a2
        }//_P_
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a2
      // @W1(a3,_P_,P) <<= +1 @W0(a3,a2,a0,_P_) @X(-2)(P,a0,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P];
          }//P
        }//_P_
      }//a3
    }
    orz::DTensor W2(_Nact,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Ndocc);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          tmp1.cptr()[(a3)*tmp1.shape(1)+(c0)] = CFip.cptr()[c0*CFip.shape(1)+a3];
        }//c0
      }//a3
      orz::DTensor tmp2(_Ndocc,_Nvirt*_Nvirt*Nrhom1);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(c0)*tmp2.shape(1)+(a*_Nvirt*Nrhom1+b*Nrhom1+_P_)] = T_m1_.cptr()[a*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+b*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+c0];
            }//_P_
          }//b
        }//a
      }//c0
      // @W2(a3,a,b,_P_) <<= +1 @Fcore(c0,a3) @T(-1)(a,b,_P_,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[a3*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+_P_] += tmp3.cptr()[a3*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+_P_];
            }//_P_
          }//b
        }//a
      }//a3
    }
    {
      orz::DTensor tmp1(Nrhom2,_Nact*Nrhom1);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a3*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)+_P_*W1.shape(2)+P];
          }//_P_
        }//a3
      }//P
      orz::DTensor tmp2(_Nact*Nrhom1,_Nvirt*_Nvirt);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(a3*Nrhom1+_P_)*tmp2.shape(1)+(a*_Nvirt+b)] = W2.cptr()[a3*W2.shape(1)*W2.shape(2)*W2.shape(3)+a*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+_P_];
            }//b
          }//a
        }//_P_
      }//a3
      // @tRFiC(a,b,P) <<= +1 _X_ @W1(a3,_P_,P) @W2(a3,a,b,_P_)
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/2 @W1(a3,_P_,P) @W2(a3,b,a,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a0
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a3,a2,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(-)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_];
            }//_P_
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*Nrhom1,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*Nrhom1+_P_)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a0
          }//a2
        }//_P_
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a2
      // @W1(a3,_P_,P) <<= +1 @W0(a3,a2,a0,_P_) @X(-2)(P,a0,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P];
          }//P
        }//_P_
      }//a3
    }
    orz::DTensor W2(_Nact,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Ndocc);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          tmp1.cptr()[(a3)*tmp1.shape(1)+(c0)] = CFip.cptr()[c0*CFip.shape(1)+a3];
        }//c0
      }//a3
      orz::DTensor tmp2(_Ndocc,_Nvirt*_Nvirt*Nrhom1);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(c0)*tmp2.shape(1)+(b*_Nvirt*Nrhom1+a*Nrhom1+_P_)] = T_m1_.cptr()[b*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+a*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+c0];
            }//_P_
          }//a
        }//b
      }//c0
      // @W2(a3,b,a,_P_) <<= +1 @Fcore(c0,a3) @T(-1)(b,a,_P_,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[a3*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+_P_] += tmp3.cptr()[a3*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+_P_];
            }//_P_
          }//a
        }//b
      }//a3
    }
    {
      orz::DTensor tmp1(Nrhom2,_Nact*Nrhom1);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a3*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)+_P_*W1.shape(2)+P];
          }//_P_
        }//a3
      }//P
      orz::DTensor tmp2(_Nact*Nrhom1,_Nvirt*_Nvirt);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(a3*Nrhom1+_P_)*tmp2.shape(1)+(b*_Nvirt+a)] = W2.cptr()[a3*W2.shape(1)*W2.shape(2)*W2.shape(3)+b*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+_P_];
            }//a
          }//b
        }//_P_
      }//a3
      // @tRFiC(a,b,P) <<= +1 _X_ @W1(a3,_P_,P) @W2(a3,b,a,_P_)
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/2 @W1(a3,_P_,P) @W2(a3,b,a,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a0
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a3,a2,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_];
            }//_P_
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*Nrhom1,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*Nrhom1+_P_)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a0
          }//a2
        }//_P_
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a2*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a2
      // @W1(a3,_P_,P) <<= +1 @W0(a3,a2,a0,_P_) @X(-2)(P,a2,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P];
          }//P
        }//_P_
      }//a3
    }
    orz::DTensor W2(_Nact,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Ndocc);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          tmp1.cptr()[(a3)*tmp1.shape(1)+(c0)] = CFip.cptr()[c0*CFip.shape(1)+a3];
        }//c0
      }//a3
      orz::DTensor tmp2(_Ndocc,_Nvirt*_Nvirt*Nrhom1);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(c0)*tmp2.shape(1)+(b*_Nvirt*Nrhom1+a*Nrhom1+_P_)] = T_m1_.cptr()[b*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+a*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+c0];
            }//_P_
          }//a
        }//b
      }//c0
      // @W2(a3,b,a,_P_) <<= +1 @Fcore(c0,a3) @T(-1)(b,a,_P_,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[a3*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+_P_] += tmp3.cptr()[a3*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+_P_];
            }//_P_
          }//a
        }//b
      }//a3
    }
    {
      orz::DTensor tmp1(Nrhom2,_Nact*Nrhom1);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a3*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)+_P_*W1.shape(2)+P];
          }//_P_
        }//a3
      }//P
      orz::DTensor tmp2(_Nact*Nrhom1,_Nvirt*_Nvirt);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(a3*Nrhom1+_P_)*tmp2.shape(1)+(b*_Nvirt+a)] = W2.cptr()[a3*W2.shape(1)*W2.shape(2)*W2.shape(3)+b*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+_P_];
            }//a
          }//b
        }//_P_
      }//a3
      // @tRFiC(a,b,P) <<= +1 _X_ @W1(a3,_P_,P) @W2(a3,b,a,_P_)
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/2 @W1(a3,_P_,P) @W2(a3,a,b,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a0
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a3,a2,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(-)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_];
            }//_P_
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*Nrhom1,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*Nrhom1+_P_)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a0
          }//a2
        }//_P_
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a2*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a2
      // @W1(a3,_P_,P) <<= +1 @W0(a3,a2,a0,_P_) @X(-2)(P,a2,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P];
          }//P
        }//_P_
      }//a3
    }
    orz::DTensor W2(_Nact,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Ndocc);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          tmp1.cptr()[(a3)*tmp1.shape(1)+(c0)] = CFip.cptr()[c0*CFip.shape(1)+a3];
        }//c0
      }//a3
      orz::DTensor tmp2(_Ndocc,_Nvirt*_Nvirt*Nrhom1);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(c0)*tmp2.shape(1)+(a*_Nvirt*Nrhom1+b*Nrhom1+_P_)] = T_m1_.cptr()[a*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+b*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+c0];
            }//_P_
          }//b
        }//a
      }//c0
      // @W2(a3,a,b,_P_) <<= +1 @Fcore(c0,a3) @T(-1)(a,b,_P_,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[a3*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+_P_] += tmp3.cptr()[a3*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+_P_];
            }//_P_
          }//b
        }//a
      }//a3
    }
    {
      orz::DTensor tmp1(Nrhom2,_Nact*Nrhom1);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a3*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)+_P_*W1.shape(2)+P];
          }//_P_
        }//a3
      }//P
      orz::DTensor tmp2(_Nact*Nrhom1,_Nvirt*_Nvirt);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(a3*Nrhom1+_P_)*tmp2.shape(1)+(a*_Nvirt+b)] = W2.cptr()[a3*W2.shape(1)*W2.shape(2)*W2.shape(3)+a*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+_P_];
            }//b
          }//a
        }//_P_
      }//a3
      // @tRFiC(a,b,P) <<= +1 _X_ @W1(a3,_P_,P) @W2(a3,a,b,_P_)
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/2 @T(-1)(a,b,_P_,c0) @W2(c0,_P_,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*_Nact*_Nact,_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&a0 : irange(0L, _Nact)){
                for(auto &&a1 : irange(0L, _Nact)){
                  tmp1.cptr()[(a5*_Nact*_Nact*_Nact*_Nact+a4*_Nact*_Nact*_Nact+a3*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a1)] = d3.cptr()[a5*d3.shape(1)*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a4*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a3*d3.shape(3)*d3.shape(4)*d3.shape(5)+a2*d3.shape(4)*d3.shape(5)+a1*d3.shape(5)+a0];
                }//a1
              }//a0
            }//a2
          }//a3
        }//a4
      }//a5
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a5,a4,a3,a2,a0,_P_) <<= +1 @D3(a5,a4,a3,a2,a1,a0) @X(-1)(+)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&a0 : irange(0L, _Nact)){
                for(auto &&_P_ : irange(0L, Nrhom1)){
                  W0.cptr()[a5*_Nact*_Nact*_Nact*_Nact*Nrhom1+a4*_Nact*_Nact*_Nact*Nrhom1+a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a5*_Nact*_Nact*_Nact*_Nact*Nrhom1+a4*_Nact*_Nact*_Nact*Nrhom1+a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_];
                }//_P_
              }//a0
            }//a2
          }//a3
        }//a4
      }//a5
    }
    orz::DTensor W1(_Nact,_Nact,_Nact,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*Nrhom1,_Nact*_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              for(auto &&a4 : irange(0L, _Nact)){
                for(auto &&a0 : irange(0L, _Nact)){
                  tmp1.cptr()[(a5*_Nact*_Nact*Nrhom1+a3*_Nact*Nrhom1+a2*Nrhom1+_P_)*tmp1.shape(1)+(a4*_Nact+a0)] = W0.cptr()[a5*W0.shape(1)*W0.shape(2)*W0.shape(3)*W0.shape(4)*W0.shape(5)+a4*W0.shape(2)*W0.shape(3)*W0.shape(4)*W0.shape(5)+a3*W0.shape(3)*W0.shape(4)*W0.shape(5)+a2*W0.shape(4)*W0.shape(5)+a0*W0.shape(5)+_P_];
                }//a0
              }//a4
            }//_P_
          }//a2
        }//a3
      }//a5
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a4*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a4*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a4
      // @W1(a5,a3,a2,_P_,P) <<= +1 @W0(a5,a4,a3,a2,a0,_P_) @X(-2)(P,a0,a4)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              for(auto &&P : irange(0L, Nrhom2)){
                W1.cptr()[a5*_Nact*_Nact*Nrhom1*Nrhom2+a3*_Nact*Nrhom1*Nrhom2+a2*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[a5*_Nact*_Nact*Nrhom1*Nrhom2+a3*_Nact*Nrhom1*Nrhom2+a2*Nrhom1*Nrhom2+_P_*Nrhom2+P];
              }//P
            }//_P_
          }//a2
        }//a3
      }//a5
    }
    orz::DTensor W2(_Ndocc,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Ndocc,_Nact*_Nact*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              tmp1.cptr()[(c0)*tmp1.shape(1)+(a2*_Nact*_Nact+a3*_Nact+a5)] = V2_FiC_aaac.cptr()[a2*V2_FiC_aaac.shape(1)*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a3*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a5*V2_FiC_aaac.shape(3)+c0];
            }//a5
          }//a3
        }//a2
      }//c0
      orz::DTensor tmp2(_Nact*_Nact*_Nact,Nrhom1*Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a5 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              for(auto &&P : irange(0L, Nrhom2)){
                tmp2.cptr()[(a2*_Nact*_Nact+a3*_Nact+a5)*tmp2.shape(1)+(_P_*Nrhom2+P)] = W1.cptr()[a5*W1.shape(1)*W1.shape(2)*W1.shape(3)*W1.shape(4)+a3*W1.shape(2)*W1.shape(3)*W1.shape(4)+a2*W1.shape(3)*W1.shape(4)+_P_*W1.shape(4)+P];
              }//P
            }//_P_
          }//a5
        }//a3
      }//a2
      // @W2(c0,_P_,P) <<= +1 @V2_FiC_aaac(a2,a3,a5,c0) @W1(a5,a3,a2,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&P : irange(0L, Nrhom2)){
            W2.cptr()[c0*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[c0*Nrhom1*Nrhom2+_P_*Nrhom2+P];
          }//P
        }//_P_
      }//c0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom1*_Ndocc);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(a*_Nvirt+b)*tmp1.shape(1)+(_P_*_Ndocc+c0)] = T_m1_.cptr()[a*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+b*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+c0];
            }//c0
          }//_P_
        }//b
      }//a
      orz::DTensor tmp2(Nrhom1*_Ndocc,Nrhom2);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(_P_*_Ndocc+c0)*tmp2.shape(1)+(P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)+_P_*W2.shape(2)+P];
          }//P
        }//c0
      }//_P_
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-1)(a,b,_P_,c0) @W2(c0,_P_,P)
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/2 @T(-1)(b,a,_P_,c0) @W2(c0,_P_,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*_Nact*_Nact,_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&a0 : irange(0L, _Nact)){
                for(auto &&a1 : irange(0L, _Nact)){
                  tmp1.cptr()[(a5*_Nact*_Nact*_Nact*_Nact+a4*_Nact*_Nact*_Nact+a3*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a1)] = d3.cptr()[a5*d3.shape(1)*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a4*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a3*d3.shape(3)*d3.shape(4)*d3.shape(5)+a2*d3.shape(4)*d3.shape(5)+a1*d3.shape(5)+a0];
                }//a1
              }//a0
            }//a2
          }//a3
        }//a4
      }//a5
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a5,a4,a3,a2,a0,_P_) <<= +1 @D3(a5,a4,a3,a2,a1,a0) @X(-1)(-)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&a0 : irange(0L, _Nact)){
                for(auto &&_P_ : irange(0L, Nrhom1)){
                  W0.cptr()[a5*_Nact*_Nact*_Nact*_Nact*Nrhom1+a4*_Nact*_Nact*_Nact*Nrhom1+a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a5*_Nact*_Nact*_Nact*_Nact*Nrhom1+a4*_Nact*_Nact*_Nact*Nrhom1+a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_];
                }//_P_
              }//a0
            }//a2
          }//a3
        }//a4
      }//a5
    }
    orz::DTensor W1(_Nact,_Nact,_Nact,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*Nrhom1,_Nact*_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              for(auto &&a4 : irange(0L, _Nact)){
                for(auto &&a0 : irange(0L, _Nact)){
                  tmp1.cptr()[(a5*_Nact*_Nact*Nrhom1+a3*_Nact*Nrhom1+a2*Nrhom1+_P_)*tmp1.shape(1)+(a4*_Nact+a0)] = W0.cptr()[a5*W0.shape(1)*W0.shape(2)*W0.shape(3)*W0.shape(4)*W0.shape(5)+a4*W0.shape(2)*W0.shape(3)*W0.shape(4)*W0.shape(5)+a3*W0.shape(3)*W0.shape(4)*W0.shape(5)+a2*W0.shape(4)*W0.shape(5)+a0*W0.shape(5)+_P_];
                }//a0
              }//a4
            }//_P_
          }//a2
        }//a3
      }//a5
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a4*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a4*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a4
      // @W1(a5,a3,a2,_P_,P) <<= +1 @W0(a5,a4,a3,a2,a0,_P_) @X(-2)(P,a0,a4)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              for(auto &&P : irange(0L, Nrhom2)){
                W1.cptr()[a5*_Nact*_Nact*Nrhom1*Nrhom2+a3*_Nact*Nrhom1*Nrhom2+a2*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[a5*_Nact*_Nact*Nrhom1*Nrhom2+a3*_Nact*Nrhom1*Nrhom2+a2*Nrhom1*Nrhom2+_P_*Nrhom2+P];
              }//P
            }//_P_
          }//a2
        }//a3
      }//a5
    }
    orz::DTensor W2(_Ndocc,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Ndocc,_Nact*_Nact*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              tmp1.cptr()[(c0)*tmp1.shape(1)+(a2*_Nact*_Nact+a3*_Nact+a5)] = V2_FiC_aaac.cptr()[a2*V2_FiC_aaac.shape(1)*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a3*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a5*V2_FiC_aaac.shape(3)+c0];
            }//a5
          }//a3
        }//a2
      }//c0
      orz::DTensor tmp2(_Nact*_Nact*_Nact,Nrhom1*Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a5 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              for(auto &&P : irange(0L, Nrhom2)){
                tmp2.cptr()[(a2*_Nact*_Nact+a3*_Nact+a5)*tmp2.shape(1)+(_P_*Nrhom2+P)] = W1.cptr()[a5*W1.shape(1)*W1.shape(2)*W1.shape(3)*W1.shape(4)+a3*W1.shape(2)*W1.shape(3)*W1.shape(4)+a2*W1.shape(3)*W1.shape(4)+_P_*W1.shape(4)+P];
              }//P
            }//_P_
          }//a5
        }//a3
      }//a2
      // @W2(c0,_P_,P) <<= +1 @V2_FiC_aaac(a2,a3,a5,c0) @W1(a5,a3,a2,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&P : irange(0L, Nrhom2)){
            W2.cptr()[c0*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[c0*Nrhom1*Nrhom2+_P_*Nrhom2+P];
          }//P
        }//_P_
      }//c0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom1*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*_Nvirt+a)*tmp1.shape(1)+(_P_*_Ndocc+c0)] = T_m1_.cptr()[b*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+a*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+c0];
            }//c0
          }//_P_
        }//a
      }//b
      orz::DTensor tmp2(Nrhom1*_Ndocc,Nrhom2);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(_P_*_Ndocc+c0)*tmp2.shape(1)+(P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)+_P_*W2.shape(2)+P];
          }//P
        }//c0
      }//_P_
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-1)(b,a,_P_,c0) @W2(c0,_P_,P)
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/2 @T(-1)(b,a,_P_,c0) @W2(c0,_P_,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*_Nact*_Nact,_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&a0 : irange(0L, _Nact)){
                for(auto &&a1 : irange(0L, _Nact)){
                  tmp1.cptr()[(a5*_Nact*_Nact*_Nact*_Nact+a4*_Nact*_Nact*_Nact+a3*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a1)] = d3.cptr()[a5*d3.shape(1)*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a4*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a3*d3.shape(3)*d3.shape(4)*d3.shape(5)+a2*d3.shape(4)*d3.shape(5)+a1*d3.shape(5)+a0];
                }//a1
              }//a0
            }//a2
          }//a3
        }//a4
      }//a5
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a5,a4,a3,a2,a0,_P_) <<= +1 @D3(a5,a4,a3,a2,a1,a0) @X(-1)(+)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&a0 : irange(0L, _Nact)){
                for(auto &&_P_ : irange(0L, Nrhom1)){
                  W0.cptr()[a5*_Nact*_Nact*_Nact*_Nact*Nrhom1+a4*_Nact*_Nact*_Nact*Nrhom1+a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a5*_Nact*_Nact*_Nact*_Nact*Nrhom1+a4*_Nact*_Nact*_Nact*Nrhom1+a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_];
                }//_P_
              }//a0
            }//a2
          }//a3
        }//a4
      }//a5
    }
    orz::DTensor W1(_Nact,_Nact,_Nact,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*Nrhom1,_Nact*_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              for(auto &&a4 : irange(0L, _Nact)){
                for(auto &&a0 : irange(0L, _Nact)){
                  tmp1.cptr()[(a5*_Nact*_Nact*Nrhom1+a3*_Nact*Nrhom1+a2*Nrhom1+_P_)*tmp1.shape(1)+(a4*_Nact+a0)] = W0.cptr()[a5*W0.shape(1)*W0.shape(2)*W0.shape(3)*W0.shape(4)*W0.shape(5)+a4*W0.shape(2)*W0.shape(3)*W0.shape(4)*W0.shape(5)+a3*W0.shape(3)*W0.shape(4)*W0.shape(5)+a2*W0.shape(4)*W0.shape(5)+a0*W0.shape(5)+_P_];
                }//a0
              }//a4
            }//_P_
          }//a2
        }//a3
      }//a5
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a4*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a4*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a4
      // @W1(a5,a3,a2,_P_,P) <<= +1 @W0(a5,a4,a3,a2,a0,_P_) @X(-2)(P,a4,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              for(auto &&P : irange(0L, Nrhom2)){
                W1.cptr()[a5*_Nact*_Nact*Nrhom1*Nrhom2+a3*_Nact*Nrhom1*Nrhom2+a2*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[a5*_Nact*_Nact*Nrhom1*Nrhom2+a3*_Nact*Nrhom1*Nrhom2+a2*Nrhom1*Nrhom2+_P_*Nrhom2+P];
              }//P
            }//_P_
          }//a2
        }//a3
      }//a5
    }
    orz::DTensor W2(_Ndocc,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Ndocc,_Nact*_Nact*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              tmp1.cptr()[(c0)*tmp1.shape(1)+(a2*_Nact*_Nact+a3*_Nact+a5)] = V2_FiC_aaac.cptr()[a2*V2_FiC_aaac.shape(1)*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a3*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a5*V2_FiC_aaac.shape(3)+c0];
            }//a5
          }//a3
        }//a2
      }//c0
      orz::DTensor tmp2(_Nact*_Nact*_Nact,Nrhom1*Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a5 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              for(auto &&P : irange(0L, Nrhom2)){
                tmp2.cptr()[(a2*_Nact*_Nact+a3*_Nact+a5)*tmp2.shape(1)+(_P_*Nrhom2+P)] = W1.cptr()[a5*W1.shape(1)*W1.shape(2)*W1.shape(3)*W1.shape(4)+a3*W1.shape(2)*W1.shape(3)*W1.shape(4)+a2*W1.shape(3)*W1.shape(4)+_P_*W1.shape(4)+P];
              }//P
            }//_P_
          }//a5
        }//a3
      }//a2
      // @W2(c0,_P_,P) <<= +1 @V2_FiC_aaac(a2,a3,a5,c0) @W1(a5,a3,a2,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&P : irange(0L, Nrhom2)){
            W2.cptr()[c0*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[c0*Nrhom1*Nrhom2+_P_*Nrhom2+P];
          }//P
        }//_P_
      }//c0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom1*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*_Nvirt+a)*tmp1.shape(1)+(_P_*_Ndocc+c0)] = T_m1_.cptr()[b*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+a*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+c0];
            }//c0
          }//_P_
        }//a
      }//b
      orz::DTensor tmp2(Nrhom1*_Ndocc,Nrhom2);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(_P_*_Ndocc+c0)*tmp2.shape(1)+(P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)+_P_*W2.shape(2)+P];
          }//P
        }//c0
      }//_P_
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-1)(b,a,_P_,c0) @W2(c0,_P_,P)
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/2 @T(-1)(a,b,_P_,c0) @W2(c0,_P_,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*_Nact*_Nact,_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&a0 : irange(0L, _Nact)){
                for(auto &&a1 : irange(0L, _Nact)){
                  tmp1.cptr()[(a5*_Nact*_Nact*_Nact*_Nact+a4*_Nact*_Nact*_Nact+a3*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a1)] = d3.cptr()[a5*d3.shape(1)*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a4*d3.shape(2)*d3.shape(3)*d3.shape(4)*d3.shape(5)+a3*d3.shape(3)*d3.shape(4)*d3.shape(5)+a2*d3.shape(4)*d3.shape(5)+a1*d3.shape(5)+a0];
                }//a1
              }//a0
            }//a2
          }//a3
        }//a4
      }//a5
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a5,a4,a3,a2,a0,_P_) <<= +1 @D3(a5,a4,a3,a2,a1,a0) @X(-1)(-)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&a0 : irange(0L, _Nact)){
                for(auto &&_P_ : irange(0L, Nrhom1)){
                  W0.cptr()[a5*_Nact*_Nact*_Nact*_Nact*Nrhom1+a4*_Nact*_Nact*_Nact*Nrhom1+a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a5*_Nact*_Nact*_Nact*_Nact*Nrhom1+a4*_Nact*_Nact*_Nact*Nrhom1+a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_];
                }//_P_
              }//a0
            }//a2
          }//a3
        }//a4
      }//a5
    }
    orz::DTensor W1(_Nact,_Nact,_Nact,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*Nrhom1,_Nact*_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              for(auto &&a4 : irange(0L, _Nact)){
                for(auto &&a0 : irange(0L, _Nact)){
                  tmp1.cptr()[(a5*_Nact*_Nact*Nrhom1+a3*_Nact*Nrhom1+a2*Nrhom1+_P_)*tmp1.shape(1)+(a4*_Nact+a0)] = W0.cptr()[a5*W0.shape(1)*W0.shape(2)*W0.shape(3)*W0.shape(4)*W0.shape(5)+a4*W0.shape(2)*W0.shape(3)*W0.shape(4)*W0.shape(5)+a3*W0.shape(3)*W0.shape(4)*W0.shape(5)+a2*W0.shape(4)*W0.shape(5)+a0*W0.shape(5)+_P_];
                }//a0
              }//a4
            }//_P_
          }//a2
        }//a3
      }//a5
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a4 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a4*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a4*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a4
      // @W1(a5,a3,a2,_P_,P) <<= +1 @W0(a5,a4,a3,a2,a0,_P_) @X(-2)(P,a4,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              for(auto &&P : irange(0L, Nrhom2)){
                W1.cptr()[a5*_Nact*_Nact*Nrhom1*Nrhom2+a3*_Nact*Nrhom1*Nrhom2+a2*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[a5*_Nact*_Nact*Nrhom1*Nrhom2+a3*_Nact*Nrhom1*Nrhom2+a2*Nrhom1*Nrhom2+_P_*Nrhom2+P];
              }//P
            }//_P_
          }//a2
        }//a3
      }//a5
    }
    orz::DTensor W2(_Ndocc,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Ndocc,_Nact*_Nact*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              tmp1.cptr()[(c0)*tmp1.shape(1)+(a2*_Nact*_Nact+a3*_Nact+a5)] = V2_FiC_aaac.cptr()[a2*V2_FiC_aaac.shape(1)*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a3*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a5*V2_FiC_aaac.shape(3)+c0];
            }//a5
          }//a3
        }//a2
      }//c0
      orz::DTensor tmp2(_Nact*_Nact*_Nact,Nrhom1*Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a5 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              for(auto &&P : irange(0L, Nrhom2)){
                tmp2.cptr()[(a2*_Nact*_Nact+a3*_Nact+a5)*tmp2.shape(1)+(_P_*Nrhom2+P)] = W1.cptr()[a5*W1.shape(1)*W1.shape(2)*W1.shape(3)*W1.shape(4)+a3*W1.shape(2)*W1.shape(3)*W1.shape(4)+a2*W1.shape(3)*W1.shape(4)+_P_*W1.shape(4)+P];
              }//P
            }//_P_
          }//a5
        }//a3
      }//a2
      // @W2(c0,_P_,P) <<= +1 @V2_FiC_aaac(a2,a3,a5,c0) @W1(a5,a3,a2,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&P : irange(0L, Nrhom2)){
            W2.cptr()[c0*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[c0*Nrhom1*Nrhom2+_P_*Nrhom2+P];
          }//P
        }//_P_
      }//c0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom1*_Ndocc);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(a*_Nvirt+b)*tmp1.shape(1)+(_P_*_Ndocc+c0)] = T_m1_.cptr()[a*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+b*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+c0];
            }//c0
          }//_P_
        }//b
      }//a
      orz::DTensor tmp2(Nrhom1*_Ndocc,Nrhom2);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(_P_*_Ndocc+c0)*tmp2.shape(1)+(P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)+_P_*W2.shape(2)+P];
          }//P
        }//c0
      }//_P_
      // @tRFiC(a,b,P) <<= +1 _X_ @T(-1)(a,b,_P_,c0) @W2(c0,_P_,P)
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/2 @W1(a3,_P_,P) @W2(a3,a,b,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a0
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a3,a2,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_];
            }//_P_
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*Nrhom1,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*Nrhom1+_P_)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a0
          }//a2
        }//_P_
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a2*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a2
      // @W1(a3,_P_,P) <<= +1 @W0(a3,a2,a0,_P_) @X(-2)(P,a2,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P];
          }//P
        }//_P_
      }//a3
    }
    orz::DTensor W2(_Nact,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*Nrhom1,_Nvirt*_Ndocc);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(a*Nrhom1+_P_)*tmp1.shape(1)+(v0*_Ndocc+c0)] = T_m1_.cptr()[v0*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+a*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+c0];
            }//c0
          }//v0
        }//_P_
      }//a
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(b*_Nact+a3)] = V2_FiC_vvac.cptr()[b*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a3*V2_FiC_vvac.shape(3)+c0];
            }//a3
          }//b
        }//c0
      }//v0
      // @W2(a3,a,b,_P_) <<= +1 @T(-1)(v0,a,_P_,c0) @V2_FiC_vvac(b,v0,a3,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              W2.cptr()[a3*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+_P_] += tmp3.cptr()[a*Nrhom1*_Nvirt*_Nact+_P_*_Nvirt*_Nact+b*_Nact+a3];
            }//a3
          }//b
        }//_P_
      }//a
    }
    {
      orz::DTensor tmp1(Nrhom2,_Nact*Nrhom1);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a3*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)+_P_*W1.shape(2)+P];
          }//_P_
        }//a3
      }//P
      orz::DTensor tmp2(_Nact*Nrhom1,_Nvirt*_Nvirt);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(a3*Nrhom1+_P_)*tmp2.shape(1)+(a*_Nvirt+b)] = W2.cptr()[a3*W2.shape(1)*W2.shape(2)*W2.shape(3)+a*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+_P_];
            }//b
          }//a
        }//_P_
      }//a3
      // @tRFiC(a,b,P) <<= +1 _X_ @W1(a3,_P_,P) @W2(a3,a,b,_P_)
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/2 @W1(a3,_P_,P) @W2(a3,a,b,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a0
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a3,a2,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(-)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_];
            }//_P_
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*Nrhom1,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*Nrhom1+_P_)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a0
          }//a2
        }//_P_
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a2*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a2
      // @W1(a3,_P_,P) <<= +1 @W0(a3,a2,a0,_P_) @X(-2)(P,a2,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P];
          }//P
        }//_P_
      }//a3
    }
    orz::DTensor W2(_Nact,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*Nrhom1,_Nvirt*_Ndocc);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(a*Nrhom1+_P_)*tmp1.shape(1)+(v0*_Ndocc+c0)] = T_m1_.cptr()[a*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+v0*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+c0];
            }//c0
          }//v0
        }//_P_
      }//a
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(b*_Nact+a3)] = V2_FiC_vvac.cptr()[b*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a3*V2_FiC_vvac.shape(3)+c0];
            }//a3
          }//b
        }//c0
      }//v0
      // @W2(a3,a,b,_P_) <<= +1 @T(-1)(a,v0,_P_,c0) @V2_FiC_vvac(b,v0,a3,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              W2.cptr()[a3*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+_P_] += tmp3.cptr()[a*Nrhom1*_Nvirt*_Nact+_P_*_Nvirt*_Nact+b*_Nact+a3];
            }//a3
          }//b
        }//_P_
      }//a
    }
    {
      orz::DTensor tmp1(Nrhom2,_Nact*Nrhom1);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a3*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)+_P_*W1.shape(2)+P];
          }//_P_
        }//a3
      }//P
      orz::DTensor tmp2(_Nact*Nrhom1,_Nvirt*_Nvirt);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(a3*Nrhom1+_P_)*tmp2.shape(1)+(a*_Nvirt+b)] = W2.cptr()[a3*W2.shape(1)*W2.shape(2)*W2.shape(3)+a*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+_P_];
            }//b
          }//a
        }//_P_
      }//a3
      // @tRFiC(a,b,P) <<= +1 _X_ @W1(a3,_P_,P) @W2(a3,a,b,_P_)
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/2 @W1(a3,_P_,P) @W2(a3,b,a,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a0
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a3,a2,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_];
            }//_P_
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*Nrhom1,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*Nrhom1+_P_)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a0
          }//a2
        }//_P_
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a2
      // @W1(a3,_P_,P) <<= +1 @W0(a3,a2,a0,_P_) @X(-2)(P,a0,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P];
          }//P
        }//_P_
      }//a3
    }
    orz::DTensor W2(_Nact,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*Nrhom1,_Nvirt*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*Nrhom1+_P_)*tmp1.shape(1)+(v0*_Ndocc+c0)] = T_m1_.cptr()[v0*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+b*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+c0];
            }//c0
          }//v0
        }//_P_
      }//b
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(a*_Nact+a3)] = V2_FiC_vvac.cptr()[a*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a3*V2_FiC_vvac.shape(3)+c0];
            }//a3
          }//a
        }//c0
      }//v0
      // @W2(a3,b,a,_P_) <<= +1 @T(-1)(v0,b,_P_,c0) @V2_FiC_vvac(a,v0,a3,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              W2.cptr()[a3*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+_P_] += tmp3.cptr()[b*Nrhom1*_Nvirt*_Nact+_P_*_Nvirt*_Nact+a*_Nact+a3];
            }//a3
          }//a
        }//_P_
      }//b
    }
    {
      orz::DTensor tmp1(Nrhom2,_Nact*Nrhom1);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a3*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)+_P_*W1.shape(2)+P];
          }//_P_
        }//a3
      }//P
      orz::DTensor tmp2(_Nact*Nrhom1,_Nvirt*_Nvirt);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(a3*Nrhom1+_P_)*tmp2.shape(1)+(b*_Nvirt+a)] = W2.cptr()[a3*W2.shape(1)*W2.shape(2)*W2.shape(3)+b*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+_P_];
            }//a
          }//b
        }//_P_
      }//a3
      // @tRFiC(a,b,P) <<= +1 _X_ @W1(a3,_P_,P) @W2(a3,b,a,_P_)
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/2 @W1(a3,_P_,P) @W2(a3,b,a,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a0
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a3,a2,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(-)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_];
            }//_P_
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*Nrhom1,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*Nrhom1+_P_)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a0
          }//a2
        }//_P_
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a2
      // @W1(a3,_P_,P) <<= +1 @W0(a3,a2,a0,_P_) @X(-2)(P,a0,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P];
          }//P
        }//_P_
      }//a3
    }
    orz::DTensor W2(_Nact,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*Nrhom1,_Nvirt*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*Nrhom1+_P_)*tmp1.shape(1)+(v0*_Ndocc+c0)] = T_m1_.cptr()[b*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+v0*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+c0];
            }//c0
          }//v0
        }//_P_
      }//b
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(a*_Nact+a3)] = V2_FiC_vvac.cptr()[a*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a3*V2_FiC_vvac.shape(3)+c0];
            }//a3
          }//a
        }//c0
      }//v0
      // @W2(a3,b,a,_P_) <<= +1 @T(-1)(b,v0,_P_,c0) @V2_FiC_vvac(a,v0,a3,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              W2.cptr()[a3*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+_P_] += tmp3.cptr()[b*Nrhom1*_Nvirt*_Nact+_P_*_Nvirt*_Nact+a*_Nact+a3];
            }//a3
          }//a
        }//_P_
      }//b
    }
    {
      orz::DTensor tmp1(Nrhom2,_Nact*Nrhom1);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a3*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)+_P_*W1.shape(2)+P];
          }//_P_
        }//a3
      }//P
      orz::DTensor tmp2(_Nact*Nrhom1,_Nvirt*_Nvirt);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(a3*Nrhom1+_P_)*tmp2.shape(1)+(b*_Nvirt+a)] = W2.cptr()[a3*W2.shape(1)*W2.shape(2)*W2.shape(3)+b*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+_P_];
            }//a
          }//b
        }//_P_
      }//a3
      // @tRFiC(a,b,P) <<= +1 _X_ @W1(a3,_P_,P) @W2(a3,b,a,_P_)
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/2 @W1(a3,_P_,P) @W2(a3,a,b,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a0
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a3,a2,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_];
            }//_P_
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*Nrhom1,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*Nrhom1+_P_)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a0
          }//a2
        }//_P_
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a2
      // @W1(a3,_P_,P) <<= +1 @W0(a3,a2,a0,_P_) @X(-2)(P,a0,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P];
          }//P
        }//_P_
      }//a3
    }
    orz::DTensor W2(_Nact,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*Nrhom1,_Nvirt*_Ndocc);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(a*Nrhom1+_P_)*tmp1.shape(1)+(v0*_Ndocc+c0)] = T_m1_.cptr()[a*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+v0*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+c0];
            }//c0
          }//v0
        }//_P_
      }//a
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(b*_Nact+a3)] = V2_FiC_vvac.cptr()[b*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a3*V2_FiC_vvac.shape(3)+c0];
            }//a3
          }//b
        }//c0
      }//v0
      // @W2(a3,a,b,_P_) <<= +1 @T(-1)(a,v0,_P_,c0) @V2_FiC_vvac(b,v0,a3,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              W2.cptr()[a3*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+_P_] += tmp3.cptr()[a*Nrhom1*_Nvirt*_Nact+_P_*_Nvirt*_Nact+b*_Nact+a3];
            }//a3
          }//b
        }//_P_
      }//a
    }
    {
      orz::DTensor tmp1(Nrhom2,_Nact*Nrhom1);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a3*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)+_P_*W1.shape(2)+P];
          }//_P_
        }//a3
      }//P
      orz::DTensor tmp2(_Nact*Nrhom1,_Nvirt*_Nvirt);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(a3*Nrhom1+_P_)*tmp2.shape(1)+(a*_Nvirt+b)] = W2.cptr()[a3*W2.shape(1)*W2.shape(2)*W2.shape(3)+a*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+_P_];
            }//b
          }//a
        }//_P_
      }//a3
      // @tRFiC(a,b,P) <<= +1 _X_ @W1(a3,_P_,P) @W2(a3,a,b,_P_)
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/2 @W1(a3,_P_,P) @W2(a3,a,b,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a0
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a3,a2,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(-)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_];
            }//_P_
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*Nrhom1,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*Nrhom1+_P_)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a0
          }//a2
        }//_P_
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a2
      // @W1(a3,_P_,P) <<= +1 @W0(a3,a2,a0,_P_) @X(-2)(P,a0,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P];
          }//P
        }//_P_
      }//a3
    }
    orz::DTensor W2(_Nact,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*Nrhom1,_Nvirt*_Ndocc);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(a*Nrhom1+_P_)*tmp1.shape(1)+(v0*_Ndocc+c0)] = T_m1_.cptr()[v0*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+a*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+c0];
            }//c0
          }//v0
        }//_P_
      }//a
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(b*_Nact+a3)] = V2_FiC_vvac.cptr()[b*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a3*V2_FiC_vvac.shape(3)+c0];
            }//a3
          }//b
        }//c0
      }//v0
      // @W2(a3,a,b,_P_) <<= +1 @T(-1)(v0,a,_P_,c0) @V2_FiC_vvac(b,v0,a3,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              W2.cptr()[a3*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+_P_] += tmp3.cptr()[a*Nrhom1*_Nvirt*_Nact+_P_*_Nvirt*_Nact+b*_Nact+a3];
            }//a3
          }//b
        }//_P_
      }//a
    }
    {
      orz::DTensor tmp1(Nrhom2,_Nact*Nrhom1);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a3*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)+_P_*W1.shape(2)+P];
          }//_P_
        }//a3
      }//P
      orz::DTensor tmp2(_Nact*Nrhom1,_Nvirt*_Nvirt);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(a3*Nrhom1+_P_)*tmp2.shape(1)+(a*_Nvirt+b)] = W2.cptr()[a3*W2.shape(1)*W2.shape(2)*W2.shape(3)+a*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+_P_];
            }//b
          }//a
        }//_P_
      }//a3
      // @tRFiC(a,b,P) <<= +1 _X_ @W1(a3,_P_,P) @W2(a3,a,b,_P_)
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/2 @W1(a3,_P_,P) @W2(a3,b,a,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a0
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a3,a2,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_];
            }//_P_
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*Nrhom1,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*Nrhom1+_P_)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a0
          }//a2
        }//_P_
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a2*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a2
      // @W1(a3,_P_,P) <<= +1 @W0(a3,a2,a0,_P_) @X(-2)(P,a2,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P];
          }//P
        }//_P_
      }//a3
    }
    orz::DTensor W2(_Nact,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*Nrhom1,_Nvirt*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*Nrhom1+_P_)*tmp1.shape(1)+(v0*_Ndocc+c0)] = T_m1_.cptr()[b*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+v0*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+c0];
            }//c0
          }//v0
        }//_P_
      }//b
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(a*_Nact+a3)] = V2_FiC_vvac.cptr()[a*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a3*V2_FiC_vvac.shape(3)+c0];
            }//a3
          }//a
        }//c0
      }//v0
      // @W2(a3,b,a,_P_) <<= +1 @T(-1)(b,v0,_P_,c0) @V2_FiC_vvac(a,v0,a3,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              W2.cptr()[a3*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+_P_] += tmp3.cptr()[b*Nrhom1*_Nvirt*_Nact+_P_*_Nvirt*_Nact+a*_Nact+a3];
            }//a3
          }//a
        }//_P_
      }//b
    }
    {
      orz::DTensor tmp1(Nrhom2,_Nact*Nrhom1);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a3*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)+_P_*W1.shape(2)+P];
          }//_P_
        }//a3
      }//P
      orz::DTensor tmp2(_Nact*Nrhom1,_Nvirt*_Nvirt);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(a3*Nrhom1+_P_)*tmp2.shape(1)+(b*_Nvirt+a)] = W2.cptr()[a3*W2.shape(1)*W2.shape(2)*W2.shape(3)+b*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+_P_];
            }//a
          }//b
        }//_P_
      }//a3
      // @tRFiC(a,b,P) <<= +1 _X_ @W1(a3,_P_,P) @W2(a3,b,a,_P_)
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P) <<= -1/2 @W1(a3,_P_,P) @W2(a3,b,a,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -0.5000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a2*_Nact+a0)*tmp1.shape(1)+(a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a0
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a3,a2,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(-)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+_P_];
            }//_P_
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,Nrhom1,Nrhom2);
    {
      orz::DTensor tmp1(_Nact*Nrhom1,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*Nrhom1+_P_)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a0
          }//a2
        }//_P_
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a2*Cm2.shape(1)*Cm2.shape(2)+a0*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a2
      // @W1(a3,_P_,P) <<= +1 @W0(a3,a2,a0,_P_) @X(-2)(P,a2,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P] += tmp3.cptr()[a3*Nrhom1*Nrhom2+_P_*Nrhom2+P];
          }//P
        }//_P_
      }//a3
    }
    orz::DTensor W2(_Nact,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*Nrhom1,_Nvirt*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*Nrhom1+_P_)*tmp1.shape(1)+(v0*_Ndocc+c0)] = T_m1_.cptr()[v0*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+b*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+c0];
            }//c0
          }//v0
        }//_P_
      }//b
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(a*_Nact+a3)] = V2_FiC_vvac.cptr()[a*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a3*V2_FiC_vvac.shape(3)+c0];
            }//a3
          }//a
        }//c0
      }//v0
      // @W2(a3,b,a,_P_) <<= +1 @T(-1)(v0,b,_P_,c0) @V2_FiC_vvac(a,v0,a3,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a3 : irange(0L, _Nact)){
              W2.cptr()[a3*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+_P_] += tmp3.cptr()[b*Nrhom1*_Nvirt*_Nact+_P_*_Nvirt*_Nact+a*_Nact+a3];
            }//a3
          }//a
        }//_P_
      }//b
    }
    {
      orz::DTensor tmp1(Nrhom2,_Nact*Nrhom1);
      for(auto &&P : irange(0L, Nrhom2)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            tmp1.cptr()[(P)*tmp1.shape(1)+(a3*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)+_P_*W1.shape(2)+P];
          }//_P_
        }//a3
      }//P
      orz::DTensor tmp2(_Nact*Nrhom1,_Nvirt*_Nvirt);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(a3*Nrhom1+_P_)*tmp2.shape(1)+(b*_Nvirt+a)] = W2.cptr()[a3*W2.shape(1)*W2.shape(2)*W2.shape(3)+b*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+_P_];
            }//a
          }//b
        }//_P_
      }//a3
      // @tRFiC(a,b,P) <<= +1 _X_ @W1(a3,_P_,P) @W2(a3,b,a,_P_)
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
        RActList.GetPair(P).Transform2EXT(RP);
        orz::DTensor Rorig = R2.GetTensor(P);
        RP += Rorig;
        R2.PutTensor(P, RP);
      }
      
      return;
    }//End func
  } }
