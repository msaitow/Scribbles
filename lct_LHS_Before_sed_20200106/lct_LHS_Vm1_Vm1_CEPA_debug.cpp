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
    void FMatrixOrthogonalizer::FormLHS_Vm1_Vm1_CEd(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &c2, const orz::DTensor &d3, DensityPack &dPack,
                                                    TensorContainer &ovl, TensorContainer &PNO4, const bool do4EXTinPNO,
                                                    const double QE0, const double Ecas, const double h0_energy, TensorContainer &DebugInts, const orz::DTensor &casfock, const orz::DTensor &FV,
                                                    TensorContainer &T2, TensorContainer &R2){

      orz::ProgressTimer pt("* V(-1)/V(-1)");
      
#include "lct_preparations_orthdens.cpp"        
      // Init amplitude container                               
      for(auto &&P : irange(0L, Nrhom1)){                       
        for(auto &&i : irange(0L, _Ndocc)){
          const long I_Pi = TActList.GetIndex(P,i);
          if (I_Pi < 0) continue;                      
          const long idx_Pi = RActList.GetIndex(P,i);           
          const long Nvirt_Pi = RActList.GetPair(P,i).GetNPNO();
          orz::DTensor Rpi(Nvirt_Pi, Nvirt_Pi);                 
          R2.PutTensor(P,i,Rpi);                                
        }                                                       
      }
      
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
  orz::DTensor tRFiC(_Nvirt,_Nvirt,Nrhom1,_Ndocc);
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 QE0 @Tc(-1)(a,b,_P_,i) @W1(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000 * QE0;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Nvirt*_Ndocc+b*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//b
      }//a
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(a,b,_P_,i) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Nvirt*_Ndocc*Nrhom1+b*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 QE0 @Tc(-1)(b,a,_P_,i) @W1(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000 * QE0;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Nvirt*_Ndocc+a*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//a
      }//b
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(b,a,_P_,i) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Nvirt*_Ndocc*Nrhom1+a*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 QE0 @Tc(-1)(b,a,_P_,i) @W1(_P_,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000 * QE0;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a0,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+_P_] += tmp3.cptr()[a0*Nrhom1+_P_];
        }//_P_
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(_P_)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+_P_];
        }//a0
      }//_P_
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @W1(_P_,P) <<= +1 @W0(a0,_P_) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          W1.cptr()[_P_*Nrhom1+P] += tmp3.cptr()[_P_*Nrhom1+P];
        }//P
      }//_P_
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Nvirt*_Ndocc+a*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//a
      }//b
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[_P_*W1.shape(1)+P];
        }//P
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(b,a,_P_,i) @W1(_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Nvirt*_Ndocc*Nrhom1+a*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 QE0 @Tc(-1)(a,b,_P_,i) @W1(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000 * QE0;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Nvirt*_Ndocc+b*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//b
      }//a
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(a,b,_P_,i) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Nvirt*_Ndocc*Nrhom1+b*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
        }//b
      }//a
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra

  if(!do4EXTinPNO){  
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @V2_FiC_vvvv(b,v0,a,v1) @W2(i,v1,v0,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&v1 : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(v1*_Nvirt*_Ndocc+v0*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[v1*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//v0
      }//v1
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @W2(i,v1,v0,P) <<= +1 @Tc(-1)(v1,v0,_P_,i) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v1 : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[i*_Nvirt*_Nvirt*Nrhom1+v1*_Nvirt*Nrhom1+v0*Nrhom1+P] += tmp3.cptr()[v1*_Nvirt*_Ndocc*Nrhom1+v0*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
        }//v0
      }//v1
    }
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
      orz::DTensor tmp2(_Nvirt*_Nvirt,_Ndocc*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&v1 : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*_Nvirt+v1)*tmp2.shape(1)+(i*Nrhom1+P)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+v1*W2.shape(2)*W2.shape(3)+v0*W2.shape(3)+P];
            }//P
          }//i
        }//v1
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @V2_FiC_vvvv(b,v0,a,v1) @W2(i,v1,v0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Nvirt*_Ndocc*Nrhom1+a*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @V2_FiC_vvvv(a,v0,b,v1) @W2(i,v1,v0,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&v1 : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(v1*_Nvirt*_Ndocc+v0*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[v1*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//v0
      }//v1
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @W2(i,v1,v0,P) <<= +1 @Tc(-1)(v1,v0,_P_,i) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v1 : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[i*_Nvirt*_Nvirt*Nrhom1+v1*_Nvirt*Nrhom1+v0*Nrhom1+P] += tmp3.cptr()[v1*_Nvirt*_Ndocc*Nrhom1+v0*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
        }//v0
      }//v1
    }
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
      orz::DTensor tmp2(_Nvirt*_Nvirt,_Ndocc*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&v1 : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*_Nvirt+v1)*tmp2.shape(1)+(i*Nrhom1+P)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+v1*W2.shape(2)*W2.shape(3)+v0*W2.shape(3)+P];
            }//P
          }//i
        }//v1
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @V2_FiC_vvvv(a,v0,b,v1) @W2(i,v1,v0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Nvirt*_Ndocc*Nrhom1+b*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @V2_FiC_vvvv(a,v0,b,v1) @W2(i,v1,v0,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a0,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+_P_] += tmp3.cptr()[a0*Nrhom1+_P_];
        }//_P_
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(_P_)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+_P_];
        }//a0
      }//_P_
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @W1(_P_,P) <<= +1 @W0(a0,_P_) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          W1.cptr()[_P_*Nrhom1+P] += tmp3.cptr()[_P_*Nrhom1+P];
        }//P
      }//_P_
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&v1 : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(v1*_Nvirt*_Ndocc+v0*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[v1*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//v0
      }//v1
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[_P_*W1.shape(1)+P];
        }//P
      }//_P_
      // @W2(i,v1,v0,P) <<= +1 @Tc(-1)(v1,v0,_P_,i) @W1(_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v1 : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[i*_Nvirt*_Nvirt*Nrhom1+v1*_Nvirt*Nrhom1+v0*Nrhom1+P] += tmp3.cptr()[v1*_Nvirt*_Ndocc*Nrhom1+v0*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
        }//v0
      }//v1
    }
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
      orz::DTensor tmp2(_Nvirt*_Nvirt,_Ndocc*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&v1 : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*_Nvirt+v1)*tmp2.shape(1)+(i*Nrhom1+P)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+v1*W2.shape(2)*W2.shape(3)+v0*W2.shape(3)+P];
            }//P
          }//i
        }//v1
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @V2_FiC_vvvv(a,v0,b,v1) @W2(i,v1,v0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Nvirt*_Ndocc*Nrhom1+b*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @V2_FiC_vvvv(b,v0,a,v1) @W2(i,v1,v0,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&v1 : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(v1*_Nvirt*_Ndocc+v0*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[v1*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//v0
      }//v1
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @W2(i,v1,v0,P) <<= +1 @Tc(-1)(v1,v0,_P_,i) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v1 : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[i*_Nvirt*_Nvirt*Nrhom1+v1*_Nvirt*Nrhom1+v0*Nrhom1+P] += tmp3.cptr()[v1*_Nvirt*_Ndocc*Nrhom1+v0*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
        }//v0
      }//v1
    }
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
      orz::DTensor tmp2(_Nvirt*_Nvirt,_Ndocc*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&v1 : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*_Nvirt+v1)*tmp2.shape(1)+(i*Nrhom1+P)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+v1*W2.shape(2)*W2.shape(3)+v0*W2.shape(3)+P];
            }//P
          }//i
        }//v1
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @V2_FiC_vvvv(b,v0,a,v1) @W2(i,v1,v0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Nvirt*_Ndocc*Nrhom1+a*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
        }//a
      }//b
    }

#ifdef _PROFILE_LHS
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif

  }//End Contra
  }
  if(do4EXTinPNO){
    orz::ProgressTimer pt("* V(-1)/V(-1) - 4ext (PNO) ");
    orz::DTensor X1(Nrhom1,Nrhom1); // (P,_P_)
    {      
      orz::DTensor X0; // (p,_P_)
      A_x_B(X0,d1,Cm1_plus,false,true); 
      orz::DTensor tmp(Cm1_plus);
      tmp.gaxpy(+2.0,Cm1_minus,-1.0);
      X1 = (tmp*X0);
    }
    orz::DTensor X2(Nrhom1,Nrhom1); // (P,_P_)
    {      
      orz::DTensor X0; // (p,_P_)
      A_x_B(X0,d1,Cm1_plus,false,true); 
      orz::DTensor tmp(Cm1_minus);
      tmp.gaxpy(+2.0,Cm1_plus,-1.0);
      X2 = (tmp*X0);
    }
    orz::DTensor X3(Nrhom1,Nrhom1); // (P,_P_)
    {      
      orz::DTensor X0; // (p,_P_)
      A_x_B(X0,d1,Cm1_minus,false,true); 
      orz::DTensor tmp(Cm1_plus);
      tmp.gaxpy(+2.0,Cm1_minus,-1.0);
      X3 = (tmp*X0);
    }
    orz::DTensor X4(Nrhom1,Nrhom1); // (P,_P_)
    {      
      orz::DTensor X0; // (p,_P_)
      A_x_B(X0,d1,Cm1_minus,false,true); 
      orz::DTensor tmp(Cm1_minus);
      tmp.gaxpy(+2.0,Cm1_plus,-1.0);
      X4 = (tmp*X0);
    }    
    for(auto &&P : irange(0L, Nrhom1)){
      for(auto &&i : irange(0L, _Ndocc)){
        const long I_Pi = TActList.GetIndex(P, i);
        if (I_Pi < 0) continue;
        PairData P_Pi = TActList.GetPair(P,i);
        const long Nvir = TActList.GetNPNO(P, i);
        orz::DTensor Up(Nvir,Nvir), Uq(Nvir,Nvir), Ur(Nvir,Nvir), Us(Nvir,Nvir);
        for(auto &&P0 : irange(0L, Nrhom1)){
          const long I_P0i = TActList.GetIndex(P0, i);
          if (I_P0i < 0) continue;
          orz::DTensor U2 = T2.GetTensor(P0, i);
          {
            orz::DTensor S = ovl.GetTensor(std::max(I_Pi,I_P0i),std::min(I_Pi,I_P0i));
            if(I_P0i>I_Pi) orz::tensor::transpose(S);
            orz::DTensor tmp;
            A_x_B_x_At(tmp, S, U2);
            U2 = tmp.copy();
          }
          Up.gaxpy(+1.0, U2, X1(P,P0));
          Uq.gaxpy(+1.0, U2, X2(P,P0));
          Ur.gaxpy(+1.0, U2, X3(P,P0));
          Us.gaxpy(+1.0, U2, X4(P,P0));
        }//_P_
        orz::DTensor I4;
        PNO4.GetTensor(P,i,I4);
        // Sort Utilde
        orz::DTensor tmp2(1,Nvir*Nvir);
        for(auto &&a : irange(0L, Nvir))        
          for(auto &&b : irange(0L, Nvir))
            tmp2.cptr()[a*Nvir+b]  = Up.cptr()[a*Up.shape(1)+b];
        orz::DTensor tmp2a(1,Nvir*Nvir);
        for(auto &&a : irange(0L, Nvir))        
          for(auto &&b : irange(0L, Nvir))
            tmp2a.cptr()[a*Nvir+b] = Uq.cptr()[a*Up.shape(1)+b];
        orz::DTensor tmp2b(1,Nvir*Nvir);
        for(auto &&a : irange(0L, Nvir))        
          for(auto &&b : irange(0L, Nvir))
            tmp2b.cptr()[a*Nvir+b] = Ur.cptr()[a*Ur.shape(1)+b];
        orz::DTensor tmp2c(1,Nvir*Nvir);
        for(auto &&a : irange(0L, Nvir))        
          for(auto &&b : irange(0L, Nvir))
            tmp2c.cptr()[a*Nvir+b] = Us.cptr()[a*Us.shape(1)+b];                        
        orz::DTensor tmp1(Nvir*Nvir,Nvir*Nvir);
        
        // term p1 + term m1
        {
          // Sort I4        
//          for(auto &&a : irange(0L, Nvir))        
//            for(auto &&b : irange(0L, Nvir))
//              for(auto &&c : irange(0L, Nvir))
//                for(auto &&d : irange(0L, Nvir))
//                  tmp1.cptr()[b*Nvir*Nvir*Nvir+d*Nvir*Nvir+a*Nvir+c] = I4.cptr()[a*Nvir*Nvir*Nvir+b*Nvir*Nvir+c*Nvir+d];
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
          
//BT          orz::DTensor tmp3(tmp2*tmp1);
//BT          tmp3.reshape_inplace(Nvir,Nvir);
//BT          P_Pi.BackTransform2EXT(tmp3);
//BT
//BT          orz::DTensor tmp3a(tmp2a*tmp1);
//BT          tmp3a.reshape_inplace(Nvir,Nvir);
//BT          P_Pi.BackTransform2EXT(tmp3a);
//BT
//BT          orz::DTensor tmp3b(tmp2b*tmp1);
//BT          tmp3b.reshape_inplace(Nvir,Nvir);
//BT          P_Pi.BackTransform2EXT(tmp3b);
//BT
//BT          orz::DTensor tmp3c(tmp2c*tmp1);
//BT          tmp3c.reshape_inplace(Nvir,Nvir);
//BT          P_Pi.BackTransform2EXT(tmp3c);
//BT          
//BT          for(auto &&a : irange(0L, _Nvirt))
//BT            for(auto &&b : irange(0L, _Nvirt))    
//BT              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3 .cptr()[a*tmp3 .shape(1)+b]
//BT                                                                               + tmp3a.cptr()[b*tmp3a.shape(1)+a]
//BT                                                                               + tmp3b.cptr()[b*tmp3b.shape(1)+a]
//BT                                                                               + tmp3c.cptr()[a*tmp3c.shape(1)+b];

          orz::DTensor tmp3(tmp2*tmp1);
          tmp3.reshape_inplace(Nvir,Nvir);

          orz::DTensor tmp3a(tmp2a*tmp1);
          tmp3a.reshape_inplace(Nvir,Nvir);

          orz::DTensor tmp3b(tmp2b*tmp1);
          tmp3b.reshape_inplace(Nvir,Nvir);

          orz::DTensor tmp3c(tmp2c*tmp1);
          tmp3c.reshape_inplace(Nvir,Nvir);

          orz::DTensor S2 = R2.GetTensor(P,i);
          for(auto &&a : irange(0L, Nvir))
            for(auto &&b : irange(0L, Nvir))    
              S2.cptr()[a*Nvir+b] += tmp3 .cptr()[a*tmp3 .shape(1)+b]
                                  +  tmp3a.cptr()[b*tmp3a.shape(1)+a]
                                  +  tmp3b.cptr()[b*tmp3b.shape(1)+a]
                                  +  tmp3c.cptr()[a*tmp3c.shape(1)+b];
          
          R2.PutTensor(P,i,S2);
        }
//work        // term p2 + term m2
//work        {
//work          // Sort I4        
//work          for(auto &&a : irange(0L, Nvir))        
//work            for(auto &&b : irange(0L, Nvir))
//work              for(auto &&c : irange(0L, Nvir))
//work                for(auto &&d : irange(0L, Nvir))
//work                  tmp1.cptr()[b*Nvir*Nvir*Nvir+d*Nvir*Nvir+c*Nvir+a] = I4.cptr()[a*Nvir*Nvir*Nvir+b*Nvir*Nvir+c*Nvir+d];
//work          orz::DTensor tmp3(tmp2a*tmp1);
//work          tmp3.reshape_inplace(Nvir,Nvir);
//work          P_Pi.BackTransform2EXT(tmp3);
//work          for(auto &&a : irange(0L, _Nvirt))
//work            for(auto &&b : irange(0L, _Nvirt))    
//work              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*tmp3.shape(1)+b];
//work        }
        
      }//i
    }//P
    
    // term1p + term1m
//work    orz::DTensor X1(Nrhom1,Nrhom1);
//work    for(auto &&P : irange(0L, Nrhom1))
//work      for(auto &&_P_ : irange(0L, Nrhom1))
//work        for(auto &&p : irange(0L, _Nact))
//work          for(auto &&a0 : irange(0L, _Nact))
//work          X1(P,_P_) += (2*Cm1_plus.cptr()[P*Cm1_plus.shape(1)+p]-Cm1_minus.cptr()[P*Cm1_plus.shape(1)+p]) * d1(p,a0) *
//work            Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
//work    orz::DTensor Y1(_Nvirt,_Nvirt,Nrhom1,_Ndocc);
//work    for(auto &&P : irange(0L, Nrhom1))
//work      for(auto &&i : irange(0L, _Ndocc))
//work        for(auto &&_P_ : irange(0L, Nrhom1))
//work          for(auto &&v1 : irange(0L, _Nvirt))
//work            for(auto &&v0 : irange(0L, _Nvirt))
//work              Y1(v1,v0,P,i) += X1(P,_P_) * T_m1_.cptr()[v1*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+v0*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+i];
//work    for(auto &&a : irange(0L, _Nvirt))
//work      for(auto &&b : irange(0L, _Nvirt))    
//work        for(auto &&v1 : irange(0L, _Nvirt))
//work          for(auto &&v0 : irange(0L, _Nvirt))
//work            for(auto &&P : irange(0L, Nrhom1))
//work              for(auto &&i : irange(0L, _Ndocc))
//work                tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += V2_FiC_vvvv.cptr()[b*V2_FiC_vvvv.shape(1)*V2_FiC_vvvv.shape(2)*V2_FiC_vvvv.shape(3)+v0*V2_FiC_vvvv.shape(2)*V2_FiC_vvvv.shape(3)+a*V2_FiC_vvvv.shape(3)+v1] * Y1(v1,v0,P,i);
  }
  
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @Tc(-1)(v0,b,_P_,i) @W2(a,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a0
          }//P
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a3,a2,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Nvirt+v0)*tmp1.shape(1)+(a2*_Nact+a3)] = V2_FiC_vvaa.cptr()[a*V2_FiC_vvaa.shape(1)*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+v0*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+a2*V2_FiC_vvaa.shape(3)+a3];
            }//a3
          }//a2
        }//v0
      }//a
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a2*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a2
      // @W2(a,v0,P,_P_) <<= +1 @V2_FiC_vvaa(a,v0,a2,a3) @W1(a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[a*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*Nrhom1+_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*Nrhom1,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*Nrhom1+_P_)*tmp2.shape(1)+(a*Nrhom1+P)] = W2.cptr()[a*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//a
        }//_P_
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(v0,b,_P_,i) @W2(a,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @Tc(-1)(v0,a,_P_,i) @W2(b,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a0
          }//P
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a3,a2,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Nvirt+v0)*tmp1.shape(1)+(a2*_Nact+a3)] = V2_FiC_vvaa.cptr()[b*V2_FiC_vvaa.shape(1)*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+v0*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+a2*V2_FiC_vvaa.shape(3)+a3];
            }//a3
          }//a2
        }//v0
      }//b
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a2*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a2
      // @W2(b,v0,P,_P_) <<= +1 @V2_FiC_vvaa(b,v0,a2,a3) @W1(a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[b*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[b*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*Nrhom1+_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*Nrhom1,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*Nrhom1+_P_)*tmp2.shape(1)+(b*Nrhom1+P)] = W2.cptr()[b*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//b
        }//_P_
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(v0,a,_P_,i) @W2(b,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @Tc(-1)(b,v0,_P_,i) @W2(a,v0,_P_,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+_P_)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a0
          }//_P_
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @W1(a3,a2,_P_,P) <<= +1 @W0(a3,a2,a0,_P_) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+_P_*Nrhom1+P] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+_P_*Nrhom1+P];
            }//P
          }//_P_
        }//a2
      }//a3
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Nvirt+v0)*tmp1.shape(1)+(a2*_Nact+a3)] = V2_FiC_vvaa.cptr()[a*V2_FiC_vvaa.shape(1)*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+v0*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+a2*V2_FiC_vvaa.shape(3)+a3];
            }//a3
          }//a2
        }//v0
      }//a
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(a2*_Nact+a3)*tmp2.shape(1)+(_P_*Nrhom1+P)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+_P_*W1.shape(3)+P];
            }//P
          }//_P_
        }//a3
      }//a2
      // @W2(a,v0,_P_,P) <<= +1 @V2_FiC_vvaa(a,v0,a2,a3) @W1(a3,a2,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[a*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+_P_*Nrhom1+P] += tmp3.cptr()[a*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+_P_*Nrhom1+P];
            }//P
          }//_P_
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*Nrhom1+_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*Nrhom1,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*Nrhom1+_P_)*tmp2.shape(1)+(a*Nrhom1+P)] = W2.cptr()[a*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+_P_*W2.shape(3)+P];
            }//P
          }//a
        }//_P_
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(b,v0,_P_,i) @W2(a,v0,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @Tc(-1)(a,v0,_P_,i) @W2(b,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a0
          }//P
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a3,a2,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Nvirt+v0)*tmp1.shape(1)+(a2*_Nact+a3)] = V2_FiC_vvaa.cptr()[b*V2_FiC_vvaa.shape(1)*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+v0*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+a2*V2_FiC_vvaa.shape(3)+a3];
            }//a3
          }//a2
        }//v0
      }//b
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a2*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a2
      // @W2(b,v0,P,_P_) <<= +1 @V2_FiC_vvaa(b,v0,a2,a3) @W1(a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[b*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[b*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*Nrhom1+_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*Nrhom1,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*Nrhom1+_P_)*tmp2.shape(1)+(b*Nrhom1+P)] = W2.cptr()[b*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//b
        }//_P_
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(a,v0,_P_,i) @W2(b,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @V2_FiC_vvcc(a,v0,c0,i) @W2(c0,v0,b,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(v0*_Nvirt*_Ndocc+b*_Ndocc+c0)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//c0
        }//b
      }//v0
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @W2(c0,v0,b,P) <<= +1 @Tc(-1)(v0,b,_P_,c0) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Nvirt*_Nvirt*Nrhom1+v0*_Nvirt*Nrhom1+b*Nrhom1+P] += tmp3.cptr()[v0*_Nvirt*_Ndocc*Nrhom1+b*_Ndocc*Nrhom1+c0*Nrhom1+P];
            }//P
          }//c0
        }//b
      }//v0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Ndocc);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*_Ndocc+c0)] = V2_FiC_vvcc.cptr()[a*V2_FiC_vvcc.shape(1)*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+v0*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+c0*V2_FiC_vvcc.shape(3)+i];
            }//c0
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(b*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+P];
            }//P
          }//b
        }//c0
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @V2_FiC_vvcc(a,v0,c0,i) @W2(c0,v0,b,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @V2_FiC_vvcc(b,v0,c0,i) @W2(c0,v0,a,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(v0*_Nvirt*_Ndocc+a*_Ndocc+c0)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//c0
        }//a
      }//v0
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @W2(c0,v0,a,P) <<= +1 @Tc(-1)(v0,a,_P_,c0) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Nvirt*_Nvirt*Nrhom1+v0*_Nvirt*Nrhom1+a*Nrhom1+P] += tmp3.cptr()[v0*_Nvirt*_Ndocc*Nrhom1+a*_Ndocc*Nrhom1+c0*Nrhom1+P];
            }//P
          }//c0
        }//a
      }//v0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*_Ndocc+c0)] = V2_FiC_vvcc.cptr()[b*V2_FiC_vvcc.shape(1)*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+v0*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+c0*V2_FiC_vvcc.shape(3)+i];
            }//c0
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(a*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+P];
            }//P
          }//a
        }//c0
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @V2_FiC_vvcc(b,v0,c0,i) @W2(c0,v0,a,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @V2_FiC_vvcc(a,v0,c0,i) @W2(c0,b,v0,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a0,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+_P_] += tmp3.cptr()[a0*Nrhom1+_P_];
        }//_P_
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(_P_)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+_P_];
        }//a0
      }//_P_
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @W1(_P_,P) <<= +1 @W0(a0,_P_) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          W1.cptr()[_P_*Nrhom1+P] += tmp3.cptr()[_P_*Nrhom1+P];
        }//P
      }//_P_
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Nvirt*_Ndocc+v0*_Ndocc+c0)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//c0
        }//v0
      }//b
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[_P_*W1.shape(1)+P];
        }//P
      }//_P_
      // @W2(c0,b,v0,P) <<= +1 @Tc(-1)(b,v0,_P_,c0) @W1(_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+v0*Nrhom1+P] += tmp3.cptr()[b*_Nvirt*_Ndocc*Nrhom1+v0*_Ndocc*Nrhom1+c0*Nrhom1+P];
            }//P
          }//c0
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Ndocc);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*_Ndocc+c0)] = V2_FiC_vvcc.cptr()[a*V2_FiC_vvcc.shape(1)*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+v0*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+c0*V2_FiC_vvcc.shape(3)+i];
            }//c0
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(b*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+b*W2.shape(2)*W2.shape(3)+v0*W2.shape(3)+P];
            }//P
          }//b
        }//c0
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @V2_FiC_vvcc(a,v0,c0,i) @W2(c0,b,v0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @V2_FiC_vvcc(b,v0,c0,i) @W2(c0,a,v0,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Nvirt*_Ndocc+v0*_Ndocc+c0)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//c0
        }//v0
      }//a
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @W2(c0,a,v0,P) <<= +1 @Tc(-1)(a,v0,_P_,c0) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+v0*Nrhom1+P] += tmp3.cptr()[a*_Nvirt*_Ndocc*Nrhom1+v0*_Ndocc*Nrhom1+c0*Nrhom1+P];
            }//P
          }//c0
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*_Ndocc+c0)] = V2_FiC_vvcc.cptr()[b*V2_FiC_vvcc.shape(1)*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+v0*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+c0*V2_FiC_vvcc.shape(3)+i];
            }//c0
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(a*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+a*W2.shape(2)*W2.shape(3)+v0*W2.shape(3)+P];
            }//P
          }//a
        }//c0
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @V2_FiC_vvcc(b,v0,c0,i) @W2(c0,a,v0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +2 @V2_FiC_vcvc(a,i,v0,c0) @W2(c0,v0,b,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   2.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(v0*_Nvirt*_Ndocc+b*_Ndocc+c0)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//c0
        }//b
      }//v0
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @W2(c0,v0,b,P) <<= +1 @Tc(-1)(v0,b,_P_,c0) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Nvirt*_Nvirt*Nrhom1+v0*_Nvirt*Nrhom1+b*Nrhom1+P] += tmp3.cptr()[v0*_Nvirt*_Ndocc*Nrhom1+b*_Ndocc*Nrhom1+c0*Nrhom1+P];
            }//P
          }//c0
        }//b
      }//v0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Ndocc);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*_Ndocc+c0)] = V2_FiC_vcvc.cptr()[a*V2_FiC_vcvc.shape(1)*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+i*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+v0*V2_FiC_vcvc.shape(3)+c0];
            }//c0
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(b*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+P];
            }//P
          }//b
        }//c0
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @V2_FiC_vcvc(a,i,v0,c0) @W2(c0,v0,b,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +2 @V2_FiC_vcvc(b,i,v0,c0) @W2(c0,v0,a,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   2.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(v0*_Nvirt*_Ndocc+a*_Ndocc+c0)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//c0
        }//a
      }//v0
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @W2(c0,v0,a,P) <<= +1 @Tc(-1)(v0,a,_P_,c0) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Nvirt*_Nvirt*Nrhom1+v0*_Nvirt*Nrhom1+a*Nrhom1+P] += tmp3.cptr()[v0*_Nvirt*_Ndocc*Nrhom1+a*_Ndocc*Nrhom1+c0*Nrhom1+P];
            }//P
          }//c0
        }//a
      }//v0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*_Ndocc+c0)] = V2_FiC_vcvc.cptr()[b*V2_FiC_vcvc.shape(1)*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+i*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+v0*V2_FiC_vcvc.shape(3)+c0];
            }//c0
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(a*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+P];
            }//P
          }//a
        }//c0
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @V2_FiC_vcvc(b,i,v0,c0) @W2(c0,v0,a,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +2 @V2_FiC_vcvc(a,i,v0,c0) @W2(c0,b,v0,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   2.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a0,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+_P_] += tmp3.cptr()[a0*Nrhom1+_P_];
        }//_P_
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(_P_)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+_P_];
        }//a0
      }//_P_
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @W1(_P_,P) <<= +1 @W0(a0,_P_) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          W1.cptr()[_P_*Nrhom1+P] += tmp3.cptr()[_P_*Nrhom1+P];
        }//P
      }//_P_
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Nvirt*_Ndocc+v0*_Ndocc+c0)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//c0
        }//v0
      }//b
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[_P_*W1.shape(1)+P];
        }//P
      }//_P_
      // @W2(c0,b,v0,P) <<= +1 @Tc(-1)(b,v0,_P_,c0) @W1(_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+v0*Nrhom1+P] += tmp3.cptr()[b*_Nvirt*_Ndocc*Nrhom1+v0*_Ndocc*Nrhom1+c0*Nrhom1+P];
            }//P
          }//c0
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Ndocc);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*_Ndocc+c0)] = V2_FiC_vcvc.cptr()[a*V2_FiC_vcvc.shape(1)*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+i*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+v0*V2_FiC_vcvc.shape(3)+c0];
            }//c0
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(b*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+b*W2.shape(2)*W2.shape(3)+v0*W2.shape(3)+P];
            }//P
          }//b
        }//c0
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @V2_FiC_vcvc(a,i,v0,c0) @W2(c0,b,v0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +2 @V2_FiC_vcvc(b,i,v0,c0) @W2(c0,a,v0,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   2.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Nvirt*_Ndocc+v0*_Ndocc+c0)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//c0
        }//v0
      }//a
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @W2(c0,a,v0,P) <<= +1 @Tc(-1)(a,v0,_P_,c0) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+v0*Nrhom1+P] += tmp3.cptr()[a*_Nvirt*_Ndocc*Nrhom1+v0*_Ndocc*Nrhom1+c0*Nrhom1+P];
            }//P
          }//c0
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*_Ndocc+c0)] = V2_FiC_vcvc.cptr()[b*V2_FiC_vcvc.shape(1)*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+i*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+v0*V2_FiC_vcvc.shape(3)+c0];
            }//c0
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(a*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+a*W2.shape(2)*W2.shape(3)+v0*W2.shape(3)+P];
            }//P
          }//a
        }//c0
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @V2_FiC_vcvc(b,i,v0,c0) @W2(c0,a,v0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @V2_FiC_vcvc(b,i,v0,c0) @W2(c0,v0,a,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(v0*_Nvirt*_Ndocc+a*_Ndocc+c0)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//c0
        }//a
      }//v0
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @W2(c0,v0,a,P) <<= +1 @Tc(-1)(v0,a,_P_,c0) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Nvirt*_Nvirt*Nrhom1+v0*_Nvirt*Nrhom1+a*Nrhom1+P] += tmp3.cptr()[v0*_Nvirt*_Ndocc*Nrhom1+a*_Ndocc*Nrhom1+c0*Nrhom1+P];
            }//P
          }//c0
        }//a
      }//v0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*_Ndocc+c0)] = V2_FiC_vcvc.cptr()[b*V2_FiC_vcvc.shape(1)*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+i*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+v0*V2_FiC_vcvc.shape(3)+c0];
            }//c0
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(a*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+P];
            }//P
          }//a
        }//c0
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @V2_FiC_vcvc(b,i,v0,c0) @W2(c0,v0,a,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @V2_FiC_vcvc(a,i,v0,c0) @W2(c0,v0,b,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(v0*_Nvirt*_Ndocc+b*_Ndocc+c0)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//c0
        }//b
      }//v0
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @W2(c0,v0,b,P) <<= +1 @Tc(-1)(v0,b,_P_,c0) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Nvirt*_Nvirt*Nrhom1+v0*_Nvirt*Nrhom1+b*Nrhom1+P] += tmp3.cptr()[v0*_Nvirt*_Ndocc*Nrhom1+b*_Ndocc*Nrhom1+c0*Nrhom1+P];
            }//P
          }//c0
        }//b
      }//v0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Ndocc);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*_Ndocc+c0)] = V2_FiC_vcvc.cptr()[a*V2_FiC_vcvc.shape(1)*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+i*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+v0*V2_FiC_vcvc.shape(3)+c0];
            }//c0
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(b*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+P];
            }//P
          }//b
        }//c0
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @V2_FiC_vcvc(a,i,v0,c0) @W2(c0,v0,b,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @V2_FiC_vcvc(b,i,v0,c0) @W2(c0,a,v0,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a0,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+_P_] += tmp3.cptr()[a0*Nrhom1+_P_];
        }//_P_
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(_P_)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+_P_];
        }//a0
      }//_P_
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @W1(_P_,P) <<= +1 @W0(a0,_P_) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          W1.cptr()[_P_*Nrhom1+P] += tmp3.cptr()[_P_*Nrhom1+P];
        }//P
      }//_P_
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Nvirt*_Ndocc+v0*_Ndocc+c0)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//c0
        }//v0
      }//a
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[_P_*W1.shape(1)+P];
        }//P
      }//_P_
      // @W2(c0,a,v0,P) <<= +1 @Tc(-1)(a,v0,_P_,c0) @W1(_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+v0*Nrhom1+P] += tmp3.cptr()[a*_Nvirt*_Ndocc*Nrhom1+v0*_Ndocc*Nrhom1+c0*Nrhom1+P];
            }//P
          }//c0
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*_Ndocc+c0)] = V2_FiC_vcvc.cptr()[b*V2_FiC_vcvc.shape(1)*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+i*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+v0*V2_FiC_vcvc.shape(3)+c0];
            }//c0
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(a*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+a*W2.shape(2)*W2.shape(3)+v0*W2.shape(3)+P];
            }//P
          }//a
        }//c0
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @V2_FiC_vcvc(b,i,v0,c0) @W2(c0,a,v0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @V2_FiC_vcvc(a,i,v0,c0) @W2(c0,b,v0,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Nvirt*_Ndocc+v0*_Ndocc+c0)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//c0
        }//v0
      }//b
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @W2(c0,b,v0,P) <<= +1 @Tc(-1)(b,v0,_P_,c0) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+v0*Nrhom1+P] += tmp3.cptr()[b*_Nvirt*_Ndocc*Nrhom1+v0*_Ndocc*Nrhom1+c0*Nrhom1+P];
            }//P
          }//c0
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Ndocc);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*_Ndocc+c0)] = V2_FiC_vcvc.cptr()[a*V2_FiC_vcvc.shape(1)*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+i*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+v0*V2_FiC_vcvc.shape(3)+c0];
            }//c0
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(b*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+b*W2.shape(2)*W2.shape(3)+v0*W2.shape(3)+P];
            }//P
          }//b
        }//c0
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @V2_FiC_vcvc(a,i,v0,c0) @W2(c0,b,v0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @Tc(-1)(a,v0,_P_,i) @W2(b,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a0
          }//P
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a3,a2,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Nvirt+v0)*tmp1.shape(1)+(a2*_Nact+a3)] = V2_FiC_vvaa.cptr()[b*V2_FiC_vvaa.shape(1)*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+v0*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+a2*V2_FiC_vvaa.shape(3)+a3];
            }//a3
          }//a2
        }//v0
      }//b
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a2*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a2
      // @W2(b,v0,P,_P_) <<= +1 @V2_FiC_vvaa(b,v0,a2,a3) @W1(a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[b*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[b*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*Nrhom1+_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*Nrhom1,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*Nrhom1+_P_)*tmp2.shape(1)+(b*Nrhom1+P)] = W2.cptr()[b*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//b
        }//_P_
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(a,v0,_P_,i) @W2(b,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @Tc(-1)(b,v0,_P_,i) @W2(a,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a0
          }//P
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a3,a2,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Nvirt+v0)*tmp1.shape(1)+(a2*_Nact+a3)] = V2_FiC_vvaa.cptr()[a*V2_FiC_vvaa.shape(1)*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+v0*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+a2*V2_FiC_vvaa.shape(3)+a3];
            }//a3
          }//a2
        }//v0
      }//a
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a2*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a2
      // @W2(a,v0,P,_P_) <<= +1 @V2_FiC_vvaa(a,v0,a2,a3) @W1(a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[a*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*Nrhom1+_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*Nrhom1,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*Nrhom1+_P_)*tmp2.shape(1)+(a*Nrhom1+P)] = W2.cptr()[a*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//a
        }//_P_
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(b,v0,_P_,i) @W2(a,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @Tc(-1)(v0,a,_P_,i) @W2(b,v0,_P_,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+_P_)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a0
          }//_P_
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @W1(a3,a2,_P_,P) <<= +1 @W0(a3,a2,a0,_P_) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+_P_*Nrhom1+P] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+_P_*Nrhom1+P];
            }//P
          }//_P_
        }//a2
      }//a3
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Nvirt+v0)*tmp1.shape(1)+(a2*_Nact+a3)] = V2_FiC_vvaa.cptr()[b*V2_FiC_vvaa.shape(1)*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+v0*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+a2*V2_FiC_vvaa.shape(3)+a3];
            }//a3
          }//a2
        }//v0
      }//b
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(a2*_Nact+a3)*tmp2.shape(1)+(_P_*Nrhom1+P)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+_P_*W1.shape(3)+P];
            }//P
          }//_P_
        }//a3
      }//a2
      // @W2(b,v0,_P_,P) <<= +1 @V2_FiC_vvaa(b,v0,a2,a3) @W1(a3,a2,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[b*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+_P_*Nrhom1+P] += tmp3.cptr()[b*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+_P_*Nrhom1+P];
            }//P
          }//_P_
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*Nrhom1+_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*Nrhom1,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*Nrhom1+_P_)*tmp2.shape(1)+(b*Nrhom1+P)] = W2.cptr()[b*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+_P_*W2.shape(3)+P];
            }//P
          }//b
        }//_P_
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(v0,a,_P_,i) @W2(b,v0,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @Tc(-1)(v0,b,_P_,i) @W2(a,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a0
          }//P
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a3,a2,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Nvirt+v0)*tmp1.shape(1)+(a2*_Nact+a3)] = V2_FiC_vvaa.cptr()[a*V2_FiC_vvaa.shape(1)*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+v0*V2_FiC_vvaa.shape(2)*V2_FiC_vvaa.shape(3)+a2*V2_FiC_vvaa.shape(3)+a3];
            }//a3
          }//a2
        }//v0
      }//a
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a2*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a2
      // @W2(a,v0,P,_P_) <<= +1 @V2_FiC_vvaa(a,v0,a2,a3) @W1(a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[a*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*Nrhom1+_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*Nrhom1,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*Nrhom1+_P_)*tmp2.shape(1)+(a*Nrhom1+P)] = W2.cptr()[a*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//a
        }//_P_
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(v0,b,_P_,i) @W2(a,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @V2_FiC_vvcc(b,v0,c0,i) @W2(c0,a,v0,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Nvirt*_Ndocc+v0*_Ndocc+c0)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//c0
        }//v0
      }//a
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @W2(c0,a,v0,P) <<= +1 @Tc(-1)(a,v0,_P_,c0) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+v0*Nrhom1+P] += tmp3.cptr()[a*_Nvirt*_Ndocc*Nrhom1+v0*_Ndocc*Nrhom1+c0*Nrhom1+P];
            }//P
          }//c0
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*_Ndocc+c0)] = V2_FiC_vvcc.cptr()[b*V2_FiC_vvcc.shape(1)*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+v0*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+c0*V2_FiC_vvcc.shape(3)+i];
            }//c0
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(a*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+a*W2.shape(2)*W2.shape(3)+v0*W2.shape(3)+P];
            }//P
          }//a
        }//c0
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @V2_FiC_vvcc(b,v0,c0,i) @W2(c0,a,v0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @V2_FiC_vvcc(a,v0,c0,i) @W2(c0,b,v0,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Nvirt*_Ndocc+v0*_Ndocc+c0)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//c0
        }//v0
      }//b
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @W2(c0,b,v0,P) <<= +1 @Tc(-1)(b,v0,_P_,c0) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+v0*Nrhom1+P] += tmp3.cptr()[b*_Nvirt*_Ndocc*Nrhom1+v0*_Ndocc*Nrhom1+c0*Nrhom1+P];
            }//P
          }//c0
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Ndocc);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*_Ndocc+c0)] = V2_FiC_vvcc.cptr()[a*V2_FiC_vvcc.shape(1)*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+v0*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+c0*V2_FiC_vvcc.shape(3)+i];
            }//c0
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(b*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+b*W2.shape(2)*W2.shape(3)+v0*W2.shape(3)+P];
            }//P
          }//b
        }//c0
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @V2_FiC_vvcc(a,v0,c0,i) @W2(c0,b,v0,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @V2_FiC_vvcc(b,v0,c0,i) @W2(c0,v0,a,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a0,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+_P_] += tmp3.cptr()[a0*Nrhom1+_P_];
        }//_P_
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(_P_)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+_P_];
        }//a0
      }//_P_
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @W1(_P_,P) <<= +1 @W0(a0,_P_) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          W1.cptr()[_P_*Nrhom1+P] += tmp3.cptr()[_P_*Nrhom1+P];
        }//P
      }//_P_
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(v0*_Nvirt*_Ndocc+a*_Ndocc+c0)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//c0
        }//a
      }//v0
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[_P_*W1.shape(1)+P];
        }//P
      }//_P_
      // @W2(c0,v0,a,P) <<= +1 @Tc(-1)(v0,a,_P_,c0) @W1(_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Nvirt*_Nvirt*Nrhom1+v0*_Nvirt*Nrhom1+a*Nrhom1+P] += tmp3.cptr()[v0*_Nvirt*_Ndocc*Nrhom1+a*_Ndocc*Nrhom1+c0*Nrhom1+P];
            }//P
          }//c0
        }//a
      }//v0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*_Ndocc+c0)] = V2_FiC_vvcc.cptr()[b*V2_FiC_vvcc.shape(1)*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+v0*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+c0*V2_FiC_vvcc.shape(3)+i];
            }//c0
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(a*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+P];
            }//P
          }//a
        }//c0
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @V2_FiC_vvcc(b,v0,c0,i) @W2(c0,v0,a,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @V2_FiC_vvcc(a,v0,c0,i) @W2(c0,v0,b,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(v0*_Nvirt*_Ndocc+b*_Ndocc+c0)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//c0
        }//b
      }//v0
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @W2(c0,v0,b,P) <<= +1 @Tc(-1)(v0,b,_P_,c0) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Nvirt*_Nvirt*Nrhom1+v0*_Nvirt*Nrhom1+b*Nrhom1+P] += tmp3.cptr()[v0*_Nvirt*_Ndocc*Nrhom1+b*_Ndocc*Nrhom1+c0*Nrhom1+P];
            }//P
          }//c0
        }//b
      }//v0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Ndocc);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*_Ndocc+c0)] = V2_FiC_vvcc.cptr()[a*V2_FiC_vvcc.shape(1)*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+v0*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+c0*V2_FiC_vvcc.shape(3)+i];
            }//c0
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*_Ndocc+c0)*tmp2.shape(1)+(b*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+P];
            }//P
          }//b
        }//c0
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @V2_FiC_vvcc(a,v0,c0,i) @W2(c0,v0,b,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @Tc(-1)(a,v0,_P_,i) @W2(b,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a0*Nrhom1+P)*tmp1.shape(1)+(a2)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a2
          }//P
        }//a0
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a2)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a2];
        }//_P_
      }//a2
      // @W1(a3,a0,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(-)(_P_,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a0
      }//a3
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Nvirt+v0)*tmp1.shape(1)+(a0*_Nact+a3)] = V2_FiC_vava.cptr()[b*V2_FiC_vava.shape(1)*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+a0*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+v0*V2_FiC_vava.shape(3)+a3];
            }//a3
          }//a0
        }//v0
      }//b
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a0*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a0
      // @W2(b,v0,P,_P_) <<= +1 @V2_FiC_vava(b,a0,v0,a3) @W1(a3,a0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[b*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[b*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*Nrhom1+_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*Nrhom1,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*Nrhom1+_P_)*tmp2.shape(1)+(b*Nrhom1+P)] = W2.cptr()[b*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//b
        }//_P_
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(a,v0,_P_,i) @W2(b,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @Tc(-1)(b,v0,_P_,i) @W2(a,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a0*Nrhom1+P)*tmp1.shape(1)+(a2)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a2
          }//P
        }//a0
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a2)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a2];
        }//_P_
      }//a2
      // @W1(a3,a0,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(-)(_P_,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a0
      }//a3
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Nvirt+v0)*tmp1.shape(1)+(a0*_Nact+a3)] = V2_FiC_vava.cptr()[a*V2_FiC_vava.shape(1)*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+a0*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+v0*V2_FiC_vava.shape(3)+a3];
            }//a3
          }//a0
        }//v0
      }//a
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a0*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a0
      // @W2(a,v0,P,_P_) <<= +1 @V2_FiC_vava(a,a0,v0,a3) @W1(a3,a0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[a*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*Nrhom1+_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*Nrhom1,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*Nrhom1+_P_)*tmp2.shape(1)+(a*Nrhom1+P)] = W2.cptr()[a*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//a
        }//_P_
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(b,v0,_P_,i) @W2(a,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @Tc(-1)(v0,a,_P_,i) @W2(b,v0,_P_,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a1*_Nact+a0)*tmp1.shape(1)+(a2)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a2
          }//a0
        }//a1
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a2)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a2];
        }//_P_
      }//a2
      // @W0(a3,a1,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(_P_,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a1*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a1*_Nact*Nrhom1+a0*Nrhom1+_P_];
            }//_P_
          }//a0
        }//a1
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a0*Nrhom1+_P_)*tmp1.shape(1)+(a1)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a1
          }//_P_
        }//a0
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W1(a3,a0,_P_,P) <<= +1 @W0(a3,a1,a0,_P_) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+_P_*Nrhom1+P] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+_P_*Nrhom1+P];
            }//P
          }//_P_
        }//a0
      }//a3
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Nvirt+v0)*tmp1.shape(1)+(a0*_Nact+a3)] = V2_FiC_vava.cptr()[b*V2_FiC_vava.shape(1)*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+a0*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+v0*V2_FiC_vava.shape(3)+a3];
            }//a3
          }//a0
        }//v0
      }//b
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(a0*_Nact+a3)*tmp2.shape(1)+(_P_*Nrhom1+P)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+_P_*W1.shape(3)+P];
            }//P
          }//_P_
        }//a3
      }//a0
      // @W2(b,v0,_P_,P) <<= +1 @V2_FiC_vava(b,a0,v0,a3) @W1(a3,a0,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[b*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+_P_*Nrhom1+P] += tmp3.cptr()[b*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+_P_*Nrhom1+P];
            }//P
          }//_P_
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*Nrhom1+_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*Nrhom1,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*Nrhom1+_P_)*tmp2.shape(1)+(b*Nrhom1+P)] = W2.cptr()[b*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+_P_*W2.shape(3)+P];
            }//P
          }//b
        }//_P_
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(v0,a,_P_,i) @W2(b,v0,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @Tc(-1)(v0,b,_P_,i) @W2(a,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a0*Nrhom1+P)*tmp1.shape(1)+(a2)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a2
          }//P
        }//a0
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a2)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a2];
        }//_P_
      }//a2
      // @W1(a3,a0,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(+)(_P_,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a0
      }//a3
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Nvirt+v0)*tmp1.shape(1)+(a0*_Nact+a3)] = V2_FiC_vava.cptr()[a*V2_FiC_vava.shape(1)*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+a0*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+v0*V2_FiC_vava.shape(3)+a3];
            }//a3
          }//a0
        }//v0
      }//a
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a0*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a0
      // @W2(a,v0,P,_P_) <<= +1 @V2_FiC_vava(a,a0,v0,a3) @W1(a3,a0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[a*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*Nrhom1+_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*Nrhom1,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*Nrhom1+_P_)*tmp2.shape(1)+(a*Nrhom1+P)] = W2.cptr()[a*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//a
        }//_P_
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(v0,b,_P_,i) @W2(a,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1/2 @Tc(-1)(a,b,_P_,i) @W2(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.5000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a5,a4,a3,a2,a0,P) <<= +1 @D3(a5,a4,a3,a2,a1,a0) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&a0 : irange(0L, _Nact)){
                for(auto &&P : irange(0L, Nrhom1)){
                  W0.cptr()[a5*_Nact*_Nact*_Nact*_Nact*Nrhom1+a4*_Nact*_Nact*_Nact*Nrhom1+a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a5*_Nact*_Nact*_Nact*_Nact*Nrhom1+a4*_Nact*_Nact*_Nact*Nrhom1+a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
                }//P
              }//a0
            }//a2
          }//a3
        }//a4
      }//a5
    }
    orz::DTensor W1(_Nact,_Nact,_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom1)){
                for(auto &&a0 : irange(0L, _Nact)){
                  tmp1.cptr()[(a5*_Nact*_Nact*_Nact*Nrhom1+a4*_Nact*_Nact*Nrhom1+a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a5*W0.shape(1)*W0.shape(2)*W0.shape(3)*W0.shape(4)*W0.shape(5)+a4*W0.shape(2)*W0.shape(3)*W0.shape(4)*W0.shape(5)+a3*W0.shape(3)*W0.shape(4)*W0.shape(5)+a2*W0.shape(4)*W0.shape(5)+a0*W0.shape(5)+P];
                }//a0
              }//P
            }//a2
          }//a3
        }//a4
      }//a5
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a5,a4,a3,a2,P,_P_) <<= +1 @W0(a5,a4,a3,a2,a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom1)){
                for(auto &&_P_ : irange(0L, Nrhom1)){
                  W1.cptr()[a5*_Nact*_Nact*_Nact*Nrhom1*Nrhom1+a4*_Nact*_Nact*Nrhom1*Nrhom1+a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a5*_Nact*_Nact*_Nact*Nrhom1*Nrhom1+a4*_Nact*_Nact*Nrhom1*Nrhom1+a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
                }//_P_
              }//P
            }//a2
          }//a3
        }//a4
      }//a5
    }
    orz::DTensor W2(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(1,_Nact*_Nact*_Nact*_Nact);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              tmp1.cptr()[(0)*tmp1.shape(1)+(a2*_Nact*_Nact*_Nact+a3*_Nact*_Nact+a4*_Nact+a5)] = V2_FiC_aaaa.cptr()[a2*V2_FiC_aaaa.shape(1)*V2_FiC_aaaa.shape(2)*V2_FiC_aaaa.shape(3)+a3*V2_FiC_aaaa.shape(2)*V2_FiC_aaaa.shape(3)+a4*V2_FiC_aaaa.shape(3)+a5];
            }//a5
          }//a4
        }//a3
      }//a2
      orz::DTensor tmp2(_Nact*_Nact*_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom1)){
                for(auto &&_P_ : irange(0L, Nrhom1)){
                  tmp2.cptr()[(a2*_Nact*_Nact*_Nact+a3*_Nact*_Nact+a4*_Nact+a5)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a5*W1.shape(1)*W1.shape(2)*W1.shape(3)*W1.shape(4)*W1.shape(5)+a4*W1.shape(2)*W1.shape(3)*W1.shape(4)*W1.shape(5)+a3*W1.shape(3)*W1.shape(4)*W1.shape(5)+a2*W1.shape(4)*W1.shape(5)+P*W1.shape(5)+_P_];
                }//_P_
              }//P
            }//a5
          }//a4
        }//a3
      }//a2
      // @W2(P,_P_) <<= +1 @V2_FiC_aaaa(a2,a3,a4,a5) @W1(a5,a4,a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W2.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Nvirt*_Ndocc+b*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//b
      }//a
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W2.cptr()[P*W2.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(a,b,_P_,i) @W2(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Nvirt*_Ndocc*Nrhom1+b*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1/2 @Tc(-1)(b,a,_P_,i) @W2(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.5000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a5,a4,a3,a2,a0,P) <<= +1 @D3(a5,a4,a3,a2,a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&a0 : irange(0L, _Nact)){
                for(auto &&P : irange(0L, Nrhom1)){
                  W0.cptr()[a5*_Nact*_Nact*_Nact*_Nact*Nrhom1+a4*_Nact*_Nact*_Nact*Nrhom1+a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a5*_Nact*_Nact*_Nact*_Nact*Nrhom1+a4*_Nact*_Nact*_Nact*Nrhom1+a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
                }//P
              }//a0
            }//a2
          }//a3
        }//a4
      }//a5
    }
    orz::DTensor W1(_Nact,_Nact,_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom1)){
                for(auto &&a0 : irange(0L, _Nact)){
                  tmp1.cptr()[(a5*_Nact*_Nact*_Nact*Nrhom1+a4*_Nact*_Nact*Nrhom1+a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a5*W0.shape(1)*W0.shape(2)*W0.shape(3)*W0.shape(4)*W0.shape(5)+a4*W0.shape(2)*W0.shape(3)*W0.shape(4)*W0.shape(5)+a3*W0.shape(3)*W0.shape(4)*W0.shape(5)+a2*W0.shape(4)*W0.shape(5)+a0*W0.shape(5)+P];
                }//a0
              }//P
            }//a2
          }//a3
        }//a4
      }//a5
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a5,a4,a3,a2,P,_P_) <<= +1 @W0(a5,a4,a3,a2,a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom1)){
                for(auto &&_P_ : irange(0L, Nrhom1)){
                  W1.cptr()[a5*_Nact*_Nact*_Nact*Nrhom1*Nrhom1+a4*_Nact*_Nact*Nrhom1*Nrhom1+a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a5*_Nact*_Nact*_Nact*Nrhom1*Nrhom1+a4*_Nact*_Nact*Nrhom1*Nrhom1+a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
                }//_P_
              }//P
            }//a2
          }//a3
        }//a4
      }//a5
    }
    orz::DTensor W2(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(1,_Nact*_Nact*_Nact*_Nact);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              tmp1.cptr()[(0)*tmp1.shape(1)+(a2*_Nact*_Nact*_Nact+a3*_Nact*_Nact+a4*_Nact+a5)] = V2_FiC_aaaa.cptr()[a2*V2_FiC_aaaa.shape(1)*V2_FiC_aaaa.shape(2)*V2_FiC_aaaa.shape(3)+a3*V2_FiC_aaaa.shape(2)*V2_FiC_aaaa.shape(3)+a4*V2_FiC_aaaa.shape(3)+a5];
            }//a5
          }//a4
        }//a3
      }//a2
      orz::DTensor tmp2(_Nact*_Nact*_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom1)){
                for(auto &&_P_ : irange(0L, Nrhom1)){
                  tmp2.cptr()[(a2*_Nact*_Nact*_Nact+a3*_Nact*_Nact+a4*_Nact+a5)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a5*W1.shape(1)*W1.shape(2)*W1.shape(3)*W1.shape(4)*W1.shape(5)+a4*W1.shape(2)*W1.shape(3)*W1.shape(4)*W1.shape(5)+a3*W1.shape(3)*W1.shape(4)*W1.shape(5)+a2*W1.shape(4)*W1.shape(5)+P*W1.shape(5)+_P_];
                }//_P_
              }//P
            }//a5
          }//a4
        }//a3
      }//a2
      // @W2(P,_P_) <<= +1 @V2_FiC_aaaa(a2,a3,a4,a5) @W1(a5,a4,a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W2.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Nvirt*_Ndocc+a*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//a
      }//b
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W2.cptr()[P*W2.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(b,a,_P_,i) @W2(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Nvirt*_Ndocc*Nrhom1+a*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1/2 @Tc(-1)(b,a,_P_,i) @W2(_P_,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.5000000000000000;
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
    orz::DTensor W1(_Nact,_Nact,_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&_P_ : irange(0L, Nrhom1)){
                for(auto &&a0 : irange(0L, _Nact)){
                  tmp1.cptr()[(a5*_Nact*_Nact*_Nact*Nrhom1+a4*_Nact*_Nact*Nrhom1+a3*_Nact*Nrhom1+a2*Nrhom1+_P_)*tmp1.shape(1)+(a0)] = W0.cptr()[a5*W0.shape(1)*W0.shape(2)*W0.shape(3)*W0.shape(4)*W0.shape(5)+a4*W0.shape(2)*W0.shape(3)*W0.shape(4)*W0.shape(5)+a3*W0.shape(3)*W0.shape(4)*W0.shape(5)+a2*W0.shape(4)*W0.shape(5)+a0*W0.shape(5)+_P_];
                }//a0
              }//_P_
            }//a2
          }//a3
        }//a4
      }//a5
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @W1(a5,a4,a3,a2,_P_,P) <<= +1 @W0(a5,a4,a3,a2,a0,_P_) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&_P_ : irange(0L, Nrhom1)){
                for(auto &&P : irange(0L, Nrhom1)){
                  W1.cptr()[a5*_Nact*_Nact*_Nact*Nrhom1*Nrhom1+a4*_Nact*_Nact*Nrhom1*Nrhom1+a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+_P_*Nrhom1+P] += tmp3.cptr()[a5*_Nact*_Nact*_Nact*Nrhom1*Nrhom1+a4*_Nact*_Nact*Nrhom1*Nrhom1+a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+_P_*Nrhom1+P];
                }//P
              }//_P_
            }//a2
          }//a3
        }//a4
      }//a5
    }
    orz::DTensor W2(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(1,_Nact*_Nact*_Nact*_Nact);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              tmp1.cptr()[(0)*tmp1.shape(1)+(a2*_Nact*_Nact*_Nact+a3*_Nact*_Nact+a4*_Nact+a5)] = V2_FiC_aaaa.cptr()[a2*V2_FiC_aaaa.shape(1)*V2_FiC_aaaa.shape(2)*V2_FiC_aaaa.shape(3)+a3*V2_FiC_aaaa.shape(2)*V2_FiC_aaaa.shape(3)+a4*V2_FiC_aaaa.shape(3)+a5];
            }//a5
          }//a4
        }//a3
      }//a2
      orz::DTensor tmp2(_Nact*_Nact*_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              for(auto &&_P_ : irange(0L, Nrhom1)){
                for(auto &&P : irange(0L, Nrhom1)){
                  tmp2.cptr()[(a2*_Nact*_Nact*_Nact+a3*_Nact*_Nact+a4*_Nact+a5)*tmp2.shape(1)+(_P_*Nrhom1+P)] = W1.cptr()[a5*W1.shape(1)*W1.shape(2)*W1.shape(3)*W1.shape(4)*W1.shape(5)+a4*W1.shape(2)*W1.shape(3)*W1.shape(4)*W1.shape(5)+a3*W1.shape(3)*W1.shape(4)*W1.shape(5)+a2*W1.shape(4)*W1.shape(5)+_P_*W1.shape(5)+P];
                }//P
              }//_P_
            }//a5
          }//a4
        }//a3
      }//a2
      // @W2(_P_,P) <<= +1 @V2_FiC_aaaa(a2,a3,a4,a5) @W1(a5,a4,a3,a2,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          W2.cptr()[_P_*Nrhom1+P] += tmp3.cptr()[_P_*Nrhom1+P];
        }//P
      }//_P_
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Nvirt*_Ndocc+a*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//a
      }//b
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W2.cptr()[_P_*W2.shape(1)+P];
        }//P
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(b,a,_P_,i) @W2(_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Nvirt*_Ndocc*Nrhom1+a*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1/2 @Tc(-1)(a,b,_P_,i) @W2(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   0.5000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a5,a4,a3,a2,a0,P) <<= +1 @D3(a5,a4,a3,a2,a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&a0 : irange(0L, _Nact)){
                for(auto &&P : irange(0L, Nrhom1)){
                  W0.cptr()[a5*_Nact*_Nact*_Nact*_Nact*Nrhom1+a4*_Nact*_Nact*_Nact*Nrhom1+a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a5*_Nact*_Nact*_Nact*_Nact*Nrhom1+a4*_Nact*_Nact*_Nact*Nrhom1+a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
                }//P
              }//a0
            }//a2
          }//a3
        }//a4
      }//a5
    }
    orz::DTensor W1(_Nact,_Nact,_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom1)){
                for(auto &&a0 : irange(0L, _Nact)){
                  tmp1.cptr()[(a5*_Nact*_Nact*_Nact*Nrhom1+a4*_Nact*_Nact*Nrhom1+a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a5*W0.shape(1)*W0.shape(2)*W0.shape(3)*W0.shape(4)*W0.shape(5)+a4*W0.shape(2)*W0.shape(3)*W0.shape(4)*W0.shape(5)+a3*W0.shape(3)*W0.shape(4)*W0.shape(5)+a2*W0.shape(4)*W0.shape(5)+a0*W0.shape(5)+P];
                }//a0
              }//P
            }//a2
          }//a3
        }//a4
      }//a5
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a5,a4,a3,a2,P,_P_) <<= +1 @W0(a5,a4,a3,a2,a0,P) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a5 : irange(0L, _Nact)){
        for(auto &&a4 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom1)){
                for(auto &&_P_ : irange(0L, Nrhom1)){
                  W1.cptr()[a5*_Nact*_Nact*_Nact*Nrhom1*Nrhom1+a4*_Nact*_Nact*Nrhom1*Nrhom1+a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a5*_Nact*_Nact*_Nact*Nrhom1*Nrhom1+a4*_Nact*_Nact*Nrhom1*Nrhom1+a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
                }//_P_
              }//P
            }//a2
          }//a3
        }//a4
      }//a5
    }
    orz::DTensor W2(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(1,_Nact*_Nact*_Nact*_Nact);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              tmp1.cptr()[(0)*tmp1.shape(1)+(a2*_Nact*_Nact*_Nact+a3*_Nact*_Nact+a4*_Nact+a5)] = V2_FiC_aaaa.cptr()[a2*V2_FiC_aaaa.shape(1)*V2_FiC_aaaa.shape(2)*V2_FiC_aaaa.shape(3)+a3*V2_FiC_aaaa.shape(2)*V2_FiC_aaaa.shape(3)+a4*V2_FiC_aaaa.shape(3)+a5];
            }//a5
          }//a4
        }//a3
      }//a2
      orz::DTensor tmp2(_Nact*_Nact*_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a4 : irange(0L, _Nact)){
            for(auto &&a5 : irange(0L, _Nact)){
              for(auto &&P : irange(0L, Nrhom1)){
                for(auto &&_P_ : irange(0L, Nrhom1)){
                  tmp2.cptr()[(a2*_Nact*_Nact*_Nact+a3*_Nact*_Nact+a4*_Nact+a5)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a5*W1.shape(1)*W1.shape(2)*W1.shape(3)*W1.shape(4)*W1.shape(5)+a4*W1.shape(2)*W1.shape(3)*W1.shape(4)*W1.shape(5)+a3*W1.shape(3)*W1.shape(4)*W1.shape(5)+a2*W1.shape(4)*W1.shape(5)+P*W1.shape(5)+_P_];
                }//_P_
              }//P
            }//a5
          }//a4
        }//a3
      }//a2
      // @W2(P,_P_) <<= +1 @V2_FiC_aaaa(a2,a3,a4,a5) @W1(a5,a4,a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W2.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Nvirt*_Ndocc+b*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//b
      }//a
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W2.cptr()[P*W2.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(a,b,_P_,i) @W2(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Nvirt*_Ndocc*Nrhom1+b*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @Tc(-1)(a,b,_P_,c0) @W2(c0,i,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a0
          }//P
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a3,a2,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
    }
    orz::DTensor W2(_Ndocc,_Ndocc,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc,_Nact*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(c0*_Ndocc+i)*tmp1.shape(1)+(a2*_Nact+a3)] = V2_FiC_aacc.cptr()[a2*V2_FiC_aacc.shape(1)*V2_FiC_aacc.shape(2)*V2_FiC_aacc.shape(3)+a3*V2_FiC_aacc.shape(2)*V2_FiC_aacc.shape(3)+c0*V2_FiC_aacc.shape(3)+i];
            }//a3
          }//a2
        }//i
      }//c0
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a2*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a2
      // @W2(c0,i,P,_P_) <<= +1 @V2_FiC_aacc(a2,a3,c0,i) @W1(a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//i
      }//c0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom1*_Ndocc);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(a*_Nvirt+b)*tmp1.shape(1)+(_P_*_Ndocc+c0)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//c0
          }//_P_
        }//b
      }//a
      orz::DTensor tmp2(Nrhom1*_Ndocc,_Ndocc*Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(_P_*_Ndocc+c0)*tmp2.shape(1)+(i*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+i*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//i
        }//c0
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(a,b,_P_,c0) @W2(c0,i,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Nvirt*_Ndocc*Nrhom1+b*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @Tc(-1)(b,a,_P_,c0) @W2(c0,i,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a0
          }//P
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a3,a2,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
    }
    orz::DTensor W2(_Ndocc,_Ndocc,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc,_Nact*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(c0*_Ndocc+i)*tmp1.shape(1)+(a2*_Nact+a3)] = V2_FiC_aacc.cptr()[a2*V2_FiC_aacc.shape(1)*V2_FiC_aacc.shape(2)*V2_FiC_aacc.shape(3)+a3*V2_FiC_aacc.shape(2)*V2_FiC_aacc.shape(3)+c0*V2_FiC_aacc.shape(3)+i];
            }//a3
          }//a2
        }//i
      }//c0
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a2*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a2
      // @W2(c0,i,P,_P_) <<= +1 @V2_FiC_aacc(a2,a3,c0,i) @W1(a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//i
      }//c0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom1*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*_Nvirt+a)*tmp1.shape(1)+(_P_*_Ndocc+c0)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//c0
          }//_P_
        }//a
      }//b
      orz::DTensor tmp2(Nrhom1*_Ndocc,_Ndocc*Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(_P_*_Ndocc+c0)*tmp2.shape(1)+(i*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+i*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//i
        }//c0
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(b,a,_P_,c0) @W2(c0,i,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Nvirt*_Ndocc*Nrhom1+a*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @Tc(-1)(b,a,_P_,c0) @W2(c0,i,_P_,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
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
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+_P_)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a0
          }//_P_
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @W1(a3,a2,_P_,P) <<= +1 @W0(a3,a2,a0,_P_) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+_P_*Nrhom1+P] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+_P_*Nrhom1+P];
            }//P
          }//_P_
        }//a2
      }//a3
    }
    orz::DTensor W2(_Ndocc,_Ndocc,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc,_Nact*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(c0*_Ndocc+i)*tmp1.shape(1)+(a2*_Nact+a3)] = V2_FiC_aacc.cptr()[a2*V2_FiC_aacc.shape(1)*V2_FiC_aacc.shape(2)*V2_FiC_aacc.shape(3)+a3*V2_FiC_aacc.shape(2)*V2_FiC_aacc.shape(3)+c0*V2_FiC_aacc.shape(3)+i];
            }//a3
          }//a2
        }//i
      }//c0
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(a2*_Nact+a3)*tmp2.shape(1)+(_P_*Nrhom1+P)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+_P_*W1.shape(3)+P];
            }//P
          }//_P_
        }//a3
      }//a2
      // @W2(c0,i,_P_,P) <<= +1 @V2_FiC_aacc(a2,a3,c0,i) @W1(a3,a2,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+_P_*Nrhom1+P] += tmp3.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+_P_*Nrhom1+P];
            }//P
          }//_P_
        }//i
      }//c0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom1*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*_Nvirt+a)*tmp1.shape(1)+(_P_*_Ndocc+c0)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//c0
          }//_P_
        }//a
      }//b
      orz::DTensor tmp2(Nrhom1*_Ndocc,_Ndocc*Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(_P_*_Ndocc+c0)*tmp2.shape(1)+(i*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+i*W2.shape(2)*W2.shape(3)+_P_*W2.shape(3)+P];
            }//P
          }//i
        }//c0
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(b,a,_P_,c0) @W2(c0,i,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Nvirt*_Ndocc*Nrhom1+a*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @Tc(-1)(a,b,_P_,c0) @W2(c0,i,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a0
          }//P
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a3,a2,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
    }
    orz::DTensor W2(_Ndocc,_Ndocc,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc,_Nact*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(c0*_Ndocc+i)*tmp1.shape(1)+(a2*_Nact+a3)] = V2_FiC_aacc.cptr()[a2*V2_FiC_aacc.shape(1)*V2_FiC_aacc.shape(2)*V2_FiC_aacc.shape(3)+a3*V2_FiC_aacc.shape(2)*V2_FiC_aacc.shape(3)+c0*V2_FiC_aacc.shape(3)+i];
            }//a3
          }//a2
        }//i
      }//c0
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a2*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a2
      // @W2(c0,i,P,_P_) <<= +1 @V2_FiC_aacc(a2,a3,c0,i) @W1(a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//i
      }//c0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom1*_Ndocc);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(a*_Nvirt+b)*tmp1.shape(1)+(_P_*_Ndocc+c0)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//c0
          }//_P_
        }//b
      }//a
      orz::DTensor tmp2(Nrhom1*_Ndocc,_Ndocc*Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(_P_*_Ndocc+c0)*tmp2.shape(1)+(i*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+i*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//i
        }//c0
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(a,b,_P_,c0) @W2(c0,i,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Nvirt*_Ndocc*Nrhom1+b*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @W1(P,_P_) @W2(i,a,b,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt,_Nvirt);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          tmp1.cptr()[(a)*tmp1.shape(1)+(v0)] = CFab.cptr()[v0*CFab.shape(1)+a];
        }//v0
      }//a
      orz::DTensor tmp2(_Nvirt,_Nvirt*Nrhom1*_Ndocc);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&i : irange(0L, _Ndocc)){
              tmp2.cptr()[(v0)*tmp2.shape(1)+(b*Nrhom1*_Ndocc+_P_*_Ndocc+i)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//i
          }//_P_
        }//b
      }//v0
      // @W2(i,a,b,_P_) <<= +1 @Fcore(v0,a) @Tc(-1)(v0,b,_P_,i)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&i : irange(0L, _Ndocc)){
              W2.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+_P_] += tmp3.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+_P_*_Ndocc+i];
            }//i
          }//_P_
        }//b
      }//a
    }
    {
      orz::DTensor tmp1(Nrhom1,Nrhom1);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(_P_)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//_P_
      }//P
      orz::DTensor tmp2(Nrhom1,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+a*_Nvirt+b)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+a*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+_P_];
            }//b
          }//a
        }//i
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(P,_P_) @W2(i,a,b,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[P*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+a*_Nvirt+b];
            }//b
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @W1(P,_P_) @W2(i,b,a,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt,_Nvirt);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          tmp1.cptr()[(b)*tmp1.shape(1)+(v0)] = CFab.cptr()[v0*CFab.shape(1)+b];
        }//v0
      }//b
      orz::DTensor tmp2(_Nvirt,_Nvirt*Nrhom1*_Ndocc);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&i : irange(0L, _Ndocc)){
              tmp2.cptr()[(v0)*tmp2.shape(1)+(a*Nrhom1*_Ndocc+_P_*_Ndocc+i)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//i
          }//_P_
        }//a
      }//v0
      // @W2(i,b,a,_P_) <<= +1 @Fcore(v0,b) @Tc(-1)(v0,a,_P_,i)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&i : irange(0L, _Ndocc)){
              W2.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+_P_] += tmp3.cptr()[b*_Nvirt*Nrhom1*_Ndocc+a*Nrhom1*_Ndocc+_P_*_Ndocc+i];
            }//i
          }//_P_
        }//a
      }//b
    }
    {
      orz::DTensor tmp1(Nrhom1,Nrhom1);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(_P_)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//_P_
      }//P
      orz::DTensor tmp2(Nrhom1,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+b*_Nvirt+a)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+b*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+_P_];
            }//a
          }//b
        }//i
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(P,_P_) @W2(i,b,a,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[P*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+b*_Nvirt+a];
            }//a
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @W1(_P_,P) @W2(i,a,b,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a0,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+_P_] += tmp3.cptr()[a0*Nrhom1+_P_];
        }//_P_
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(_P_)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+_P_];
        }//a0
      }//_P_
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @W1(_P_,P) <<= +1 @W0(a0,_P_) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          W1.cptr()[_P_*Nrhom1+P] += tmp3.cptr()[_P_*Nrhom1+P];
        }//P
      }//_P_
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt,_Nvirt);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          tmp1.cptr()[(a)*tmp1.shape(1)+(v0)] = CFab.cptr()[v0*CFab.shape(1)+a];
        }//v0
      }//a
      orz::DTensor tmp2(_Nvirt,_Nvirt*Nrhom1*_Ndocc);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&i : irange(0L, _Ndocc)){
              tmp2.cptr()[(v0)*tmp2.shape(1)+(b*Nrhom1*_Ndocc+_P_*_Ndocc+i)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//i
          }//_P_
        }//b
      }//v0
      // @W2(i,a,b,_P_) <<= +1 @Fcore(v0,a) @Tc(-1)(b,v0,_P_,i)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&i : irange(0L, _Ndocc)){
              W2.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+_P_] += tmp3.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+_P_*_Ndocc+i];
            }//i
          }//_P_
        }//b
      }//a
    }
    {
      orz::DTensor tmp1(Nrhom1,Nrhom1);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(_P_)] = W1.cptr()[_P_*W1.shape(1)+P];
        }//_P_
      }//P
      orz::DTensor tmp2(Nrhom1,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+a*_Nvirt+b)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+a*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+_P_];
            }//b
          }//a
        }//i
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(_P_,P) @W2(i,a,b,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[P*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+a*_Nvirt+b];
            }//b
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @W1(P,_P_) @W2(i,b,a,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt,_Nvirt);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          tmp1.cptr()[(b)*tmp1.shape(1)+(v0)] = CFab.cptr()[v0*CFab.shape(1)+b];
        }//v0
      }//b
      orz::DTensor tmp2(_Nvirt,_Nvirt*Nrhom1*_Ndocc);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&i : irange(0L, _Ndocc)){
              tmp2.cptr()[(v0)*tmp2.shape(1)+(a*Nrhom1*_Ndocc+_P_*_Ndocc+i)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//i
          }//_P_
        }//a
      }//v0
      // @W2(i,b,a,_P_) <<= +1 @Fcore(v0,b) @Tc(-1)(a,v0,_P_,i)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&i : irange(0L, _Ndocc)){
              W2.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+_P_] += tmp3.cptr()[b*_Nvirt*Nrhom1*_Ndocc+a*Nrhom1*_Ndocc+_P_*_Ndocc+i];
            }//i
          }//_P_
        }//a
      }//b
    }
    {
      orz::DTensor tmp1(Nrhom1,Nrhom1);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(_P_)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//_P_
      }//P
      orz::DTensor tmp2(Nrhom1,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+b*_Nvirt+a)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+b*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+_P_];
            }//a
          }//b
        }//i
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(P,_P_) @W2(i,b,a,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[P*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+b*_Nvirt+a];
            }//a
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @W1(P,_P_) @W2(i,b,a,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt,_Nvirt);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          tmp1.cptr()[(b)*tmp1.shape(1)+(v0)] = CFab.cptr()[v0*CFab.shape(1)+b];
        }//v0
      }//b
      orz::DTensor tmp2(_Nvirt,_Nvirt*Nrhom1*_Ndocc);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&i : irange(0L, _Ndocc)){
              tmp2.cptr()[(v0)*tmp2.shape(1)+(a*Nrhom1*_Ndocc+_P_*_Ndocc+i)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//i
          }//_P_
        }//a
      }//v0
      // @W2(i,b,a,_P_) <<= +1 @Fcore(v0,b) @Tc(-1)(a,v0,_P_,i)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&i : irange(0L, _Ndocc)){
              W2.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+_P_] += tmp3.cptr()[b*_Nvirt*Nrhom1*_Ndocc+a*Nrhom1*_Ndocc+_P_*_Ndocc+i];
            }//i
          }//_P_
        }//a
      }//b
    }
    {
      orz::DTensor tmp1(Nrhom1,Nrhom1);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(_P_)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//_P_
      }//P
      orz::DTensor tmp2(Nrhom1,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+b*_Nvirt+a)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+b*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+_P_];
            }//a
          }//b
        }//i
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(P,_P_) @W2(i,b,a,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[P*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+b*_Nvirt+a];
            }//a
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @W1(P,_P_) @W2(i,a,b,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt,_Nvirt);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          tmp1.cptr()[(a)*tmp1.shape(1)+(v0)] = CFab.cptr()[v0*CFab.shape(1)+a];
        }//v0
      }//a
      orz::DTensor tmp2(_Nvirt,_Nvirt*Nrhom1*_Ndocc);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&i : irange(0L, _Ndocc)){
              tmp2.cptr()[(v0)*tmp2.shape(1)+(b*Nrhom1*_Ndocc+_P_*_Ndocc+i)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//i
          }//_P_
        }//b
      }//v0
      // @W2(i,a,b,_P_) <<= +1 @Fcore(v0,a) @Tc(-1)(b,v0,_P_,i)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&i : irange(0L, _Ndocc)){
              W2.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+_P_] += tmp3.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+_P_*_Ndocc+i];
            }//i
          }//_P_
        }//b
      }//a
    }
    {
      orz::DTensor tmp1(Nrhom1,Nrhom1);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(_P_)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//_P_
      }//P
      orz::DTensor tmp2(Nrhom1,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+a*_Nvirt+b)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+a*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+_P_];
            }//b
          }//a
        }//i
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(P,_P_) @W2(i,a,b,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[P*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+a*_Nvirt+b];
            }//b
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @W1(_P_,P) @W2(i,b,a,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a0,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+_P_] += tmp3.cptr()[a0*Nrhom1+_P_];
        }//_P_
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(_P_)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+_P_];
        }//a0
      }//_P_
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @W1(_P_,P) <<= +1 @W0(a0,_P_) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          W1.cptr()[_P_*Nrhom1+P] += tmp3.cptr()[_P_*Nrhom1+P];
        }//P
      }//_P_
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt,_Nvirt);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          tmp1.cptr()[(b)*tmp1.shape(1)+(v0)] = CFab.cptr()[v0*CFab.shape(1)+b];
        }//v0
      }//b
      orz::DTensor tmp2(_Nvirt,_Nvirt*Nrhom1*_Ndocc);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&i : irange(0L, _Ndocc)){
              tmp2.cptr()[(v0)*tmp2.shape(1)+(a*Nrhom1*_Ndocc+_P_*_Ndocc+i)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//i
          }//_P_
        }//a
      }//v0
      // @W2(i,b,a,_P_) <<= +1 @Fcore(v0,b) @Tc(-1)(v0,a,_P_,i)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&i : irange(0L, _Ndocc)){
              W2.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+_P_] += tmp3.cptr()[b*_Nvirt*Nrhom1*_Ndocc+a*Nrhom1*_Ndocc+_P_*_Ndocc+i];
            }//i
          }//_P_
        }//a
      }//b
    }
    {
      orz::DTensor tmp1(Nrhom1,Nrhom1);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(_P_)] = W1.cptr()[_P_*W1.shape(1)+P];
        }//_P_
      }//P
      orz::DTensor tmp2(Nrhom1,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+b*_Nvirt+a)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+b*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+_P_];
            }//a
          }//b
        }//i
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(_P_,P) @W2(i,b,a,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[P*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+b*_Nvirt+a];
            }//a
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @W1(P,_P_) @W2(i,a,b,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt,_Nvirt);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          tmp1.cptr()[(a)*tmp1.shape(1)+(v0)] = CFab.cptr()[v0*CFab.shape(1)+a];
        }//v0
      }//a
      orz::DTensor tmp2(_Nvirt,_Nvirt*Nrhom1*_Ndocc);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&i : irange(0L, _Ndocc)){
              tmp2.cptr()[(v0)*tmp2.shape(1)+(b*Nrhom1*_Ndocc+_P_*_Ndocc+i)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//i
          }//_P_
        }//b
      }//v0
      // @W2(i,a,b,_P_) <<= +1 @Fcore(v0,a) @Tc(-1)(v0,b,_P_,i)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&i : irange(0L, _Ndocc)){
              W2.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+_P_] += tmp3.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+_P_*_Ndocc+i];
            }//i
          }//_P_
        }//b
      }//a
    }
    {
      orz::DTensor tmp1(Nrhom1,Nrhom1);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(_P_)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//_P_
      }//P
      orz::DTensor tmp2(Nrhom1,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+a*_Nvirt+b)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+a*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+_P_];
            }//b
          }//a
        }//i
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(P,_P_) @W2(i,a,b,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[P*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+a*_Nvirt+b];
            }//b
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @Tc(-1)(a,b,_P_,i) @W2(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a0
          }//P
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a3,a2,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
    }
    orz::DTensor W2(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(1,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          tmp1.cptr()[(0)*tmp1.shape(1)+(a3*_Nact+a2)] = CFpq.cptr()[a3*CFpq.shape(1)+a2];
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a3*_Nact+a2)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
      // @W2(P,_P_) <<= +1 @Fcore(a3,a2) @W1(a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W2.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Nvirt*_Ndocc+b*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//b
      }//a
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W2.cptr()[P*W2.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(a,b,_P_,i) @W2(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Nvirt*_Ndocc*Nrhom1+b*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @Tc(-1)(b,a,_P_,i) @W2(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a0
          }//P
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a3,a2,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
    }
    orz::DTensor W2(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(1,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          tmp1.cptr()[(0)*tmp1.shape(1)+(a3*_Nact+a2)] = CFpq.cptr()[a3*CFpq.shape(1)+a2];
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a3*_Nact+a2)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
      // @W2(P,_P_) <<= +1 @Fcore(a3,a2) @W1(a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W2.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Nvirt*_Ndocc+a*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//a
      }//b
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W2.cptr()[P*W2.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(b,a,_P_,i) @W2(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Nvirt*_Ndocc*Nrhom1+a*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @Tc(-1)(b,a,_P_,i) @W2(_P_,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+_P_)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a0
          }//_P_
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @W1(a3,a2,_P_,P) <<= +1 @W0(a3,a2,a0,_P_) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+_P_*Nrhom1+P] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+_P_*Nrhom1+P];
            }//P
          }//_P_
        }//a2
      }//a3
    }
    orz::DTensor W2(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(1,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          tmp1.cptr()[(0)*tmp1.shape(1)+(a3*_Nact+a2)] = CFpq.cptr()[a3*CFpq.shape(1)+a2];
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(a3*_Nact+a2)*tmp2.shape(1)+(_P_*Nrhom1+P)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+_P_*W1.shape(3)+P];
            }//P
          }//_P_
        }//a2
      }//a3
      // @W2(_P_,P) <<= +1 @Fcore(a3,a2) @W1(a3,a2,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          W2.cptr()[_P_*Nrhom1+P] += tmp3.cptr()[_P_*Nrhom1+P];
        }//P
      }//_P_
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Nvirt*_Ndocc+a*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//a
      }//b
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W2.cptr()[_P_*W2.shape(1)+P];
        }//P
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(b,a,_P_,i) @W2(_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Nvirt*_Ndocc*Nrhom1+a*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @Tc(-1)(a,b,_P_,i) @W2(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a0
          }//P
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a3,a2,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
    }
    orz::DTensor W2(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(1,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          tmp1.cptr()[(0)*tmp1.shape(1)+(a3*_Nact+a2)] = CFpq.cptr()[a3*CFpq.shape(1)+a2];
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a3*_Nact+a2)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
      // @W2(P,_P_) <<= +1 @Fcore(a3,a2) @W1(a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W2.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Nvirt*_Ndocc+b*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//b
      }//a
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W2.cptr()[P*W2.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(a,b,_P_,i) @W2(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Nvirt*_Ndocc*Nrhom1+b*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @W1(P,_P_) @W2(i,a,b,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc,_Ndocc);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          tmp1.cptr()[(i)*tmp1.shape(1)+(c0)] = CFij.cptr()[i*CFij.shape(1)+c0];
        }//c0
      }//i
      orz::DTensor tmp2(_Ndocc,_Nvirt*_Nvirt*Nrhom1);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(c0)*tmp2.shape(1)+(a*_Nvirt*Nrhom1+b*Nrhom1+_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//b
        }//a
      }//c0
      // @W2(i,a,b,_P_) <<= +1 @Fcore(i,c0) @Tc(-1)(a,b,_P_,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+_P_] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+_P_];
            }//_P_
          }//b
        }//a
      }//i
    }
    {
      orz::DTensor tmp1(Nrhom1,Nrhom1);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(_P_)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//_P_
      }//P
      orz::DTensor tmp2(Nrhom1,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+a*_Nvirt+b)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+a*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+_P_];
            }//b
          }//a
        }//i
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(P,_P_) @W2(i,a,b,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[P*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+a*_Nvirt+b];
            }//b
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @W1(P,_P_) @W2(i,b,a,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc,_Ndocc);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          tmp1.cptr()[(i)*tmp1.shape(1)+(c0)] = CFij.cptr()[i*CFij.shape(1)+c0];
        }//c0
      }//i
      orz::DTensor tmp2(_Ndocc,_Nvirt*_Nvirt*Nrhom1);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(c0)*tmp2.shape(1)+(b*_Nvirt*Nrhom1+a*Nrhom1+_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//a
        }//b
      }//c0
      // @W2(i,b,a,_P_) <<= +1 @Fcore(i,c0) @Tc(-1)(b,a,_P_,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+_P_] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+_P_];
            }//_P_
          }//a
        }//b
      }//i
    }
    {
      orz::DTensor tmp1(Nrhom1,Nrhom1);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(_P_)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//_P_
      }//P
      orz::DTensor tmp2(Nrhom1,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+b*_Nvirt+a)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+b*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+_P_];
            }//a
          }//b
        }//i
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(P,_P_) @W2(i,b,a,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[P*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+b*_Nvirt+a];
            }//a
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @W1(_P_,P) @W2(i,b,a,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a0,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+_P_] += tmp3.cptr()[a0*Nrhom1+_P_];
        }//_P_
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(_P_)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+_P_];
        }//a0
      }//_P_
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @W1(_P_,P) <<= +1 @W0(a0,_P_) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          W1.cptr()[_P_*Nrhom1+P] += tmp3.cptr()[_P_*Nrhom1+P];
        }//P
      }//_P_
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc,_Ndocc);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          tmp1.cptr()[(i)*tmp1.shape(1)+(c0)] = CFij.cptr()[i*CFij.shape(1)+c0];
        }//c0
      }//i
      orz::DTensor tmp2(_Ndocc,_Nvirt*_Nvirt*Nrhom1);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(c0)*tmp2.shape(1)+(b*_Nvirt*Nrhom1+a*Nrhom1+_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//a
        }//b
      }//c0
      // @W2(i,b,a,_P_) <<= +1 @Fcore(i,c0) @Tc(-1)(b,a,_P_,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+_P_] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+_P_];
            }//_P_
          }//a
        }//b
      }//i
    }
    {
      orz::DTensor tmp1(Nrhom1,Nrhom1);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(_P_)] = W1.cptr()[_P_*W1.shape(1)+P];
        }//_P_
      }//P
      orz::DTensor tmp2(Nrhom1,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+b*_Nvirt+a)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+b*W2.shape(2)*W2.shape(3)+a*W2.shape(3)+_P_];
            }//a
          }//b
        }//i
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(_P_,P) @W2(i,b,a,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[P*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+b*_Nvirt+a];
            }//a
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @W1(P,_P_) @W2(i,a,b,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    orz::DTensor W2(_Ndocc,_Nvirt,_Nvirt,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc,_Ndocc);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          tmp1.cptr()[(i)*tmp1.shape(1)+(c0)] = CFij.cptr()[i*CFij.shape(1)+c0];
        }//c0
      }//i
      orz::DTensor tmp2(_Ndocc,_Nvirt*_Nvirt*Nrhom1);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(c0)*tmp2.shape(1)+(a*_Nvirt*Nrhom1+b*Nrhom1+_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//b
        }//a
      }//c0
      // @W2(i,a,b,_P_) <<= +1 @Fcore(i,c0) @Tc(-1)(a,b,_P_,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+_P_] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+_P_];
            }//_P_
          }//b
        }//a
      }//i
    }
    {
      orz::DTensor tmp1(Nrhom1,Nrhom1);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(_P_)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//_P_
      }//P
      orz::DTensor tmp2(Nrhom1,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+a*_Nvirt+b)] = W2.cptr()[i*W2.shape(1)*W2.shape(2)*W2.shape(3)+a*W2.shape(2)*W2.shape(3)+b*W2.shape(3)+_P_];
            }//b
          }//a
        }//i
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(P,_P_) @W2(i,a,b,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[P*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+a*_Nvirt+b];
            }//b
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 Ecore @Tc(-1)(a,b,_P_,i) @W1(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000 * Ecore;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Nvirt*_Ndocc+b*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//b
      }//a
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(a,b,_P_,i) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Nvirt*_Ndocc*Nrhom1+b*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 Ecore @Tc(-1)(b,a,_P_,i) @W1(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000 * Ecore;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Nvirt*_Ndocc+a*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//a
      }//b
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(b,a,_P_,i) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Nvirt*_Ndocc*Nrhom1+a*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 Ecore @Tc(-1)(b,a,_P_,i) @W1(_P_,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000 * Ecore;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a1];
        }//_P_
      }//a1
      // @W0(a0,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+_P_] += tmp3.cptr()[a0*Nrhom1+_P_];
        }//_P_
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(_P_)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+_P_];
        }//a0
      }//_P_
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @W1(_P_,P) <<= +1 @W0(a0,_P_) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          W1.cptr()[_P_*Nrhom1+P] += tmp3.cptr()[_P_*Nrhom1+P];
        }//P
      }//_P_
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Nvirt*_Ndocc+a*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//a
      }//b
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[_P_*W1.shape(1)+P];
        }//P
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(b,a,_P_,i) @W1(_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Nvirt*_Ndocc*Nrhom1+a*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 Ecore @Tc(-1)(a,b,_P_,i) @W1(P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000 * Ecore;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a0,P) <<= +1 @D1(a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          W0.cptr()[a0*Nrhom1+P] += tmp3.cptr()[a0*Nrhom1+P];
        }//P
      }//a0
    }
    orz::DTensor W1(Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(Nrhom1,_Nact);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(P)*tmp1.shape(1)+(a0)] = W0.cptr()[a0*W0.shape(1)+P];
        }//a0
      }//P
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(P,_P_) <<= +1 @W0(a0,P) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[P*Nrhom1+_P_] += tmp3.cptr()[P*Nrhom1+_P_];
        }//_P_
      }//P
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Nvirt*_Ndocc+b*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//b
      }//a
      orz::DTensor tmp2(Nrhom1,Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(P)] = W1.cptr()[P*W1.shape(1)+_P_];
        }//P
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @Tc(-1)(a,b,_P_,i) @W1(P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Nvirt*_Ndocc*Nrhom1+b*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @T(-1)(a,b,_P_,c0) @W2(c0,i,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a0
          }//P
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a3,a2,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
    }
    orz::DTensor W2(_Ndocc,_Ndocc,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc,_Nact*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(c0*_Ndocc+i)*tmp1.shape(1)+(a2*_Nact+a3)] = V2_FiC_acac.cptr()[a2*V2_FiC_acac.shape(1)*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+c0*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+a3*V2_FiC_acac.shape(3)+i];
            }//a3
          }//a2
        }//i
      }//c0
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a2*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a2
      // @W2(c0,i,P,_P_) <<= +1 @V2_FiC_acac(a2,c0,a3,i) @W1(a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//i
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
      orz::DTensor tmp2(Nrhom1*_Ndocc,_Ndocc*Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(_P_*_Ndocc+c0)*tmp2.shape(1)+(i*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+i*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//i
        }//c0
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @T(-1)(a,b,_P_,c0) @W2(c0,i,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Nvirt*_Ndocc*Nrhom1+b*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @T(-1)(b,a,_P_,c0) @W2(c0,i,_P_,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a2*_Nact+a1)*tmp1.shape(1)+(a0)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a3,a2,a1,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a1*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a1*Nrhom1+_P_];
            }//_P_
          }//a1
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+_P_)*tmp1.shape(1)+(a1)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a1*W0.shape(3)+_P_];
            }//a1
          }//_P_
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W1(a3,a2,_P_,P) <<= +1 @W0(a3,a2,a1,_P_) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+_P_*Nrhom1+P] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+_P_*Nrhom1+P];
            }//P
          }//_P_
        }//a2
      }//a3
    }
    orz::DTensor W2(_Ndocc,_Ndocc,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc,_Nact*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(c0*_Ndocc+i)*tmp1.shape(1)+(a2*_Nact+a3)] = V2_FiC_acac.cptr()[a2*V2_FiC_acac.shape(1)*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+c0*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+a3*V2_FiC_acac.shape(3)+i];
            }//a3
          }//a2
        }//i
      }//c0
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(a2*_Nact+a3)*tmp2.shape(1)+(_P_*Nrhom1+P)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+_P_*W1.shape(3)+P];
            }//P
          }//_P_
        }//a3
      }//a2
      // @W2(c0,i,_P_,P) <<= +1 @V2_FiC_acac(a2,c0,a3,i) @W1(a3,a2,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+_P_*Nrhom1+P] += tmp3.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+_P_*Nrhom1+P];
            }//P
          }//_P_
        }//i
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
      orz::DTensor tmp2(Nrhom1*_Ndocc,_Ndocc*Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(_P_*_Ndocc+c0)*tmp2.shape(1)+(i*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+i*W2.shape(2)*W2.shape(3)+_P_*W2.shape(3)+P];
            }//P
          }//i
        }//c0
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @T(-1)(b,a,_P_,c0) @W2(c0,i,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Nvirt*_Ndocc*Nrhom1+a*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @T(-1)(b,a,_P_,c0) @W2(c0,i,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a0
          }//P
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a3,a2,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
    }
    orz::DTensor W2(_Ndocc,_Ndocc,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc,_Nact*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(c0*_Ndocc+i)*tmp1.shape(1)+(a2*_Nact+a3)] = V2_FiC_acac.cptr()[a2*V2_FiC_acac.shape(1)*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+c0*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+a3*V2_FiC_acac.shape(3)+i];
            }//a3
          }//a2
        }//i
      }//c0
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a2*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a2
      // @W2(c0,i,P,_P_) <<= +1 @V2_FiC_acac(a2,c0,a3,i) @W1(a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//i
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
      orz::DTensor tmp2(Nrhom1*_Ndocc,_Ndocc*Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(_P_*_Ndocc+c0)*tmp2.shape(1)+(i*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+i*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//i
        }//c0
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @T(-1)(b,a,_P_,c0) @W2(c0,i,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Nvirt*_Ndocc*Nrhom1+a*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @T(-1)(a,b,_P_,c0) @W2(c0,i,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a0
          }//P
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a3,a2,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
    }
    orz::DTensor W2(_Ndocc,_Ndocc,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc,_Nact*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(c0*_Ndocc+i)*tmp1.shape(1)+(a2*_Nact+a3)] = V2_FiC_acac.cptr()[a2*V2_FiC_acac.shape(1)*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+c0*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+a3*V2_FiC_acac.shape(3)+i];
            }//a3
          }//a2
        }//i
      }//c0
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a2*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a2
      // @W2(c0,i,P,_P_) <<= +1 @V2_FiC_acac(a2,c0,a3,i) @W1(a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//i
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
      orz::DTensor tmp2(Nrhom1*_Ndocc,_Ndocc*Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(_P_*_Ndocc+c0)*tmp2.shape(1)+(i*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+i*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//i
        }//c0
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @T(-1)(a,b,_P_,c0) @W2(c0,i,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Nvirt*_Ndocc*Nrhom1+b*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @T(-1)(b,a,_P_,c0) @W2(c0,i,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a0*Nrhom1+P)*tmp1.shape(1)+(a2)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a2
          }//P
        }//a0
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a2)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a2];
        }//_P_
      }//a2
      // @W1(a3,a0,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(+)(_P_,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a0
      }//a3
    }
    orz::DTensor W2(_Ndocc,_Ndocc,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc,_Nact*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(c0*_Ndocc+i)*tmp1.shape(1)+(a0*_Nact+a3)] = V2_FiC_acac.cptr()[a0*V2_FiC_acac.shape(1)*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+c0*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+a3*V2_FiC_acac.shape(3)+i];
            }//a3
          }//a0
        }//i
      }//c0
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a0*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a0
      // @W2(c0,i,P,_P_) <<= +1 @V2_FiC_acac(a0,c0,a3,i) @W1(a3,a0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//i
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
      orz::DTensor tmp2(Nrhom1*_Ndocc,_Ndocc*Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(_P_*_Ndocc+c0)*tmp2.shape(1)+(i*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+i*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//i
        }//c0
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @T(-1)(b,a,_P_,c0) @W2(c0,i,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Nvirt*_Ndocc*Nrhom1+a*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @T(-1)(a,b,_P_,c0) @W2(c0,i,_P_,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a1*_Nact+a0)*tmp1.shape(1)+(a2)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a2
          }//a0
        }//a1
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a2)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a2];
        }//_P_
      }//a2
      // @W0(a3,a1,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(_P_,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a1*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a1*_Nact*Nrhom1+a0*Nrhom1+_P_];
            }//_P_
          }//a0
        }//a1
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a0*Nrhom1+_P_)*tmp1.shape(1)+(a1)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a1
          }//_P_
        }//a0
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W1(a3,a0,_P_,P) <<= +1 @W0(a3,a1,a0,_P_) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+_P_*Nrhom1+P] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+_P_*Nrhom1+P];
            }//P
          }//_P_
        }//a0
      }//a3
    }
    orz::DTensor W2(_Ndocc,_Ndocc,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc,_Nact*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(c0*_Ndocc+i)*tmp1.shape(1)+(a0*_Nact+a3)] = V2_FiC_acac.cptr()[a0*V2_FiC_acac.shape(1)*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+c0*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+a3*V2_FiC_acac.shape(3)+i];
            }//a3
          }//a0
        }//i
      }//c0
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(a0*_Nact+a3)*tmp2.shape(1)+(_P_*Nrhom1+P)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+_P_*W1.shape(3)+P];
            }//P
          }//_P_
        }//a3
      }//a0
      // @W2(c0,i,_P_,P) <<= +1 @V2_FiC_acac(a0,c0,a3,i) @W1(a3,a0,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+_P_*Nrhom1+P] += tmp3.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+_P_*Nrhom1+P];
            }//P
          }//_P_
        }//i
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
      orz::DTensor tmp2(Nrhom1*_Ndocc,_Ndocc*Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(_P_*_Ndocc+c0)*tmp2.shape(1)+(i*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+i*W2.shape(2)*W2.shape(3)+_P_*W2.shape(3)+P];
            }//P
          }//i
        }//c0
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @T(-1)(a,b,_P_,c0) @W2(c0,i,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Nvirt*_Ndocc*Nrhom1+b*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @T(-1)(a,b,_P_,c0) @W2(c0,i,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a0*Nrhom1+P)*tmp1.shape(1)+(a2)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a2
          }//P
        }//a0
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a2)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a2];
        }//_P_
      }//a2
      // @W1(a3,a0,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(-)(_P_,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a0
      }//a3
    }
    orz::DTensor W2(_Ndocc,_Ndocc,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc,_Nact*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(c0*_Ndocc+i)*tmp1.shape(1)+(a0*_Nact+a3)] = V2_FiC_acac.cptr()[a0*V2_FiC_acac.shape(1)*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+c0*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+a3*V2_FiC_acac.shape(3)+i];
            }//a3
          }//a0
        }//i
      }//c0
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a0*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a0
      // @W2(c0,i,P,_P_) <<= +1 @V2_FiC_acac(a0,c0,a3,i) @W1(a3,a0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//i
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
      orz::DTensor tmp2(Nrhom1*_Ndocc,_Ndocc*Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(_P_*_Ndocc+c0)*tmp2.shape(1)+(i*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+i*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//i
        }//c0
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @T(-1)(a,b,_P_,c0) @W2(c0,i,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Nvirt*_Ndocc*Nrhom1+b*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= +1 @T(-1)(b,a,_P_,c0) @W2(c0,i,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a0*Nrhom1+P)*tmp1.shape(1)+(a2)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a2
          }//P
        }//a0
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a2)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a2];
        }//_P_
      }//a2
      // @W1(a3,a0,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(-)(_P_,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a0
      }//a3
    }
    orz::DTensor W2(_Ndocc,_Ndocc,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc,_Nact*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(c0*_Ndocc+i)*tmp1.shape(1)+(a0*_Nact+a3)] = V2_FiC_acac.cptr()[a0*V2_FiC_acac.shape(1)*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+c0*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+a3*V2_FiC_acac.shape(3)+i];
            }//a3
          }//a0
        }//i
      }//c0
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a0*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a0
      // @W2(c0,i,P,_P_) <<= +1 @V2_FiC_acac(a0,c0,a3,i) @W1(a3,a0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[c0*_Ndocc*Nrhom1*Nrhom1+i*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//i
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
      orz::DTensor tmp2(Nrhom1*_Ndocc,_Ndocc*Nrhom1);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(_P_*_Ndocc+c0)*tmp2.shape(1)+(i*Nrhom1+P)] = W2.cptr()[c0*W2.shape(1)*W2.shape(2)*W2.shape(3)+i*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//i
        }//c0
      }//_P_
      // @tRFiC(a,b,P,i) <<= +1 _X_ @T(-1)(b,a,_P_,c0) @W2(c0,i,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Nvirt*_Ndocc*Nrhom1+a*_Ndocc*Nrhom1+i*Nrhom1+P];
            }//P
          }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @T(-1)(v0,a,_P_,i) @W2(b,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a0*Nrhom1+P)*tmp1.shape(1)+(a2)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a2
          }//P
        }//a0
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a2)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a2];
        }//_P_
      }//a2
      // @W1(a3,a0,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(+)(_P_,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a0
      }//a3
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Nvirt+v0)*tmp1.shape(1)+(a0*_Nact+a3)] = V2_FiC_vava.cptr()[b*V2_FiC_vava.shape(1)*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+a0*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+v0*V2_FiC_vava.shape(3)+a3];
            }//a3
          }//a0
        }//v0
      }//b
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a0*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a0
      // @W2(b,v0,P,_P_) <<= +1 @V2_FiC_vava(b,a0,v0,a3) @W1(a3,a0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[b*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[b*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*Nrhom1+_P_)] = T_m1_.cptr()[v0*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+a*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+i];
            }//_P_
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*Nrhom1,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*Nrhom1+_P_)*tmp2.shape(1)+(b*Nrhom1+P)] = W2.cptr()[b*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//b
        }//_P_
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @T(-1)(v0,a,_P_,i) @W2(b,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @T(-1)(v0,b,_P_,i) @W2(a,v0,_P_,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a1*_Nact+a0)*tmp1.shape(1)+(a2)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a2
          }//a0
        }//a1
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a2)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a2];
        }//_P_
      }//a2
      // @W0(a3,a1,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(_P_,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a1*_Nact*Nrhom1+a0*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a1*_Nact*Nrhom1+a0*Nrhom1+_P_];
            }//_P_
          }//a0
        }//a1
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a0*Nrhom1+_P_)*tmp1.shape(1)+(a1)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//a1
          }//_P_
        }//a0
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W1(a3,a0,_P_,P) <<= +1 @W0(a3,a1,a0,_P_) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+_P_*Nrhom1+P] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+_P_*Nrhom1+P];
            }//P
          }//_P_
        }//a0
      }//a3
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Nvirt+v0)*tmp1.shape(1)+(a0*_Nact+a3)] = V2_FiC_vava.cptr()[a*V2_FiC_vava.shape(1)*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+a0*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+v0*V2_FiC_vava.shape(3)+a3];
            }//a3
          }//a0
        }//v0
      }//a
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(a0*_Nact+a3)*tmp2.shape(1)+(_P_*Nrhom1+P)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+_P_*W1.shape(3)+P];
            }//P
          }//_P_
        }//a3
      }//a0
      // @W2(a,v0,_P_,P) <<= +1 @V2_FiC_vava(a,a0,v0,a3) @W1(a3,a0,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[a*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+_P_*Nrhom1+P] += tmp3.cptr()[a*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+_P_*Nrhom1+P];
            }//P
          }//_P_
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*Nrhom1+_P_)] = T_m1_.cptr()[v0*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+b*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+i];
            }//_P_
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*Nrhom1,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*Nrhom1+_P_)*tmp2.shape(1)+(a*Nrhom1+P)] = W2.cptr()[a*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+_P_*W2.shape(3)+P];
            }//P
          }//a
        }//_P_
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @T(-1)(v0,b,_P_,i) @W2(a,v0,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @T(-1)(a,v0,_P_,i) @W2(b,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a0*Nrhom1+P)*tmp1.shape(1)+(a2)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a2
          }//P
        }//a0
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a2)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a2];
        }//_P_
      }//a2
      // @W1(a3,a0,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(-)(_P_,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a0
      }//a3
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Nvirt+v0)*tmp1.shape(1)+(a0*_Nact+a3)] = V2_FiC_vava.cptr()[b*V2_FiC_vava.shape(1)*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+a0*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+v0*V2_FiC_vava.shape(3)+a3];
            }//a3
          }//a0
        }//v0
      }//b
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a0*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a0
      // @W2(b,v0,P,_P_) <<= +1 @V2_FiC_vava(b,a0,v0,a3) @W1(a3,a0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[b*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[b*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*Nrhom1+_P_)] = T_m1_.cptr()[a*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+v0*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+i];
            }//_P_
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*Nrhom1,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*Nrhom1+_P_)*tmp2.shape(1)+(b*Nrhom1+P)] = W2.cptr()[b*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//b
        }//_P_
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @T(-1)(a,v0,_P_,i) @W2(b,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @T(-1)(b,v0,_P_,i) @W2(a,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a0*Nrhom1+P)*tmp1.shape(1)+(a2)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a2
          }//P
        }//a0
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a2)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a2];
        }//_P_
      }//a2
      // @W1(a3,a0,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(-)(_P_,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a0
      }//a3
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Nvirt+v0)*tmp1.shape(1)+(a0*_Nact+a3)] = V2_FiC_vava.cptr()[a*V2_FiC_vava.shape(1)*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+a0*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+v0*V2_FiC_vava.shape(3)+a3];
            }//a3
          }//a0
        }//v0
      }//a
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a0*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a0
      // @W2(a,v0,P,_P_) <<= +1 @V2_FiC_vava(a,a0,v0,a3) @W1(a3,a0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[a*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*Nrhom1+_P_)] = T_m1_.cptr()[b*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+v0*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+i];
            }//_P_
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*Nrhom1,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*Nrhom1+_P_)*tmp2.shape(1)+(a*Nrhom1+P)] = W2.cptr()[a*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//a
        }//_P_
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @T(-1)(b,v0,_P_,i) @W2(a,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @T(-1)(a,v0,_P_,i) @W2(b,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a0
          }//P
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a3,a2,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Nvirt+v0)*tmp1.shape(1)+(a2*_Nact+a3)] = V2_FiC_vava.cptr()[b*V2_FiC_vava.shape(1)*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+a2*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+v0*V2_FiC_vava.shape(3)+a3];
            }//a3
          }//a2
        }//v0
      }//b
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a2*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a2
      // @W2(b,v0,P,_P_) <<= +1 @V2_FiC_vava(b,a2,v0,a3) @W1(a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[b*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[b*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*Nrhom1+_P_)] = T_m1_.cptr()[a*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+v0*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+i];
            }//_P_
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*Nrhom1,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*Nrhom1+_P_)*tmp2.shape(1)+(b*Nrhom1+P)] = W2.cptr()[b*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//b
        }//_P_
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @T(-1)(a,v0,_P_,i) @W2(b,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @T(-1)(b,v0,_P_,i) @W2(a,v0,_P_,P)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,_Nact,_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*_Nact,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*_Nact+a2*_Nact+a1)*tmp1.shape(1)+(a0)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a3,a2,a1,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a1*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a1*Nrhom1+_P_];
            }//_P_
          }//a1
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+_P_)*tmp1.shape(1)+(a1)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a1*W0.shape(3)+_P_];
            }//a1
          }//_P_
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W1(a3,a2,_P_,P) <<= +1 @W0(a3,a2,a1,_P_) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+_P_*Nrhom1+P] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+_P_*Nrhom1+P];
            }//P
          }//_P_
        }//a2
      }//a3
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Nvirt+v0)*tmp1.shape(1)+(a2*_Nact+a3)] = V2_FiC_vava.cptr()[a*V2_FiC_vava.shape(1)*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+a2*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+v0*V2_FiC_vava.shape(3)+a3];
            }//a3
          }//a2
        }//v0
      }//a
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(a2*_Nact+a3)*tmp2.shape(1)+(_P_*Nrhom1+P)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+_P_*W1.shape(3)+P];
            }//P
          }//_P_
        }//a3
      }//a2
      // @W2(a,v0,_P_,P) <<= +1 @V2_FiC_vava(a,a2,v0,a3) @W1(a3,a2,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom1)){
            for(auto &&P : irange(0L, Nrhom1)){
              W2.cptr()[a*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+_P_*Nrhom1+P] += tmp3.cptr()[a*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+_P_*Nrhom1+P];
            }//P
          }//_P_
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*Nrhom1+_P_)] = T_m1_.cptr()[b*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+v0*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+i];
            }//_P_
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*Nrhom1,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*Nrhom1+_P_)*tmp2.shape(1)+(a*Nrhom1+P)] = W2.cptr()[a*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+_P_*W2.shape(3)+P];
            }//P
          }//a
        }//_P_
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @T(-1)(b,v0,_P_,i) @W2(a,v0,_P_,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @T(-1)(v0,a,_P_,i) @W2(b,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(+)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a0
          }//P
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a3,a2,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Nvirt+v0)*tmp1.shape(1)+(a2*_Nact+a3)] = V2_FiC_vava.cptr()[b*V2_FiC_vava.shape(1)*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+a2*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+v0*V2_FiC_vava.shape(3)+a3];
            }//a3
          }//a2
        }//v0
      }//b
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a2*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a2
      // @W2(b,v0,P,_P_) <<= +1 @V2_FiC_vava(b,a2,v0,a3) @W1(a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[b*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[b*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*Nrhom1+_P_)] = T_m1_.cptr()[v0*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+a*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+i];
            }//_P_
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*Nrhom1,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*Nrhom1+_P_)*tmp2.shape(1)+(b*Nrhom1+P)] = W2.cptr()[b*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//b
        }//_P_
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @T(-1)(v0,a,_P_,i) @W2(b,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[a*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//i
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
std::cout << boost::format("%80s") % "@tRFiC(a,b,P,i) <<= -1 @T(-1)(v0,b,_P_,i) @W2(a,v0,P,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
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
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a1];
        }//P
      }//a1
      // @W0(a3,a2,a0,P) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(-)(P,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&P : irange(0L, Nrhom1)){
              W0.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P] += tmp3.cptr()[a3*_Nact*_Nact*Nrhom1+a2*_Nact*Nrhom1+a0*Nrhom1+P];
            }//P
          }//a0
        }//a2
      }//a3
    }
    orz::DTensor W1(_Nact,_Nact,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nact*_Nact*Nrhom1,_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact*Nrhom1+a2*Nrhom1+P)*tmp1.shape(1)+(a0)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a2*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+P];
            }//a0
          }//P
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W1(a3,a2,P,_P_) <<= +1 @W0(a3,a2,a0,P) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W1.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a3*_Nact*Nrhom1*Nrhom1+a2*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//a2
      }//a3
    }
    orz::DTensor W2(_Nvirt,_Nvirt,Nrhom1,Nrhom1);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Nact*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Nvirt+v0)*tmp1.shape(1)+(a2*_Nact+a3)] = V2_FiC_vava.cptr()[a*V2_FiC_vava.shape(1)*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+a2*V2_FiC_vava.shape(2)*V2_FiC_vava.shape(3)+v0*V2_FiC_vava.shape(3)+a3];
            }//a3
          }//a2
        }//v0
      }//a
      orz::DTensor tmp2(_Nact*_Nact,Nrhom1*Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a2*_Nact+a3)*tmp2.shape(1)+(P*Nrhom1+_P_)] = W1.cptr()[a3*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+P*W1.shape(3)+_P_];
            }//_P_
          }//P
        }//a3
      }//a2
      // @W2(a,v0,P,_P_) <<= +1 @V2_FiC_vava(a,a2,v0,a3) @W1(a3,a2,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom1)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              W2.cptr()[a*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_] += tmp3.cptr()[a*_Nvirt*Nrhom1*Nrhom1+v0*Nrhom1*Nrhom1+P*Nrhom1+_P_];
            }//_P_
          }//P
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*Nrhom1+_P_)] = T_m1_.cptr()[v0*T_m1_.shape(1)*T_m1_.shape(2)*T_m1_.shape(3)+b*T_m1_.shape(2)*T_m1_.shape(3)+_P_*T_m1_.shape(3)+i];
            }//_P_
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*Nrhom1,_Nvirt*Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tmp2.cptr()[(v0*Nrhom1+_P_)*tmp2.shape(1)+(a*Nrhom1+P)] = W2.cptr()[a*W2.shape(1)*W2.shape(2)*W2.shape(3)+v0*W2.shape(2)*W2.shape(3)+P*W2.shape(3)+_P_];
            }//P
          }//a
        }//_P_
      }//v0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @T(-1)(v0,b,_P_,i) @W2(a,v0,P,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[b*_Ndocc*_Nvirt*Nrhom1+i*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//i
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
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          const long I_Pi = RActList.GetIndex(P,i);
          if (I_Pi < 0) continue;
          orz::DTensor RP(_Nvirt, _Nvirt);
          for(auto &&a : irange(0L, _Nvirt))
            for(auto &&b : irange(0L, _Nvirt))
              RP(a,b) = tRFiC(a,b,P,i);
          orz::DTensor Rorig = R2.GetTensor(P,i);
          RActList.GetPair(P,i).Transform2EXT(RP);
          RP += Rorig;
          R2.PutTensor(P,i, RP);
        }
      }
      
      return;
    }//End func
  } }
