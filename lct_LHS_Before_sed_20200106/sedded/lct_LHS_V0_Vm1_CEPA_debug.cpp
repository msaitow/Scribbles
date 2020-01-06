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
    void FMatrixOrthogonalizer::FormLHS_V0_Vm1_CEd(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &c2, const orz::DTensor &d3, DensityPack &dPack,
                                                   TensorContainer &ovl, TensorContainer &PNO4, const bool do4EXTinPNO,
                                                   const double QE0, const double Ecas, const double h0_energy, TensorContainer &DebugInts, const orz::DTensor &casfock, const orz::DTensor &FV,
                                                   TensorContainer &T2, TensorContainer &R2){
      //#include "lct_preparations.cpp"
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
  orz::DTensor tRFiC(_Ndocc,_Ndocc,_Nvirt,_Nvirt);
  {

#ifdef _PROFILE_LHS
double t_trans = 0;
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= -1 @V2_FiC_vvac(a,v0,a1,j) @W1(i,a1,v0,b)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(v0*_Nvirt*_Ndocc+b*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//b
      }//v0
      orz::DTensor tmp2(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//a1
      }//_P_
      // @W1(i,a1,v0,b) <<= +1 @Tc(-1)(v0,b,_P_,i) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&a1 : irange(0L, _Nact)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+v0*_Nvirt+b] += tmp3.cptr()[v0*_Nvirt*_Ndocc*_Nact+b*_Ndocc*_Nact+i*_Nact+a1];
            }//a1
          }//i
        }//b
      }//v0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Ndocc+j)*tmp1.shape(1)+(v0*_Nact+a1)] = V2_FiC_vvac.cptr()[a*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a1*V2_FiC_vvac.shape(3)+j];
            }//a1
          }//v0
        }//j
      }//a
      orz::DTensor tmp2(_Nvirt*_Nact,_Ndocc*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a1)*tmp2.shape(1)+(i*_Nvirt+b)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a1*W1.shape(2)*W1.shape(3)+v0*W1.shape(3)+b];
            }//b
          }//i
        }//a1
      }//v0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @V2_FiC_vvac(a,v0,a1,j) @W1(i,a1,v0,b)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&b : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Ndocc*_Ndocc*_Nvirt+j*_Ndocc*_Nvirt+i*_Nvirt+b];
            }//b
          }//i
        }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= -1 @V2_FiC_vvac(a,v0,a1,j) @W1(i,a1,b,v0)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Nvirt*_Ndocc+v0*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//v0
      }//b
      orz::DTensor tmp2(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//a1
      }//_P_
      // @W1(i,a1,b,v0) <<= +1 @Tc(-1)(b,v0,_P_,i) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&a1 : irange(0L, _Nact)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+b*_Nvirt+v0] += tmp3.cptr()[b*_Nvirt*_Ndocc*_Nact+v0*_Ndocc*_Nact+i*_Nact+a1];
            }//a1
          }//i
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Ndocc+j)*tmp1.shape(1)+(v0*_Nact+a1)] = V2_FiC_vvac.cptr()[a*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a1*V2_FiC_vvac.shape(3)+j];
            }//a1
          }//v0
        }//j
      }//a
      orz::DTensor tmp2(_Nvirt*_Nact,_Ndocc*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a1)*tmp2.shape(1)+(i*_Nvirt+b)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a1*W1.shape(2)*W1.shape(3)+b*W1.shape(3)+v0];
            }//b
          }//i
        }//a1
      }//v0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @V2_FiC_vvac(a,v0,a1,j) @W1(i,a1,b,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&b : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Ndocc*_Ndocc*_Nvirt+j*_Ndocc*_Nvirt+i*_Nvirt+b];
            }//b
          }//i
        }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= -1 @V2_FiC_vvac(b,v0,a1,i) @W1(j,a1,v0,a)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(v0*_Nvirt*_Ndocc+a*_Ndocc+j)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+j];
            }//_P_
          }//j
        }//a
      }//v0
      orz::DTensor tmp2(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//a1
      }//_P_
      // @W1(j,a1,v0,a) <<= +1 @Tc(-1)(v0,a,_P_,j) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&a1 : irange(0L, _Nact)){
              W1.cptr()[j*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+v0*_Nvirt+a] += tmp3.cptr()[v0*_Nvirt*_Ndocc*_Nact+a*_Ndocc*_Nact+j*_Nact+a1];
            }//a1
          }//j
        }//a
      }//v0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a1)] = V2_FiC_vvac.cptr()[b*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a1*V2_FiC_vvac.shape(3)+i];
            }//a1
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*_Nact,_Ndocc*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a1)*tmp2.shape(1)+(j*_Nvirt+a)] = W1.cptr()[j*W1.shape(1)*W1.shape(2)*W1.shape(3)+a1*W1.shape(2)*W1.shape(3)+v0*W1.shape(3)+a];
            }//a
          }//j
        }//a1
      }//v0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @V2_FiC_vvac(b,v0,a1,i) @W1(j,a1,v0,a)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&a : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[b*_Ndocc*_Ndocc*_Nvirt+i*_Ndocc*_Nvirt+j*_Nvirt+a];
            }//a
          }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= -1 @V2_FiC_vvac(b,v0,a1,i) @W1(j,a1,a,v0)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Nvirt*_Ndocc+v0*_Ndocc+j)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+j];
            }//_P_
          }//j
        }//v0
      }//a
      orz::DTensor tmp2(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//a1
      }//_P_
      // @W1(j,a1,a,v0) <<= +1 @Tc(-1)(a,v0,_P_,j) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&a1 : irange(0L, _Nact)){
              W1.cptr()[j*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+a*_Nvirt+v0] += tmp3.cptr()[a*_Nvirt*_Ndocc*_Nact+v0*_Ndocc*_Nact+j*_Nact+a1];
            }//a1
          }//j
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a1)] = V2_FiC_vvac.cptr()[b*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a1*V2_FiC_vvac.shape(3)+i];
            }//a1
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*_Nact,_Ndocc*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a1)*tmp2.shape(1)+(j*_Nvirt+a)] = W1.cptr()[j*W1.shape(1)*W1.shape(2)*W1.shape(3)+a1*W1.shape(2)*W1.shape(3)+a*W1.shape(3)+v0];
            }//a
          }//j
        }//a1
      }//v0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @V2_FiC_vvac(b,v0,a1,i) @W1(j,a1,a,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&a : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[b*_Ndocc*_Ndocc*_Nvirt+i*_Ndocc*_Nvirt+j*_Nvirt+a];
            }//a
          }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= -1 @V2_FiC_vvac(a,v0,a1,i) @W1(j,a1,b,v0)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Nvirt*_Ndocc+v0*_Ndocc+j)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+j];
            }//_P_
          }//j
        }//v0
      }//b
      orz::DTensor tmp2(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//a1
      }//_P_
      // @W1(j,a1,b,v0) <<= +1 @Tc(-1)(b,v0,_P_,j) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&a1 : irange(0L, _Nact)){
              W1.cptr()[j*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+b*_Nvirt+v0] += tmp3.cptr()[b*_Nvirt*_Ndocc*_Nact+v0*_Ndocc*_Nact+j*_Nact+a1];
            }//a1
          }//j
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a1)] = V2_FiC_vvac.cptr()[a*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a1*V2_FiC_vvac.shape(3)+i];
            }//a1
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*_Nact,_Ndocc*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a1)*tmp2.shape(1)+(j*_Nvirt+b)] = W1.cptr()[j*W1.shape(1)*W1.shape(2)*W1.shape(3)+a1*W1.shape(2)*W1.shape(3)+b*W1.shape(3)+v0];
            }//b
          }//j
        }//a1
      }//v0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @V2_FiC_vvac(a,v0,a1,i) @W1(j,a1,b,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&b : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Ndocc*_Ndocc*_Nvirt+i*_Ndocc*_Nvirt+j*_Nvirt+b];
            }//b
          }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= -1 @V2_FiC_vvac(a,v0,a1,i) @W1(j,a1,v0,b)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(v0*_Nvirt*_Ndocc+b*_Ndocc+j)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+j];
            }//_P_
          }//j
        }//b
      }//v0
      orz::DTensor tmp2(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//a1
      }//_P_
      // @W1(j,a1,v0,b) <<= +1 @Tc(-1)(v0,b,_P_,j) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&a1 : irange(0L, _Nact)){
              W1.cptr()[j*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+v0*_Nvirt+b] += tmp3.cptr()[v0*_Nvirt*_Ndocc*_Nact+b*_Ndocc*_Nact+j*_Nact+a1];
            }//a1
          }//j
        }//b
      }//v0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a1)] = V2_FiC_vvac.cptr()[a*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a1*V2_FiC_vvac.shape(3)+i];
            }//a1
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*_Nact,_Ndocc*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a1)*tmp2.shape(1)+(j*_Nvirt+b)] = W1.cptr()[j*W1.shape(1)*W1.shape(2)*W1.shape(3)+a1*W1.shape(2)*W1.shape(3)+v0*W1.shape(3)+b];
            }//b
          }//j
        }//a1
      }//v0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @V2_FiC_vvac(a,v0,a1,i) @W1(j,a1,v0,b)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&b : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Ndocc*_Ndocc*_Nvirt+i*_Ndocc*_Nvirt+j*_Nvirt+b];
            }//b
          }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= -1 @V2_FiC_vvac(b,v0,a1,j) @W1(i,a1,a,v0)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Nvirt*_Ndocc+v0*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//v0
      }//a
      orz::DTensor tmp2(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//a1
      }//_P_
      // @W1(i,a1,a,v0) <<= +1 @Tc(-1)(a,v0,_P_,i) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&a1 : irange(0L, _Nact)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+a*_Nvirt+v0] += tmp3.cptr()[a*_Nvirt*_Ndocc*_Nact+v0*_Ndocc*_Nact+i*_Nact+a1];
            }//a1
          }//i
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Ndocc+j)*tmp1.shape(1)+(v0*_Nact+a1)] = V2_FiC_vvac.cptr()[b*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a1*V2_FiC_vvac.shape(3)+j];
            }//a1
          }//v0
        }//j
      }//b
      orz::DTensor tmp2(_Nvirt*_Nact,_Ndocc*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a1)*tmp2.shape(1)+(i*_Nvirt+a)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a1*W1.shape(2)*W1.shape(3)+a*W1.shape(3)+v0];
            }//a
          }//i
        }//a1
      }//v0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @V2_FiC_vvac(b,v0,a1,j) @W1(i,a1,a,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&a : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[b*_Ndocc*_Ndocc*_Nvirt+j*_Ndocc*_Nvirt+i*_Nvirt+a];
            }//a
          }//i
        }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= -1 @V2_FiC_vvac(b,v0,a1,j) @W1(i,a1,v0,a)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(v0*_Nvirt*_Ndocc+a*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//a
      }//v0
      orz::DTensor tmp2(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//a1
      }//_P_
      // @W1(i,a1,v0,a) <<= +1 @Tc(-1)(v0,a,_P_,i) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&a1 : irange(0L, _Nact)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+v0*_Nvirt+a] += tmp3.cptr()[v0*_Nvirt*_Ndocc*_Nact+a*_Ndocc*_Nact+i*_Nact+a1];
            }//a1
          }//i
        }//a
      }//v0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Ndocc+j)*tmp1.shape(1)+(v0*_Nact+a1)] = V2_FiC_vvac.cptr()[b*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a1*V2_FiC_vvac.shape(3)+j];
            }//a1
          }//v0
        }//j
      }//b
      orz::DTensor tmp2(_Nvirt*_Nact,_Ndocc*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a1)*tmp2.shape(1)+(i*_Nvirt+a)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a1*W1.shape(2)*W1.shape(3)+v0*W1.shape(3)+a];
            }//a
          }//i
        }//a1
      }//v0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @V2_FiC_vvac(b,v0,a1,j) @W1(i,a1,v0,a)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&a : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[b*_Ndocc*_Ndocc*_Nvirt+j*_Ndocc*_Nvirt+i*_Nvirt+a];
            }//a
          }//i
        }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= -1 @V2_FiC_vavc(v0,a1,b,i) @W1(j,a1,a,v0)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Nvirt*_Ndocc+v0*_Ndocc+j)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+j];
            }//_P_
          }//j
        }//v0
      }//a
      orz::DTensor tmp2(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//a1
      }//_P_
      // @W1(j,a1,a,v0) <<= +1 @Tc(-1)(a,v0,_P_,j) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&a1 : irange(0L, _Nact)){
              W1.cptr()[j*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+a*_Nvirt+v0] += tmp3.cptr()[a*_Nvirt*_Ndocc*_Nact+v0*_Ndocc*_Nact+j*_Nact+a1];
            }//a1
          }//j
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a1)] = V2_FiC_vavc.cptr()[v0*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a1*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+b*V2_FiC_vavc.shape(3)+i];
            }//a1
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*_Nact,_Ndocc*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a1)*tmp2.shape(1)+(j*_Nvirt+a)] = W1.cptr()[j*W1.shape(1)*W1.shape(2)*W1.shape(3)+a1*W1.shape(2)*W1.shape(3)+a*W1.shape(3)+v0];
            }//a
          }//j
        }//a1
      }//v0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @V2_FiC_vavc(v0,a1,b,i) @W1(j,a1,a,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&a : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[b*_Ndocc*_Ndocc*_Nvirt+i*_Ndocc*_Nvirt+j*_Nvirt+a];
            }//a
          }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= -1 @V2_FiC_vavc(v0,a1,b,i) @W1(j,a1,v0,a)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(v0*_Nvirt*_Ndocc+a*_Ndocc+j)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+j];
            }//_P_
          }//j
        }//a
      }//v0
      orz::DTensor tmp2(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//a1
      }//_P_
      // @W1(j,a1,v0,a) <<= +1 @Tc(-1)(v0,a,_P_,j) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&a1 : irange(0L, _Nact)){
              W1.cptr()[j*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+v0*_Nvirt+a] += tmp3.cptr()[v0*_Nvirt*_Ndocc*_Nact+a*_Ndocc*_Nact+j*_Nact+a1];
            }//a1
          }//j
        }//a
      }//v0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a1)] = V2_FiC_vavc.cptr()[v0*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a1*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+b*V2_FiC_vavc.shape(3)+i];
            }//a1
          }//v0
        }//i
      }//b
      orz::DTensor tmp2(_Nvirt*_Nact,_Ndocc*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a1)*tmp2.shape(1)+(j*_Nvirt+a)] = W1.cptr()[j*W1.shape(1)*W1.shape(2)*W1.shape(3)+a1*W1.shape(2)*W1.shape(3)+v0*W1.shape(3)+a];
            }//a
          }//j
        }//a1
      }//v0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @V2_FiC_vavc(v0,a1,b,i) @W1(j,a1,v0,a)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&a : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[b*_Ndocc*_Ndocc*_Nvirt+i*_Ndocc*_Nvirt+j*_Nvirt+a];
            }//a
          }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= +2 @V2_FiC_vavc(v0,a1,b,j) @W1(i,a1,a,v0)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   2.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Nvirt*_Ndocc+v0*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//v0
      }//a
      orz::DTensor tmp2(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//a1
      }//_P_
      // @W1(i,a1,a,v0) <<= +1 @Tc(-1)(a,v0,_P_,i) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&a1 : irange(0L, _Nact)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+a*_Nvirt+v0] += tmp3.cptr()[a*_Nvirt*_Ndocc*_Nact+v0*_Ndocc*_Nact+i*_Nact+a1];
            }//a1
          }//i
        }//v0
      }//a
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Ndocc+j)*tmp1.shape(1)+(v0*_Nact+a1)] = V2_FiC_vavc.cptr()[v0*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a1*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+b*V2_FiC_vavc.shape(3)+j];
            }//a1
          }//v0
        }//j
      }//b
      orz::DTensor tmp2(_Nvirt*_Nact,_Ndocc*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a1)*tmp2.shape(1)+(i*_Nvirt+a)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a1*W1.shape(2)*W1.shape(3)+a*W1.shape(3)+v0];
            }//a
          }//i
        }//a1
      }//v0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @V2_FiC_vavc(v0,a1,b,j) @W1(i,a1,a,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&a : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[b*_Ndocc*_Ndocc*_Nvirt+j*_Ndocc*_Nvirt+i*_Nvirt+a];
            }//a
          }//i
        }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= +2 @V2_FiC_vavc(v0,a1,b,j) @W1(i,a1,v0,a)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   2.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(v0*_Nvirt*_Ndocc+a*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//a
      }//v0
      orz::DTensor tmp2(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//a1
      }//_P_
      // @W1(i,a1,v0,a) <<= +1 @Tc(-1)(v0,a,_P_,i) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&a1 : irange(0L, _Nact)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+v0*_Nvirt+a] += tmp3.cptr()[v0*_Nvirt*_Ndocc*_Nact+a*_Ndocc*_Nact+i*_Nact+a1];
            }//a1
          }//i
        }//a
      }//v0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(b*_Ndocc+j)*tmp1.shape(1)+(v0*_Nact+a1)] = V2_FiC_vavc.cptr()[v0*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a1*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+b*V2_FiC_vavc.shape(3)+j];
            }//a1
          }//v0
        }//j
      }//b
      orz::DTensor tmp2(_Nvirt*_Nact,_Ndocc*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a1)*tmp2.shape(1)+(i*_Nvirt+a)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a1*W1.shape(2)*W1.shape(3)+v0*W1.shape(3)+a];
            }//a
          }//i
        }//a1
      }//v0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @V2_FiC_vavc(v0,a1,b,j) @W1(i,a1,v0,a)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&a : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[b*_Ndocc*_Ndocc*_Nvirt+j*_Ndocc*_Nvirt+i*_Nvirt+a];
            }//a
          }//i
        }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= -1 @V2_FiC_vavc(v0,a1,a,j) @W1(i,a1,b,v0)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Nvirt*_Ndocc+v0*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//v0
      }//b
      orz::DTensor tmp2(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//a1
      }//_P_
      // @W1(i,a1,b,v0) <<= +1 @Tc(-1)(b,v0,_P_,i) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&a1 : irange(0L, _Nact)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+b*_Nvirt+v0] += tmp3.cptr()[b*_Nvirt*_Ndocc*_Nact+v0*_Ndocc*_Nact+i*_Nact+a1];
            }//a1
          }//i
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Ndocc+j)*tmp1.shape(1)+(v0*_Nact+a1)] = V2_FiC_vavc.cptr()[v0*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a1*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a*V2_FiC_vavc.shape(3)+j];
            }//a1
          }//v0
        }//j
      }//a
      orz::DTensor tmp2(_Nvirt*_Nact,_Ndocc*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a1)*tmp2.shape(1)+(i*_Nvirt+b)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a1*W1.shape(2)*W1.shape(3)+b*W1.shape(3)+v0];
            }//b
          }//i
        }//a1
      }//v0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @V2_FiC_vavc(v0,a1,a,j) @W1(i,a1,b,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&b : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Ndocc*_Ndocc*_Nvirt+j*_Ndocc*_Nvirt+i*_Nvirt+b];
            }//b
          }//i
        }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= -1 @V2_FiC_vavc(v0,a1,a,j) @W1(i,a1,v0,b)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(v0*_Nvirt*_Ndocc+b*_Ndocc+i)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+i];
            }//_P_
          }//i
        }//b
      }//v0
      orz::DTensor tmp2(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//a1
      }//_P_
      // @W1(i,a1,v0,b) <<= +1 @Tc(-1)(v0,b,_P_,i) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&a1 : irange(0L, _Nact)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+v0*_Nvirt+b] += tmp3.cptr()[v0*_Nvirt*_Ndocc*_Nact+b*_Ndocc*_Nact+i*_Nact+a1];
            }//a1
          }//i
        }//b
      }//v0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Ndocc+j)*tmp1.shape(1)+(v0*_Nact+a1)] = V2_FiC_vavc.cptr()[v0*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a1*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a*V2_FiC_vavc.shape(3)+j];
            }//a1
          }//v0
        }//j
      }//a
      orz::DTensor tmp2(_Nvirt*_Nact,_Ndocc*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a1)*tmp2.shape(1)+(i*_Nvirt+b)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a1*W1.shape(2)*W1.shape(3)+v0*W1.shape(3)+b];
            }//b
          }//i
        }//a1
      }//v0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @V2_FiC_vavc(v0,a1,a,j) @W1(i,a1,v0,b)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&b : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Ndocc*_Ndocc*_Nvirt+j*_Ndocc*_Nvirt+i*_Nvirt+b];
            }//b
          }//i
        }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= +2 @V2_FiC_vavc(v0,a1,a,i) @W1(j,a1,b,v0)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   2.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Nvirt*_Ndocc+v0*_Ndocc+j)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+v0*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+j];
            }//_P_
          }//j
        }//v0
      }//b
      orz::DTensor tmp2(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//a1
      }//_P_
      // @W1(j,a1,b,v0) <<= +1 @Tc(-1)(b,v0,_P_,j) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&a1 : irange(0L, _Nact)){
              W1.cptr()[j*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+b*_Nvirt+v0] += tmp3.cptr()[b*_Nvirt*_Ndocc*_Nact+v0*_Ndocc*_Nact+j*_Nact+a1];
            }//a1
          }//j
        }//v0
      }//b
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a1)] = V2_FiC_vavc.cptr()[v0*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a1*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a*V2_FiC_vavc.shape(3)+i];
            }//a1
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*_Nact,_Ndocc*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a1)*tmp2.shape(1)+(j*_Nvirt+b)] = W1.cptr()[j*W1.shape(1)*W1.shape(2)*W1.shape(3)+a1*W1.shape(2)*W1.shape(3)+b*W1.shape(3)+v0];
            }//b
          }//j
        }//a1
      }//v0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @V2_FiC_vavc(v0,a1,a,i) @W1(j,a1,b,v0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&b : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Ndocc*_Ndocc*_Nvirt+i*_Ndocc*_Nvirt+j*_Nvirt+b];
            }//b
          }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= +2 @V2_FiC_vavc(v0,a1,a,i) @W1(j,a1,v0,b)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   2.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(v0*_Nvirt*_Ndocc+b*_Ndocc+j)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[v0*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+j];
            }//_P_
          }//j
        }//b
      }//v0
      orz::DTensor tmp2(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//a1
      }//_P_
      // @W1(j,a1,v0,b) <<= +1 @Tc(-1)(v0,b,_P_,j) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&a1 : irange(0L, _Nact)){
              W1.cptr()[j*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+v0*_Nvirt+b] += tmp3.cptr()[v0*_Nvirt*_Ndocc*_Nact+b*_Ndocc*_Nact+j*_Nact+a1];
            }//a1
          }//j
        }//b
      }//v0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Ndocc,_Nvirt*_Nact);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v0 : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a*_Ndocc+i)*tmp1.shape(1)+(v0*_Nact+a1)] = V2_FiC_vavc.cptr()[v0*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a1*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a*V2_FiC_vavc.shape(3)+i];
            }//a1
          }//v0
        }//i
      }//a
      orz::DTensor tmp2(_Nvirt*_Nact,_Ndocc*_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(v0*_Nact+a1)*tmp2.shape(1)+(j*_Nvirt+b)] = W1.cptr()[j*W1.shape(1)*W1.shape(2)*W1.shape(3)+a1*W1.shape(2)*W1.shape(3)+v0*W1.shape(3)+b];
            }//b
          }//j
        }//a1
      }//v0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @V2_FiC_vavc(v0,a1,a,i) @W1(j,a1,v0,b)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&b : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Ndocc*_Ndocc*_Nvirt+i*_Ndocc*_Nvirt+j*_Nvirt+b];
            }//b
          }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= -1 @Tc(-1)(a,b,_P_,i) @W1(j,_P_)";
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
          tmp2.cptr()[(a2)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a2];
        }//_P_
      }//a2
      // @W0(a3,a1,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(-)(_P_,a2)
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
    orz::DTensor W1(_Ndocc,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc,_Nact*_Nact*_Nact);
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(j)*tmp1.shape(1)+(a0*_Nact*_Nact+a1*_Nact+a3)] = V2_FiC_aaac.cptr()[a0*V2_FiC_aaac.shape(1)*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a1*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a3*V2_FiC_aaac.shape(3)+j];
            }//a3
          }//a1
        }//a0
      }//j
      orz::DTensor tmp2(_Nact*_Nact*_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a0*_Nact*_Nact+a1*_Nact+a3)*tmp2.shape(1)+(_P_)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//_P_
          }//a3
        }//a1
      }//a0
      // @W1(j,_P_) <<= +1 @V2_FiC_aaac(a0,a1,a3,j) @W0(a3,a1,a0,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[j*Nrhom1+_P_] += tmp3.cptr()[j*Nrhom1+_P_];
        }//_P_
      }//j
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
      orz::DTensor tmp2(Nrhom1,_Ndocc);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&j : irange(0L, _Ndocc)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(j)] = W1.cptr()[j*W1.shape(1)+_P_];
        }//j
      }//_P_
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(-1)(a,b,_P_,i) @W1(j,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&j : irange(0L, _Ndocc)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Nvirt*_Ndocc*_Ndocc+b*_Ndocc*_Ndocc+i*_Ndocc+j];
            }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= -1 @Tc(-1)(b,a,_P_,i) @W1(j,_P_)";
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
    orz::DTensor W1(_Ndocc,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc,_Nact*_Nact*_Nact);
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(j)*tmp1.shape(1)+(a0*_Nact*_Nact+a1*_Nact+a3)] = V2_FiC_aaac.cptr()[a0*V2_FiC_aaac.shape(1)*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a1*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a3*V2_FiC_aaac.shape(3)+j];
            }//a3
          }//a1
        }//a0
      }//j
      orz::DTensor tmp2(_Nact*_Nact*_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a0*_Nact*_Nact+a1*_Nact+a3)*tmp2.shape(1)+(_P_)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//_P_
          }//a3
        }//a1
      }//a0
      // @W1(j,_P_) <<= +1 @V2_FiC_aaac(a0,a1,a3,j) @W0(a3,a1,a0,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[j*Nrhom1+_P_] += tmp3.cptr()[j*Nrhom1+_P_];
        }//_P_
      }//j
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
      orz::DTensor tmp2(Nrhom1,_Ndocc);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&j : irange(0L, _Ndocc)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(j)] = W1.cptr()[j*W1.shape(1)+_P_];
        }//j
      }//_P_
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(-1)(b,a,_P_,i) @W1(j,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&j : irange(0L, _Ndocc)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[b*_Nvirt*_Ndocc*_Ndocc+a*_Ndocc*_Ndocc+i*_Ndocc+j];
            }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= +1 @V2_FiC_accc(a1,j,c0,i) @W1(c0,a1,a,b)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Nvirt*_Ndocc+b*_Ndocc+c0)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//c0
        }//b
      }//a
      orz::DTensor tmp2(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//a1
      }//_P_
      // @W1(c0,a1,a,b) <<= +1 @Tc(-1)(a,b,_P_,c0) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&a1 : irange(0L, _Nact)){
              W1.cptr()[c0*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Nvirt*_Ndocc*_Nact+b*_Ndocc*_Nact+c0*_Nact+a1];
            }//a1
          }//c0
        }//b
      }//a
    }
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc,_Nact*_Ndocc);
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(j*_Ndocc+i)*tmp1.shape(1)+(a1*_Ndocc+c0)] = V2_FiC_accc.cptr()[a1*V2_FiC_accc.shape(1)*V2_FiC_accc.shape(2)*V2_FiC_accc.shape(3)+j*V2_FiC_accc.shape(2)*V2_FiC_accc.shape(3)+c0*V2_FiC_accc.shape(3)+i];
            }//c0
          }//a1
        }//i
      }//j
      orz::DTensor tmp2(_Nact*_Ndocc,_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1*_Ndocc+c0)*tmp2.shape(1)+(a*_Nvirt+b)] = W1.cptr()[c0*W1.shape(1)*W1.shape(2)*W1.shape(3)+a1*W1.shape(2)*W1.shape(3)+a*W1.shape(3)+b];
            }//b
          }//a
        }//c0
      }//a1
      // @tRFiC(i,j,a,b) <<= +1 _X_ @V2_FiC_accc(a1,j,c0,i) @W1(c0,a1,a,b)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[j*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+a*_Nvirt+b];
            }//b
          }//a
        }//i
      }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= +1 @V2_FiC_accc(a1,j,c0,i) @W1(c0,a1,b,a)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Nvirt*_Ndocc+a*_Ndocc+c0)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//c0
        }//a
      }//b
      orz::DTensor tmp2(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//a1
      }//_P_
      // @W1(c0,a1,b,a) <<= +1 @Tc(-1)(b,a,_P_,c0) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&a1 : irange(0L, _Nact)){
              W1.cptr()[c0*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[b*_Nvirt*_Ndocc*_Nact+a*_Ndocc*_Nact+c0*_Nact+a1];
            }//a1
          }//c0
        }//a
      }//b
    }
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc,_Nact*_Ndocc);
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(j*_Ndocc+i)*tmp1.shape(1)+(a1*_Ndocc+c0)] = V2_FiC_accc.cptr()[a1*V2_FiC_accc.shape(1)*V2_FiC_accc.shape(2)*V2_FiC_accc.shape(3)+j*V2_FiC_accc.shape(2)*V2_FiC_accc.shape(3)+c0*V2_FiC_accc.shape(3)+i];
            }//c0
          }//a1
        }//i
      }//j
      orz::DTensor tmp2(_Nact*_Ndocc,_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1*_Ndocc+c0)*tmp2.shape(1)+(b*_Nvirt+a)] = W1.cptr()[c0*W1.shape(1)*W1.shape(2)*W1.shape(3)+a1*W1.shape(2)*W1.shape(3)+b*W1.shape(3)+a];
            }//a
          }//b
        }//c0
      }//a1
      // @tRFiC(i,j,a,b) <<= +1 _X_ @V2_FiC_accc(a1,j,c0,i) @W1(c0,a1,b,a)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[j*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+b*_Nvirt+a];
            }//a
          }//b
        }//i
      }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= -1 @Tc(-1)(b,a,_P_,j) @W1(i,_P_)";
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
          tmp2.cptr()[(a2)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a2];
        }//_P_
      }//a2
      // @W0(a3,a1,a0,_P_) <<= +1 @D2(a3,a2,a1,a0) @X(-1)(-)(_P_,a2)
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
    orz::DTensor W1(_Ndocc,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc,_Nact*_Nact*_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(i)*tmp1.shape(1)+(a0*_Nact*_Nact+a1*_Nact+a3)] = V2_FiC_aaac.cptr()[a0*V2_FiC_aaac.shape(1)*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a1*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a3*V2_FiC_aaac.shape(3)+i];
            }//a3
          }//a1
        }//a0
      }//i
      orz::DTensor tmp2(_Nact*_Nact*_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a0*_Nact*_Nact+a1*_Nact+a3)*tmp2.shape(1)+(_P_)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//_P_
          }//a3
        }//a1
      }//a0
      // @W1(i,_P_) <<= +1 @V2_FiC_aaac(a0,a1,a3,i) @W0(a3,a1,a0,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[i*Nrhom1+_P_] += tmp3.cptr()[i*Nrhom1+_P_];
        }//_P_
      }//i
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Nvirt*_Ndocc+a*_Ndocc+j)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+j];
            }//_P_
          }//j
        }//a
      }//b
      orz::DTensor tmp2(Nrhom1,_Ndocc);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(i)] = W1.cptr()[i*W1.shape(1)+_P_];
        }//i
      }//_P_
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(-1)(b,a,_P_,j) @W1(i,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&i : irange(0L, _Ndocc)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[b*_Nvirt*_Ndocc*_Ndocc+a*_Ndocc*_Ndocc+j*_Ndocc+i];
            }//i
          }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= -1 @Tc(-1)(a,b,_P_,j) @W1(i,_P_)";
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
    orz::DTensor W1(_Ndocc,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc,_Nact*_Nact*_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a3 : irange(0L, _Nact)){
              tmp1.cptr()[(i)*tmp1.shape(1)+(a0*_Nact*_Nact+a1*_Nact+a3)] = V2_FiC_aaac.cptr()[a0*V2_FiC_aaac.shape(1)*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a1*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a3*V2_FiC_aaac.shape(3)+i];
            }//a3
          }//a1
        }//a0
      }//i
      orz::DTensor tmp2(_Nact*_Nact*_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp2.cptr()[(a0*_Nact*_Nact+a1*_Nact+a3)*tmp2.shape(1)+(_P_)] = W0.cptr()[a3*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+a0*W0.shape(3)+_P_];
            }//_P_
          }//a3
        }//a1
      }//a0
      // @W1(i,_P_) <<= +1 @V2_FiC_aaac(a0,a1,a3,i) @W0(a3,a1,a0,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[i*Nrhom1+_P_] += tmp3.cptr()[i*Nrhom1+_P_];
        }//_P_
      }//i
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Nvirt*_Ndocc+b*_Ndocc+j)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+j];
            }//_P_
          }//j
        }//b
      }//a
      orz::DTensor tmp2(Nrhom1,_Ndocc);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(i)] = W1.cptr()[i*W1.shape(1)+_P_];
        }//i
      }//_P_
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(-1)(a,b,_P_,j) @W1(i,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&i : irange(0L, _Ndocc)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Nvirt*_Ndocc*_Ndocc+b*_Ndocc*_Ndocc+j*_Ndocc+i];
            }//i
          }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= +1 @V2_FiC_accc(a1,i,c0,j) @W1(c0,a1,b,a)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Nvirt*_Ndocc+a*_Ndocc+c0)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//c0
        }//a
      }//b
      orz::DTensor tmp2(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//a1
      }//_P_
      // @W1(c0,a1,b,a) <<= +1 @Tc(-1)(b,a,_P_,c0) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&a1 : irange(0L, _Nact)){
              W1.cptr()[c0*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[b*_Nvirt*_Ndocc*_Nact+a*_Ndocc*_Nact+c0*_Nact+a1];
            }//a1
          }//c0
        }//a
      }//b
    }
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc,_Nact*_Ndocc);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(i*_Ndocc+j)*tmp1.shape(1)+(a1*_Ndocc+c0)] = V2_FiC_accc.cptr()[a1*V2_FiC_accc.shape(1)*V2_FiC_accc.shape(2)*V2_FiC_accc.shape(3)+i*V2_FiC_accc.shape(2)*V2_FiC_accc.shape(3)+c0*V2_FiC_accc.shape(3)+j];
            }//c0
          }//a1
        }//j
      }//i
      orz::DTensor tmp2(_Nact*_Ndocc,_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1*_Ndocc+c0)*tmp2.shape(1)+(b*_Nvirt+a)] = W1.cptr()[c0*W1.shape(1)*W1.shape(2)*W1.shape(3)+a1*W1.shape(2)*W1.shape(3)+b*W1.shape(3)+a];
            }//a
          }//b
        }//c0
      }//a1
      // @tRFiC(i,j,a,b) <<= +1 _X_ @V2_FiC_accc(a1,i,c0,j) @W1(c0,a1,b,a)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+b*_Nvirt+a];
            }//a
          }//b
        }//j
      }//i
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= +1 @V2_FiC_accc(a1,i,c0,j) @W1(c0,a1,a,b)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Nvirt*_Ndocc+b*_Ndocc+c0)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+c0];
            }//_P_
          }//c0
        }//b
      }//a
      orz::DTensor tmp2(Nrhom1,_Nact);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(a1)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//a1
      }//_P_
      // @W1(c0,a1,a,b) <<= +1 @Tc(-1)(a,b,_P_,c0) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&a1 : irange(0L, _Nact)){
              W1.cptr()[c0*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Nvirt*_Ndocc*_Nact+b*_Ndocc*_Nact+c0*_Nact+a1];
            }//a1
          }//c0
        }//b
      }//a
    }
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc,_Nact*_Ndocc);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(i*_Ndocc+j)*tmp1.shape(1)+(a1*_Ndocc+c0)] = V2_FiC_accc.cptr()[a1*V2_FiC_accc.shape(1)*V2_FiC_accc.shape(2)*V2_FiC_accc.shape(3)+i*V2_FiC_accc.shape(2)*V2_FiC_accc.shape(3)+c0*V2_FiC_accc.shape(3)+j];
            }//c0
          }//a1
        }//j
      }//i
      orz::DTensor tmp2(_Nact*_Ndocc,_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1*_Ndocc+c0)*tmp2.shape(1)+(a*_Nvirt+b)] = W1.cptr()[c0*W1.shape(1)*W1.shape(2)*W1.shape(3)+a1*W1.shape(2)*W1.shape(3)+a*W1.shape(3)+b];
            }//b
          }//a
        }//c0
      }//a1
      // @tRFiC(i,j,a,b) <<= +1 _X_ @V2_FiC_accc(a1,i,c0,j) @W1(c0,a1,a,b)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b];
            }//b
          }//a
        }//j
      }//i
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= -1 @Tc(-1)(a,b,_P_,i) @W1(j,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc,_Nact);
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(j)*tmp1.shape(1)+(a1)] = CFip.cptr()[j*CFip.shape(1)+a1];
        }//a1
      }//j
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//_P_
      }//a1
      // @W1(j,_P_) <<= +1 @Fcore(j,a1) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[j*Nrhom1+_P_] += tmp3.cptr()[j*Nrhom1+_P_];
        }//_P_
      }//j
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
      orz::DTensor tmp2(Nrhom1,_Ndocc);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&j : irange(0L, _Ndocc)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(j)] = W1.cptr()[j*W1.shape(1)+_P_];
        }//j
      }//_P_
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(-1)(a,b,_P_,i) @W1(j,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&j : irange(0L, _Ndocc)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Nvirt*_Ndocc*_Ndocc+b*_Ndocc*_Ndocc+i*_Ndocc+j];
            }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= -1 @Tc(-1)(b,a,_P_,i) @W1(j,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc,_Nact);
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(j)*tmp1.shape(1)+(a1)] = CFip.cptr()[j*CFip.shape(1)+a1];
        }//a1
      }//j
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//_P_
      }//a1
      // @W1(j,_P_) <<= +1 @Fcore(j,a1) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[j*Nrhom1+_P_] += tmp3.cptr()[j*Nrhom1+_P_];
        }//_P_
      }//j
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
      orz::DTensor tmp2(Nrhom1,_Ndocc);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&j : irange(0L, _Ndocc)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(j)] = W1.cptr()[j*W1.shape(1)+_P_];
        }//j
      }//_P_
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(-1)(b,a,_P_,i) @W1(j,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&j : irange(0L, _Ndocc)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[b*_Nvirt*_Ndocc*_Ndocc+a*_Ndocc*_Ndocc+i*_Ndocc+j];
            }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= -1 @Tc(-1)(b,a,_P_,j) @W1(i,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_minus.cptr()[_P_*Cm1_minus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(-)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(i)*tmp1.shape(1)+(a1)] = CFip.cptr()[i*CFip.shape(1)+a1];
        }//a1
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//_P_
      }//a1
      // @W1(i,_P_) <<= +1 @Fcore(i,a1) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[i*Nrhom1+_P_] += tmp3.cptr()[i*Nrhom1+_P_];
        }//_P_
      }//i
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(b*_Nvirt*_Ndocc+a*_Ndocc+j)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[b*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+a*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+j];
            }//_P_
          }//j
        }//a
      }//b
      orz::DTensor tmp2(Nrhom1,_Ndocc);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(i)] = W1.cptr()[i*W1.shape(1)+_P_];
        }//i
      }//_P_
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(-1)(b,a,_P_,j) @W1(i,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&i : irange(0L, _Ndocc)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[b*_Nvirt*_Ndocc*_Ndocc+a*_Ndocc*_Ndocc+j*_Ndocc+i];
            }//i
          }//j
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
std::cout << boost::format("%80s") % "@tRFiC(i,j,a,b) <<= -1 @Tc(-1)(a,b,_P_,j) @W1(i,_P_)";
high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    const double _X_ =  -1.0000000000000000;
    orz::DTensor W0(_Nact,Nrhom1);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp1.cptr()[(a1)*tmp1.shape(1)+(a0)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a0
      }//a1
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(_P_)] = Cm1_plus.cptr()[_P_*Cm1_plus.shape(1)+a0];
        }//_P_
      }//a0
      // @W0(a1,_P_) <<= +1 @D1(a1,a0) @X(-1)(+)(_P_,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W0.cptr()[a1*Nrhom1+_P_] += tmp3.cptr()[a1*Nrhom1+_P_];
        }//_P_
      }//a1
    }
    orz::DTensor W1(_Ndocc,Nrhom1);
    {
      orz::DTensor tmp1(_Ndocc,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(i)*tmp1.shape(1)+(a1)] = CFip.cptr()[i*CFip.shape(1)+a1];
        }//a1
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(_P_)] = W0.cptr()[a1*W0.shape(1)+_P_];
        }//_P_
      }//a1
      // @W1(i,_P_) <<= +1 @Fcore(i,a1) @W0(a1,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&_P_ : irange(0L, Nrhom1)){
          W1.cptr()[i*Nrhom1+_P_] += tmp3.cptr()[i*Nrhom1+_P_];
        }//_P_
      }//i
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt*_Ndocc,Nrhom1);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&_P_ : irange(0L, Nrhom1)){
              tmp1.cptr()[(a*_Nvirt*_Ndocc+b*_Ndocc+j)*tmp1.shape(1)+(_P_)] = Tc_m1_.cptr()[a*Tc_m1_.shape(1)*Tc_m1_.shape(2)*Tc_m1_.shape(3)+b*Tc_m1_.shape(2)*Tc_m1_.shape(3)+_P_*Tc_m1_.shape(3)+j];
            }//_P_
          }//j
        }//b
      }//a
      orz::DTensor tmp2(Nrhom1,_Ndocc);
      for(auto &&_P_ : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          tmp2.cptr()[(_P_)*tmp2.shape(1)+(i)] = W1.cptr()[i*W1.shape(1)+_P_];
        }//i
      }//_P_
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(-1)(a,b,_P_,j) @W1(i,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&j : irange(0L, _Ndocc)){
            for(auto &&i : irange(0L, _Ndocc)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Nvirt*_Ndocc*_Ndocc+b*_Ndocc*_Ndocc+j*_Ndocc+i];
            }//i
          }//j
        }//b
      }//a
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
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&j : irange(0L, _Ndocc)){
          const long I_ij = RActList.GetIndex(i, j);
          if (I_ij < 0) continue;
          orz::DTensor Rij(_Nvirt, _Nvirt);
          for(auto &&a : irange(0L, _Nvirt))
            for(auto &&b : irange(0L, _Nvirt))
              Rij(a,b) = tRFiC(i,j,a,b);
          RActList.GetPair(i,j).Transform2EXT(Rij);
          orz::DTensor Rorig = R2.GetTensor(i,j);
          Rij += Rorig;
          R2.PutTensor(i, j, Rij);
        }
      }
        
      return;
    }//End func
  } }
