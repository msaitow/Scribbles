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
void FMatrixOrthogonalizer::FormLHS_Vm1_V0p(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &d3, DensityPack &dPack,
                                            TensorContainer &ovl, const double &Ecas, const double &h0_energy, const orz::DTensor &casfock, const orz::DTensor &FV,
                                            TensorContainer &T2, TensorContainer &R2){
#include "lct_preparations.cpp"
orz::DTensor &Jint_cavv_caav_0_ = (*dPack["Jint_cavv_caav_0_"]);
orz::DTensor &Jint_cavv_caav_1_ = (*dPack["Jint_cavv_caav_1_"]);

#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "JF0 and JF1";
#endif
high_resolution_clock::time_point t_jf1 = high_resolution_clock::now();
TensorContainer JF0(Nrhom1), JF1(Nrhom1);
for(auto &&P : irange(0L, Nrhom1)){
  orz::DTensor X1;
  {
    orz::DTensor tmp1(Nrho0p, _Nact);
    for(auto &&P0 : irange(0L, Nrho0p)){
      for(auto &&a1 : irange(0L, _Nact)){
        tmp1.cptr()[P0*tmp1.shape(1)+a1] = Jint_cavv_caav_0_.cptr()[P*Jint_cavv_caav_0_.shape(1)*Jint_cavv_caav_0_.shape(2)+P0*Jint_cavv_caav_0_.shape(2)+a1];
      }//a1
    }//P0
    X1 = (tmp1*Fpa);
  }
  orz::DTensor X2;
  {
    orz::DTensor tmp2(Nrho0p, _Nact);
    for(auto &&P0 : irange(0L, Nrho0p)){
      for(auto &&a1 : irange(0L, _Nact)){
        tmp2.cptr()[P0*tmp2.shape(1)+a1] = Jint_cavv_caav_1_.cptr()[P*Jint_cavv_caav_1_.shape(1)*Jint_cavv_caav_1_.shape(2)+P0*Jint_cavv_caav_1_.shape(2)+a1];
      }//a1
    }//P0
    X2 = (tmp2*Fpa);
  }
  JF0.PutTensor(P,X1);
  JF1.PutTensor(P,X2);  
}//P
high_resolution_clock::time_point t_jf2 = high_resolution_clock::now();
duration<double> time_span_jf = duration_cast<duration<double>>(t_jf2 - t_jf1);
#ifdef _PROFILE_LHS
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span_jf.count() % 0 << endl;
#endif
 
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(P,i,a,b) <<= +1 @Jint_cavv_caav_0_(P,P0,a1) @cF(a1,a) @T(0')(i,b,P0)" << std::endl;
std::cout << boost::format("%80s") % "@S2(P,i,a,b) <<= +1 @Jint_cavv_caav_1_(P,P0,a1) @cF(a1,b) @T(0')(i,a,P0)"; 
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(P,i,a,b) <<= +1 @Jint_cavv_caav_0_(P,P0,a1) @cF(a1,a) @T(0')(i,b,P0)
for(auto &&P : irange(0L, Nrhom1)){
  orz::DTensor X0 = JF0.GetTensor(P);
  orz::DTensor X1 = JF1.GetTensor(P);  
  for(auto &&i : irange(0L, _Ndocc)){
    const long I_Pi = RActList.GetIndex(P, i);
    if (I_Pi < 0) continue;
    orz::DTensor S2 = R2.GetTensor(P, i);
    const long _Nvirt_Pi = RActList.GetNPNO(P, i);
    if (I_Pi < 0) continue;
    orz::DTensor X(X0.copy());
    {
      RActList.GetPair(P, i).Transform1EXT_R(X);
    }
    orz::DTensor Y(X1.copy());
    {
      RActList.GetPair(P, i).Transform1EXT_R(Y);
    }    
    for(auto &&P0 : irange(0L, Nrho0p)){
      const long I_P0 = TActList.GetIndex(P0);
      if (I_P0 < 0) continue;
      orz::DTensor U2 = T2.GetTensor(P0);
      {
        high_resolution_clock::time_point o1 = high_resolution_clock::now();
        orz::DTensor tmp(U2.copy());
        RActList.GetPair(P, i).Transform1EXT_R(tmp);
        U2 = tmp.copy();
        high_resolution_clock::time_point o2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
        t_trans += time_span.count();
      }
      for(auto &&a : irange(0L, _Nvirt_Pi)){
        for(auto &&b : irange(0L, _Nvirt_Pi)){
          S2.cptr()[a*S2.shape(1)+b] += X.cptr()[P0*X.shape(1)+a] * U2.cptr()[i*U2.shape(1)+b];
        }//b
      }//a
      for(auto &&b : irange(0L, _Nvirt_Pi)){
        for(auto &&a : irange(0L, _Nvirt_Pi)){
          S2.cptr()[a*S2.shape(1)+b] += Y.cptr()[P0*Y.shape(1)+b] * U2.cptr()[i*U2.shape(1)+a];
        }//a
      }//b
    }//P0
    R2.PutTensor(P, i, S2);
  }//i
}//P

high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
#ifdef _PROFILE_LHS
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif
}
//{
//double t_trans=0.0;
//#ifdef _PROFILE_LHS
//std::cout << boost::format("%80s") % "@S2(P,i,a,b) <<= +1 @Jint_cavv_caav_1_(P,P0,a1) @cF(a1,b) @T(0')(i,a,P0)";
//#endif
//high_resolution_clock::time_point t1 = high_resolution_clock::now();
//// @S2(P,i,a,b) <<= +1 @Jint_cavv_caav_1_(P,P0,a1) @cF(a1,b) @T(0')(i,a,P0)
//for(auto &&P : irange(0L, Nrhom1)){
//  orz::DTensor X1 = JF1.GetTensor(P);  
//  for(auto &&i : irange(0L, _Ndocc)){
//    const long I_Pi = RActList.GetIndex(P, i);
//    if (I_Pi < 0) continue;
//    orz::DTensor S2 = R2.GetTensor(P, i);
//    const long _Nvirt_Pi = RActList.GetNPNO(P, i);
//    if (I_Pi < 0) continue;
//    for(auto &&P0 : irange(0L, Nrho0p)){
//      const long I_P0 = TActList.GetIndex(P0);
//      if (I_P0 < 0) continue;
//      orz::DTensor U2 = T2.GetTensor(P0);
//      {
//        high_resolution_clock::time_point o1 = high_resolution_clock::now();
//        orz::DTensor tmp(U2.copy());
//        RActList.GetPair(P, i).Transform1EXT_R(tmp);
//        U2 = tmp.copy();
//        high_resolution_clock::time_point o2 = high_resolution_clock::now();
//        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
//        t_trans += time_span.count();
//      }
//      orz::DTensor X(X1.copy());
//      {
//        RActList.GetPair(P, i).Transform1EXT_R(X);
//      }
//      for(auto &&b : irange(0L, _Nvirt_Pi)){
//        for(auto &&a : irange(0L, _Nvirt_Pi)){
//          S2.cptr()[a*S2.shape(1)+b] += X.cptr()[P0*X.shape(1)+a] * U2.cptr()[i*U2.shape(1)+a];
//        }//a
//      }//b
//    }//P0
//    R2.PutTensor(P, i, S2);
//  }//i
//}//P
//
//high_resolution_clock::time_point t2 = high_resolution_clock::now();
//duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
//#ifdef _PROFILE_LHS
//cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
//#endif
//}
      return;
    }//End func
  } }
