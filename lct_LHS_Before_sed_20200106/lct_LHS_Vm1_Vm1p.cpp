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
void FMatrixOrthogonalizer::FormLHS_Vm1_Vm1p(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &d3, DensityPack &dPack,
                                            TensorContainer &ovl, const double &Ecas, const double &h0_energy, const orz::DTensor &casfock, const orz::DTensor &FV,
                                            TensorContainer &T2, TensorContainer &R2){
#include "lct_preparations.cpp"
orz::DTensor &Jint_cavv_aaav_0_ = (*dPack["Jint_cavv_aaav_0_"]);
orz::DTensor &Jint_cavv_aaav_1_ = (*dPack["Jint_cavv_aaav_1_"]);

orz::DTensor JT0(Nrhom1, _Nvirt), JT1(Nrhom1, _Nvirt);
for(auto &&P0 : irange(0L, Nrhom1p)){
  const long I_P0 = TActList.GetIndex(P0);
  if (I_P0 < 0) continue;
  orz::DTensor U2 = T2.GetTensor(P0);
  for(auto &&P : irange(0L, Nrhom1)){
    orz::DTensor X0(_Nvirt);
    X0.gaxpy(+1.0, U2, Jint_cavv_aaav_0_.cptr()[P*Jint_cavv_aaav_0_.shape(1)+P0]);
    for(auto &&a : irange(0L, _Nvirt)) JT0.cptr()[P*JT0.shape(1)+a] += X0.cptr()[a];    
    orz::DTensor X1(_Nvirt);
    X1.gaxpy(+1.0, U2, Jint_cavv_aaav_1_.cptr()[P*Jint_cavv_aaav_1_.shape(1)+P0]);
    for(auto &&a : irange(0L, _Nvirt)) JT1.cptr()[P*JT1.shape(1)+a] += X1.cptr()[a];    
  }
} 

{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(P,i,a,b) <<= +1 @Jint_cavv_aaav_0_(P,P0) @cF(i,a) @T(-1')(b,P0)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(P,i,a,b) <<= +1 @Jint_cavv_aaav_0_(P,P0) @cF(i,a) @T(-1')(b,P0)
for(auto &&P : irange(0L, Nrhom1)){
  for(auto &&i : irange(0L, _Ndocc)){
    const long I_Pi = RActList.GetIndex(P, i);
    if (I_Pi < 0) continue;
    orz::DTensor S2 = R2.GetTensor(P, i);
    const long _Nvirt_Pi = RActList.GetNPNO(P, i);
    if (I_Pi < 0) continue;
    orz::DTensor U2(JT0);
    {
      high_resolution_clock::time_point o1 = high_resolution_clock::now();
      orz::DTensor tmp(U2.copy());
      RActList.GetPair(P, i).Transform1EXT_R(tmp);
      U2 = tmp.copy();
      high_resolution_clock::time_point o2 = high_resolution_clock::now();
      duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
      t_trans += time_span.count();
    }
    orz::DTensor Fia(_Fia_.copy());
    {
      RActList.GetPair(P, i).Transform1EXT_R(Fia);
    }
    for(auto &&a : irange(0L, _Nvirt_Pi)){
      for(auto &&b : irange(0L, _Nvirt_Pi)){
        S2.cptr()[a*S2.shape(1)+b] += U2.cptr()[P*U2.shape(1)+b] * Fia.cptr()[i*Fia.shape(1)+a];
      }//b
    }//a
    R2.PutTensor(P, i, S2);
  }//i
}//P

high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
#ifdef _PROFILE_LHS
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif
}

{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(P,i,a,b) <<= +1 @Jint_cavv_aaav_1_(P,P0) @cF(i,b) @T(-1')(a,P0)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(P,i,a,b) <<= +1 @Jint_cavv_aaav_1_(P,P0) @cF(i,b) @T(-1')(a,P0)
for(auto &&P : irange(0L, Nrhom1)){
  for(auto &&i : irange(0L, _Ndocc)){
    const long I_Pi = RActList.GetIndex(P, i);
    if (I_Pi < 0) continue;
    orz::DTensor S2 = R2.GetTensor(P, i);
    const long _Nvirt_Pi = RActList.GetNPNO(P, i);
    if (I_Pi < 0) continue;
    orz::DTensor U2(JT1);
    {
      high_resolution_clock::time_point o1 = high_resolution_clock::now();
      orz::DTensor tmp(U2.copy());
      RActList.GetPair(P, i).Transform1EXT_R(tmp);
      U2 = tmp.copy();
      high_resolution_clock::time_point o2 = high_resolution_clock::now();
      duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
      t_trans += time_span.count();
    }
    orz::DTensor Fia(_Fia_.copy());
    {
      RActList.GetPair(P, i).Transform1EXT_R(Fia);
    }
    for(auto &&b : irange(0L, _Nvirt_Pi)){
      for(auto &&a : irange(0L, _Nvirt_Pi)){
        S2.cptr()[a*S2.shape(1)+b] += U2.cptr()[P*U2.shape(1)+a] * Fia.cptr()[i*Fia.shape(1)+b];
      }//a
    }//b
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
//std::cout << boost::format("%80s") % "@S2(P,i,a,b) <<= +1 @Jint_cavv_aaav_0_(P,P0) @cF(i,a) @T(-1')(b,P0)";
//#endif
//high_resolution_clock::time_point t1 = high_resolution_clock::now();
//// @S2(P,i,a,b) <<= +1 @Jint_cavv_aaav_0_(P,P0) @cF(i,a) @T(-1')(b,P0)
//for(auto &&P : irange(0L, Nrhom1)){
//  for(auto &&i : irange(0L, _Ndocc)){
//    const long I_Pi = RActList.GetIndex(P, i);
//    if (I_Pi < 0) continue;
//    orz::DTensor S2 = R2.GetTensor(P, i);
//    const long _Nvirt_Pi = RActList.GetNPNO(P, i);
//    if (I_Pi < 0) continue;
//    for(auto &&P0 : irange(0L, Nrhom1p)){
//      const long I_P0 = TActList.GetIndex(P0);
//      if (I_P0 < 0) continue;
//      orz::DTensor U2 = T2.GetTensor(P0);
//      {
//        high_resolution_clock::time_point o1 = high_resolution_clock::now();
//        orz::DTensor tmp(U2.copy());
//        RActList.GetPair(P, i).Transform1EXT_L(tmp);
//        U2 = tmp.copy();
//        high_resolution_clock::time_point o2 = high_resolution_clock::now();
//        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
//        t_trans += time_span.count();
//      }
//      orz::DTensor Fia(_Fia_.copy());
//      {
//        RActList.GetPair(P, i).Transform1EXT_R(Fia);
//      }
//      for(auto &&a : irange(0L, _Nvirt_Pi)){
//        for(auto &&b : irange(0L, _Nvirt_Pi)){
//          S2.cptr()[a*S2.shape(1)+b] += Jint_cavv_aaav_0_.cptr()[P*Jint_cavv_aaav_0_.shape(1)+P0] * Fia.cptr()[i*Fia.shape(1)+a] * U2.cptr()[b];
//        }//b
//      }//a
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
//{
//double t_trans=0.0;
//#ifdef _PROFILE_LHS
//std::cout << boost::format("%80s") % "@S2(P,i,a,b) <<= +1 @Jint_cavv_aaav_1_(P,P0) @cF(i,b) @T(-1')(a,P0)";
//#endif
//high_resolution_clock::time_point t1 = high_resolution_clock::now();
//// @S2(P,i,a,b) <<= +1 @Jint_cavv_aaav_1_(P,P0) @cF(i,b) @T(-1')(a,P0)
//for(auto &&P : irange(0L, Nrhom1)){
//  for(auto &&i : irange(0L, _Ndocc)){
//    const long I_Pi = RActList.GetIndex(P, i);
//    if (I_Pi < 0) continue;
//    orz::DTensor S2 = R2.GetTensor(P, i);
//    const long _Nvirt_Pi = RActList.GetNPNO(P, i);
//    if (I_Pi < 0) continue;
//    for(auto &&P0 : irange(0L, Nrhom1p)){
//      const long I_P0 = TActList.GetIndex(P0);
//      if (I_P0 < 0) continue;
//      orz::DTensor U2 = T2.GetTensor(P0);
//      {
//        high_resolution_clock::time_point o1 = high_resolution_clock::now();
//        orz::DTensor tmp(U2.copy());
//        RActList.GetPair(P, i).Transform1EXT_L(tmp);
//        U2 = tmp.copy();
//        high_resolution_clock::time_point o2 = high_resolution_clock::now();
//        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
//        t_trans += time_span.count();
//      }
//      orz::DTensor Fia(_Fia_.copy());
//      {
//        RActList.GetPair(P, i).Transform1EXT_R(Fia);
//      }
//      for(auto &&b : irange(0L, _Nvirt_Pi)){
//        for(auto &&a : irange(0L, _Nvirt_Pi)){
//          S2.cptr()[a*S2.shape(1)+b] += Jint_cavv_aaav_1_.cptr()[P*Jint_cavv_aaav_1_.shape(1)+P0] * Fia.cptr()[i*Fia.shape(1)+b] * U2.cptr()[a];
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