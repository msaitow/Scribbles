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
void FMatrixOrthogonalizer::FormLHS_Vm1p_Vm2(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &d3, DensityPack &dPack,
                                            TensorContainer &ovl, const double &Ecas, const double &h0_energy, const orz::DTensor &casfock, const orz::DTensor &FV,
                                            TensorContainer &T2, TensorContainer &R2){
#include "lct_preparations.cpp"
// @cF(a5,v0) @T(-2)(v0,a,P0)
TensorContainer FavT(Nrhom2);
for(auto &&P0 : irange(0L, Nrhom2)){
  const long I_P0 = TActList.GetIndex(P0);
  if (I_P0 < 0) continue;
  const long  _Nvirt_P0 = TActList.GetPair(P0).GetNPNO();  orz::DTensor U2 = T2.GetTensor(P0);
  {
    high_resolution_clock::time_point o1 = high_resolution_clock::now();
    orz::DTensor tmp(U2.copy());
    TActList.GetPair(0,P0).BackTransform1EXT_L(tmp);
    U2 = tmp.copy();
    high_resolution_clock::time_point o2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
    //t_trans += time_span.count();
  }
  orz::DTensor X(_Nact, _Nvirt_P0);
//  for(auto &&a5 : irange(0L, _Nact)){
//    for(auto &&v0 : irange(0L, _Nvirt)){
//      for(auto &&a : irange(0L, _Nvirt_P0)){
//        X.cptr()[a5*X.shape(1)+a] += Fpa.cptr()[a5*Fpa.shape(1)+v0] * U2.cptr()[v0*U2.shape(1)+a];
//      }//a
//    }//v0
//  }//a5
  A_x_B(X, Fpa, U2, false, false);
  FavT.PutTensor(P0,X);
}//P0
// @cF(a5,v0) @T(-2)(a,v0,P0)
TensorContainer TFav(Nrhom2);
for(auto &&P0 : irange(0L, Nrhom2)){
  const long I_P0 = TActList.GetIndex(P0);
  if (I_P0 < 0) continue;
  const long  _Nvirt_P0 = TActList.GetPair(P0).GetNPNO();  orz::DTensor U2 = T2.GetTensor(P0);
  {
    high_resolution_clock::time_point o1 = high_resolution_clock::now();
    orz::DTensor tmp(U2.copy());
    TActList.GetPair(0,P0).BackTransform1EXT_R(tmp);
    U2 = tmp.copy();
    high_resolution_clock::time_point o2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
    //t_trans += time_span.count();
  }
  orz::DTensor X(_Nvirt_P0, _Nact);
//  for(auto &&a5 : irange(0L, _Nact)){
//    for(auto &&v0 : irange(0L, _Nvirt)){
//      for(auto &&a : irange(0L, _Nvirt_P0)){
//        X.cptr()[a*X.shape(1)+a5] += Fpa.cptr()[a5*Fpa.shape(1)+v0] * U2.cptr()[a*U2.shape(1)+v0];
//      }//a
//    }//v0
//  }//a5
  A_x_B(X, U2, Fpa, false, true);
  TFav.PutTensor(P0,X);
}//P0
orz::DTensor &Jint_aaav_aavv_0_ = (*dPack["Jint_aaav_aavv_0_"]);
orz::DTensor &Jint_aaav_aavv_1_ = (*dPack["Jint_aaav_aavv_1_"]);
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(P,a) <<= +1 @FavT(a5,a,P0) @Jint_aaav_aavv_0_(P,P0,a5)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(P,a) <<= +1 @FavT(a5,a,P0) @Jint_aaav_aavv_0_(P,P0,a5)
for(auto &&P : irange(0L, Nrhom1p)){
  const long I_P = RActList.GetIndex(P);
  if (I_P < 0) continue;
  orz::DTensor S2 = R2.GetTensor(P);
  for(auto &&P0 : irange(0L, Nrhom2)){
    const long I_P0 = TActList.GetIndex(P0);
    if (I_P0 < 0) continue;
    orz::DTensor U2 = FavT.GetTensor(P0);
    {
      high_resolution_clock::time_point o1 = high_resolution_clock::now();
      orz::DTensor tmp(U2.copy());
      TActList.GetPair(0,P0).BackTransform1EXT_R(tmp);
      U2 = tmp.copy();
      high_resolution_clock::time_point o2 = high_resolution_clock::now();
      duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
      t_trans += time_span.count();
    }
    for(auto &&a5 : irange(0L, _Nact)){
      for(auto &&a : irange(0L, _Nvirt)){
        S2.cptr()[a] += U2.cptr()[a5*U2.shape(1)+a] * Jint_aaav_aavv_0_.cptr()[P*Jint_aaav_aavv_0_.shape(1)*Jint_aaav_aavv_0_.shape(2)+P0*Jint_aaav_aavv_0_.shape(2)+a5];
      }//a
    }//a5
  }//P0
  R2.PutTensor(P, S2);
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
std::cout << boost::format("%80s") % "@S2(P,a) <<= +1 @TFav(a,a5,P0) @Jint_aaav_aavv_1_(P,P0,a5)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(P,a) <<= +1 @TFav(a,a5,P0) @Jint_aaav_aavv_1_(P,P0,a5)
for(auto &&P : irange(0L, Nrhom1p)){
  const long I_P = RActList.GetIndex(P);
  if (I_P < 0) continue;
  orz::DTensor S2 = R2.GetTensor(P);
  for(auto &&P0 : irange(0L, Nrhom2)){
    const long I_P0 = TActList.GetIndex(P0);
    if (I_P0 < 0) continue;
    orz::DTensor U2 = TFav.GetTensor(P0);
    {
      high_resolution_clock::time_point o1 = high_resolution_clock::now();
      orz::DTensor tmp(U2.copy());
      TActList.GetPair(0,P0).BackTransform1EXT_L(tmp);
      U2 = tmp.copy();
      high_resolution_clock::time_point o2 = high_resolution_clock::now();
      duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
      t_trans += time_span.count();
    }
    for(auto &&a5 : irange(0L, _Nact)){
      for(auto &&a : irange(0L, _Nvirt)){
        S2.cptr()[a] += U2.cptr()[a*U2.shape(1)+a5] * Jint_aaav_aavv_1_.cptr()[P*Jint_aaav_aavv_1_.shape(1)*Jint_aaav_aavv_1_.shape(2)+P0*Jint_aaav_aavv_1_.shape(2)+a5];
      }//a
    }//a5
  }//P0
  R2.PutTensor(P, S2);
}//P

high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
#ifdef _PROFILE_LHS
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif
}
      return;
    }//End func
  } }
