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
void FMatrixOrthogonalizer::FormLHS_Vp1p_Vp1p(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &d3, DensityPack &dPack,
                                            TensorContainer &ovl, const double &Ecas, const double &h0_energy, const orz::DTensor &casfock, const orz::DTensor &FV,
                                            TensorContainer &T2, TensorContainer &R2){
#include "lct_preparations.cpp"
for(auto &&rho : irange(0L, Nrhop1p)){
  orz::DTensor tmp(_Ndocc);           
  R2.PutTensor(rho,tmp);              
}                                     
// @cF(i,c0) @T(+1')(c0,P0)
TensorContainer FccT(Nrhop1p);
for(auto &&P0 : irange(0L, Nrhop1p)){
  const long I_P0 = TActList.GetIndex(P0);
  if (I_P0 < 0) continue;
  orz::DTensor U2 = T2.GetTensor(P0);
  orz::DTensor X(_Ndocc);
//  for(auto &&i : irange(0L, _Ndocc)){
//    for(auto &&c0 : irange(0L, _Ndocc)){
//      X.cptr()[i] += Fij.cptr()[i*Fij.shape(1)+c0] * U2.cptr()[c0];
//    }//c0
//  }//i
  A_x_B(X, Fij, U2, false, false);
  FccT.PutTensor(P0,X);
}//P0
orz::DTensor &Jint_caaa_caaa_0_ = (*dPack["Jint_caaa_caaa_0_"]);
orz::DTensor &Jint_caaa_caaa_1_ = (*dPack["Jint_caaa_caaa_1_"]);
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(P,i) <<= +1 @FccT(i,P0) @Jint_caaa_caaa_0_(P,P0)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(P,i) <<= +1 @FccT(i,P0) @Jint_caaa_caaa_0_(P,P0)
for(auto &&P : irange(0L, Nrhop1p)){
  const long I_P = RActList.GetIndex(P);
  if (I_P < 0) continue;
  orz::DTensor S2 = R2.GetTensor(P);
  for(auto &&P0 : irange(0L, Nrhop1p)){
    const long I_P0 = TActList.GetIndex(P0);
    if (I_P0 < 0) continue;
    orz::DTensor U2 = FccT.GetTensor(P0);
    for(auto &&i : irange(0L, _Ndocc)){
      S2.cptr()[i] += U2.cptr()[i] * Jint_caaa_caaa_0_.cptr()[P*Jint_caaa_caaa_0_.shape(1)+P0];
    }//i
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
std::cout << boost::format("%80s") % "@S2(P,i) <<= +1 @Jint_caaa_caaa_1_(P,P0) @T(+1')(i,P0)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(P,i) <<= +1 @Jint_caaa_caaa_1_(P,P0) @T(+1')(i,P0)
for(auto &&P : irange(0L, Nrhop1p)){
  const long I_P = RActList.GetIndex(P);
  if (I_P < 0) continue;
  orz::DTensor S2 = R2.GetTensor(P);
  for(auto &&P0 : irange(0L, Nrhop1p)){
    const long I_P0 = TActList.GetIndex(P0);
    if (I_P0 < 0) continue;
    orz::DTensor U2 = T2.GetTensor(P0);
    for(auto &&i : irange(0L, _Ndocc)){
      S2.cptr()[i] += Jint_caaa_caaa_1_.cptr()[P*Jint_caaa_caaa_1_.shape(1)+P0] * U2.cptr()[i];
    }//i
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
