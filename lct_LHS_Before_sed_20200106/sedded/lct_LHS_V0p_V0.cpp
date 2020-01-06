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
void FMatrixOrthogonalizer::FormLHS_V0p_V0(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &d3, DensityPack &dPack,
                                            TensorContainer &ovl, const double &Ecas, const double &h0_energy, const orz::DTensor &casfock, const orz::DTensor &FV,
                                            TensorContainer &T2, TensorContainer &R2){
#include "lct_preparations.cpp"
orz::DTensor FTE(_Nvirt,_Ndocc);
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&c0 : irange(0L, _Ndocc)){
    const long I_c0i = TActList.GetIndex(c0, i);
    if (I_c0i < 0) continue;
    orz::DTensor U2 = T2.GetTensor(i, c0);
    {
      high_resolution_clock::time_point o1 = high_resolution_clock::now();
      orz::DTensor tmp(U2.copy());
      TActList.GetPair(c0, i).BackTransform2EXT(tmp);
      U2 = tmp.copy();
      high_resolution_clock::time_point o2 = high_resolution_clock::now();
      duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
      //t_trans += time_span.count();
    }
//    for(auto &&v0 : irange(0L, _Nvirt)){
//      for(auto &&a : irange(0L, _Nvirt)){
//        FTE.cptr()[a*FTE.shape(1)+i] += U2.cptr()[v0*U2.shape(1)+a] * Fia.cptr()[c0*Fia.shape(1)+v0];
//      }//a
//    }//v0
    orz::DTensor Fia_c0(_Nvirt);
    for(auto &&v0 : irange(0L, _Nvirt)) Fia_c0.cptr()[v0] = Fia.cptr()[c0*Fia.shape(1)+v0];
    orz::DTensor X;
    A_x_B(X, U2, Fia_c0, true, false);
    for(auto &&a : irange(0L, _Nvirt)) FTE.cptr()[a*FTE.shape(1)+i] += X.cptr()[a];
  }//c0
}//i
orz::DTensor FTD(_Nvirt,_Ndocc);
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&c0 : irange(0L, _Ndocc)){
    const long I_c0i = TActList.GetIndex(c0, i);
    if (I_c0i < 0) continue;
    orz::DTensor U2 = T2.GetTensor(c0, i);
    {
      high_resolution_clock::time_point o1 = high_resolution_clock::now();
      orz::DTensor tmp(U2.copy());
      TActList.GetPair(c0, i).BackTransform2EXT(tmp);
      U2 = tmp.copy();
      high_resolution_clock::time_point o2 = high_resolution_clock::now();
      duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
      //t_trans += time_span.count();
    }
//    for(auto &&v0 : irange(0L, _Nvirt)){
//      for(auto &&a : irange(0L, _Nvirt)){
//        FTD.cptr()[a*FTD.shape(1)+i] += U2.cptr()[v0*U2.shape(1)+a] * Fia.cptr()[c0*Fia.shape(1)+v0];
//      }//a
//    }//v0
    orz::DTensor Fia_c0(_Nvirt);
    for(auto &&v0 : irange(0L, _Nvirt)) Fia_c0.cptr()[v0] = Fia.cptr()[c0*Fia.shape(1)+v0];
    orz::DTensor X;
    A_x_B(X, U2, Fia_c0, true, false);
    for(auto &&a : irange(0L, _Nvirt)) FTD.cptr()[a*FTD.shape(1)+i] += X.cptr()[a];
  }//c0
}//i
orz::DTensor &Jint_caav_ccvv_0_ = (*dPack["Jint_caav_ccvv_0_"]);
orz::DTensor &Jint_caav_ccvv_1_ = (*dPack["Jint_caav_ccvv_1_"]);
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(P,i,a) <<= +1 @FTE(a,i) @Jint_caav_ccvv_0_(P)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(P,i,a) <<= +1 @FTE(a,i) @Jint_caav_ccvv_0_(P)
for(auto &&P : irange(0L, Nrho0p)){
  const long I_P = RActList.GetIndex(P);
  if (I_P < 0) continue;
  orz::DTensor S2 = R2.GetTensor(P);
  for(auto &&i : irange(0L, _Ndocc)){
    for(auto &&a : irange(0L, _Nvirt)){
      S2.cptr()[i*S2.shape(1)+a] += FTE.cptr()[a*FTE.shape(1)+i] * Jint_caav_ccvv_0_.cptr()[P];
    }//a
  }//i
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
std::cout << boost::format("%80s") % "@S2(P,i,a) <<= +1 @FTD(a,i) @Jint_caav_ccvv_1_(P)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(P,i,a) <<= +1 @FTD(a,i) @Jint_caav_ccvv_1_(P)
for(auto &&P : irange(0L, Nrho0p)){
  const long I_P = RActList.GetIndex(P);
  if (I_P < 0) continue;
  orz::DTensor S2 = R2.GetTensor(P);
  for(auto &&i : irange(0L, _Ndocc)){
    for(auto &&a : irange(0L, _Nvirt)){
      S2.cptr()[i*S2.shape(1)+a] += FTD.cptr()[a*FTD.shape(1)+i] * Jint_caav_ccvv_1_.cptr()[P];
    }//a
  }//i
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
