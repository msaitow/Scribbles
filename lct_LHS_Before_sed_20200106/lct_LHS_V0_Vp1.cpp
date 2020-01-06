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
void FMatrixOrthogonalizer::FormLHS_V0_Vp1(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &d3, DensityPack &dPack,
                                            TensorContainer &ovl, const double &Ecas, const double &h0_energy, const orz::DTensor &casfock, const orz::DTensor &FV,
                                            TensorContainer &T2, TensorContainer &R2){
#include "lct_preparations.cpp"
orz::DTensor &Jint_ccvv_ccav_0_ = (*dPack["Jint_ccvv_ccav_0_"]);
orz::DTensor &Jint_ccvv_ccav_1_ = (*dPack["Jint_ccvv_ccav_1_"]);
orz::DTensor &Jint_ccvv_ccav_2_ = (*dPack["Jint_ccvv_ccav_2_"]);
orz::DTensor &Jint_ccvv_ccav_3_ = (*dPack["Jint_ccvv_ccav_3_"]);

orz::DTensor JF0(Jint_ccvv_ccav_0_*Fpa);
orz::DTensor JF1(Jint_ccvv_ccav_1_*Fpa);
orz::DTensor JF2(Jint_ccvv_ccav_2_*Fpa);
orz::DTensor JF3(Jint_ccvv_ccav_3_*Fpa); 
 
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(i,j,a,b) <<= +1 @Jint_ccvv_ccav_0_(P0,a1) @cF(a1,a) @T(+1)(b,P0,j,i)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,b) <<= +1 @Jint_ccvv_ccav_0_(P0,a1) @cF(a1,a) @T(+1)(b,P0,j,i)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = T2.GetTensor(j, i);
    orz::DTensor Fpa(JF0.copy());
    {
      RActList.GetPair(i, j).Transform1EXT_R(Fpa);
    }
    Add_A_x_B(S2, Fpa, U2, true, true);
    R2.PutTensor(i, j, S2);
  }//j
}//i

high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
#ifdef _PROFILE_LHS
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif
}
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(i,j,a,b) <<= +1 @Jint_ccvv_ccav_1_(P0,a1) @cF(a1,a) @T(+1)(b,P0,i,j)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,b) <<= +1 @Jint_ccvv_ccav_1_(P0,a1) @cF(a1,a) @T(+1)(b,P0,i,j)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = T2.GetTensor(i, j);
    orz::DTensor Fpa(JF1.copy());
    {
      RActList.GetPair(i, j).Transform1EXT_R(Fpa);
    }
    Add_A_x_B(S2, Fpa, U2, true, true);
    R2.PutTensor(i, j, S2);
  }//j
}//i

high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
#ifdef _PROFILE_LHS
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif
}
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(i,j,a,b) <<= +1 @Jint_ccvv_ccav_2_(P0,a1) @cF(a1,b) @T(+1)(a,P0,j,i)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,b) <<= +1 @Jint_ccvv_ccav_2_(P0,a1) @cF(a1,b) @T(+1)(a,P0,j,i)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = T2.GetTensor(j, i);
    orz::DTensor Fpa(JF2.copy());
    {
      RActList.GetPair(i, j).Transform1EXT_R(Fpa);
    }
    Add_A_x_B(S2, U2, Fpa, false, false);
    R2.PutTensor(i, j, S2);
  }//j
}//i

high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
#ifdef _PROFILE_LHS
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif
}
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(i,j,a,b) <<= +1 @Jint_ccvv_ccav_3_(P0,a1) @cF(a1,b) @T(+1)(a,P0,i,j)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,b) <<= +1 @Jint_ccvv_ccav_3_(P0,a1) @cF(a1,b) @T(+1)(a,P0,i,j)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = T2.GetTensor(i, j);
    orz::DTensor Fpa(JF3.copy());
    {
      RActList.GetPair(i, j).Transform1EXT_R(Fpa);
    }
    Add_A_x_B(S2, U2, Fpa, false, false);
    R2.PutTensor(i, j, S2);
  }//j
}//i

high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
#ifdef _PROFILE_LHS
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif
}
      return;
    }//End func
  } }
