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
void FMatrixOrthogonalizer::FormLHS_Vp1_V0(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &d3, DensityPack &dPack,
                                            TensorContainer &ovl, const double &Ecas, const double &h0_energy, const orz::DTensor &casfock, const orz::DTensor &FV,
                                            TensorContainer &T2, TensorContainer &R2){
#include "lct_preparations.cpp"
TensorContainer FavT(_Ndocc,_Ndocc);
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = T2.GetTensor(i, j);
    orz::DTensor Fpa(_Fpa_.copy());
    {
      RActList.GetPair(i, j).Transform1EXT_R(Fpa);
    }
    orz::DTensor X(_Nact, _Nvirt_ij);
//    for(auto &&a1 : irange(0L, _Nact)){
//      for(auto &&v0 : irange(0L, _Nvirt_ij)){
//        for(auto &&a : irange(0L, _Nvirt_ij)){
//          X.cptr()[a1*X.shape(1)+a] += U2.cptr()[v0*U2.shape(1)+a] * Fpa.cptr()[a1*Fpa.shape(1)+v0];
//        }//a
//      }//v0
//    }//a1
    A_x_B(X, Fpa, U2, false, false);
    FavT.PutTensor(i,j,X);
  }//j
}//i
orz::DTensor &Jint_ccav_ccvv_0_ = (*dPack["Jint_ccav_ccvv_0_"]);
orz::DTensor &Jint_ccav_ccvv_1_ = (*dPack["Jint_ccav_ccvv_1_"]);
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(i,j,a,P) <<= +1 @FavT(a1,a,j,i) @Jint_ccav_ccvv_0_(P,a1)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,P) <<= +1 @FavT(a1,a,j,i) @Jint_ccav_ccvv_0_(P,a1)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = FavT.GetTensor(j, i);
    for(auto &&a1 : irange(0L, _Nact)){
      for(auto &&P : irange(0L, Nrhop1)){
        for(auto &&a : irange(0L, _Nvirt_ij)){
          S2.cptr()[a*S2.shape(1)+P] += U2.cptr()[a1*U2.shape(1)+a] * Jint_ccav_ccvv_0_.cptr()[P*Jint_ccav_ccvv_0_.shape(1)+a1];
        }//a
      }//P
    }//a1
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
std::cout << boost::format("%80s") % "@S2(i,j,a,P) <<= +1 @FavT(a1,a,i,j) @Jint_ccav_ccvv_1_(P,a1)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,P) <<= +1 @FavT(a1,a,i,j) @Jint_ccav_ccvv_1_(P,a1)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = FavT.GetTensor(i, j);
    for(auto &&a1 : irange(0L, _Nact)){
      for(auto &&P : irange(0L, Nrhop1)){
        for(auto &&a : irange(0L, _Nvirt_ij)){
          S2.cptr()[a*S2.shape(1)+P] += U2.cptr()[a1*U2.shape(1)+a] * Jint_ccav_ccvv_1_.cptr()[P*Jint_ccav_ccvv_1_.shape(1)+a1];
        }//a
      }//P
    }//a1
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
