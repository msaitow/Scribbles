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
void FMatrixOrthogonalizer::FormLHS_Vp1_Vp2(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &d3, DensityPack &dPack,
                                            TensorContainer &ovl, const double &Ecas, const double &h0_energy, const orz::DTensor &casfock, const orz::DTensor &FV,
                                            TensorContainer &T2, TensorContainer &R2){
#include "lct_preparations.cpp"
orz::DTensor &Jint_ccav_ccaa_0_ = (*dPack["Jint_ccav_ccaa_0_"]);
orz::DTensor &Jint_ccav_ccaa_1_ = (*dPack["Jint_ccav_ccaa_1_"]);

TensorContainer JF0(Nrhop1), JF1(Nrhop1);
 {
   orz::DTensor fpa(_Fpa_.copy());
   for(auto &&P : irange(0L, Nrhop1)){
     orz::DTensor J0_P(Nrhop2,_Nact), J1_P(Nrhop2,_Nact);
     for(auto &&P0 : irange(0L, Nrhop2)){
       for(auto &&a2 : irange(0L, _Nact)){
         J0_P.cptr()[P0*J0_P.shape(1)+a2] = Jint_ccav_ccaa_0_.cptr()[P*Jint_ccav_ccaa_0_.shape(1)*Jint_ccav_ccaa_0_.shape(2)+P0*Jint_ccav_ccaa_0_.shape(2)+a2];
         J1_P.cptr()[P0*J1_P.shape(1)+a2] = Jint_ccav_ccaa_1_.cptr()[P*Jint_ccav_ccaa_1_.shape(1)*Jint_ccav_ccaa_1_.shape(2)+P0*Jint_ccav_ccaa_1_.shape(2)+a2];
       }
     }
     orz::DTensor JF0_P(J0_P*fpa);
     JF0.PutTensor(P, JF0_P);
     orz::DTensor JF1_P(J1_P*fpa);
     JF1.PutTensor(P, JF1_P);
   }
 }
 
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(i,j,a,P) <<= +1 @Jint_ccav_ccaa_0_(P,P0,a2) @cF(a2,a) @T(+2)(j,i,P0)" << std::endl;
std::cout << boost::format("%80s") % "@S2(i,j,a,P) <<= +1 @Jint_ccav_ccaa_1_(P,P0,a2) @cF(a2,a) @T(+2)(i,j,P0)"; 
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,P) <<= +1 @Jint_ccav_ccaa_0_(P,P0,a2) @cF(a2,a) @T(+2)(j,i,P0)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor Fpa(_Fpa_.copy());
    {
      RActList.GetPair(i, j).Transform1EXT_R(Fpa);
    }    
    for(auto &&P : irange(0L, Nrhop1)){
      orz::DTensor JF0_ij = JF0.GetTensor(P);
      {
        RActList.GetPair(i, j).Transform1EXT_R(JF0_ij);
      }
      orz::DTensor JF1_ij = JF1.GetTensor(P);
      {
        RActList.GetPair(i, j).Transform1EXT_R(JF1_ij);
      }                
      for(auto &&P0 : irange(0L, Nrhop2)){
        const long I_P0 = TActList.GetIndex(P0);
        if (I_P0 < 0) continue;
        orz::DTensor U2 = T2.GetTensor(P0);
        for(auto &&a : irange(0L, _Nvirt_ij)) S2.cptr()[a*S2.shape(1)+P] += JF0_ij.cptr()[P0*JF0_ij.shape(1)+a] * U2.cptr()[j*U2.shape(1)+i];
        for(auto &&a : irange(0L, _Nvirt_ij)) S2.cptr()[a*S2.shape(1)+P] += JF1_ij.cptr()[P0*JF0_ij.shape(1)+a] * U2.cptr()[i*U2.shape(1)+j];        
      }//P0
    }//P
    R2.PutTensor(i, j, S2);
  }//j
}//i

high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
#ifdef _PROFILE_LHS
cout << boost::format(" ... %8.5e sec. (t_ovl ... %8.5e sec. )") % time_span.count() % t_trans << endl;
#endif
}
//{
//double t_trans=0.0;
//#ifdef _PROFILE_LHS
//std::cout << boost::format("%80s") % "@S2(i,j,a,P) <<= +1 @Jint_ccav_ccaa_1_(P,P0,a2) @cF(a2,a) @T(+2)(i,j,P0)";
//#endif
//high_resolution_clock::time_point t1 = high_resolution_clock::now();
//// @S2(i,j,a,P) <<= +1 @Jint_ccav_ccaa_1_(P,P0,a2) @cF(a2,a) @T(+2)(i,j,P0)
//for(auto &&i : irange(0L, _Ndocc)){
//  for(auto &&j : irange(0L, _Ndocc)){
//    const long I_ij = RActList.GetIndex(i, j);
//    if (I_ij < 0) continue;
//    orz::DTensor S2 = R2.GetTensor(i, j);
//    const long _Nvirt_ij = RActList.GetNPNO(i, j);
//    if (I_ij < 0) continue;
//    for(auto &&P0 : irange(0L, Nrhop2)){
//      const long I_P0 = TActList.GetIndex(P0);
//      if (I_P0 < 0) continue;
//      orz::DTensor U2 = T2.GetTensor(P0);
//      orz::DTensor Fpa(_Fpa_.copy());
//      {
//        RActList.GetPair(i, j).Transform1EXT_R(Fpa);
//      }
//      for(auto &&P : irange(0L, Nrhop1)){
//        for(auto &&a2 : irange(0L, _Nact)){
//          for(auto &&a : irange(0L, _Nvirt_ij)){
//            S2.cptr()[a*S2.shape(1)+P] += Jint_ccav_ccaa_1_.cptr()[P*Jint_ccav_ccaa_1_.shape(1)*Jint_ccav_ccaa_1_.shape(2)+P0*Jint_ccav_ccaa_1_.shape(2)+a2] * Fpa.cptr()[a2*Fpa.shape(1)+a] * U2.cptr()[i*U2.shape(1)+j];
//          }//a
//        }//a2
//      }//P
//    }//P0
//    R2.PutTensor(i, j, S2);
//  }//j
//}//i
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
