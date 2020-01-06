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
void FMatrixOrthogonalizer::FormLHS_V0_Vm1(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &d3, DensityPack &dPack,
                                            TensorContainer &ovl, const double &Ecas, const double &h0_energy, const orz::DTensor &casfock, const orz::DTensor &FV,
                                            TensorContainer &T2, TensorContainer &R2){
#include "lct_preparations.cpp"
orz::DTensor &Jint_ccvv_cavv_0_ = (*dPack["Jint_ccvv_cavv_0_"]);
orz::DTensor &Jint_ccvv_cavv_1_ = (*dPack["Jint_ccvv_cavv_1_"]);
orz::DTensor &Jint_ccvv_cavv_2_ = (*dPack["Jint_ccvv_cavv_2_"]);
orz::DTensor &Jint_ccvv_cavv_3_ = (*dPack["Jint_ccvv_cavv_3_"]);
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(i,j,a,b) <<= +1 @Jint_ccvv_cavv_0_(P0,i) @T(-1)(a,b,P0,j)" << std::endl;
std::cout << boost::format("%80s") % "@S2(i,j,a,b) <<= +1 @Jint_ccvv_cavv_1_(P0,i) @T(-1)(b,a,P0,j)"; 
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,b) <<= +1 @Jint_ccvv_cavv_0_(P0,i) @T(-1)(a,b,P0,j)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    for(auto &&P0 : irange(0L, Nrhom1)){
      const long I_P0j = TActList.GetIndex(P0, j);
      if (I_P0j < 0) continue;
      orz::DTensor U2 = T2.GetTensor(P0, j);
      {
        high_resolution_clock::time_point o1 = high_resolution_clock::now();
        orz::DTensor S = ovl.GetTensor(I_P0j,I_ij);
        orz::DTensor tmp;
        At_x_B_x_A(tmp, S, U2);
        U2 = tmp.copy();
        high_resolution_clock::time_point o2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
        t_trans += time_span.count();
      }
//      for(auto &&a : irange(0L, _Nvirt_ij)){
//        for(auto &&b : irange(0L, _Nvirt_ij)){
//          S2.cptr()[a*S2.shape(1)+b] += Jint_ccvv_cavv_0_.cptr()[P0*Jint_ccvv_cavv_0_.shape(1)+i] * U2.cptr()[a*U2.shape(1)+b];
//        }//b
//      }//a
//      for(auto &&b : irange(0L, _Nvirt_ij)){
//        for(auto &&a : irange(0L, _Nvirt_ij)){
//          S2.cptr()[a*S2.shape(1)+b] += Jint_ccvv_cavv_1_.cptr()[P0*Jint_ccvv_cavv_1_.shape(1)+i] * U2.cptr()[b*U2.shape(1)+a];
//        }//a
//      }//b
      S2.gaxpy(+1.0, U2             , Jint_ccvv_cavv_0_.cptr()[P0*Jint_ccvv_cavv_0_.shape(1)+i]);
      S2.gaxpy(+1.0, U2.swapdim(1,0), Jint_ccvv_cavv_1_.cptr()[P0*Jint_ccvv_cavv_1_.shape(1)+i]);
    }//P0
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
//std::cout << boost::format("%80s") % "@S2(i,j,a,b) <<= +1 @Jint_ccvv_cavv_1_(P0,i) @T(-1)(b,a,P0,j)";
//#endif
//high_resolution_clock::time_point t1 = high_resolution_clock::now();
//// @S2(i,j,a,b) <<= +1 @Jint_ccvv_cavv_1_(P0,i) @T(-1)(b,a,P0,j)
//for(auto &&i : irange(0L, _Ndocc)){
//  for(auto &&j : irange(0L, _Ndocc)){
//    const long I_ij = RActList.GetIndex(i, j);
//    if (I_ij < 0) continue;
//    orz::DTensor S2 = R2.GetTensor(i, j);
//    const long _Nvirt_ij = RActList.GetNPNO(i, j);
//    if (I_ij < 0) continue;
//    for(auto &&P0 : irange(0L, Nrhom1)){
//      const long I_P0j = TActList.GetIndex(P0, j);
//      if (I_P0j < 0) continue;
//      orz::DTensor U2 = T2.GetTensor(P0, j);
//      {
//        high_resolution_clock::time_point o1 = high_resolution_clock::now();
//        orz::DTensor S = ovl.GetTensor(I_P0j,I_ij);
//        orz::DTensor tmp;
//        At_x_B_x_A(tmp, S, U2);
//        U2 = tmp.copy();
//        high_resolution_clock::time_point o2 = high_resolution_clock::now();
//        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
//        t_trans += time_span.count();
//      }
//      for(auto &&b : irange(0L, _Nvirt_ij)){
//        for(auto &&a : irange(0L, _Nvirt_ij)){
//          S2.cptr()[a*S2.shape(1)+b] += Jint_ccvv_cavv_1_.cptr()[P0*Jint_ccvv_cavv_1_.shape(1)+i] * U2.cptr()[b*U2.shape(1)+a];
//        }//a
//      }//b
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
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(i,j,a,b) <<= +1 @Jint_ccvv_cavv_2_(P0,j) @T(-1)(a,b,P0,i)" << std::endl;
std::cout << boost::format("%80s") % "@S2(i,j,a,b) <<= +1 @Jint_ccvv_cavv_3_(P0,j) @T(-1)(b,a,P0,i)"; 
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,b) <<= +1 @Jint_ccvv_cavv_2_(P0,j) @T(-1)(a,b,P0,i)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    for(auto &&P0 : irange(0L, Nrhom1)){
      const long I_P0i = TActList.GetIndex(P0, i);
      if (I_P0i < 0) continue;
      orz::DTensor U2 = T2.GetTensor(P0, i);
      {
        high_resolution_clock::time_point o1 = high_resolution_clock::now();
        orz::DTensor S = ovl.GetTensor(I_P0i,I_ij);
        orz::DTensor tmp;
        At_x_B_x_A(tmp, S, U2);
        U2 = tmp.copy();
        high_resolution_clock::time_point o2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
        t_trans += time_span.count();
      }
//      for(auto &&a : irange(0L, _Nvirt_ij)){
//        for(auto &&b : irange(0L, _Nvirt_ij)){
//          S2.cptr()[a*S2.shape(1)+b] += Jint_ccvv_cavv_2_.cptr()[P0*Jint_ccvv_cavv_2_.shape(1)+j] * U2.cptr()[a*U2.shape(1)+b];
//        }//b
//      }//a
//      for(auto &&b : irange(0L, _Nvirt_ij)){
//        for(auto &&a : irange(0L, _Nvirt_ij)){
//          S2.cptr()[a*S2.shape(1)+b] += Jint_ccvv_cavv_3_.cptr()[P0*Jint_ccvv_cavv_3_.shape(1)+j] * U2.cptr()[b*U2.shape(1)+a];
//        }//a
//      }//b
      S2.gaxpy(+1.0, U2             , Jint_ccvv_cavv_2_.cptr()[P0*Jint_ccvv_cavv_2_.shape(1)+j]);
      S2.gaxpy(+1.0, U2.swapdim(1,0), Jint_ccvv_cavv_3_.cptr()[P0*Jint_ccvv_cavv_3_.shape(1)+j]);      
    }//P0
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
//std::cout << boost::format("%80s") % "@S2(i,j,a,b) <<= +1 @Jint_ccvv_cavv_3_(P0,j) @T(-1)(b,a,P0,i)";
//#endif
//high_resolution_clock::time_point t1 = high_resolution_clock::now();
//// @S2(i,j,a,b) <<= +1 @Jint_ccvv_cavv_3_(P0,j) @T(-1)(b,a,P0,i)
//for(auto &&i : irange(0L, _Ndocc)){
//  for(auto &&j : irange(0L, _Ndocc)){
//    const long I_ij = RActList.GetIndex(i, j);
//    if (I_ij < 0) continue;
//    orz::DTensor S2 = R2.GetTensor(i, j);
//    const long _Nvirt_ij = RActList.GetNPNO(i, j);
//    if (I_ij < 0) continue;
//    for(auto &&P0 : irange(0L, Nrhom1)){
//      const long I_P0i = TActList.GetIndex(P0, i);
//      if (I_P0i < 0) continue;
//      orz::DTensor U2 = T2.GetTensor(P0, i);
//      {
//        high_resolution_clock::time_point o1 = high_resolution_clock::now();
//        orz::DTensor S = ovl.GetTensor(I_P0i,I_ij);
//        orz::DTensor tmp;
//        At_x_B_x_A(tmp, S, U2);
//        U2 = tmp.copy();
//        high_resolution_clock::time_point o2 = high_resolution_clock::now();
//        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
//        t_trans += time_span.count();
//      }
//      for(auto &&b : irange(0L, _Nvirt_ij)){
//        for(auto &&a : irange(0L, _Nvirt_ij)){
//          S2.cptr()[a*S2.shape(1)+b] += Jint_ccvv_cavv_3_.cptr()[P0*Jint_ccvv_cavv_3_.shape(1)+j] * U2.cptr()[b*U2.shape(1)+a];
//        }//a
//      }//b
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
