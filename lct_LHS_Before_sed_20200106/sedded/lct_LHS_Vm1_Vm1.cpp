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
void FMatrixOrthogonalizer::FormLHS_Vm1_Vm1(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &d3, DensityPack &dPack,
                                            TensorContainer &ovl, const double &Ecas, const double &h0_energy, const orz::DTensor &casfock, const orz::DTensor &FV,
                                            TensorContainer &T2, TensorContainer &R2){
#include "lct_preparations.cpp"
// Init amplitude container                               
for(auto &&P : irange(0L, Nrhom1)){                       
  for(auto &&i : irange(0L, _Ndocc)){                     
    const long idx_Pi = RActList.GetIndex(P,i);
    if (idx_Pi < 0) continue;
    const long Nvirt_Pi = RActList.GetPair(P,i).GetNPNO();
    orz::DTensor Rpi(Nvirt_Pi, Nvirt_Pi);                 
    R2.PutTensor(P,i,Rpi);                                
  }                                                       
}                                                         
TensorContainer TFcc(Nrhom1,_Ndocc);
for(auto &&P0 : irange(0L, Nrhom1)){
  for(auto &&i : irange(0L, _Ndocc)){
    const long I_P0i = RActList.GetIndex(P0, i);
    if (I_P0i < 0) continue;
    const long _Nvirt_P0i = TActList.GetNPNO(P0, i);
    orz::DTensor X(_Nvirt_P0i,_Nvirt_P0i);
    for(auto &&c0 : irange(0L, _Ndocc)){
      const long I_P0c0 = RActList.GetIndex(P0, c0);
      if (I_P0c0 < 0) continue;
      orz::DTensor U2 = T2.GetTensor(P0, c0);
      {
        high_resolution_clock::time_point o1 = high_resolution_clock::now();
        orz::DTensor S = ovl.CopyTensor(std::max(I_P0i,I_P0c0),std::min(I_P0i,I_P0c0));
        if(I_P0c0>I_P0i) orz::tensor::transpose(S);
        orz::DTensor tmp;
        A_x_B_x_At(tmp, S, U2);
        U2 = tmp.copy();
        high_resolution_clock::time_point o2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
        //t_trans += time_span.count();
      }
      //X += Fij.cptr()[i*Fij.shape(1)+c0] * U2;
      X.gaxpy(+1.0, U2, Fij.cptr()[i*Fij.shape(1)+c0]);
    }//c0
    TFcc.PutTensor(P0,i,X);
  }//i
}//P0
 TensorContainer FvvT(Nrhom1,_Ndocc), TFvv(Nrhom1,_Ndocc);
for(auto &&P : irange(0L, Nrhom1)){
  for(auto &&i : irange(0L, _Ndocc)){
    const long I_Pi = RActList.GetIndex(P, i);
    if (I_Pi < 0) continue;
    const long _Nvirt_Pi = RActList.GetNPNO(P, i);
    orz::DTensor U2 = T2.GetTensor(P, i);
    orz::DTensor Fab(_Fab_.copy());
    RActList.GetPair(P, i).Transform2EXT(Fab);
    {
      orz::DTensor X(_Nvirt_Pi,_Nvirt_Pi);
//      for(auto &&b : irange(0L, _Nvirt_Pi)){
//        for(auto &&a : irange(0L, _Nvirt_Pi)){
//          for(auto &&v0 : irange(0L, _Nvirt_Pi)){
//            X.cptr()[a*X.shape(1)+b] += Fab.cptr()[v0*Fab.shape(1)+a] * U2.cptr()[v0*U2.shape(1)+b];
//          }
//        }//a
//      }//b
      A_x_B(X, Fab, U2, false, false);
      FvvT.PutTensor(P,i,X);
    }
    {
      orz::DTensor X(_Nvirt_Pi,_Nvirt_Pi);
//      for(auto &&b : irange(0L, _Nvirt_Pi)){
//        for(auto &&a : irange(0L, _Nvirt_Pi)){
//          for(auto &&v0 : irange(0L, _Nvirt_Pi)){
//            X.cptr()[a*X.shape(1)+b] += U2.cptr()[a*U2.shape(1)+v0] * Fab.cptr()[v0*Fab.shape(1)+b];
//          }
//        }//a
//      }//b
      A_x_B(X, U2, Fab, false, false);
      TFvv.PutTensor(P,i,X);
    }    
  }
}
//TensorContainer TFvv(Nrhom1,_Ndocc);
//for(auto &&P : irange(0L, Nrhom1)){
//  for(auto &&i : irange(0L, _Ndocc)){
//    const long I_Pi = RActList.GetIndex(P, i);
//    if (I_Pi < 0) continue;
//    const long _Nvirt_Pi = RActList.GetNPNO(P, i);
//    orz::DTensor U2 = T2.GetTensor(P, i);
//    orz::DTensor Fab(_Fab_.copy());
//    RActList.GetPair(P, i).Transform2EXT(Fab);
//    {
//      orz::DTensor X(_Nvirt_Pi,_Nvirt_Pi);
////      for(auto &&b : irange(0L, _Nvirt_Pi)){
////        for(auto &&a : irange(0L, _Nvirt_Pi)){
////          for(auto &&v0 : irange(0L, _Nvirt_Pi)){
////            X.cptr()[a*X.shape(1)+b] += U2.cptr()[a*U2.shape(1)+v0] * Fab.cptr()[v0*Fab.shape(1)+b];
////          }
////        }//a
////      }//b
//      A_x_B(X, U2, Fab, false, false);
//      TFvv.PutTensor(P,i,X);
//    }
//  }
//}
orz::DTensor &Jint_cavv_cavv_0_ = (*dPack["Jint_cavv_cavv_0_"]);
orz::DTensor &Jint_cavv_cavv_1_ = (*dPack["Jint_cavv_cavv_1_"]);
orz::DTensor &Jint_cavv_cavv_2_ = (*dPack["Jint_cavv_cavv_2_"]);
orz::DTensor &Jint_cavv_cavv_3_ = (*dPack["Jint_cavv_cavv_3_"]);
orz::DTensor &Jint_cavv_cavv_4_ = (*dPack["Jint_cavv_cavv_4_"]);
orz::DTensor &Jint_cavv_cavv_5_ = (*dPack["Jint_cavv_cavv_5_"]);
orz::DTensor &Jint_cavv_cavv_6_ = (*dPack["Jint_cavv_cavv_6_"]);
orz::DTensor &Jint_cavv_cavv_7_ = (*dPack["Jint_cavv_cavv_7_"]);
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(P,i,a,b) <<= +1 @Jint_cavv_cavv_0_(P,P0) @T(-1)(a,b,P0,i)" << std::endl;
std::cout << boost::format("%80s") % "@S2(P,i,a,b) <<= +1 @Jint_cavv_cavv_1_(P,P0) @T(-1)(b,a,P0,i)"; 
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(P,i,a,b) <<= +1 @Jint_cavv_cavv_0_(P,P0) @T(-1)(a,b,P0,i)
for(auto &&P : irange(0L, Nrhom1)){
  for(auto &&i : irange(0L, _Ndocc)){
    const long I_Pi = RActList.GetIndex(P, i);
    if (I_Pi < 0) continue;
    orz::DTensor S2 = R2.GetTensor(P, i);
    const long _Nvirt_Pi = RActList.GetNPNO(P, i);
    if (I_Pi < 0) continue;
    for(auto &&P0 : irange(0L, Nrhom1)){
      const long I_P0i = TActList.GetIndex(P0, i);
      if (I_P0i < 0) continue;
      orz::DTensor U2 = T2.GetTensor(P0, i);
      {
        high_resolution_clock::time_point o1 = high_resolution_clock::now();
        orz::DTensor S = ovl.CopyTensor(std::max(I_Pi,I_P0i),std::min(I_Pi,I_P0i));
        if(I_P0i>I_Pi) orz::tensor::transpose(S);
        orz::DTensor tmp;
        A_x_B_x_At(tmp, S, U2);
        U2 = tmp.copy();
        high_resolution_clock::time_point o2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
        t_trans += time_span.count();
      }
//      for(auto &&a : irange(0L, _Nvirt_Pi)){
//        for(auto &&b : irange(0L, _Nvirt_Pi)){
//          S2.cptr()[a*S2.shape(1)+b] += Jint_cavv_cavv_0_.cptr()[P*Jint_cavv_cavv_0_.shape(1)+P0] * U2.cptr()[a*U2.shape(1)+b];
//        }//b
//      }//a
//      for(auto &&b : irange(0L, _Nvirt_Pi)){
//        for(auto &&a : irange(0L, _Nvirt_Pi)){
//          S2.cptr()[a*S2.shape(1)+b] += Jint_cavv_cavv_1_.cptr()[P*Jint_cavv_cavv_1_.shape(1)+P0] * U2.cptr()[b*U2.shape(1)+a];
//        }//a
//      }//b
      S2.gaxpy(+1.0, U2             , Jint_cavv_cavv_0_.cptr()[P*Jint_cavv_cavv_0_.shape(1)+P0]);
      S2.gaxpy(+1.0, U2.swapdim(1,0), Jint_cavv_cavv_1_.cptr()[P*Jint_cavv_cavv_1_.shape(1)+P0]);
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
//std::cout << boost::format("%80s") % "@S2(P,i,a,b) <<= +1 @Jint_cavv_cavv_1_(P,P0) @T(-1)(b,a,P0,i)";
//#endif
//high_resolution_clock::time_point t1 = high_resolution_clock::now();
//// @S2(P,i,a,b) <<= +1 @Jint_cavv_cavv_1_(P,P0) @T(-1)(b,a,P0,i)
//for(auto &&P : irange(0L, Nrhom1)){
//  for(auto &&i : irange(0L, _Ndocc)){
//    const long I_Pi = RActList.GetIndex(P, i);
//    if (I_Pi < 0) continue;
//    orz::DTensor S2 = R2.GetTensor(P, i);
//    const long _Nvirt_Pi = RActList.GetNPNO(P, i);
//    if (I_Pi < 0) continue;
//    for(auto &&P0 : irange(0L, Nrhom1)){
//      const long I_P0i = TActList.GetIndex(P0, i);
//      if (I_P0i < 0) continue;
//      orz::DTensor U2 = T2.GetTensor(P0, i);
//      {
//        high_resolution_clock::time_point o1 = high_resolution_clock::now();
//        orz::DTensor S = ovl.CopyTensor(std::max(I_Pi,I_P0i),std::min(I_Pi,I_P0i));
//        if(I_P0i>I_Pi) orz::tensor::transpose(S);
//        orz::DTensor tmp;
//        A_x_B_x_At(tmp, S, U2);
//        U2 = tmp.copy();
//        high_resolution_clock::time_point o2 = high_resolution_clock::now();
//        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
//        t_trans += time_span.count();
//      }
//      for(auto &&b : irange(0L, _Nvirt_Pi)){
//        for(auto &&a : irange(0L, _Nvirt_Pi)){
//          S2.cptr()[a*S2.shape(1)+b] += Jint_cavv_cavv_1_.cptr()[P*Jint_cavv_cavv_1_.shape(1)+P0] * U2.cptr()[b*U2.shape(1)+a];
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
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(P,i,a,b) <<= +1 @TFcc(a,b,P0,i) @Jint_cavv_cavv_2_(P,P0)" << std::endl;
std::cout << boost::format("%80s") % "@S2(P,i,a,b) <<= +1 @TFcc(b,a,P0,i) @Jint_cavv_cavv_3_(P,P0)"; 
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(P,i,a,b) <<= +1 @TFcc(a,b,P0,i) @Jint_cavv_cavv_2_(P,P0)
for(auto &&P : irange(0L, Nrhom1)){
  for(auto &&i : irange(0L, _Ndocc)){
    const long I_Pi = RActList.GetIndex(P, i);
    if (I_Pi < 0) continue;
    orz::DTensor S2 = R2.GetTensor(P, i);
    const long _Nvirt_Pi = RActList.GetNPNO(P, i);
    if (I_Pi < 0) continue;
    for(auto &&P0 : irange(0L, Nrhom1)){
      const long I_P0i = TActList.GetIndex(P0, i);
      if (I_P0i < 0) continue;
      orz::DTensor U2 = TFcc.GetTensor(P0, i);
      {
        high_resolution_clock::time_point o1 = high_resolution_clock::now();
        orz::DTensor S = ovl.CopyTensor(std::max(I_Pi,I_P0i),std::min(I_Pi,I_P0i));
        if(I_P0i>I_Pi) orz::tensor::transpose(S);
        orz::DTensor tmp;
        A_x_B_x_At(tmp, S, U2);
        U2 = tmp.copy();
        high_resolution_clock::time_point o2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
        t_trans += time_span.count();
      }
//      for(auto &&a : irange(0L, _Nvirt_Pi)){
//        for(auto &&b : irange(0L, _Nvirt_Pi)){
//          S2.cptr()[a*S2.shape(1)+b] += U2.cptr()[a*U2.shape(1)+b] * Jint_cavv_cavv_2_.cptr()[P*Jint_cavv_cavv_2_.shape(1)+P0];
//        }//b
//      }//a
//      for(auto &&b : irange(0L, _Nvirt_Pi)){
//        for(auto &&a : irange(0L, _Nvirt_Pi)){
//          S2.cptr()[a*S2.shape(1)+b] += U2.cptr()[b*U2.shape(1)+a] * Jint_cavv_cavv_3_.cptr()[P*Jint_cavv_cavv_3_.shape(1)+P0];
//        }//a
//      }//b
      S2.gaxpy(+1.0, U2             , Jint_cavv_cavv_2_.cptr()[P*Jint_cavv_cavv_2_.shape(1)+P0]);
      S2.gaxpy(+1.0, U2.swapdim(1,0), Jint_cavv_cavv_3_.cptr()[P*Jint_cavv_cavv_3_.shape(1)+P0]);
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
//std::cout << boost::format("%80s") % "@S2(P,i,a,b) <<= +1 @TFcc(b,a,P0,i) @Jint_cavv_cavv_3_(P,P0)";
//#endif
//high_resolution_clock::time_point t1 = high_resolution_clock::now();
//// @S2(P,i,a,b) <<= +1 @TFcc(b,a,P0,i) @Jint_cavv_cavv_3_(P,P0)
//for(auto &&P : irange(0L, Nrhom1)){
//  for(auto &&i : irange(0L, _Ndocc)){
//    const long I_Pi = RActList.GetIndex(P, i);
//    if (I_Pi < 0) continue;
//    orz::DTensor S2 = R2.GetTensor(P, i);
//    const long _Nvirt_Pi = RActList.GetNPNO(P, i);
//    if (I_Pi < 0) continue;
//    for(auto &&P0 : irange(0L, Nrhom1)){
//      const long I_P0i = TActList.GetIndex(P0, i);
//      if (I_P0i < 0) continue;
//      orz::DTensor U2 = TFcc.GetTensor(P0, i);
//      {
//        high_resolution_clock::time_point o1 = high_resolution_clock::now();
//        orz::DTensor S = ovl.CopyTensor(std::max(I_Pi,I_P0i),std::min(I_Pi,I_P0i));
//        if(I_P0i>I_Pi) orz::tensor::transpose(S);
//        orz::DTensor tmp;
//        A_x_B_x_At(tmp, S, U2);
//        U2 = tmp.copy();
//        high_resolution_clock::time_point o2 = high_resolution_clock::now();
//        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
//        t_trans += time_span.count();
//      }
//      for(auto &&b : irange(0L, _Nvirt_Pi)){
//        for(auto &&a : irange(0L, _Nvirt_Pi)){
//          S2.cptr()[a*S2.shape(1)+b] += U2.cptr()[b*U2.shape(1)+a] * Jint_cavv_cavv_3_.cptr()[P*Jint_cavv_cavv_3_.shape(1)+P0];
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
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(P,i,a,b) <<= +1 @FvvT(b,a,P0,i) @Jint_cavv_cavv_4_(P,P0)" << std::endl;
std::cout << boost::format("%80s") % "@S2(P,i,a,b) <<= +1 @FvvT(a,b,P0,i) @Jint_cavv_cavv_6_(P,P0)"; 
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(P,i,a,b) <<= +1 @FvvT(b,a,P0,i) @Jint_cavv_cavv_4_(P,P0)
for(auto &&P : irange(0L, Nrhom1)){
  for(auto &&i : irange(0L, _Ndocc)){
    const long I_Pi = RActList.GetIndex(P, i);
    if (I_Pi < 0) continue;
    orz::DTensor S2 = R2.GetTensor(P, i);
    const long _Nvirt_Pi = RActList.GetNPNO(P, i);
    if (I_Pi < 0) continue;
    for(auto &&P0 : irange(0L, Nrhom1)){
      const long I_P0i = TActList.GetIndex(P0, i);
      if (I_P0i < 0) continue;
      orz::DTensor U2 = FvvT.GetTensor(P0, i);
      {
        high_resolution_clock::time_point o1 = high_resolution_clock::now();
        orz::DTensor S = ovl.CopyTensor(std::max(I_Pi,I_P0i),std::min(I_Pi,I_P0i));
        if(I_P0i>I_Pi) orz::tensor::transpose(S);
        orz::DTensor tmp;
        A_x_B_x_At(tmp, S, U2);
        U2 = tmp.copy();
        high_resolution_clock::time_point o2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
        t_trans += time_span.count();
      }
//      for(auto &&b : irange(0L, _Nvirt_Pi)){
//        for(auto &&a : irange(0L, _Nvirt_Pi)){
//          S2.cptr()[a*S2.shape(1)+b] += U2.cptr()[b*U2.shape(1)+a] * Jint_cavv_cavv_4_.cptr()[P*Jint_cavv_cavv_4_.shape(1)+P0];
//        }//a
//      }//b
//      for(auto &&a : irange(0L, _Nvirt_Pi)){
//        for(auto &&b : irange(0L, _Nvirt_Pi)){
//          S2.cptr()[a*S2.shape(1)+b] += U2.cptr()[a*U2.shape(1)+b] * Jint_cavv_cavv_6_.cptr()[P*Jint_cavv_cavv_6_.shape(1)+P0];
//        }//b
//      }//a
      S2.gaxpy(+1.0, U2.swapdim(1,0), Jint_cavv_cavv_4_.cptr()[P*Jint_cavv_cavv_4_.shape(1)+P0]);
      S2.gaxpy(+1.0, U2             , Jint_cavv_cavv_6_.cptr()[P*Jint_cavv_cavv_6_.shape(1)+P0]);      
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
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(P,i,a,b) <<= +1 @TFvv(a,b,P0,i) @Jint_cavv_cavv_5_(P,P0)" << std::endl;
std::cout << boost::format("%80s") % "@S2(P,i,a,b) <<= +1 @TFvv(b,a,P0,i) @Jint_cavv_cavv_7_(P,P0)"; 
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(P,i,a,b) <<= +1 @TFvv(a,b,P0,i) @Jint_cavv_cavv_5_(P,P0)
for(auto &&P : irange(0L, Nrhom1)){
  for(auto &&i : irange(0L, _Ndocc)){
    const long I_Pi = RActList.GetIndex(P, i);
    if (I_Pi < 0) continue;
    orz::DTensor S2 = R2.GetTensor(P, i);
    const long _Nvirt_Pi = RActList.GetNPNO(P, i);
    if (I_Pi < 0) continue;
    for(auto &&P0 : irange(0L, Nrhom1)){
      const long I_P0i = TActList.GetIndex(P0, i);
      if (I_P0i < 0) continue;
      orz::DTensor U2 = TFvv.GetTensor(P0, i);
      {
        high_resolution_clock::time_point o1 = high_resolution_clock::now();
        orz::DTensor S = ovl.CopyTensor(std::max(I_Pi,I_P0i),std::min(I_Pi,I_P0i));
        if(I_P0i>I_Pi) orz::tensor::transpose(S);
        orz::DTensor tmp;
        A_x_B_x_At(tmp, S, U2);
        U2 = tmp.copy();
        high_resolution_clock::time_point o2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
        t_trans += time_span.count();
      }
//      for(auto &&a : irange(0L, _Nvirt_Pi)){
//        for(auto &&b : irange(0L, _Nvirt_Pi)){
//          S2.cptr()[a*S2.shape(1)+b] += U2.cptr()[a*U2.shape(1)+b] * Jint_cavv_cavv_5_.cptr()[P*Jint_cavv_cavv_5_.shape(1)+P0];
//        }//b
//      }//a
//      for(auto &&b : irange(0L, _Nvirt_Pi)){
//        for(auto &&a : irange(0L, _Nvirt_Pi)){
//          S2.cptr()[a*S2.shape(1)+b] += U2.cptr()[b*U2.shape(1)+a] * Jint_cavv_cavv_7_.cptr()[P*Jint_cavv_cavv_7_.shape(1)+P0];
//        }//a
//      }//b
      S2.gaxpy(+1.0, U2             , Jint_cavv_cavv_5_.cptr()[P*Jint_cavv_cavv_5_.shape(1)+P0]);
      S2.gaxpy(+1.0, U2.swapdim(1,0), Jint_cavv_cavv_7_.cptr()[P*Jint_cavv_cavv_7_.shape(1)+P0]);
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
//std::cout << boost::format("%80s") % "@S2(P,i,a,b) <<= +1 @FvvT(a,b,P0,i) @Jint_cavv_cavv_6_(P,P0)";
//#endif
//high_resolution_clock::time_point t1 = high_resolution_clock::now();
//// @S2(P,i,a,b) <<= +1 @FvvT(a,b,P0,i) @Jint_cavv_cavv_6_(P,P0)
//for(auto &&P : irange(0L, Nrhom1)){
//  for(auto &&i : irange(0L, _Ndocc)){
//    const long I_Pi = RActList.GetIndex(P, i);
//    if (I_Pi < 0) continue;
//    orz::DTensor S2 = R2.GetTensor(P, i);
//    const long _Nvirt_Pi = RActList.GetNPNO(P, i);
//    if (I_Pi < 0) continue;
//    for(auto &&P0 : irange(0L, Nrhom1)){
//      const long I_P0i = TActList.GetIndex(P0, i);
//      if (I_P0i < 0) continue;
//      orz::DTensor U2 = FvvT.GetTensor(P0, i);
//      {
//        high_resolution_clock::time_point o1 = high_resolution_clock::now();
//        orz::DTensor S = ovl.CopyTensor(std::max(I_Pi,I_P0i),std::min(I_Pi,I_P0i));
//        if(I_P0i>I_Pi) orz::tensor::transpose(S);
//        orz::DTensor tmp;
//        A_x_B_x_At(tmp, S, U2);
//        U2 = tmp.copy();
//        high_resolution_clock::time_point o2 = high_resolution_clock::now();
//        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
//        t_trans += time_span.count();
//      }
//      for(auto &&a : irange(0L, _Nvirt_Pi)){
//        for(auto &&b : irange(0L, _Nvirt_Pi)){
//          S2.cptr()[a*S2.shape(1)+b] += U2.cptr()[a*U2.shape(1)+b] * Jint_cavv_cavv_6_.cptr()[P*Jint_cavv_cavv_6_.shape(1)+P0];
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
//std::cout << boost::format("%80s") % "@S2(P,i,a,b) <<= +1 @TFvv(b,a,P0,i) @Jint_cavv_cavv_7_(P,P0)";
//#endif
//high_resolution_clock::time_point t1 = high_resolution_clock::now();
//// @S2(P,i,a,b) <<= +1 @TFvv(b,a,P0,i) @Jint_cavv_cavv_7_(P,P0)
//for(auto &&P : irange(0L, Nrhom1)){
//  for(auto &&i : irange(0L, _Ndocc)){
//    const long I_Pi = RActList.GetIndex(P, i);
//    if (I_Pi < 0) continue;
//    orz::DTensor S2 = R2.GetTensor(P, i);
//    const long _Nvirt_Pi = RActList.GetNPNO(P, i);
//    if (I_Pi < 0) continue;
//    for(auto &&P0 : irange(0L, Nrhom1)){
//      const long I_P0i = TActList.GetIndex(P0, i);
//      if (I_P0i < 0) continue;
//      orz::DTensor U2 = TFvv.GetTensor(P0, i);
//      {
//        high_resolution_clock::time_point o1 = high_resolution_clock::now();
//        orz::DTensor S = ovl.CopyTensor(std::max(I_Pi,I_P0i),std::min(I_Pi,I_P0i));
//        if(I_P0i>I_Pi) orz::tensor::transpose(S);
//        orz::DTensor tmp;
//        A_x_B_x_At(tmp, S, U2);
//        U2 = tmp.copy();
//        high_resolution_clock::time_point o2 = high_resolution_clock::now();
//        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
//        t_trans += time_span.count();
//      }
//      for(auto &&b : irange(0L, _Nvirt_Pi)){
//        for(auto &&a : irange(0L, _Nvirt_Pi)){
//          S2.cptr()[a*S2.shape(1)+b] += U2.cptr()[b*U2.shape(1)+a] * Jint_cavv_cavv_7_.cptr()[P*Jint_cavv_cavv_7_.shape(1)+P0];
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
