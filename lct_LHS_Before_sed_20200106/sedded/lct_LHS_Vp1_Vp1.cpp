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
void FMatrixOrthogonalizer::FormLHS_Vp1_Vp1(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &d3, DensityPack &dPack,
                                            TensorContainer &ovl, const double &Ecas, const double &h0_energy, const orz::DTensor &casfock, const orz::DTensor &FV,
                                            TensorContainer &T2, TensorContainer &R2){
#include "lct_preparations.cpp"
for(auto &&i : irange(0L, _Ndocc)){          
  for(auto &&j : irange(0L, _Ndocc)){        
    const long I_ij = RActList.GetIndex(i,j);
    if (I_ij < 0) continue;                  
    PairData &P_ij = RActList.GetPair(i,j);  
    const int Nvirt = P_ij.GetNPNO();        
    orz::DTensor tmp(Nvirt,Nrhop1);          
    R2.PutTensor(i,j,tmp);                   
  }                                          
}                                            
TensorContainer TFcc(_Ndocc,_Ndocc);
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor X(_Nvirt_ij, Nrhop1);
    for(auto &&c0 : irange(0L, _Ndocc)){
      const long I_c0i = TActList.GetIndex(c0, i);
      if (I_c0i < 0) continue;
      orz::DTensor U2 = T2.GetTensor(i, c0);
      {
        high_resolution_clock::time_point o1 = high_resolution_clock::now();
        orz::DTensor S = ovl.CopyTensor(std::max(I_ij,I_c0i),std::min(I_ij,I_c0i));
        if(I_c0i>I_ij) orz::tensor::transpose(S);
        orz::DTensor tmp(S * U2);
        U2 = tmp.copy();
        high_resolution_clock::time_point o2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
        //t_trans += time_span.count();
      }
//      for(auto &&P0 : irange(0L, Nrhop1)){
//        for(auto &&a : irange(0L, _Nvirt_ij)){
//          X.cptr()[a*X.shape(1)+P0] += Fij.cptr()[j*Fij.shape(1)+c0] * U2.cptr()[a*U2.shape(1)+P0];
//        }//a
//      }//P0
      X.gaxpy(+1.0, U2, Fij.cptr()[j*Fij.shape(1)+c0]);
    }//c0
    TFcc.PutTensor(i,j,X);
  }//j
}//i
TensorContainer FccT(_Ndocc,_Ndocc);
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor X(_Nvirt_ij, Nrhop1);
    for(auto &&c0 : irange(0L, _Ndocc)){
      const long I_c0j = TActList.GetIndex(c0, j);
      if (I_c0j < 0) continue;
      orz::DTensor U2 = T2.GetTensor(c0,j);
      {
        high_resolution_clock::time_point o1 = high_resolution_clock::now();
        orz::DTensor S = ovl.CopyTensor(std::max(I_ij,I_c0j),std::min(I_ij,I_c0j));
        if(I_c0j>I_ij) orz::tensor::transpose(S);
        orz::DTensor tmp(S * U2);
        U2 = tmp.copy();
        high_resolution_clock::time_point o2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
        //t_trans += time_span.count();
      }
//      for(auto &&P0 : irange(0L, Nrhop1)){
//        for(auto &&a : irange(0L, _Nvirt_ij)){
//          X.cptr()[a*X.shape(1)+P0] += Fij.cptr()[i*Fij.shape(1)+c0] * U2.cptr()[a*U2.shape(1)+P0];
//        }//a
//      }//P0
      X.gaxpy(+1.0, U2, Fij.cptr()[i*Fij.shape(1)+c0]);
    }//c0
    FccT.PutTensor(i,j,X);
  }//j
}//i
TensorContainer FvvT(_Ndocc,_Ndocc);
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = T2.GetTensor(i, j);
    orz::DTensor Fab(_Fab_.copy());
    {
      RActList.GetPair(i, j).Transform2EXT(Fab);
    }
    orz::DTensor X(_Nvirt_ij,Nrhop1);
//    for(auto &&P0 : irange(0L, Nrhop1)){
//      for(auto &&v0 : irange(0L, _Nvirt_ij)){
//        for(auto &&a : irange(0L, _Nvirt_ij)){
//          X.cptr()[a*X.shape(1)+P0] += Fab.cptr()[v0*Fab.shape(1)+a] * U2.cptr()[v0*U2.shape(1)+P0];
//        }//a
//      }//v0
//    }//P0
    A_x_B(X, Fab, U2, true, false);
    FvvT.PutTensor(i, j, X);
  }//j
}//i
orz::DTensor &Jint_ccav_ccav_0_ = (*dPack["Jint_ccav_ccav_0_"]);
orz::DTensor &Jint_ccav_ccav_1_ = (*dPack["Jint_ccav_ccav_1_"]);
orz::DTensor &Jint_ccav_ccav_2_ = (*dPack["Jint_ccav_ccav_2_"]);
orz::DTensor &Jint_ccav_ccav_3_ = (*dPack["Jint_ccav_ccav_3_"]);
orz::DTensor &Jint_ccav_ccav_4_ = (*dPack["Jint_ccav_ccav_4_"]);
orz::DTensor &Jint_ccav_ccav_5_ = (*dPack["Jint_ccav_ccav_5_"]);
orz::DTensor &Jint_ccav_ccav_6_ = (*dPack["Jint_ccav_ccav_6_"]);
orz::DTensor &Jint_ccav_ccav_7_ = (*dPack["Jint_ccav_ccav_7_"]);
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(i,j,a,P) <<= +1 @TFcc(a,P0,i,j) @Jint_ccav_ccav_0_(P,P0)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,P) <<= +1 @TFcc(a,P0,i,j) @Jint_ccav_ccav_0_(P,P0)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = TFcc.GetTensor(i, j);
    for(auto &&P0 : irange(0L, Nrhop1)){
      for(auto &&P : irange(0L, Nrhop1)){
        for(auto &&a : irange(0L, _Nvirt_ij)){
          S2.cptr()[a*S2.shape(1)+P] += U2.cptr()[a*U2.shape(1)+P0] * Jint_ccav_ccav_0_.cptr()[P*Jint_ccav_ccav_0_.shape(1)+P0];
        }//a
      }//P
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
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(i,j,a,P) <<= +1 @FccT(a,P0,j,i) @Jint_ccav_ccav_1_(P,P0)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,P) <<= +1 @FccT(a,P0,j,i) @Jint_ccav_ccav_1_(P,P0)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = FccT.GetTensor(j, i);
    for(auto &&P0 : irange(0L, Nrhop1)){
      for(auto &&P : irange(0L, Nrhop1)){
        for(auto &&a : irange(0L, _Nvirt_ij)){
          S2.cptr()[a*S2.shape(1)+P] += U2.cptr()[a*U2.shape(1)+P0] * Jint_ccav_ccav_1_.cptr()[P*Jint_ccav_ccav_1_.shape(1)+P0];
        }//a
      }//P
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
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(i,j,a,P) <<= +1 @TFcc(a,P0,j,i) @Jint_ccav_ccav_2_(P,P0)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,P) <<= +1 @TFcc(a,P0,j,i) @Jint_ccav_ccav_2_(P,P0)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = TFcc.GetTensor(j, i);
    for(auto &&P0 : irange(0L, Nrhop1)){
      for(auto &&P : irange(0L, Nrhop1)){
        for(auto &&a : irange(0L, _Nvirt_ij)){
          S2.cptr()[a*S2.shape(1)+P] += U2.cptr()[a*U2.shape(1)+P0] * Jint_ccav_ccav_2_.cptr()[P*Jint_ccav_ccav_2_.shape(1)+P0];
        }//a
      }//P
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
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(i,j,a,P) <<= +1 @FccT(a,P0,i,j) @Jint_ccav_ccav_3_(P,P0)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,P) <<= +1 @FccT(a,P0,i,j) @Jint_ccav_ccav_3_(P,P0)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = FccT.GetTensor(i, j);
    for(auto &&P0 : irange(0L, Nrhop1)){
      for(auto &&P : irange(0L, Nrhop1)){
        for(auto &&a : irange(0L, _Nvirt_ij)){
          S2.cptr()[a*S2.shape(1)+P] += U2.cptr()[a*U2.shape(1)+P0] * Jint_ccav_ccav_3_.cptr()[P*Jint_ccav_ccav_3_.shape(1)+P0];
        }//a
      }//P
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
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(i,j,a,P) <<= +1 @Jint_ccav_ccav_4_(P,P0) @T(+1)(a,P0,i,j)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,P) <<= +1 @Jint_ccav_ccav_4_(P,P0) @T(+1)(a,P0,i,j)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = T2.GetTensor(i, j);
    for(auto &&P : irange(0L, Nrhop1)){
      for(auto &&P0 : irange(0L, Nrhop1)){
        for(auto &&a : irange(0L, _Nvirt_ij)){
          S2.cptr()[a*S2.shape(1)+P] += Jint_ccav_ccav_4_.cptr()[P*Jint_ccav_ccav_4_.shape(1)+P0] * U2.cptr()[a*U2.shape(1)+P0];
        }//a
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
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(i,j,a,P) <<= +1 @Jint_ccav_ccav_5_(P,P0) @T(+1)(a,P0,j,i)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,P) <<= +1 @Jint_ccav_ccav_5_(P,P0) @T(+1)(a,P0,j,i)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = T2.GetTensor(j, i);
    for(auto &&P : irange(0L, Nrhop1)){
      for(auto &&P0 : irange(0L, Nrhop1)){
        for(auto &&a : irange(0L, _Nvirt_ij)){
          S2.cptr()[a*S2.shape(1)+P] += Jint_ccav_ccav_5_.cptr()[P*Jint_ccav_ccav_5_.shape(1)+P0] * U2.cptr()[a*U2.shape(1)+P0];
        }//a
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
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(i,j,a,P) <<= +1 @FvvT(a,P0,i,j) @Jint_ccav_ccav_6_(P,P0)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,P) <<= +1 @FvvT(a,P0,i,j) @Jint_ccav_ccav_6_(P,P0)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = FvvT.GetTensor(i, j);
    for(auto &&P0 : irange(0L, Nrhop1)){
      for(auto &&P : irange(0L, Nrhop1)){
        for(auto &&a : irange(0L, _Nvirt_ij)){
          S2.cptr()[a*S2.shape(1)+P] += U2.cptr()[a*U2.shape(1)+P0] * Jint_ccav_ccav_6_.cptr()[P*Jint_ccav_ccav_6_.shape(1)+P0];
        }//a
      }//P
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
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(i,j,a,P) <<= +1 @FvvT(a,P0,j,i) @Jint_ccav_ccav_7_(P,P0)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,P) <<= +1 @FvvT(a,P0,j,i) @Jint_ccav_ccav_7_(P,P0)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = FvvT.GetTensor(j, i);
    for(auto &&P0 : irange(0L, Nrhop1)){
      for(auto &&P : irange(0L, Nrhop1)){
        for(auto &&a : irange(0L, _Nvirt_ij)){
          S2.cptr()[a*S2.shape(1)+P] += U2.cptr()[a*U2.shape(1)+P0] * Jint_ccav_ccav_7_.cptr()[P*Jint_ccav_ccav_7_.shape(1)+P0];
        }//a
      }//P
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
      return;
    }//End func
  } }
