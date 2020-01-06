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
void FMatrixOrthogonalizer::FormLHS_V0_V0(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &d3, DensityPack &dPack,
                                            TensorContainer &ovl, const double &Ecas, const double &h0_energy, const orz::DTensor &casfock, const orz::DTensor &FV,
                                            TensorContainer &T2, TensorContainer &R2){
#include "lct_preparations.cpp"
R2.initContainer(_Ndocc, _Ndocc);             
// Init amplitude container                   
for(auto &&p : irange(0L, _Ndocc)){           
  for(auto &&q : irange(0L, _Ndocc)){         
    const long I_pq = RActList.GetIndex(p, q);
    if (I_pq < 0) continue;                   
    PairData &P_pq = RActList.GetPair(p,q);   
    const int Nvirt = P_pq.GetNPNO();         
    orz::DTensor Rpq(Nvirt, Nvirt);           
    R2.PutTensor(p, q, Rpq);                  
  }                                           
}                                             
TensorContainer TFcc(_Ndocc,_Ndocc);
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = TActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    const long _Nvirt_ij = TActList.GetNPNO(i, j);
    orz::DTensor X(_Nvirt_ij,_Nvirt_ij);
    for(auto &&c0 : irange(0L, _Ndocc)){
      const long I_ic0 = RActList.GetIndex(i, c0);
      if (I_ic0 < 0) continue;
      orz::DTensor U2 = T2.GetTensor(i, c0);
      {
        high_resolution_clock::time_point o1 = high_resolution_clock::now();
        orz::DTensor S = ovl.CopyTensor(std::max(I_ij,I_ic0),std::min(I_ij,I_ic0));
        if(I_ic0>I_ij) orz::tensor::transpose(S);
        orz::DTensor tmp;
        A_x_B_x_At(tmp, S, U2);
        U2 = tmp.copy();
        high_resolution_clock::time_point o2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
        //t_trans += time_span.count();
      }
      X.gaxpy(+1.0, U2, Fij.cptr()[j*Fij.shape(1)+c0]);
    }//c0
    TFcc.PutTensor(i,j,X);
  }//j
}//i
TensorContainer FccT(_Ndocc,_Ndocc);
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = TActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    const long _Nvirt_ij = TActList.GetNPNO(i, j);
    orz::DTensor X(_Nvirt_ij,_Nvirt_ij);
    for(auto &&c0 : irange(0L, _Ndocc)){
      const long I_c0j = RActList.GetIndex(c0, j);
      if (I_c0j < 0) continue;
      orz::DTensor U2 = T2.GetTensor(c0, j);
      {
        high_resolution_clock::time_point o1 = high_resolution_clock::now();
        orz::DTensor S = ovl.CopyTensor(std::max(I_ij,I_c0j),std::min(I_ij,I_c0j));
        if(I_c0j>I_ij) orz::tensor::transpose(S);
        orz::DTensor tmp;
        A_x_B_x_At(tmp, S, U2);
        U2 = tmp.copy();
        high_resolution_clock::time_point o2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(o2 - o1);
        //t_trans += time_span.count();
      }
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
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    orz::DTensor U2 = T2.GetTensor(i, j);
    orz::DTensor Fab(_Fab_.copy());
    RActList.GetPair(i, j).Transform2EXT(Fab);
    orz::DTensor X(_Nvirt_ij,_Nvirt_ij);
//    for(auto &&v0 : irange(0L, _Nvirt_ij)){
//      for(auto &&b : irange(0L, _Nvirt_ij)){
//        for(auto &&a : irange(0L, _Nvirt_ij)){
//          X.cptr()[a*X.shape(1)+b] += Fab.cptr()[v0*Fab.shape(1)+a] * U2.cptr()[v0*U2.shape(1)+b];
//        }//a
//      }//b
//    }//v0
    A_x_B(X, Fab, U2, false, false);
    FvvT.PutTensor(i, j, X);
  }//j
}//i
orz::DTensor &Jint_ccvv_ccvv_0_ = (*dPack["Jint_ccvv_ccvv_0_"]);
orz::DTensor &Jint_ccvv_ccvv_1_ = (*dPack["Jint_ccvv_ccvv_1_"]);
orz::DTensor &Jint_ccvv_ccvv_2_ = (*dPack["Jint_ccvv_ccvv_2_"]);
orz::DTensor &Jint_ccvv_ccvv_3_ = (*dPack["Jint_ccvv_ccvv_3_"]);
orz::DTensor &Jint_ccvv_ccvv_4_ = (*dPack["Jint_ccvv_ccvv_4_"]);
orz::DTensor &Jint_ccvv_ccvv_5_ = (*dPack["Jint_ccvv_ccvv_5_"]);
orz::DTensor &Jint_ccvv_ccvv_6_ = (*dPack["Jint_ccvv_ccvv_6_"]);
orz::DTensor &Jint_ccvv_ccvv_7_ = (*dPack["Jint_ccvv_ccvv_7_"]);
orz::DTensor &Jint_ccvv_ccvv_8_ = (*dPack["Jint_ccvv_ccvv_8_"]);
orz::DTensor &Jint_ccvv_ccvv_9_ = (*dPack["Jint_ccvv_ccvv_9_"]);
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(i,j,a,b) <<= +1 @Jint_ccvv_ccvv_0_() @T(0)(b,a,j,i)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,b) <<= +1 @Jint_ccvv_ccvv_0_() @T(0)(b,a,j,i)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = T2.GetTensor(j, i);
    for(auto &&b : irange(0L, _Nvirt_ij)){
      for(auto &&a : irange(0L, _Nvirt_ij)){
        S2.cptr()[a*S2.shape(1)+b] += Jint_ccvv_ccvv_0_.cptr()[0] * U2.cptr()[b*U2.shape(1)+a];
      }//a
    }//b
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
std::cout << boost::format("%80s") % "@S2(i,j,a,b) <<= +1 @TFcc(b,a,j,i) @Jint_ccvv_ccvv_1_()";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,b) <<= +1 @TFcc(b,a,j,i) @Jint_ccvv_ccvv_1_()
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = TFcc.GetTensor(j, i);
    for(auto &&b : irange(0L, _Nvirt_ij)){
      for(auto &&a : irange(0L, _Nvirt_ij)){
        S2.cptr()[a*S2.shape(1)+b] += U2.cptr()[b*U2.shape(1)+a] * Jint_ccvv_ccvv_1_.cptr()[0];
      }//a
    }//b
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
std::cout << boost::format("%80s") % "@S2(i,j,a,b) <<= +1 @FccT(b,a,j,i) @Jint_ccvv_ccvv_2_()";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,b) <<= +1 @FccT(b,a,j,i) @Jint_ccvv_ccvv_2_()
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = FccT.GetTensor(j, i);
    for(auto &&b : irange(0L, _Nvirt_ij)){
      for(auto &&a : irange(0L, _Nvirt_ij)){
        S2.cptr()[a*S2.shape(1)+b] += U2.cptr()[b*U2.shape(1)+a] * Jint_ccvv_ccvv_2_.cptr()[0];
      }//a
    }//b
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
std::cout << boost::format("%80s") % "@S2(i,j,a,b) <<= +1 @TFcc(b,a,i,j) @Jint_ccvv_ccvv_3_()";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,b) <<= +1 @TFcc(b,a,i,j) @Jint_ccvv_ccvv_3_()
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = TFcc.GetTensor(i, j);
    for(auto &&b : irange(0L, _Nvirt_ij)){
      for(auto &&a : irange(0L, _Nvirt_ij)){
        S2.cptr()[a*S2.shape(1)+b] += U2.cptr()[b*U2.shape(1)+a] * Jint_ccvv_ccvv_3_.cptr()[0];
      }//a
    }//b
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
std::cout << boost::format("%80s") % "@S2(i,j,a,b) <<= +1 @FccT(b,a,i,j) @Jint_ccvv_ccvv_4_()";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,b) <<= +1 @FccT(b,a,i,j) @Jint_ccvv_ccvv_4_()
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = FccT.GetTensor(i, j);
    for(auto &&b : irange(0L, _Nvirt_ij)){
      for(auto &&a : irange(0L, _Nvirt_ij)){
        S2.cptr()[a*S2.shape(1)+b] += U2.cptr()[b*U2.shape(1)+a] * Jint_ccvv_ccvv_4_.cptr()[0];
      }//a
    }//b
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
std::cout << boost::format("%80s") % "@S2(i,j,a,b) <<= +1 @Jint_ccvv_ccvv_5_() @T(0)(b,a,i,j)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,b) <<= +1 @Jint_ccvv_ccvv_5_() @T(0)(b,a,i,j)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = T2.GetTensor(i, j);
    for(auto &&b : irange(0L, _Nvirt_ij)){
      for(auto &&a : irange(0L, _Nvirt_ij)){
        S2.cptr()[a*S2.shape(1)+b] += Jint_ccvv_ccvv_5_.cptr()[0] * U2.cptr()[b*U2.shape(1)+a];
      }//a
    }//b
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
std::cout << boost::format("%80s") % "@S2(i,j,a,b) <<= +1 @FvvT(b,a,j,i) @Jint_ccvv_ccvv_6_()";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,b) <<= +1 @FvvT(b,a,j,i) @Jint_ccvv_ccvv_6_()
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = FvvT.GetTensor(j, i);
    for(auto &&b : irange(0L, _Nvirt_ij)){
      for(auto &&a : irange(0L, _Nvirt_ij)){
        S2.cptr()[a*S2.shape(1)+b] += U2.cptr()[b*U2.shape(1)+a] * Jint_ccvv_ccvv_6_.cptr()[0];
      }//a
    }//b
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
std::cout << boost::format("%80s") % "@S2(i,j,a,b) <<= +1 @FvvT(b,a,i,j) @Jint_ccvv_ccvv_7_()";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,b) <<= +1 @FvvT(b,a,i,j) @Jint_ccvv_ccvv_7_()
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = FvvT.GetTensor(i, j);
    for(auto &&b : irange(0L, _Nvirt_ij)){
      for(auto &&a : irange(0L, _Nvirt_ij)){
        S2.cptr()[a*S2.shape(1)+b] += U2.cptr()[b*U2.shape(1)+a] * Jint_ccvv_ccvv_7_.cptr()[0];
      }//a
    }//b
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
std::cout << boost::format("%80s") % "@S2(i,j,a,b) <<= +1 @FvvT(a,b,i,j) @Jint_ccvv_ccvv_8_()";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,b) <<= +1 @FvvT(a,b,i,j) @Jint_ccvv_ccvv_8_()
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = FvvT.GetTensor(i, j);
    for(auto &&a : irange(0L, _Nvirt_ij)){
      for(auto &&b : irange(0L, _Nvirt_ij)){
        S2.cptr()[a*S2.shape(1)+b] += U2.cptr()[a*U2.shape(1)+b] * Jint_ccvv_ccvv_8_.cptr()[0];
      }//b
    }//a
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
std::cout << boost::format("%80s") % "@S2(i,j,a,b) <<= +1 @FvvT(a,b,j,i) @Jint_ccvv_ccvv_9_()";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(i,j,a,b) <<= +1 @FvvT(a,b,j,i) @Jint_ccvv_ccvv_9_()
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = FvvT.GetTensor(j, i);
    for(auto &&a : irange(0L, _Nvirt_ij)){
      for(auto &&b : irange(0L, _Nvirt_ij)){
        S2.cptr()[a*S2.shape(1)+b] += U2.cptr()[a*U2.shape(1)+b] * Jint_ccvv_ccvv_9_.cptr()[0];
      }//b
    }//a
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
