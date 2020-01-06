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
namespace orz { namespace lct {
void FMatrixOrthogonalizer::FormLHS_Vm1_Vp1(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &d3, DensityPack &dPack,
                                            TensorContainer &ovl, const double &Ecas, const double &h0_energy, const orz::DTensor &casfock, const orz::DTensor &FV,
                                            TensorContainer &T2, TensorContainer &R2){
#include "lct_preparations.cpp"
orz::DTensor &Jint_cavv_ccav_0_ = (*dPack["Jint_cavv_ccav_0_"]);
orz::DTensor &Jint_cavv_ccav_1_ = (*dPack["Jint_cavv_ccav_1_"]);
orz::DTensor &Jint_cavv_ccav_2_ = (*dPack["Jint_cavv_ccav_2_"]);
orz::DTensor &Jint_cavv_ccav_3_ = (*dPack["Jint_cavv_ccav_3_"]);
orz::DTensor &Jint_cavv_ccav_4_ = (*dPack["Jint_cavv_ccav_4_"]);
orz::DTensor &Jint_cavv_ccav_5_ = (*dPack["Jint_cavv_ccav_5_"]);
orz::DTensor &Jint_cavv_ccav_6_ = (*dPack["Jint_cavv_ccav_6_"]);
orz::DTensor &Jint_cavv_ccav_7_ = (*dPack["Jint_cavv_ccav_7_"]);
std::cout << "@S2(i,j,a,P) <<= +1 @Jint_cavv_ccav_0_(P,P0) @cF(j,c0) @T(+1)(a,P0,i,c0)" << std::endl;
// @S2(i,j,a,P) <<= +1 @Jint_cavv_ccav_0_(P,P0) @cF(j,c0) @T(+1)(a,P0,i,c0)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    for(auto &&c0 : irange(0L, _Ndocc)){
      const long I_c0i = TActList.GetIndex(c0, i);
      if (I_c0i < 0) continue;
      orz::DTensor U2 = T2.GetTensor(i, c0);
      {
        orz::DTensor tmp(U2.copy());
        RActList.GetPair(i, j).Transform1EXT_L(tmp);
        U2 = tmp.copy();
      }
      for(auto &&P : irange(0L, Nrhop1)){
        for(auto &&P0 : irange(0L, Nrhop1)){
          for(auto &&a : irange(0L, _Nvirt_ij)){
            S2.cptr()[a*S2.shape(1)+P] += Jint_cavv_ccav_0_.cptr()[P*Jint_cavv_ccav_0_.shape(1)+P0] * Fij.cptr()[j*Fij.shape(1)+c0] * U2.cptr()[a*U2.shape(1)+P0];
          }//a
        }//P0
      }//P
    }//c0
    R2.PutTensor(i, j, S2);
  }//j
}//i
std::cout << "@S2(i,j,a,P) <<= +1 @Jint_cavv_ccav_1_(P,P0) @cF(j,c0) @T(+1)(a,P0,c0,i)" << std::endl;
// @S2(i,j,a,P) <<= +1 @Jint_cavv_ccav_1_(P,P0) @cF(j,c0) @T(+1)(a,P0,c0,i)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    for(auto &&c0 : irange(0L, _Ndocc)){
      const long I_ic0 = TActList.GetIndex(i, c0);
      if (I_ic0 < 0) continue;
      orz::DTensor U2 = T2.GetTensor(c0, i);
      {
        orz::DTensor tmp(U2.copy());
        RActList.GetPair(i, j).Transform1EXT_L(tmp);
        U2 = tmp.copy();
      }
      for(auto &&P : irange(0L, Nrhop1)){
        for(auto &&P0 : irange(0L, Nrhop1)){
          for(auto &&a : irange(0L, _Nvirt_ij)){
            S2.cptr()[a*S2.shape(1)+P] += Jint_cavv_ccav_1_.cptr()[P*Jint_cavv_ccav_1_.shape(1)+P0] * Fij.cptr()[j*Fij.shape(1)+c0] * U2.cptr()[a*U2.shape(1)+P0];
          }//a
        }//P0
      }//P
    }//c0
    R2.PutTensor(i, j, S2);
  }//j
}//i
std::cout << "@S2(i,j,a,P) <<= +1 @Jint_cavv_ccav_2_(P,P0) @cF(i,c0) @T(+1)(a,P0,j,c0)" << std::endl;
// @S2(i,j,a,P) <<= +1 @Jint_cavv_ccav_2_(P,P0) @cF(i,c0) @T(+1)(a,P0,j,c0)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    for(auto &&c0 : irange(0L, _Ndocc)){
      const long I_c0j = TActList.GetIndex(c0, j);
      if (I_c0j < 0) continue;
      orz::DTensor U2 = T2.GetTensor(j, c0);
      {
        orz::DTensor tmp(U2.copy());
        RActList.GetPair(i, j).Transform1EXT_L(tmp);
        U2 = tmp.copy();
      }
      for(auto &&P : irange(0L, Nrhop1)){
        for(auto &&P0 : irange(0L, Nrhop1)){
          for(auto &&a : irange(0L, _Nvirt_ij)){
            S2.cptr()[a*S2.shape(1)+P] += Jint_cavv_ccav_2_.cptr()[P*Jint_cavv_ccav_2_.shape(1)+P0] * Fij.cptr()[i*Fij.shape(1)+c0] * U2.cptr()[a*U2.shape(1)+P0];
          }//a
        }//P0
      }//P
    }//c0
    R2.PutTensor(i, j, S2);
  }//j
}//i
std::cout << "@S2(i,j,a,P) <<= +1 @Jint_cavv_ccav_3_(P,P0) @cF(i,c0) @T(+1)(a,P0,c0,j)" << std::endl;
// @S2(i,j,a,P) <<= +1 @Jint_cavv_ccav_3_(P,P0) @cF(i,c0) @T(+1)(a,P0,c0,j)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    for(auto &&c0 : irange(0L, _Ndocc)){
      const long I_jc0 = TActList.GetIndex(j, c0);
      if (I_jc0 < 0) continue;
      orz::DTensor U2 = T2.GetTensor(c0, j);
      {
        orz::DTensor tmp(U2.copy());
        RActList.GetPair(i, j).Transform1EXT_L(tmp);
        U2 = tmp.copy();
      }
      for(auto &&P : irange(0L, Nrhop1)){
        for(auto &&P0 : irange(0L, Nrhop1)){
          for(auto &&a : irange(0L, _Nvirt_ij)){
            S2.cptr()[a*S2.shape(1)+P] += Jint_cavv_ccav_3_.cptr()[P*Jint_cavv_ccav_3_.shape(1)+P0] * Fij.cptr()[i*Fij.shape(1)+c0] * U2.cptr()[a*U2.shape(1)+P0];
          }//a
        }//P0
      }//P
    }//c0
    R2.PutTensor(i, j, S2);
  }//j
}//i
std::cout << "@S2(i,j,a,P) <<= +1 @Jint_cavv_ccav_4_(P,P0) @T(+1)(a,P0,i,j)" << std::endl;
// @S2(i,j,a,P) <<= +1 @Jint_cavv_ccav_4_(P,P0) @T(+1)(a,P0,i,j)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = T2.GetTensor(i, j);
    {
      orz::DTensor tmp(U2.copy());
      RActList.GetPair(i, j).Transform1EXT_L(tmp);
      U2 = tmp.copy();
    }
    for(auto &&P : irange(0L, Nrhop1)){
      for(auto &&P0 : irange(0L, Nrhop1)){
        for(auto &&a : irange(0L, _Nvirt_ij)){
          S2.cptr()[a*S2.shape(1)+P] += Jint_cavv_ccav_4_.cptr()[P*Jint_cavv_ccav_4_.shape(1)+P0] * U2.cptr()[a*U2.shape(1)+P0];
        }//a
      }//P0
    }//P
    R2.PutTensor(i, j, S2);
  }//j
}//i
std::cout << "@S2(i,j,a,P) <<= +1 @Jint_cavv_ccav_5_(P,P0) @T(+1)(a,P0,j,i)" << std::endl;
// @S2(i,j,a,P) <<= +1 @Jint_cavv_ccav_5_(P,P0) @T(+1)(a,P0,j,i)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = T2.GetTensor(j, i);
    {
      orz::DTensor tmp(U2.copy());
      RActList.GetPair(i, j).Transform1EXT_L(tmp);
      U2 = tmp.copy();
    }
    for(auto &&P : irange(0L, Nrhop1)){
      for(auto &&P0 : irange(0L, Nrhop1)){
        for(auto &&a : irange(0L, _Nvirt_ij)){
          S2.cptr()[a*S2.shape(1)+P] += Jint_cavv_ccav_5_.cptr()[P*Jint_cavv_ccav_5_.shape(1)+P0] * U2.cptr()[a*U2.shape(1)+P0];
        }//a
      }//P0
    }//P
    R2.PutTensor(i, j, S2);
  }//j
}//i
std::cout << "@S2(i,j,a,P) <<= +1 @Jint_cavv_ccav_6_(P,P0) @cF(v0,a) @T(+1)(v0,P0,i,j)" << std::endl;
// @S2(i,j,a,P) <<= +1 @Jint_cavv_ccav_6_(P,P0) @cF(v0,a) @T(+1)(v0,P0,i,j)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = T2.GetTensor(i, j);
    {
      orz::DTensor tmp(U2.copy());
      RActList.GetPair(i, j).Transform1EXT_L(tmp);
      U2 = tmp.copy();
    }
    orz::DTensor Fab(_Fab_.copy());
    {
      RActList.GetPair(i, j).Transform2EXT(Fab);
    }
    for(auto &&P : irange(0L, Nrhop1)){
      for(auto &&P0 : irange(0L, Nrhop1)){
        for(auto &&v0 : irange(0L, _Nvirt_ij)){
          for(auto &&a : irange(0L, _Nvirt_ij)){
            S2.cptr()[a*S2.shape(1)+P] += Jint_cavv_ccav_6_.cptr()[P*Jint_cavv_ccav_6_.shape(1)+P0] * Fab.cptr()[v0*Fab.shape(1)+a] * U2.cptr()[v0*U2.shape(1)+P0];
          }//a
        }//v0
      }//P0
    }//P
    R2.PutTensor(i, j, S2);
  }//j
}//i
std::cout << "@S2(i,j,a,P) <<= +1 @Jint_cavv_ccav_7_(P,P0) @cF(v0,a) @T(+1)(v0,P0,j,i)" << std::endl;
// @S2(i,j,a,P) <<= +1 @Jint_cavv_ccav_7_(P,P0) @cF(v0,a) @T(+1)(v0,P0,j,i)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    orz::DTensor U2 = T2.GetTensor(j, i);
    {
      orz::DTensor tmp(U2.copy());
      RActList.GetPair(i, j).Transform1EXT_L(tmp);
      U2 = tmp.copy();
    }
    orz::DTensor Fab(_Fab_.copy());
    {
      RActList.GetPair(i, j).Transform2EXT(Fab);
    }
    for(auto &&P : irange(0L, Nrhop1)){
      for(auto &&P0 : irange(0L, Nrhop1)){
        for(auto &&v0 : irange(0L, _Nvirt_ij)){
          for(auto &&a : irange(0L, _Nvirt_ij)){
            S2.cptr()[a*S2.shape(1)+P] += Jint_cavv_ccav_7_.cptr()[P*Jint_cavv_ccav_7_.shape(1)+P0] * Fab.cptr()[v0*Fab.shape(1)+a] * U2.cptr()[v0*U2.shape(1)+P0];
          }//a
        }//v0
      }//P0
    }//P
    R2.PutTensor(i, j, S2);
  }//j
}//i
      return;
    }//End func
  } }
