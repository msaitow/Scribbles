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
void FMatrixOrthogonalizer::FormLHS_Vm1_Vp2(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &d3, DensityPack &dPack,
                                            TensorContainer &ovl, const double &Ecas, const double &h0_energy, const orz::DTensor &casfock, const orz::DTensor &FV,
                                            TensorContainer &T2, TensorContainer &R2){
#include "lct_preparations.cpp"
orz::DTensor &Jint_cavv_ccaa_0_ = (*dPack["Jint_cavv_ccaa_0_"]);
orz::DTensor &Jint_cavv_ccaa_1_ = (*dPack["Jint_cavv_ccaa_1_"]);
std::cout << "@S2(i,j,a,P) <<= +1 @Jint_cavv_ccaa_0_(P,P0,a2) @cF(a2,a) @T(+2)(j,i,P0)" << std::endl;
// @S2(i,j,a,P) <<= +1 @Jint_cavv_ccaa_0_(P,P0,a2) @cF(a2,a) @T(+2)(j,i,P0)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&j : irange(0L, _Ndocc)){
    const long I_ij = RActList.GetIndex(i, j);
    if (I_ij < 0) continue;
    orz::DTensor S2 = R2.GetTensor(i, j);
    const long _Nvirt_ij = RActList.GetNPNO(i, j);
    if (I_ij < 0) continue;
    for(auto &&P0 : irange(0L, Nrhop2)){
      c