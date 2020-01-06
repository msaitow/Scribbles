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
void FMatrixOrthogonalizer::FormLHS_V0p_Vp1p(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &d3, DensityPack &dPack,
                                            TensorContainer &ovl, const double &Ecas, const double &h0_energy, const orz::DTensor &casfock, const orz::DTensor &FV,
                                            TensorContainer &T2, TensorContainer &R2){
#include "lct_preparations.cpp"
orz::DTensor &Jint_caav_caaa_0_ = (*dPack["Jint_caav_caaa_0_"]);

TensorContainer JF(Nrho0p);
for(auto &&P : irange(0L, Nrho0p)){
  const long I_P = RActList.GetIndex(P);
  if (I_P < 0) continue;  
  orz::DTensor J_P(Nrhop1p, _Nact);
  for(auto &&P0 : irange(0L, Nrhop1p)){
    for(auto &&a2 : irange(0L, _Nact)){    
      J_P.cptr()[P0*J_P.shape(1)+a2] = Jint_caav_caaa_0_.cptr()[P*Jint_caav_caaa_0_.shape(1)*Jint_caav_caaa_0_.shape(2)+P0*Jint_caav_caaa_0_.shape(2)+a2];
    }//a2
  }//P0
  orz::DTensor JF_P(J_P*Fpa);
  JF.PutTensor(P, JF_P);
}//P
 
{
double t_trans=0.0;
#ifdef _PROFILE_LHS
std::cout << boost::format("%80s") % "@S2(P,i,a) <<= +1 @Jint_caav_caaa_0_(P,P0,a2) @cF(a2,a) @T(+1')(i,P0)";
#endif
high_resolution_clock::time_point t1 = high_resolution_clock::now();
// @S2(P,i,a) <<= +1 @Jint_caav_caaa_0_(P,P0,a2) @cF(a2,a) @T(+1')(i,P0)
for(auto &&P : irange(0L, Nrho0p)){
  const long I_P = RActList.GetIndex(P);
  if (I_P < 0) continue;
  orz::DTensor JF_P = JF.GetTensor(P);
  orz::DTensor S2 = R2.GetTensor(P);
  for(auto &&P0 : irange(0L, Nrhop1p)){
    const long I_P0 = TActList.GetIndex(P0);
    if (I_P0 < 0) continue;
    orz::DTensor U2 = T2.GetTensor(P0);
    for(auto &&i : irange(0L, _Ndocc)){
      for(auto &&a : irange(0L, _Nvirt)){
        S2.cptr()[i*S2.shape(1)+a] += JF_P.cptr()[P0*JF_P.shape(1)+a] * U2.cptr()[i];
      }//a
   }//i
  }//P0
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
