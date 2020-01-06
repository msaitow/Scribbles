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
    void FMatrixOrthogonalizer::FormLHS_Vm2_V0_CEd(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &c2, const orz::DTensor &d3, DensityPack &dPack,
                                                   TensorContainer &ovl, TensorContainer &PNO4, const bool do4EXTinPNO,
                                                   const double QE0, const double Ecas, const double h0_energy, TensorContainer &DebugInts, const orz::DTensor &casfock, const orz::DTensor &FV,
                                                   TensorContainer &T2, TensorContainer &R2){

#include "lct_preparations_orthdens.cpp"      

      // Back-transform PNO ampitude into canonical basis
      //TensorContainer Ttilde(_Ndocc, _Ndocc);
      orz::DTensor Tc(_Ndocc, _Ndocc, _Nvirt, _Nvirt);
      orz::DTensor T (_Ndocc, _Ndocc, _Nvirt, _Nvirt);      
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&j : irange(0L, _Ndocc)){
          const long I_ij = TActList.GetIndex(i, j);
          if (I_ij < 0) continue;
          orz::DTensor U2 = T2.GetTensor(i, j);
          orz::DTensor U2tilde(U2.copy());
          U2tilde.gaxpy(+2.0, U2.swapdim(1,0), -1.0);
          //Ttilde.PutTensor(i, j, U2tilde);
          // Back-transform the PNO amplitude into canonical basis
          TActList.GetPair(i,j).BackTransform2EXT(U2     );
          TActList.GetPair(i,j).BackTransform2EXT(U2tilde);          
          for(auto &&a : irange(0L, _Nvirt))
            for(auto &&b : irange(0L, _Nvirt)){
              T (i,j,a,b) = U2     (a,b);
              Tc(i,j,a,b) = U2tilde(a,b);              
            }
        }
      }      
        
#include "lct_integrals4EXT.cpp"

 //Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn
  orz::DTensor tRFiC(_Nvirt,_Nvirt,Nrhom2);
  {
    const double _X_ =   1.0000000000000000;
    orz::DTensor W0(_Ndocc,_Ndocc,_Nact,_Nact);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp1.cptr()[(a2*_Nact+a0)*tmp1.shape(1)+(a3*_Nact+a1)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a1
          }//a3
        }//a0
      }//a2
      orz::DTensor tmp2(_Nact*_Nact,_Ndocc*_Ndocc);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&c1 : irange(0L, _Ndocc)){
              tmp2.cptr()[(a3*_Nact+a1)*tmp2.shape(1)+(c0*_Ndocc+c1)] = V2_FiC_acac.cptr()[a1*V2_FiC_acac.shape(1)*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+c0*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+a3*V2_FiC_acac.shape(3)+c1];
            }//c1
          }//c0
        }//a1
      }//a3
      // @W0(c0,c1,a2,a0) <<= +1 @D2(a3,a2,a1,a0) @V2_FiC_acac(a1,c0,a3,c1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&c1 : irange(0L, _Ndocc)){
              W0.cptr()[c0*_Ndocc*_Nact*_Nact+c1*_Nact*_Nact+a2*_Nact+a0] += tmp3.cptr()[a2*_Nact*_Ndocc*_Ndocc+a0*_Ndocc*_Ndocc+c0*_Ndocc+c1];
            }//c1
          }//c0
        }//a0
      }//a2
    }
    orz::DTensor W1(_Ndocc,_Ndocc,Nrhom2);
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc,_Nact*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&c1 : irange(0L, _Ndocc)){
          for(auto &&a2 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(c0*_Ndocc+c1)*tmp1.shape(1)+(a2*_Nact+a0)] = W0.cptr()[c0*W0.shape(1)*W0.shape(2)*W0.shape(3)+c1*W0.shape(2)*W0.shape(3)+a2*W0.shape(3)+a0];
            }//a0
          }//a2
        }//c1
      }//c0
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(a2*_Nact+a0)*tmp2.shape(1)+(P)] = Cm2.cptr()[a0*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+P];
          }//P
        }//a0
      }//a2
      // @W1(c0,c1,P) <<= +1 @W0(c0,c1,a2,a0) @X(-2)(P,a0,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&c1 : irange(0L, _Ndocc)){
          for(auto &&P : irange(0L, Nrhom2)){
            W1.cptr()[c0*_Ndocc*Nrhom2+c1*Nrhom2+P] += tmp3.cptr()[c0*_Ndocc*Nrhom2+c1*Nrhom2+P];
          }//P
        }//c1
      }//c0
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Ndocc*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c1 : irange(0L, _Ndocc)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*_Nvirt+a)*tmp1.shape(1)+(c1*_Ndocc+c0)] = T.cptr()[c1*T.shape(1)*T.shape(2)*T.shape(3)+c0*T.shape(2)*T.shape(3)+b*T.shape(3)+a];
            }//c0
          }//c1
        }//a
      }//b
      orz::DTensor tmp2(_Ndocc*_Ndocc,Nrhom2);
      for(auto &&c1 : irange(0L, _Ndocc)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&P : irange(0L, Nrhom2)){
            tmp2.cptr()[(c1*_Ndocc+c0)*tmp2.shape(1)+(P)] = W1.cptr()[c0*W1.shape(1)*W1.shape(2)+c1*W1.shape(2)+P];
          }//P
        }//c0
      }//c1
      // @tRFiC(a,b,P) <<= +1 _X_ @T(c1,c0,b,a) @W1(c0,c1,P)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&P : irange(0L, Nrhom2)){
            tRFiC.cptr()[a*_Nvirt*Nrhom2+b*Nrhom2+P] += tmp3.cptr()[b*_Nvirt*Nrhom2+a*Nrhom2+P];
          }//P
        }//a
      }//b
    }
  }//End Contra
  //genCode.insert(std::make_pair("tRFiC", boost::shared_ptr<orz::DTensor>( new orz::DTensor(tRFiC) )));
 //Code-End-Code-End-Code-End-Code-End-Code-End-Code-End-Code-End      
    
      // Transform tRFiC into PNO basis
      for(auto &&P : irange(0L, Nrhom2)){
        const long I_P = RActList.GetIndex(P);
        if (I_P < 0) continue;
        orz::DTensor RP(_Nvirt, _Nvirt);
        for(auto &&a : irange(0L, _Nvirt))
          for(auto &&b : irange(0L, _Nvirt))
            RP(a,b) = tRFiC(a,b,P);
        RActList.GetPair(P).Transform2EXT(RP);
        orz::DTensor Rorig = R2.GetTensor(P);
        RP += Rorig;
        R2.PutTensor(P, RP);
      }
            
      return;
    }//End func
  } }
