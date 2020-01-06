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
    void FMatrixOrthogonalizer::FormLHS_Vm1_V0_CEd(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &c2, const orz::DTensor &d3, DensityPack &dPack,
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
  orz::DTensor tRFiC(_Nvirt,_Nvirt,Nrhom1,_Ndocc);
  {
    const double _X_ =  -2.0000000000000000;
    orz::DTensor W0(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt,_Ndocc*_Nvirt);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp1.cptr()[(i*_Nvirt+a)*tmp1.shape(1)+(c0*_Nvirt+v0)] = Tc.cptr()[i*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+a*Tc.shape(3)+v0];
            }//v0
          }//c0
        }//a
      }//i
      orz::DTensor tmp2(_Ndocc*_Nvirt,_Nvirt*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp2.cptr()[(c0*_Nvirt+v0)*tmp2.shape(1)+(b*_Nact+a1)] = V2_FiC_vvac.cptr()[b*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a1*V2_FiC_vvac.shape(3)+c0];
            }//a1
          }//b
        }//v0
      }//c0
      // @W0(i,a1,a,b) <<= +1 @Tc(i,c0,a,v0) @V2_FiC_vvac(b,v0,a1,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              W0.cptr()[i*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[i*_Nvirt*_Nvirt*_Nact+a*_Nvirt*_Nact+b*_Nact+a1];
            }//a1
          }//b
        }//a
      }//i
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+a*_Nvirt+b)] = W0.cptr()[i*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+a*W0.shape(3)+b];
            }//b
          }//a
        }//i
      }//a1
      // @W1(i,a0,a,b) <<= +1 @D1(a1,a0) @W0(i,a1,a,b)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a0*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a0*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+a*_Nvirt+b];
            }//b
          }//a
        }//i
      }//a0
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+a*_Nvirt+b)*tmp1.shape(1)+(a0)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+a*W1.shape(3)+b];
            }//a0
          }//b
        }//a
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(i,a0,a,b) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -2.0000000000000000;
    orz::DTensor W0(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt,_Ndocc*_Nvirt);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp1.cptr()[(i*_Nvirt+b)*tmp1.shape(1)+(c0*_Nvirt+v0)] = Tc.cptr()[i*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+b*Tc.shape(3)+v0];
            }//v0
          }//c0
        }//b
      }//i
      orz::DTensor tmp2(_Ndocc*_Nvirt,_Nvirt*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp2.cptr()[(c0*_Nvirt+v0)*tmp2.shape(1)+(a*_Nact+a1)] = V2_FiC_vvac.cptr()[a*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a1*V2_FiC_vvac.shape(3)+c0];
            }//a1
          }//a
        }//v0
      }//c0
      // @W0(i,a1,b,a) <<= +1 @Tc(i,c0,b,v0) @V2_FiC_vvac(a,v0,a1,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              W0.cptr()[i*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[i*_Nvirt*_Nvirt*_Nact+b*_Nvirt*_Nact+a*_Nact+a1];
            }//a1
          }//a
        }//b
      }//i
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+b*_Nvirt+a)] = W0.cptr()[i*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+b*W0.shape(3)+a];
            }//a
          }//b
        }//i
      }//a1
      // @W1(i,a0,b,a) <<= +1 @D1(a1,a0) @W0(i,a1,b,a)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a0*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[a0*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+b*_Nvirt+a];
            }//a
          }//b
        }//i
      }//a0
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(a0)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+b*W1.shape(3)+a];
            }//a0
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a0];
        }//P
      }//a0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(i,a0,b,a) @X(-1)(+)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//b
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -2.0000000000000000;
    orz::DTensor W0(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt,_Ndocc*_Nvirt);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp1.cptr()[(i*_Nvirt+b)*tmp1.shape(1)+(c0*_Nvirt+v0)] = Tc.cptr()[i*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+v0*Tc.shape(3)+b];
            }//v0
          }//c0
        }//b
      }//i
      orz::DTensor tmp2(_Ndocc*_Nvirt,_Nvirt*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp2.cptr()[(c0*_Nvirt+v0)*tmp2.shape(1)+(a*_Nact+a1)] = V2_FiC_vvac.cptr()[a*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a1*V2_FiC_vvac.shape(3)+c0];
            }//a1
          }//a
        }//v0
      }//c0
      // @W0(i,a1,b,a) <<= +1 @Tc(i,c0,v0,b) @V2_FiC_vvac(a,v0,a1,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              W0.cptr()[i*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[i*_Nvirt*_Nvirt*_Nact+b*_Nvirt*_Nact+a*_Nact+a1];
            }//a1
          }//a
        }//b
      }//i
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+b*_Nvirt+a)] = W0.cptr()[i*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+b*W0.shape(3)+a];
            }//a
          }//b
        }//i
      }//a1
      // @W1(i,a0,b,a) <<= +1 @D1(a1,a0) @W0(i,a1,b,a)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a0*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[a0*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+b*_Nvirt+a];
            }//a
          }//b
        }//i
      }//a0
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(a0)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+b*W1.shape(3)+a];
            }//a0
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(i,a0,b,a) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//b
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -2.0000000000000000;
    orz::DTensor W0(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt,_Ndocc*_Nvirt);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp1.cptr()[(i*_Nvirt+a)*tmp1.shape(1)+(c0*_Nvirt+v0)] = Tc.cptr()[i*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+v0*Tc.shape(3)+a];
            }//v0
          }//c0
        }//a
      }//i
      orz::DTensor tmp2(_Ndocc*_Nvirt,_Nvirt*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp2.cptr()[(c0*_Nvirt+v0)*tmp2.shape(1)+(b*_Nact+a1)] = V2_FiC_vvac.cptr()[b*V2_FiC_vvac.shape(1)*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+v0*V2_FiC_vvac.shape(2)*V2_FiC_vvac.shape(3)+a1*V2_FiC_vvac.shape(3)+c0];
            }//a1
          }//b
        }//v0
      }//c0
      // @W0(i,a1,a,b) <<= +1 @Tc(i,c0,v0,a) @V2_FiC_vvac(b,v0,a1,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              W0.cptr()[i*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[i*_Nvirt*_Nvirt*_Nact+a*_Nvirt*_Nact+b*_Nact+a1];
            }//a1
          }//b
        }//a
      }//i
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+a*_Nvirt+b)] = W0.cptr()[i*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+a*W0.shape(3)+b];
            }//b
          }//a
        }//i
      }//a1
      // @W1(i,a0,a,b) <<= +1 @D1(a1,a0) @W0(i,a1,a,b)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a0*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a0*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+a*_Nvirt+b];
            }//b
          }//a
        }//i
      }//a0
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+a*_Nvirt+b)*tmp1.shape(1)+(a0)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+a*W1.shape(3)+b];
            }//a0
          }//b
        }//a
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a0];
        }//P
      }//a0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(i,a0,a,b) @X(-1)(+)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -2.0000000000000000;
    orz::DTensor W0(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt,_Ndocc*_Nvirt);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp1.cptr()[(i*_Nvirt+b)*tmp1.shape(1)+(c0*_Nvirt+v0)] = Tc.cptr()[i*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+b*Tc.shape(3)+v0];
            }//v0
          }//c0
        }//b
      }//i
      orz::DTensor tmp2(_Ndocc*_Nvirt,_Nvirt*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp2.cptr()[(c0*_Nvirt+v0)*tmp2.shape(1)+(a*_Nact+a1)] = V2_FiC_vavc.cptr()[a*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a1*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+v0*V2_FiC_vavc.shape(3)+c0];
            }//a1
          }//a
        }//v0
      }//c0
      // @W0(i,a1,b,a) <<= +1 @Tc(i,c0,b,v0) @V2_FiC_vavc(a,a1,v0,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              W0.cptr()[i*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[i*_Nvirt*_Nvirt*_Nact+b*_Nvirt*_Nact+a*_Nact+a1];
            }//a1
          }//a
        }//b
      }//i
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+b*_Nvirt+a)] = W0.cptr()[i*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+b*W0.shape(3)+a];
            }//a
          }//b
        }//i
      }//a1
      // @W1(i,a0,b,a) <<= +1 @D1(a1,a0) @W0(i,a1,b,a)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a0*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[a0*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+b*_Nvirt+a];
            }//a
          }//b
        }//i
      }//a0
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(a0)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+b*W1.shape(3)+a];
            }//a0
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(i,a0,b,a) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//b
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -2.0000000000000000;
    orz::DTensor W0(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt,_Ndocc*_Nvirt);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp1.cptr()[(i*_Nvirt+a)*tmp1.shape(1)+(c0*_Nvirt+v0)] = Tc.cptr()[i*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+a*Tc.shape(3)+v0];
            }//v0
          }//c0
        }//a
      }//i
      orz::DTensor tmp2(_Ndocc*_Nvirt,_Nvirt*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp2.cptr()[(c0*_Nvirt+v0)*tmp2.shape(1)+(b*_Nact+a1)] = V2_FiC_vavc.cptr()[b*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a1*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+v0*V2_FiC_vavc.shape(3)+c0];
            }//a1
          }//b
        }//v0
      }//c0
      // @W0(i,a1,a,b) <<= +1 @Tc(i,c0,a,v0) @V2_FiC_vavc(b,a1,v0,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              W0.cptr()[i*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[i*_Nvirt*_Nvirt*_Nact+a*_Nvirt*_Nact+b*_Nact+a1];
            }//a1
          }//b
        }//a
      }//i
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+a*_Nvirt+b)] = W0.cptr()[i*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+a*W0.shape(3)+b];
            }//b
          }//a
        }//i
      }//a1
      // @W1(i,a0,a,b) <<= +1 @D1(a1,a0) @W0(i,a1,a,b)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a0*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a0*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+a*_Nvirt+b];
            }//b
          }//a
        }//i
      }//a0
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+a*_Nvirt+b)*tmp1.shape(1)+(a0)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+a*W1.shape(3)+b];
            }//a0
          }//b
        }//a
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a0];
        }//P
      }//a0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(i,a0,a,b) @X(-1)(+)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =   4.0000000000000000;
    orz::DTensor W0(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt,_Ndocc*_Nvirt);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp1.cptr()[(i*_Nvirt+a)*tmp1.shape(1)+(c0*_Nvirt+v0)] = Tc.cptr()[i*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+a*Tc.shape(3)+v0];
            }//v0
          }//c0
        }//a
      }//i
      orz::DTensor tmp2(_Ndocc*_Nvirt,_Nvirt*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp2.cptr()[(c0*_Nvirt+v0)*tmp2.shape(1)+(b*_Nact+a1)] = V2_FiC_vavc.cptr()[b*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a1*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+v0*V2_FiC_vavc.shape(3)+c0];
            }//a1
          }//b
        }//v0
      }//c0
      // @W0(i,a1,a,b) <<= +1 @Tc(i,c0,a,v0) @V2_FiC_vavc(b,a1,v0,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              W0.cptr()[i*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[i*_Nvirt*_Nvirt*_Nact+a*_Nvirt*_Nact+b*_Nact+a1];
            }//a1
          }//b
        }//a
      }//i
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+a*_Nvirt+b)] = W0.cptr()[i*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+a*W0.shape(3)+b];
            }//b
          }//a
        }//i
      }//a1
      // @W1(i,a0,a,b) <<= +1 @D1(a1,a0) @W0(i,a1,a,b)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a0*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a0*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+a*_Nvirt+b];
            }//b
          }//a
        }//i
      }//a0
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+a*_Nvirt+b)*tmp1.shape(1)+(a0)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+a*W1.shape(3)+b];
            }//a0
          }//b
        }//a
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(i,a0,a,b) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =   4.0000000000000000;
    orz::DTensor W0(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt,_Ndocc*_Nvirt);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp1.cptr()[(i*_Nvirt+b)*tmp1.shape(1)+(c0*_Nvirt+v0)] = Tc.cptr()[i*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+b*Tc.shape(3)+v0];
            }//v0
          }//c0
        }//b
      }//i
      orz::DTensor tmp2(_Ndocc*_Nvirt,_Nvirt*_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              tmp2.cptr()[(c0*_Nvirt+v0)*tmp2.shape(1)+(a*_Nact+a1)] = V2_FiC_vavc.cptr()[a*V2_FiC_vavc.shape(1)*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+a1*V2_FiC_vavc.shape(2)*V2_FiC_vavc.shape(3)+v0*V2_FiC_vavc.shape(3)+c0];
            }//a1
          }//a
        }//v0
      }//c0
      // @W0(i,a1,b,a) <<= +1 @Tc(i,c0,b,v0) @V2_FiC_vavc(a,a1,v0,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a1 : irange(0L, _Nact)){
              W0.cptr()[i*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[i*_Nvirt*_Nvirt*_Nact+b*_Nvirt*_Nact+a*_Nact+a1];
            }//a1
          }//a
        }//b
      }//i
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+b*_Nvirt+a)] = W0.cptr()[i*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+b*W0.shape(3)+a];
            }//a
          }//b
        }//i
      }//a1
      // @W1(i,a0,b,a) <<= +1 @D1(a1,a0) @W0(i,a1,b,a)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a0*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[a0*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+b*_Nvirt+a];
            }//a
          }//b
        }//i
      }//a0
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(a0)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+b*W1.shape(3)+a];
            }//a0
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a0];
        }//P
      }//a0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(i,a0,b,a) @X(-1)(+)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//b
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -2.0000000000000000;
    orz::DTensor W0(_Ndocc,_Nact);
    {
      orz::DTensor tmp1(_Nact,_Nact*_Nact*_Nact);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a2)*tmp1.shape(1)+(a3*_Nact*_Nact+a1*_Nact+a0)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a3
      }//a2
      orz::DTensor tmp2(_Nact*_Nact*_Nact,_Ndocc);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp2.cptr()[(a3*_Nact*_Nact+a1*_Nact+a0)*tmp2.shape(1)+(c0)] = V2_FiC_aaac.cptr()[a0*V2_FiC_aaac.shape(1)*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a1*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a3*V2_FiC_aaac.shape(3)+c0];
            }//c0
          }//a0
        }//a1
      }//a3
      // @W0(c0,a2) <<= +1 @D2(a3,a2,a1,a0) @V2_FiC_aaac(a0,a1,a3,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          W0.cptr()[c0*_Nact+a2] += tmp3.cptr()[a2*_Ndocc+c0];
        }//c0
      }//a2
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Ndocc);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+a*_Nvirt+b)*tmp1.shape(1)+(c0)] = Tc.cptr()[i*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+a*Tc.shape(3)+b];
            }//c0
          }//b
        }//a
      }//i
      orz::DTensor tmp2(_Ndocc,_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&a2 : irange(0L, _Nact)){
          tmp2.cptr()[(c0)*tmp2.shape(1)+(a2)] = W0.cptr()[c0*W0.shape(1)+a2];
        }//a2
      }//c0
      // @W1(i,a2,a,b) <<= +1 @Tc(i,c0,a,b) @W0(c0,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a2 : irange(0L, _Nact)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[i*_Nvirt*_Nvirt*_Nact+a*_Nvirt*_Nact+b*_Nact+a2];
            }//a2
          }//b
        }//a
      }//i
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+a*_Nvirt+b)*tmp1.shape(1)+(a2)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+a*W1.shape(3)+b];
            }//a2
          }//b
        }//a
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a2)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a2];
        }//P
      }//a2
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(i,a2,a,b) @X(-1)(-)(P,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -2.0000000000000000;
    orz::DTensor W0(_Ndocc,_Nact);
    {
      orz::DTensor tmp1(_Nact,_Nact*_Nact*_Nact);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&a3 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a2)*tmp1.shape(1)+(a3*_Nact*_Nact+a1*_Nact+a0)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a2*d2.shape(2)*d2.shape(3)+a1*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a3
      }//a2
      orz::DTensor tmp2(_Nact*_Nact*_Nact,_Ndocc);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          for(auto &&a0 : irange(0L, _Nact)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp2.cptr()[(a3*_Nact*_Nact+a1*_Nact+a0)*tmp2.shape(1)+(c0)] = V2_FiC_aaac.cptr()[a0*V2_FiC_aaac.shape(1)*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a1*V2_FiC_aaac.shape(2)*V2_FiC_aaac.shape(3)+a3*V2_FiC_aaac.shape(3)+c0];
            }//c0
          }//a0
        }//a1
      }//a3
      // @W0(c0,a2) <<= +1 @D2(a3,a2,a1,a0) @V2_FiC_aaac(a0,a1,a3,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          W0.cptr()[c0*_Nact+a2] += tmp3.cptr()[a2*_Ndocc+c0];
        }//c0
      }//a2
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Ndocc);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(c0)] = Tc.cptr()[i*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+b*Tc.shape(3)+a];
            }//c0
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Ndocc,_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&a2 : irange(0L, _Nact)){
          tmp2.cptr()[(c0)*tmp2.shape(1)+(a2)] = W0.cptr()[c0*W0.shape(1)+a2];
        }//a2
      }//c0
      // @W1(i,a2,b,a) <<= +1 @Tc(i,c0,b,a) @W0(c0,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a2 : irange(0L, _Nact)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a2*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[i*_Nvirt*_Nvirt*_Nact+b*_Nvirt*_Nact+a*_Nact+a2];
            }//a2
          }//a
        }//b
      }//i
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(a2)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a2*W1.shape(2)*W1.shape(3)+b*W1.shape(3)+a];
            }//a2
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a2 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a2)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a2];
        }//P
      }//a2
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(i,a2,b,a) @X(-1)(+)(P,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//b
      }//i
    }
  }//End Contra
  {
    const double _X_ =   2.0000000000000000;
    orz::DTensor W0(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Ndocc*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c1 : irange(0L, _Ndocc)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*_Nvirt+a)*tmp1.shape(1)+(c1*_Ndocc+c0)] = Tc.cptr()[c1*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+b*Tc.shape(3)+a];
            }//c0
          }//c1
        }//a
      }//b
      orz::DTensor tmp2(_Ndocc*_Ndocc,_Nact*_Ndocc);
      for(auto &&c1 : irange(0L, _Ndocc)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&i : irange(0L, _Ndocc)){
              tmp2.cptr()[(c1*_Ndocc+c0)*tmp2.shape(1)+(a1*_Ndocc+i)] = V2_FiC_accc.cptr()[a1*V2_FiC_accc.shape(1)*V2_FiC_accc.shape(2)*V2_FiC_accc.shape(3)+c1*V2_FiC_accc.shape(2)*V2_FiC_accc.shape(3)+c0*V2_FiC_accc.shape(3)+i];
            }//i
          }//a1
        }//c0
      }//c1
      // @W0(i,a1,b,a) <<= +1 @Tc(c1,c0,b,a) @V2_FiC_accc(a1,c1,c0,i)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&i : irange(0L, _Ndocc)){
              W0.cptr()[i*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[b*_Nvirt*_Nact*_Ndocc+a*_Nact*_Ndocc+a1*_Ndocc+i];
            }//i
          }//a1
        }//a
      }//b
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+b*_Nvirt+a)] = W0.cptr()[i*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+b*W0.shape(3)+a];
            }//a
          }//b
        }//i
      }//a1
      // @W1(i,a0,b,a) <<= +1 @D1(a1,a0) @W0(i,a1,b,a)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a0*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[a0*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+b*_Nvirt+a];
            }//a
          }//b
        }//i
      }//a0
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(a0)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+b*W1.shape(3)+a];
            }//a0
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(i,a0,b,a) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//b
      }//i
    }
  }//End Contra
  {
    const double _X_ =   2.0000000000000000;
    orz::DTensor W0(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,_Ndocc*_Ndocc);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c1 : irange(0L, _Ndocc)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(b*_Nvirt+a)*tmp1.shape(1)+(c1*_Ndocc+c0)] = Tc.cptr()[c1*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+b*Tc.shape(3)+a];
            }//c0
          }//c1
        }//a
      }//b
      orz::DTensor tmp2(_Ndocc*_Ndocc,_Nact*_Ndocc);
      for(auto &&c1 : irange(0L, _Ndocc)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&i : irange(0L, _Ndocc)){
              tmp2.cptr()[(c1*_Ndocc+c0)*tmp2.shape(1)+(a1*_Ndocc+i)] = V2_FiC_accc.cptr()[a1*V2_FiC_accc.shape(1)*V2_FiC_accc.shape(2)*V2_FiC_accc.shape(3)+c0*V2_FiC_accc.shape(2)*V2_FiC_accc.shape(3)+c1*V2_FiC_accc.shape(3)+i];
            }//i
          }//a1
        }//c0
      }//c1
      // @W0(i,a1,b,a) <<= +1 @Tc(c1,c0,b,a) @V2_FiC_accc(a1,c0,c1,i)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&i : irange(0L, _Ndocc)){
              W0.cptr()[i*_Nact*_Nvirt*_Nvirt+a1*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[b*_Nvirt*_Nact*_Ndocc+a*_Nact*_Ndocc+a1*_Ndocc+i];
            }//i
          }//a1
        }//a
      }//b
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,_Ndocc*_Nvirt*_Nvirt);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tmp2.cptr()[(a1)*tmp2.shape(1)+(i*_Nvirt*_Nvirt+b*_Nvirt+a)] = W0.cptr()[i*W0.shape(1)*W0.shape(2)*W0.shape(3)+a1*W0.shape(2)*W0.shape(3)+b*W0.shape(3)+a];
            }//a
          }//b
        }//i
      }//a1
      // @W1(i,a0,b,a) <<= +1 @D1(a1,a0) @W0(i,a1,b,a)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a0*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[a0*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+b*_Nvirt+a];
            }//a
          }//b
        }//i
      }//a0
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(a0)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+b*W1.shape(3)+a];
            }//a0
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a0];
        }//P
      }//a0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(i,a0,b,a) @X(-1)(+)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//b
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -2.0000000000000000;
    orz::DTensor W0(_Ndocc,_Nact);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,_Ndocc);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(c0)] = CFip.cptr()[c0*CFip.shape(1)+a1];
        }//c0
      }//a1
      // @W0(c0,a0) <<= +1 @D1(a1,a0) @Fcore(c0,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          W0.cptr()[c0*_Nact+a0] += tmp3.cptr()[a0*_Ndocc+c0];
        }//c0
      }//a0
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Ndocc);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+a*_Nvirt+b)*tmp1.shape(1)+(c0)] = Tc.cptr()[i*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+a*Tc.shape(3)+b];
            }//c0
          }//b
        }//a
      }//i
      orz::DTensor tmp2(_Ndocc,_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp2.cptr()[(c0)*tmp2.shape(1)+(a0)] = W0.cptr()[c0*W0.shape(1)+a0];
        }//a0
      }//c0
      // @W1(i,a0,a,b) <<= +1 @Tc(i,c0,a,b) @W0(c0,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a0 : irange(0L, _Nact)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a0*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[i*_Nvirt*_Nvirt*_Nact+a*_Nvirt*_Nact+b*_Nact+a0];
            }//a0
          }//b
        }//a
      }//i
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+a*_Nvirt+b)*tmp1.shape(1)+(a0)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+a*W1.shape(3)+b];
            }//a0
          }//b
        }//a
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_minus.cptr()[P*Cm1_minus.shape(1)+a0];
        }//P
      }//a0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(i,a0,a,b) @X(-1)(-)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+a*_Nvirt*Nrhom1+b*Nrhom1+P];
            }//P
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -2.0000000000000000;
    orz::DTensor W0(_Ndocc,_Nact);
    {
      orz::DTensor tmp1(_Nact,_Nact);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&a1 : irange(0L, _Nact)){
          tmp1.cptr()[(a0)*tmp1.shape(1)+(a1)] = d1.cptr()[a1*d1.shape(1)+a0];
        }//a1
      }//a0
      orz::DTensor tmp2(_Nact,_Ndocc);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          tmp2.cptr()[(a1)*tmp2.shape(1)+(c0)] = CFip.cptr()[c0*CFip.shape(1)+a1];
        }//c0
      }//a1
      // @W0(c0,a0) <<= +1 @D1(a1,a0) @Fcore(c0,a1)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          W0.cptr()[c0*_Nact+a0] += tmp3.cptr()[a0*_Ndocc+c0];
        }//c0
      }//a0
    }
    orz::DTensor W1(_Ndocc,_Nact,_Nvirt,_Nvirt);
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Ndocc);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(c0)] = Tc.cptr()[i*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+b*Tc.shape(3)+a];
            }//c0
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Ndocc,_Nact);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&a0 : irange(0L, _Nact)){
          tmp2.cptr()[(c0)*tmp2.shape(1)+(a0)] = W0.cptr()[c0*W0.shape(1)+a0];
        }//a0
      }//c0
      // @W1(i,a0,b,a) <<= +1 @Tc(i,c0,b,a) @W0(c0,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a0 : irange(0L, _Nact)){
              W1.cptr()[i*_Nact*_Nvirt*_Nvirt+a0*_Nvirt*_Nvirt+b*_Nvirt+a] += tmp3.cptr()[i*_Nvirt*_Nvirt*_Nact+b*_Nvirt*_Nact+a*_Nact+a0];
            }//a0
          }//a
        }//b
      }//i
    }
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(a0)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)*W1.shape(3)+a0*W1.shape(2)*W1.shape(3)+b*W1.shape(3)+a];
            }//a0
          }//a
        }//b
      }//i
      orz::DTensor tmp2(_Nact,Nrhom1);
      for(auto &&a0 : irange(0L, _Nact)){
        for(auto &&P : irange(0L, Nrhom1)){
          tmp2.cptr()[(a0)*tmp2.shape(1)+(P)] = Cm1_plus.cptr()[P*Cm1_plus.shape(1)+a0];
        }//P
      }//a0
      // @tRFiC(a,b,P,i) <<= +1 _X_ @W1(i,a0,b,a) @X(-1)(+)(P,a0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&P : irange(0L, Nrhom1)){
              tRFiC.cptr()[a*_Nvirt*Nrhom1*_Ndocc+b*Nrhom1*_Ndocc+P*_Ndocc+i] += tmp3.cptr()[i*_Nvirt*_Nvirt*Nrhom1+b*_Nvirt*Nrhom1+a*Nrhom1+P];
            }//P
          }//a
        }//b
      }//i
    }
  }//End Contra
  //genCode.insert(std::make_pair("tRFiC", boost::shared_ptr<orz::DTensor>( new orz::DTensor(tRFiC) )));
 //Code-End-Code-End-Code-End-Code-End-Code-End-Code-End-Code-End
    
      // Transform tRFiC into PNO basis
      for(auto &&P : irange(0L, Nrhom1)){
        for(auto &&i : irange(0L, _Ndocc)){
          const long I_Pi = RActList.GetIndex(P,i);
          if (I_Pi < 0) continue;
          orz::DTensor RP(_Nvirt, _Nvirt);
          for(auto &&a : irange(0L, _Nvirt))
            for(auto &&b : irange(0L, _Nvirt))
              RP(a,b) = tRFiC(a,b,P,i);
          RActList.GetPair(P,i).Transform2EXT(RP);
          orz::DTensor Rorig = R2.GetTensor(P,i);
          RP += Rorig;          
          R2.PutTensor(P,i, RP);
        }
      }
      
      return;
    }//End func
  } }
