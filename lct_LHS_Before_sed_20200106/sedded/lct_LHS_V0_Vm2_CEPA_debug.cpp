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
    void FMatrixOrthogonalizer::FormLHS_V0_Vm2_CEd(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &c2, const orz::DTensor &d3, DensityPack &dPack,
                                                   TensorContainer &ovl, TensorContainer &PNO4, const bool do4EXTinPNO,
                                                   const double QE0, const double Ecas, const double h0_energy, TensorContainer &DebugInts, const orz::DTensor &casfock, const orz::DTensor &FV,
                                                   TensorContainer &T2, TensorContainer &R2){
      
#include "lct_preparations_orthdens.cpp"      
      
      // Back-transform PNO ampitude into canonical basis
      orz::DTensor T_m2_(_Nvirt,_Nvirt,Nrhom2);
      for(auto &&P : irange(0L, Nrhom2)){
        const long I_P = TActList.GetIndex(P);
        if (I_P < 0) continue;
        orz::DTensor U2 = T2.GetTensor(P);
        TActList.GetPair(P).BackTransform2EXT(U2);
        for(auto &&a : irange(0L, _Nvirt))
          for(auto &&b : irange(0L, _Nvirt))
            T_m2_(a,b,P) = U2(a,b);
      }
        
#include "lct_integrals4EXT.cpp"

 //Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn
  orz::DTensor tRFiC(_Ndocc,_Ndocc,_Nvirt,_Nvirt);
  {
    const double _X_ =   0.5000000000000000;
    orz::DTensor W0(_Ndocc,_Ndocc,_Nact,_Nact);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact+a2)*tmp1.shape(1)+(a1*_Nact+a0)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,_Ndocc*_Ndocc);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&j : irange(0L, _Ndocc)){
              tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(i*_Ndocc+j)] = V2_FiC_acac.cptr()[a1*V2_FiC_acac.shape(1)*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+i*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+a0*V2_FiC_acac.shape(3)+j];
            }//j
          }//i
        }//a0
      }//a1
      // @W0(i,j,a3,a2) <<= +1 @D2(a3,a1,a2,a0) @V2_FiC_acac(a1,i,a0,j)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&j : irange(0L, _Ndocc)){
              W0.cptr()[i*_Ndocc*_Nact*_Nact+j*_Nact*_Nact+a3*_Nact+a2] += tmp3.cptr()[a3*_Nact*_Ndocc*_Ndocc+a2*_Ndocc*_Ndocc+i*_Ndocc+j];
            }//j
          }//i
        }//a2
      }//a3
    }
    orz::DTensor W1(_Ndocc,_Ndocc,Nrhom2);
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc,_Nact*_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Ndocc+j)*tmp1.shape(1)+(a3*_Nact+a2)] = W0.cptr()[i*W0.shape(1)*W0.shape(2)*W0.shape(3)+j*W0.shape(2)*W0.shape(3)+a3*W0.shape(3)+a2];
            }//a2
          }//a3
        }//j
      }//i
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a2)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a3*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+_P_];
          }//_P_
        }//a2
      }//a3
      // @W1(i,j,_P_) <<= +1 @W0(i,j,a3,a2) @X(-2)(_P_,a3,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            W1.cptr()[i*_Ndocc*Nrhom2+j*Nrhom2+_P_] += tmp3.cptr()[i*_Ndocc*Nrhom2+j*Nrhom2+_P_];
          }//_P_
        }//j
      }//i
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(a*_Nvirt+b)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[a*T_m2_.shape(1)*T_m2_.shape(2)+b*T_m2_.shape(2)+_P_];
          }//_P_
        }//b
      }//a
      orz::DTensor tmp2(Nrhom2,_Ndocc*_Ndocc);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&j : irange(0L, _Ndocc)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Ndocc+j)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)+j*W1.shape(2)+_P_];
          }//j
        }//i
      }//_P_
      // @tRFiC(i,j,a,b) <<= +1 _X_ @T(-2)(a,b,_P_) @W1(i,j,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&j : irange(0L, _Ndocc)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[a*_Nvirt*_Ndocc*_Ndocc+b*_Ndocc*_Ndocc+i*_Ndocc+j];
            }//j
          }//i
        }//b
      }//a
    }
  }//End Contra
  {
    const double _X_ =   0.5000000000000000;
    orz::DTensor W0(_Ndocc,_Ndocc,_Nact,_Nact);
    {
      orz::DTensor tmp1(_Nact*_Nact,_Nact*_Nact);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&a1 : irange(0L, _Nact)){
            for(auto &&a0 : irange(0L, _Nact)){
              tmp1.cptr()[(a3*_Nact+a2)*tmp1.shape(1)+(a1*_Nact+a0)] = d2.cptr()[a3*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a0];
            }//a0
          }//a1
        }//a2
      }//a3
      orz::DTensor tmp2(_Nact*_Nact,_Ndocc*_Ndocc);
      for(auto &&a1 : irange(0L, _Nact)){
        for(auto &&a0 : irange(0L, _Nact)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&j : irange(0L, _Ndocc)){
              tmp2.cptr()[(a1*_Nact+a0)*tmp2.shape(1)+(i*_Ndocc+j)] = V2_FiC_acac.cptr()[a0*V2_FiC_acac.shape(1)*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+i*V2_FiC_acac.shape(2)*V2_FiC_acac.shape(3)+a1*V2_FiC_acac.shape(3)+j];
            }//j
          }//i
        }//a0
      }//a1
      // @W0(i,j,a3,a2) <<= +1 @D2(a3,a1,a2,a0) @V2_FiC_acac(a0,i,a1,j)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&j : irange(0L, _Ndocc)){
              W0.cptr()[i*_Ndocc*_Nact*_Nact+j*_Nact*_Nact+a3*_Nact+a2] += tmp3.cptr()[a3*_Nact*_Ndocc*_Ndocc+a2*_Ndocc*_Ndocc+i*_Ndocc+j];
            }//j
          }//i
        }//a2
      }//a3
    }
    orz::DTensor W1(_Ndocc,_Ndocc,Nrhom2);
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc,_Nact*_Nact);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&a3 : irange(0L, _Nact)){
            for(auto &&a2 : irange(0L, _Nact)){
              tmp1.cptr()[(i*_Ndocc+j)*tmp1.shape(1)+(a3*_Nact+a2)] = W0.cptr()[i*W0.shape(1)*W0.shape(2)*W0.shape(3)+j*W0.shape(2)*W0.shape(3)+a3*W0.shape(3)+a2];
            }//a2
          }//a3
        }//j
      }//i
      orz::DTensor tmp2(_Nact*_Nact,Nrhom2);
      for(auto &&a3 : irange(0L, _Nact)){
        for(auto &&a2 : irange(0L, _Nact)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp2.cptr()[(a3*_Nact+a2)*tmp2.shape(1)+(_P_)] = Cm2.cptr()[a3*Cm2.shape(1)*Cm2.shape(2)+a2*Cm2.shape(2)+_P_];
          }//_P_
        }//a2
      }//a3
      // @W1(i,j,_P_) <<= +1 @W0(i,j,a3,a2) @X(-2)(_P_,a3,a2)
      orz::DTensor tmp3(tmp1*tmp2);
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&j : irange(0L, _Ndocc)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            W1.cptr()[i*_Ndocc*Nrhom2+j*Nrhom2+_P_] += tmp3.cptr()[i*_Ndocc*Nrhom2+j*Nrhom2+_P_];
          }//_P_
        }//j
      }//i
    }
    {
      orz::DTensor tmp1(_Nvirt*_Nvirt,Nrhom2);
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&_P_ : irange(0L, Nrhom2)){
            tmp1.cptr()[(b*_Nvirt+a)*tmp1.shape(1)+(_P_)] = T_m2_.cptr()[b*T_m2_.shape(1)*T_m2_.shape(2)+a*T_m2_.shape(2)+_P_];
          }//_P_
        }//a
      }//b
      orz::DTensor tmp2(Nrhom2,_Ndocc*_Ndocc);
      for(auto &&_P_ : irange(0L, Nrhom2)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&j : irange(0L, _Ndocc)){
            tmp2.cptr()[(_P_)*tmp2.shape(1)+(i*_Ndocc+j)] = W1.cptr()[i*W1.shape(1)*W1.shape(2)+j*W1.shape(2)+_P_];
          }//j
        }//i
      }//_P_
      // @tRFiC(i,j,a,b) <<= +1 _X_ @T(-2)(b,a,_P_) @W1(i,j,_P_)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&j : irange(0L, _Ndocc)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[b*_Nvirt*_Ndocc*_Ndocc+a*_Ndocc*_Ndocc+i*_Ndocc+j];
            }//j
          }//i
        }//a
      }//b
    }
  }//End Contra
  //genCode.insert(std::make_pair("tRFiC", boost::shared_ptr<orz::DTensor>( new orz::DTensor(tRFiC) )));
 //Code-End-Code-End-Code-End-Code-End-Code-End-Code-End-Code-End      
    
      // Transform tRFiC into PNO basis
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&j : irange(0L, _Ndocc)){
          const long I_ij = RActList.GetIndex(i, j);
          if (I_ij < 0) continue;
          orz::DTensor Rij(_Nvirt, _Nvirt);
          for(auto &&a : irange(0L, _Nvirt))
            for(auto &&b : irange(0L, _Nvirt))
              Rij(a,b) = tRFiC(i,j,a,b);
          RActList.GetPair(i,j).Transform2EXT(Rij);
          orz::DTensor Rorig = R2.GetTensor(i,j);
          Rij += Rorig;
          R2.PutTensor(i, j, Rij);
        }
      }
            
      return;
    }//End func
  } }
