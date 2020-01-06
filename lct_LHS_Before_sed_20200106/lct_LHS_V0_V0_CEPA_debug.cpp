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
    void FMatrixOrthogonalizer::FormLHS_V0_V0_CEd(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &d3, DensityPack &dPack,
                                                  TensorContainer &ovl, TensorContainer &PNO4, const bool do4EXTinPNO,
                                                  const double QE0, const double Ecas, const double h0_energy, TensorContainer &DebugInts, const orz::DTensor &casfock, const orz::DTensor &FV,
                                                  TensorContainer &T2, TensorContainer &R2){

      orz::ProgressTimer pt("* V(0)/V(0)");
      
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
  orz::DTensor tRFiC(_Ndocc,_Ndocc,_Nvirt,_Nvirt);
  if(!do4EXTinPNO)
  {
    orz::ProgressTimer pt("* V(0)/V(0) - 4ext (canonical) ");    
    const double _X_ =   4.0000000000000000;
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc,_Nvirt*_Nvirt);
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&v1 : irange(0L, _Nvirt)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp1.cptr()[(j*_Ndocc+i)*tmp1.shape(1)+(v1*_Nvirt+v0)] = Tc.cptr()[j*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+i*Tc.shape(2)*Tc.shape(3)+v1*Tc.shape(3)+v0];
            }//v0
          }//v1
        }//i
      }//j
      orz::DTensor tmp2(_Nvirt*_Nvirt,_Nvirt*_Nvirt);
      for(auto &&v1 : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tmp2.cptr()[(v1*_Nvirt+v0)*tmp2.shape(1)+(a*_Nvirt+b)] = V2_FiC_vvvv.cptr()[a*V2_FiC_vvvv.shape(1)*V2_FiC_vvvv.shape(2)*V2_FiC_vvvv.shape(3)+v0*V2_FiC_vvvv.shape(2)*V2_FiC_vvvv.shape(3)+b*V2_FiC_vvvv.shape(3)+v1];
            }//b
          }//a
        }//v0
      }//v1
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(j,i,v1,v0) @V2_FiC_vvvv(a,v0,b,v1)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[j*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+a*_Nvirt+b];
            }//b
          }//a
        }//i
      }//j
    }
  }//End Contra
  else
  {
    // PNO-based construction of 4-external terms
    orz::ProgressTimer pt("* V(0)/V(0) - 4ext (PNO) ");
    for(auto &&j : irange(0L, _Ndocc)){
      for(auto &&i : irange(0L, _Ndocc)){
        const long I_ji = TActList.GetIndex(j, i);
        if (I_ji < 0) continue;
        PairData Pji = TActList.GetPair(j, i);
        const long Nvir = Pji.GetNPNO();
        orz::DTensor U2 = T2.GetTensor(j, i);
        orz::DTensor U2tilde(U2.copy());
        U2tilde.gaxpy(+2.0, U2.swapdim(1,0), -1.0);
        //orz::DTensor I4 = PNO4.GetTensor(std::max(i,j),std::min(i,j)).copy();
        orz::DTensor I4;
        PNO4.GetTensor(std::max(i,j),std::min(i,j),I4);
        // Sort I4        
        orz::DTensor tmp1(Nvir*Nvir,Nvir*Nvir);        
        const long Nvir2 = Nvir*(Nvir+1)/2;
        for(auto &&a : irange(0L, Nvir)){
          for(auto &&c : irange(0L, Nvir)){
            const long idx_ac = (c<a ? a*(a+1)/2+c : c*(c+1)/2+a);  
            for(auto &&b : irange(0L, Nvir)){
              for(auto &&d : irange(0L, Nvir)){
                const long idx_bd = (d<b ? b*(b+1)/2+d : d*(d+1)/2+b);
                const long idx_acbd = (idx_bd<idx_ac ? idx_ac*(idx_ac+1)/2+idx_bd : idx_bd*(idx_bd+1)/2+idx_ac);
                tmp1.cptr()[c*Nvir*Nvir*Nvir+d*Nvir*Nvir+a*Nvir+b] = I4.cptr()[idx_acbd];
              }
            }
          }
        }
        
        // Sort Utilde
        orz::DTensor tmp2(1,Nvir*Nvir);
        for(auto &&a : irange(0L, Nvir))        
          for(auto &&b : irange(0L, Nvir))
            tmp2.cptr()[a*Nvir+b] = U2tilde.cptr()[a*Nvir+b];
        orz::DTensor tmp3(tmp2*tmp1);
        tmp3 *= 4.0;
        // Sort S2
        orz::DTensor S2 = R2.GetTensor(i, j);
        for(auto &&a : irange(0L, Nvir))        
          for(auto &&b : irange(0L, Nvir))
            S2.cptr()[a*S2.shape(1)+b] += tmp3.cptr()[b*Nvir+a];
        R2.PutTensor(i,j,S2);
      }//i
    }//j
  }
  {
    const double _X_ =  -4.0000000000000000;
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
      orz::DTensor tmp2(_Ndocc*_Nvirt,_Nvirt*_Ndocc);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&j : irange(0L, _Ndocc)){
              tmp2.cptr()[(c0*_Nvirt+v0)*tmp2.shape(1)+(a*_Ndocc+j)] = V2_FiC_vvcc.cptr()[a*V2_FiC_vvcc.shape(1)*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+v0*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+c0*V2_FiC_vvcc.shape(3)+j];
            }//j
          }//a
        }//v0
      }//c0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(i,c0,v0,b) @V2_FiC_vvcc(a,v0,c0,j)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&j : irange(0L, _Ndocc)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[i*_Nvirt*_Nvirt*_Ndocc+b*_Nvirt*_Ndocc+a*_Ndocc+j];
            }//j
          }//a
        }//b
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -4.0000000000000000;
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt,_Ndocc*_Nvirt);
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp1.cptr()[(j*_Nvirt+b)*tmp1.shape(1)+(c0*_Nvirt+v0)] = Tc.cptr()[j*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+b*Tc.shape(3)+v0];
            }//v0
          }//c0
        }//b
      }//j
      orz::DTensor tmp2(_Ndocc*_Nvirt,_Nvirt*_Ndocc);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&i : irange(0L, _Ndocc)){
              tmp2.cptr()[(c0*_Nvirt+v0)*tmp2.shape(1)+(a*_Ndocc+i)] = V2_FiC_vvcc.cptr()[a*V2_FiC_vvcc.shape(1)*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+v0*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+c0*V2_FiC_vvcc.shape(3)+i];
            }//i
          }//a
        }//v0
      }//c0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(j,c0,b,v0) @V2_FiC_vvcc(a,v0,c0,i)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&i : irange(0L, _Ndocc)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[j*_Nvirt*_Nvirt*_Ndocc+b*_Nvirt*_Ndocc+a*_Ndocc+i];
            }//i
          }//a
        }//b
      }//j
    }
  }//End Contra
  {
    const double _X_ =  -4.0000000000000000;
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
      orz::DTensor tmp2(_Ndocc*_Nvirt,_Nvirt*_Ndocc);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&j : irange(0L, _Ndocc)){
              tmp2.cptr()[(c0*_Nvirt+v0)*tmp2.shape(1)+(b*_Ndocc+j)] = V2_FiC_vvcc.cptr()[b*V2_FiC_vvcc.shape(1)*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+v0*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+c0*V2_FiC_vvcc.shape(3)+j];
            }//j
          }//b
        }//v0
      }//c0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(i,c0,a,v0) @V2_FiC_vvcc(b,v0,c0,j)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&j : irange(0L, _Ndocc)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[i*_Nvirt*_Nvirt*_Ndocc+a*_Nvirt*_Ndocc+b*_Ndocc+j];
            }//j
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -4.0000000000000000;
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt,_Ndocc*_Nvirt);
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp1.cptr()[(j*_Nvirt+a)*tmp1.shape(1)+(c0*_Nvirt+v0)] = Tc.cptr()[j*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+v0*Tc.shape(3)+a];
            }//v0
          }//c0
        }//a
      }//j
      orz::DTensor tmp2(_Ndocc*_Nvirt,_Nvirt*_Ndocc);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&i : irange(0L, _Ndocc)){
              tmp2.cptr()[(c0*_Nvirt+v0)*tmp2.shape(1)+(b*_Ndocc+i)] = V2_FiC_vvcc.cptr()[b*V2_FiC_vvcc.shape(1)*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+v0*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+c0*V2_FiC_vvcc.shape(3)+i];
            }//i
          }//b
        }//v0
      }//c0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(j,c0,v0,a) @V2_FiC_vvcc(b,v0,c0,i)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&i : irange(0L, _Ndocc)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[j*_Nvirt*_Nvirt*_Ndocc+a*_Nvirt*_Ndocc+b*_Ndocc+i];
            }//i
          }//b
        }//a
      }//j
    }
  }//End Contra
  {
    const double _X_ =  -4.0000000000000000;
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt,_Ndocc*_Nvirt);
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp1.cptr()[(j*_Nvirt+a)*tmp1.shape(1)+(c0*_Nvirt+v0)] = Tc.cptr()[j*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+a*Tc.shape(3)+v0];
            }//v0
          }//c0
        }//a
      }//j
      orz::DTensor tmp2(_Ndocc*_Nvirt,_Nvirt*_Ndocc);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&i : irange(0L, _Ndocc)){
              tmp2.cptr()[(c0*_Nvirt+v0)*tmp2.shape(1)+(b*_Ndocc+i)] = V2_FiC_vcvc.cptr()[b*V2_FiC_vcvc.shape(1)*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+i*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+v0*V2_FiC_vcvc.shape(3)+c0];
            }//i
          }//b
        }//v0
      }//c0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(j,c0,a,v0) @V2_FiC_vcvc(b,i,v0,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&i : irange(0L, _Ndocc)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[j*_Nvirt*_Nvirt*_Ndocc+a*_Nvirt*_Ndocc+b*_Ndocc+i];
            }//i
          }//b
        }//a
      }//j
    }
  }//End Contra
  {
    const double _X_ =   8.0000000000000000;
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
      orz::DTensor tmp2(_Ndocc*_Nvirt,_Nvirt*_Ndocc);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&j : irange(0L, _Ndocc)){
              tmp2.cptr()[(c0*_Nvirt+v0)*tmp2.shape(1)+(b*_Ndocc+j)] = V2_FiC_vcvc.cptr()[b*V2_FiC_vcvc.shape(1)*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+j*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+v0*V2_FiC_vcvc.shape(3)+c0];
            }//j
          }//b
        }//v0
      }//c0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(i,c0,a,v0) @V2_FiC_vcvc(b,j,v0,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&j : irange(0L, _Ndocc)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[i*_Nvirt*_Nvirt*_Ndocc+a*_Nvirt*_Ndocc+b*_Ndocc+j];
            }//j
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =  -4.0000000000000000;
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
      orz::DTensor tmp2(_Ndocc*_Nvirt,_Nvirt*_Ndocc);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&j : irange(0L, _Ndocc)){
              tmp2.cptr()[(c0*_Nvirt+v0)*tmp2.shape(1)+(a*_Ndocc+j)] = V2_FiC_vcvc.cptr()[a*V2_FiC_vcvc.shape(1)*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+j*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+v0*V2_FiC_vcvc.shape(3)+c0];
            }//j
          }//a
        }//v0
      }//c0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(i,c0,b,v0) @V2_FiC_vcvc(a,j,v0,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&j : irange(0L, _Ndocc)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[i*_Nvirt*_Nvirt*_Ndocc+b*_Nvirt*_Ndocc+a*_Ndocc+j];
            }//j
          }//a
        }//b
      }//i
    }
  }//End Contra
  {
    const double _X_ =   8.0000000000000000;
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt,_Ndocc*_Nvirt);
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&c0 : irange(0L, _Ndocc)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp1.cptr()[(j*_Nvirt+b)*tmp1.shape(1)+(c0*_Nvirt+v0)] = Tc.cptr()[j*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+b*Tc.shape(3)+v0];
            }//v0
          }//c0
        }//b
      }//j
      orz::DTensor tmp2(_Ndocc*_Nvirt,_Nvirt*_Ndocc);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&i : irange(0L, _Ndocc)){
              tmp2.cptr()[(c0*_Nvirt+v0)*tmp2.shape(1)+(a*_Ndocc+i)] = V2_FiC_vcvc.cptr()[a*V2_FiC_vcvc.shape(1)*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+i*V2_FiC_vcvc.shape(2)*V2_FiC_vcvc.shape(3)+v0*V2_FiC_vcvc.shape(3)+c0];
            }//i
          }//a
        }//v0
      }//c0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(j,c0,b,v0) @V2_FiC_vcvc(a,i,v0,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&i : irange(0L, _Ndocc)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[j*_Nvirt*_Nvirt*_Ndocc+b*_Nvirt*_Ndocc+a*_Ndocc+i];
            }//i
          }//a
        }//b
      }//j
    }
  }//End Contra
  {
    const double _X_ =   4.0000000000000000;
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
      orz::DTensor tmp2(_Ndocc*_Ndocc,_Ndocc*_Ndocc);
      for(auto &&c1 : irange(0L, _Ndocc)){
        for(auto &&c0 : irange(0L, _Ndocc)){
          for(auto &&i : irange(0L, _Ndocc)){
            for(auto &&j : irange(0L, _Ndocc)){
              tmp2.cptr()[(c1*_Ndocc+c0)*tmp2.shape(1)+(i*_Ndocc+j)] = V2_FiC_cccc.cptr()[c0*V2_FiC_cccc.shape(1)*V2_FiC_cccc.shape(2)*V2_FiC_cccc.shape(3)+i*V2_FiC_cccc.shape(2)*V2_FiC_cccc.shape(3)+c1*V2_FiC_cccc.shape(3)+j];
            }//j
          }//i
        }//c0
      }//c1
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(c1,c0,b,a) @V2_FiC_cccc(c0,i,c1,j)
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
  {
    const double _X_ =  -4.0000000000000000;
    {
      orz::DTensor tmp1(_Ndocc*_Nvirt*_Nvirt,_Ndocc);
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&c0 : irange(0L, _Ndocc)){
              tmp1.cptr()[(j*_Nvirt*_Nvirt+b*_Nvirt+a)*tmp1.shape(1)+(c0)] = Tc.cptr()[j*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+b*Tc.shape(3)+a];
            }//c0
          }//a
        }//b
      }//j
      orz::DTensor tmp2(_Ndocc,_Ndocc);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          tmp2.cptr()[(c0)*tmp2.shape(1)+(i)] = Fij.cptr()[i*Fij.shape(1)+c0];
        }//i
      }//c0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(j,c0,b,a) @cF(i,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&i : irange(0L, _Ndocc)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[j*_Nvirt*_Nvirt*_Ndocc+b*_Nvirt*_Ndocc+a*_Ndocc+i];
            }//i
          }//a
        }//b
      }//j
    }
  }//End Contra
  {
    const double _X_ =  -4.0000000000000000;
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
      orz::DTensor tmp2(_Ndocc,_Ndocc);
      for(auto &&c0 : irange(0L, _Ndocc)){
        for(auto &&j : irange(0L, _Ndocc)){
          tmp2.cptr()[(c0)*tmp2.shape(1)+(j)] = Fij.cptr()[j*Fij.shape(1)+c0];
        }//j
      }//c0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(i,c0,a,b) @cF(j,c0)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&j : irange(0L, _Ndocc)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[i*_Nvirt*_Nvirt*_Ndocc+a*_Nvirt*_Ndocc+b*_Ndocc+j];
            }//j
          }//b
        }//a
      }//i
    }
  }//End Contra
  {
    const double _X_ =   4.0000000000000000;
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc*_Nvirt,_Nvirt);
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp1.cptr()[(j*_Ndocc*_Nvirt+i*_Nvirt+a)*tmp1.shape(1)+(v0)] = Tc.cptr()[j*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+i*Tc.shape(2)*Tc.shape(3)+v0*Tc.shape(3)+a];
            }//v0
          }//a
        }//i
      }//j
      orz::DTensor tmp2(_Nvirt,_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          tmp2.cptr()[(v0)*tmp2.shape(1)+(b)] = Fab.cptr()[v0*Fab.shape(1)+b];
        }//b
      }//v0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(j,i,v0,a) @cF(v0,b)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&a : irange(0L, _Nvirt)){
            for(auto &&b : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[j*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+a*_Nvirt+b];
            }//b
          }//a
        }//i
      }//j
    }
  }//End Contra
  {
    const double _X_ =   4.0000000000000000;
    {
      orz::DTensor tmp1(_Ndocc*_Ndocc*_Nvirt,_Nvirt);
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&v0 : irange(0L, _Nvirt)){
              tmp1.cptr()[(j*_Ndocc*_Nvirt+i*_Nvirt+b)*tmp1.shape(1)+(v0)] = Tc.cptr()[j*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+i*Tc.shape(2)*Tc.shape(3)+b*Tc.shape(3)+v0];
            }//v0
          }//b
        }//i
      }//j
      orz::DTensor tmp2(_Nvirt,_Nvirt);
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          tmp2.cptr()[(v0)*tmp2.shape(1)+(a)] = Fab.cptr()[v0*Fab.shape(1)+a];
        }//a
      }//v0
      // @tRFiC(i,j,a,b) <<= +1 _X_ @Tc(j,i,b,v0) @cF(v0,a)
      orz::DTensor tmp3(tmp1*tmp2);
      tmp3 *= 1.0 * _X_;
      for(auto &&j : irange(0L, _Ndocc)){
        for(auto &&i : irange(0L, _Ndocc)){
          for(auto &&b : irange(0L, _Nvirt)){
            for(auto &&a : irange(0L, _Nvirt)){
              tRFiC.cptr()[i*_Ndocc*_Nvirt*_Nvirt+j*_Nvirt*_Nvirt+a*_Nvirt+b] += tmp3.cptr()[j*_Ndocc*_Nvirt*_Nvirt+i*_Nvirt*_Nvirt+b*_Nvirt+a];
            }//a
          }//b
        }//i
      }//j
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
