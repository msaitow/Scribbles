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
                                                  const double QE0, const double Ecas, const double h0_energy, TensorContainer &DebugInts, const orz::DTensor &casfock, const orz::DTensor &FV,
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
      orz::DTensor tRFiC(_Ndocc, _Ndocc, _Nvirt, _Nvirt);
      // Fetch the integrals from the integral container
      orz::DTensor V2_FiC_ccvv = DebugInts.GetTensor( 0);
      orz::DTensor V2_FiC_vccv = DebugInts.GetTensor( 1);
      orz::DTensor V2_FiC_cccc = DebugInts.GetTensor( 2);
      orz::DTensor V2_FiC_vvvv = DebugInts.GetTensor( 3);
      orz::DTensor V2_FiC_vvcc = DebugInts.GetTensor( 4);
      orz::DTensor V2_FiC_aavv = DebugInts.GetTensor( 5);
      orz::DTensor V2_FiC_vaav = DebugInts.GetTensor( 6);
      orz::DTensor Fcore       = DebugInts.GetTensor( 7);
      orz::DTensor Ecore_      = DebugInts.GetTensor( 8);
      orz::DTensor V2_FiC_vvaa = DebugInts.GetTensor( 9);
      orz::DTensor V2_FiC_caca = DebugInts.GetTensor(10);
      orz::DTensor V2_FiC_vava = DebugInts.GetTensor(11);
      orz::DTensor V2_FiC_ccaa = DebugInts.GetTensor(12);
      orz::DTensor V2_FiC_aaaa = DebugInts.GetTensor(13);                        
      const double Ecore       = Ecore_(0); 

      assert(V2_FiC_cccc.rank() == 4);
      assert(V2_FiC_vvvv.rank() == 4);

      assert(V2_FiC_cccc.size() == _Ndocc*_Ndocc*_Ndocc*_Ndocc);
      assert(V2_FiC_vvvv.size() == _Nvirt*_Nvirt*_Nvirt*_Nvirt);
      
      assert(V2_FiC_vvvv.shape(0)==_Nvirt);
      assert(V2_FiC_vvvv.shape(1)==_Nvirt);
      assert(V2_FiC_vvvv.shape(2)==_Nvirt);
      assert(V2_FiC_vvvv.shape(3)==_Nvirt);
      
      assert(V2_FiC_cccc.shape(0)==_Ndocc);
      assert(V2_FiC_cccc.shape(1)==_Ndocc);
      assert(V2_FiC_cccc.shape(2)==_Ndocc);
      assert(V2_FiC_cccc.shape(3)==_Ndocc);      
      
//      orz::DTensor FcoreC(Fcore(orz::Slice(           0,              _Ndocc), orz::Slice(           0,              _Ndocc)));
//      orz::DTensor FcoreA(Fcore(orz::Slice(      _Ndocc,        _Ndocc+_Nact), orz::Slice(      _Ndocc,        _Ndocc+_Nact)));
//      orz::DTensor FcoreV(Fcore(orz::Slice(_Ndocc+_Nact, _Ndocc+_Nact+_Nvirt), orz::Slice(_Ndocc+_Nact, _Ndocc+_Nact+_Nvirt)));      
//
//      orz::DTensor Fpq = this->_Fock.copy();
      
      {
 //Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn-Code-Bgn
std::cout << "@tRFiC(i,j,a,b) <<= +4 @Tc(j,i,v1,v0) @V2_FiC_vvvv(v1,b,v0,a)" << std::endl;
// @tRFiC(i,j,a,b) <<= +4 @Tc(j,i,v1,v0) @V2_FiC_vvvv(v1,b,v0,a)
//orz::DTensor tRFiC(_Ndocc,_Ndocc,_Nvirt,_Nvirt);
for(auto &&j : irange(0L, _Ndocc)){
  for(auto &&i : irange(0L, _Ndocc)){
    for(auto &&v1 : irange(0L, _Nvirt)){
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            tRFiC.cptr()[i*tRFiC.shape(1)*tRFiC.shape(2)*tRFiC.shape(3)+j*tRFiC.shape(2)*tRFiC.shape(3)+a*tRFiC.shape(3)+b] +=  4.0 * Tc.cptr()[j*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+i*Tc.shape(2)*Tc.shape(3)+v1*Tc.shape(3)+v0] * V2_FiC_vvvv.cptr()[v1*V2_FiC_vvvv.shape(1)*V2_FiC_vvvv.shape(2)*V2_FiC_vvvv.shape(3)+b*V2_FiC_vvvv.shape(2)*V2_FiC_vvvv.shape(3)+v0*V2_FiC_vvvv.shape(3)+a];
          }//a
        }//b
      }//v0
    }//v1
  }//i
}//j
std::cout << "@tRFiC(i,j,a,b) <<= -4 @Tc(i,c0,v0,b) @V2_FiC_vvcc(v0,a,j,c0)" << std::endl;
// @tRFiC(i,j,a,b) <<= -4 @Tc(i,c0,v0,b) @V2_FiC_vvcc(v0,a,j,c0)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&c0 : irange(0L, _Ndocc)){
    for(auto &&j : irange(0L, _Ndocc)){
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            tRFiC.cptr()[i*tRFiC.shape(1)*tRFiC.shape(2)*tRFiC.shape(3)+j*tRFiC.shape(2)*tRFiC.shape(3)+a*tRFiC.shape(3)+b] -=  4.0 * Tc.cptr()[i*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+v0*Tc.shape(3)+b] * V2_FiC_vvcc.cptr()[v0*V2_FiC_vvcc.shape(1)*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+a*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+j*V2_FiC_vvcc.shape(3)+c0];
          }//a
        }//b
      }//v0
    }//j
  }//c0
}//i
std::cout << "@tRFiC(i,j,a,b) <<= -4 @Tc(j,c0,b,v0) @V2_FiC_vvcc(v0,a,i,c0)" << std::endl;
// @tRFiC(i,j,a,b) <<= -4 @Tc(j,c0,b,v0) @V2_FiC_vvcc(v0,a,i,c0)
for(auto &&j : irange(0L, _Ndocc)){
  for(auto &&c0 : irange(0L, _Ndocc)){
    for(auto &&i : irange(0L, _Ndocc)){
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            tRFiC.cptr()[i*tRFiC.shape(1)*tRFiC.shape(2)*tRFiC.shape(3)+j*tRFiC.shape(2)*tRFiC.shape(3)+a*tRFiC.shape(3)+b] -=  4.0 * Tc.cptr()[j*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+b*Tc.shape(3)+v0] * V2_FiC_vvcc.cptr()[v0*V2_FiC_vvcc.shape(1)*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+a*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+i*V2_FiC_vvcc.shape(3)+c0];
          }//a
        }//v0
      }//b
    }//i
  }//c0
}//j
std::cout << "@tRFiC(i,j,a,b) <<= -4 @Tc(i,c0,a,v0) @V2_FiC_vvcc(v0,b,j,c0)" << std::endl;
// @tRFiC(i,j,a,b) <<= -4 @Tc(i,c0,a,v0) @V2_FiC_vvcc(v0,b,j,c0)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&c0 : irange(0L, _Ndocc)){
    for(auto &&j : irange(0L, _Ndocc)){
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            tRFiC.cptr()[i*tRFiC.shape(1)*tRFiC.shape(2)*tRFiC.shape(3)+j*tRFiC.shape(2)*tRFiC.shape(3)+a*tRFiC.shape(3)+b] -=  4.0 * Tc.cptr()[i*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+a*Tc.shape(3)+v0] * V2_FiC_vvcc.cptr()[v0*V2_FiC_vvcc.shape(1)*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+b*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+j*V2_FiC_vvcc.shape(3)+c0];
          }//b
        }//v0
      }//a
    }//j
  }//c0
}//i
std::cout << "@tRFiC(i,j,a,b) <<= -4 @Tc(j,c0,v0,a) @V2_FiC_vvcc(v0,b,i,c0)" << std::endl;
// @tRFiC(i,j,a,b) <<= -4 @Tc(j,c0,v0,a) @V2_FiC_vvcc(v0,b,i,c0)
for(auto &&j : irange(0L, _Ndocc)){
  for(auto &&c0 : irange(0L, _Ndocc)){
    for(auto &&i : irange(0L, _Ndocc)){
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            tRFiC.cptr()[i*tRFiC.shape(1)*tRFiC.shape(2)*tRFiC.shape(3)+j*tRFiC.shape(2)*tRFiC.shape(3)+a*tRFiC.shape(3)+b] -=  4.0 * Tc.cptr()[j*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+v0*Tc.shape(3)+a] * V2_FiC_vvcc.cptr()[v0*V2_FiC_vvcc.shape(1)*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+b*V2_FiC_vvcc.shape(2)*V2_FiC_vvcc.shape(3)+i*V2_FiC_vvcc.shape(3)+c0];
          }//b
        }//a
      }//v0
    }//i
  }//c0
}//j
std::cout << "@tRFiC(i,j,a,b) <<= -4 @Tc(j,c0,a,v0) @V2_FiC_vccv(v0,c0,i,b)" << std::endl;
// @tRFiC(i,j,a,b) <<= -4 @Tc(j,c0,a,v0) @V2_FiC_vccv(v0,c0,i,b)
for(auto &&j : irange(0L, _Ndocc)){
  for(auto &&c0 : irange(0L, _Ndocc)){
    for(auto &&i : irange(0L, _Ndocc)){
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            tRFiC.cptr()[i*tRFiC.shape(1)*tRFiC.shape(2)*tRFiC.shape(3)+j*tRFiC.shape(2)*tRFiC.shape(3)+a*tRFiC.shape(3)+b] -=  4.0 * Tc.cptr()[j*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+a*Tc.shape(3)+v0] * V2_FiC_vccv.cptr()[v0*V2_FiC_vccv.shape(1)*V2_FiC_vccv.shape(2)*V2_FiC_vccv.shape(3)+c0*V2_FiC_vccv.shape(2)*V2_FiC_vccv.shape(3)+i*V2_FiC_vccv.shape(3)+b];
          }//b
        }//v0
      }//a
    }//i
  }//c0
}//j
std::cout << "@tRFiC(i,j,a,b) <<= +8 @Tc(i,c0,a,v0) @V2_FiC_vccv(v0,c0,j,b)" << std::endl;
// @tRFiC(i,j,a,b) <<= +8 @Tc(i,c0,a,v0) @V2_FiC_vccv(v0,c0,j,b)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&c0 : irange(0L, _Ndocc)){
    for(auto &&j : irange(0L, _Ndocc)){
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&b : irange(0L, _Nvirt)){
            tRFiC.cptr()[i*tRFiC.shape(1)*tRFiC.shape(2)*tRFiC.shape(3)+j*tRFiC.shape(2)*tRFiC.shape(3)+a*tRFiC.shape(3)+b] +=  8.0 * Tc.cptr()[i*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+a*Tc.shape(3)+v0] * V2_FiC_vccv.cptr()[v0*V2_FiC_vccv.shape(1)*V2_FiC_vccv.shape(2)*V2_FiC_vccv.shape(3)+c0*V2_FiC_vccv.shape(2)*V2_FiC_vccv.shape(3)+j*V2_FiC_vccv.shape(3)+b];
          }//b
        }//v0
      }//a
    }//j
  }//c0
}//i
std::cout << "@tRFiC(i,j,a,b) <<= -4 @Tc(i,c0,b,v0) @V2_FiC_vccv(v0,c0,j,a)" << std::endl;
// @tRFiC(i,j,a,b) <<= -4 @Tc(i,c0,b,v0) @V2_FiC_vccv(v0,c0,j,a)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&c0 : irange(0L, _Ndocc)){
    for(auto &&j : irange(0L, _Ndocc)){
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            tRFiC.cptr()[i*tRFiC.shape(1)*tRFiC.shape(2)*tRFiC.shape(3)+j*tRFiC.shape(2)*tRFiC.shape(3)+a*tRFiC.shape(3)+b] -=  4.0 * Tc.cptr()[i*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+b*Tc.shape(3)+v0] * V2_FiC_vccv.cptr()[v0*V2_FiC_vccv.shape(1)*V2_FiC_vccv.shape(2)*V2_FiC_vccv.shape(3)+c0*V2_FiC_vccv.shape(2)*V2_FiC_vccv.shape(3)+j*V2_FiC_vccv.shape(3)+a];
          }//a
        }//v0
      }//b
    }//j
  }//c0
}//i
std::cout << "@tRFiC(i,j,a,b) <<= +8 @Tc(j,c0,b,v0) @V2_FiC_vccv(v0,c0,i,a)" << std::endl;
// @tRFiC(i,j,a,b) <<= +8 @Tc(j,c0,b,v0) @V2_FiC_vccv(v0,c0,i,a)
for(auto &&j : irange(0L, _Ndocc)){
  for(auto &&c0 : irange(0L, _Ndocc)){
    for(auto &&i : irange(0L, _Ndocc)){
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&v0 : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            tRFiC.cptr()[i*tRFiC.shape(1)*tRFiC.shape(2)*tRFiC.shape(3)+j*tRFiC.shape(2)*tRFiC.shape(3)+a*tRFiC.shape(3)+b] +=  8.0 * Tc.cptr()[j*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+b*Tc.shape(3)+v0] * V2_FiC_vccv.cptr()[v0*V2_FiC_vccv.shape(1)*V2_FiC_vccv.shape(2)*V2_FiC_vccv.shape(3)+c0*V2_FiC_vccv.shape(2)*V2_FiC_vccv.shape(3)+i*V2_FiC_vccv.shape(3)+a];
          }//a
        }//v0
      }//b
    }//i
  }//c0
}//j
std::cout << "@tRFiC(i,j,a,b) <<= +4 @Tc(c1,c0,b,a) @V2_FiC_cccc(j,c1,i,c0)" << std::endl;
// @tRFiC(i,j,a,b) <<= +4 @Tc(c1,c0,b,a) @V2_FiC_cccc(j,c1,i,c0)
for(auto &&c1 : irange(0L, _Ndocc)){
  for(auto &&c0 : irange(0L, _Ndocc)){
    for(auto &&j : irange(0L, _Ndocc)){
      for(auto &&i : irange(0L, _Ndocc)){
        for(auto &&b : irange(0L, _Nvirt)){
          for(auto &&a : irange(0L, _Nvirt)){
            tRFiC.cptr()[i*tRFiC.shape(1)*tRFiC.shape(2)*tRFiC.shape(3)+j*tRFiC.shape(2)*tRFiC.shape(3)+a*tRFiC.shape(3)+b] +=  4.0 * Tc.cptr()[c1*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+b*Tc.shape(3)+a] * V2_FiC_cccc.cptr()[j*V2_FiC_cccc.shape(1)*V2_FiC_cccc.shape(2)*V2_FiC_cccc.shape(3)+c1*V2_FiC_cccc.shape(2)*V2_FiC_cccc.shape(3)+i*V2_FiC_cccc.shape(3)+c0];
          }//a
        }//b
      }//i
    }//j
  }//c0
}//c1
std::cout << "@tRFiC(i,j,a,b) <<= -4 @Tc(j,c0,b,a) @cF(i,c0)" << std::endl;
// @tRFiC(i,j,a,b) <<= -4 @Tc(j,c0,b,a) @cF(i,c0)
for(auto &&j : irange(0L, _Ndocc)){
  for(auto &&c0 : irange(0L, _Ndocc)){
    for(auto &&i : irange(0L, _Ndocc)){
      for(auto &&b : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          tRFiC.cptr()[i*tRFiC.shape(1)*tRFiC.shape(2)*tRFiC.shape(3)+j*tRFiC.shape(2)*tRFiC.shape(3)+a*tRFiC.shape(3)+b] -=  4.0 * Tc.cptr()[j*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+b*Tc.shape(3)+a] * Fij.cptr()[i*Fij.shape(1)+c0];
        }//a
      }//b
    }//i
  }//c0
}//j
std::cout << "@tRFiC(i,j,a,b) <<= -4 @Tc(i,c0,a,b) @cF(j,c0)" << std::endl;
// @tRFiC(i,j,a,b) <<= -4 @Tc(i,c0,a,b) @cF(j,c0)
for(auto &&i : irange(0L, _Ndocc)){
  for(auto &&c0 : irange(0L, _Ndocc)){
    for(auto &&j : irange(0L, _Ndocc)){
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          tRFiC.cptr()[i*tRFiC.shape(1)*tRFiC.shape(2)*tRFiC.shape(3)+j*tRFiC.shape(2)*tRFiC.shape(3)+a*tRFiC.shape(3)+b] -=  4.0 * Tc.cptr()[i*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+c0*Tc.shape(2)*Tc.shape(3)+a*Tc.shape(3)+b] * Fij.cptr()[j*Fij.shape(1)+c0];
        }//b
      }//a
    }//j
  }//c0
}//i
std::cout << "@tRFiC(i,j,a,b) <<= +4 @Tc(j,i,v0,a) @cF(v0,b)" << std::endl;
// @tRFiC(i,j,a,b) <<= +4 @Tc(j,i,v0,a) @cF(v0,b)
for(auto &&j : irange(0L, _Ndocc)){
  for(auto &&i : irange(0L, _Ndocc)){
    for(auto &&v0 : irange(0L, _Nvirt)){
      for(auto &&a : irange(0L, _Nvirt)){
        for(auto &&b : irange(0L, _Nvirt)){
          tRFiC.cptr()[i*tRFiC.shape(1)*tRFiC.shape(2)*tRFiC.shape(3)+j*tRFiC.shape(2)*tRFiC.shape(3)+a*tRFiC.shape(3)+b] +=  4.0 * Tc.cptr()[j*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+i*Tc.shape(2)*Tc.shape(3)+v0*Tc.shape(3)+a] * Fab.cptr()[v0*Fab.shape(1)+b];
        }//b
      }//a
    }//v0
  }//i
}//j
std::cout << "@tRFiC(i,j,a,b) <<= +4 @Tc(j,i,b,v0) @cF(v0,a)" << std::endl;
// @tRFiC(i,j,a,b) <<= +4 @Tc(j,i,b,v0) @cF(v0,a)
for(auto &&j : irange(0L, _Ndocc)){
  for(auto &&i : irange(0L, _Ndocc)){
    for(auto &&b : irange(0L, _Nvirt)){
      for(auto &&v0 : irange(0L, _Nvirt)){
        for(auto &&a : irange(0L, _Nvirt)){
          tRFiC.cptr()[i*tRFiC.shape(1)*tRFiC.shape(2)*tRFiC.shape(3)+j*tRFiC.shape(2)*tRFiC.shape(3)+a*tRFiC.shape(3)+b] +=  4.0 * Tc.cptr()[j*Tc.shape(1)*Tc.shape(2)*Tc.shape(3)+i*Tc.shape(2)*Tc.shape(3)+b*Tc.shape(3)+v0] * Fab.cptr()[v0*Fab.shape(1)+a];
        }//a
      }//v0
    }//b
  }//i
}//j
 //Code-End-Code-End-Code-End-Code-End-Code-End-Code-End-Code-End
      }
      
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
          R2.PutTensor(i, j, Rij);
        }
      }
      
      return;
    }//End func
  } }
