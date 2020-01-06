>// Copyright (C) 2005-2014 by Takeshi Yanai (yanait@gmail.com)
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

#define sqrtHalf (1/sqrt(2.0))

    void FMatrixOrthogonalizer::FormLHS_V0p_V0p(PairList &RActList, PairList &TActList, const orz::DTensor &d1, const orz::DTensor &d2, const orz::DTensor &d3, const orz::DTensor &c6,
                                                DensityPack &dPack, const double &Ecas, const double &h0_energy, const orz::DTensor &casfock, const orz::DTensor &FV,
                                                TensorContainer &T2, TensorContainer &R2){

      orz::DTensor Fcore(casfock(orz::Slice(0, _Ndocc), orz::Slice(0, _Ndocc)).copy());
      orz::DTensor Fact(casfock(orz::Slice(_Ndocc, _Ndocc+_Nact), orz::Slice(_Ndocc, _Ndocc+_Nact)).copy());
      const long Nrho0p = _Pvec["(0')"]->shape(0);
      orz::DTensor C0p1(_Cmat["(0')(1)"]->copy());
      orz::DTensor C0p2(_Cmat["(0')(2)"]->copy());
      R2.initContainer(Nrho0p);
      for(auto &&rho : irange(0L, Nrho0p)){      
        orz::DTensor tmp(_Ndocc,_Nvirt);
        R2.PutTensor(rho,tmp);
      }

      //-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN      
      {
        std::cout << "@S2(P,i,a) <<= +1 @D2(a0,a1,a2,a3) @cF(i,c0) @X(0')(1)(P,a0,a3) @T(0')(c0,a,P0) @X(0')(1)(P0,a1,a2)" << std::endl;
        // @S2(P,i,a) <<= +1 @D2(a0,a1,a2,a3) @cF(i,c0) @X(0')(1)(P,a0,a3) @T(0')(c0,a,P0) @X(0')(1)(P0,a1,a2)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&c0 : irange(0L, _Ndocc)){
                        for(auto &&a : irange(0L, _Nvirt)){
                          S2.cptr()[i*S2.shape(1)+a] += d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * Fcore.cptr()[i*Fcore.shape(1)+c0] * C0p1.cptr()[a0*C0p1.shape(1)*C0p1.shape(2)+a3*C0p1.shape(2)+P] * U2.cptr()[c0*U2.shape(1)+a] * C0p1.cptr()[a1*C0p1.shape(1)*C0p1.shape(2)+a2*C0p1.shape(2)+P0];
                        }//a
                      }//c0
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -2 Ecas @D2(a0,a1,a2,a3) @X(0')(1)(P,a0,a3) @T(0')(i,a,P0) @X(0')(1)(P0,a1,a2)" << std::endl;
        // @S2(P,i,a) <<= -2 Ecas @D2(a0,a1,a2,a3) @X(0')(1)(P,a0,a3) @T(0')(i,a,P0) @X(0')(1)(P0,a1,a2)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] -=  2.0 * Ecas * d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * C0p1.cptr()[a0*C0p1.shape(1)*C0p1.shape(2)+a3*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p1.cptr()[a1*C0p1.shape(1)*C0p1.shape(2)+a2*C0p1.shape(2)+P0];
                      }//a
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -1 @C3(a0,a1,a2,a3) @X(0')(1)(P,a3,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a2,a1)" << std::endl;
        // @S2(P,i,a) <<= -1 @C3(a0,a1,a2,a3) @X(0')(1)(P,a3,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a2,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] -= C3.cptr()[a0*C3.shape(1)*C3.shape(2)*C3.shape(3)+a1*C3.shape(2)*C3.shape(3)+a2*C3.shape(3)+a3] * C0p1.cptr()[a3*C0p1.shape(1)*C0p1.shape(2)+a0*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a1*C0p1.shape(2)+P0];
                      }//a
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(a4,a3) @X(0')(1)(P,a4,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a2,a1)" << std::endl;
        // @S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(a4,a3) @X(0')(1)(P,a4,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a2,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&a4 : irange(0L, _Nact)){
                      for(auto &&i : irange(0L, _Ndocc)){
                        for(auto &&a : irange(0L, _Nvirt)){
                          S2.cptr()[i*S2.shape(1)+a] -= d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * Fact.cptr()[a4*Fact.shape(1)+a3] * C0p1.cptr()[a4*C0p1.shape(1)*C0p1.shape(2)+a0*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a1*C0p1.shape(2)+P0];
                        }//a
                      }//i
                    }//a4
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +4 Ecas @D1(a0,a1) @X(0')(1)(P,a2,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a2,a1)" << std::endl;
        // @S2(P,i,a) <<= +4 Ecas @D1(a0,a1) @X(0')(1)(P,a2,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a2,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&i : irange(0L, _Ndocc)){
                    for(auto &&a : irange(0L, _Nvirt)){
                      S2.cptr()[i*S2.shape(1)+a] +=  4.0 * Ecas * d1.cptr()[a0*d1.shape(1)+a1] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a0*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a1*C0p1.shape(2)+P0];
                    }//a
                  }//i
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -2 @D1(a0,a1) @cF(i,c0) @X(0')(1)(P,a2,a0) @T(0')(c0,a,P0) @X(0')(1)(P0,a2,a1)" << std::endl;
        // @S2(P,i,a) <<= -2 @D1(a0,a1) @cF(i,c0) @X(0')(1)(P,a2,a0) @T(0')(c0,a,P0) @X(0')(1)(P0,a2,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&i : irange(0L, _Ndocc)){
                  for(auto &&c0 : irange(0L, _Ndocc)){
                    for(auto &&a2 : irange(0L, _Nact)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] -=  2.0 * d1.cptr()[a0*d1.shape(1)+a1] * Fcore.cptr()[i*Fcore.shape(1)+c0] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a0*C0p1.shape(2)+P] * U2.cptr()[c0*U2.shape(1)+a] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a1*C0p1.shape(2)+P0];
                      }//a
                    }//a2
                  }//c0
                }//i
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +2 @D2(a0,a1,a2,a3) @cF(a2,a3) @X(0')(1)(P,a4,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a4,a1)" << std::endl;
        // @S2(P,i,a) <<= +2 @D2(a0,a1,a2,a3) @cF(a2,a3) @X(0')(1)(P,a4,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a4,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&a4 : irange(0L, _Nact)){
                      for(auto &&i : irange(0L, _Ndocc)){
                        for(auto &&a : irange(0L, _Nvirt)){
                          S2.cptr()[i*S2.shape(1)+a] +=  2.0 * d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * Fact.cptr()[a2*Fact.shape(1)+a3] * C0p1.cptr()[a4*C0p1.shape(1)*C0p1.shape(2)+a0*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p1.cptr()[a4*C0p1.shape(1)*C0p1.shape(2)+a1*C0p1.shape(2)+P0];
                        }//a
                      }//i
                    }//a4
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(a1,a4) @X(0')(1)(P,a0,a3) @T(0')(i,a,P0) @X(0')(1)(P0,a4,a2)" << std::endl;
        // @S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(a1,a4) @X(0')(1)(P,a0,a3) @T(0')(i,a,P0) @X(0')(1)(P0,a4,a2)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&a4 : irange(0L, _Nact)){
                      for(auto &&i : irange(0L, _Ndocc)){
                        for(auto &&a : irange(0L, _Nvirt)){
                          S2.cptr()[i*S2.shape(1)+a] -= d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * Fact.cptr()[a1*Fact.shape(1)+a4] * C0p1.cptr()[a0*C0p1.shape(1)*C0p1.shape(2)+a3*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p1.cptr()[a4*C0p1.shape(1)*C0p1.shape(2)+a2*C0p1.shape(2)+P0];
                        }//a
                      }//i
                    }//a4
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +2 @D1(a0,a1) @cF(a2,a3) @X(0')(1)(P,a2,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a3,a1)" << std::endl;
        // @S2(P,i,a) <<= +2 @D1(a0,a1) @cF(a2,a3) @X(0')(1)(P,a2,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a3,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] +=  2.0 * d1.cptr()[a0*d1.shape(1)+a1] * Fact.cptr()[a2*Fact.shape(1)+a3] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a0*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p1.cptr()[a3*C0p1.shape(1)*C0p1.shape(2)+a1*C0p1.shape(2)+P0];
                      }//a
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(v0,a) @X(0')(1)(P,a0,a3) @T(0')(i,v0,P0) @X(0')(1)(P0,a1,a2)" << std::endl;
        // @S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(v0,a) @X(0')(1)(P,a0,a3) @T(0')(i,v0,P0) @X(0')(1)(P0,a1,a2)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&v0 : irange(0L, _Nvirt)){
                        for(auto &&a : irange(0L, _Nvirt)){
                          S2.cptr()[i*S2.shape(1)+a] -= d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * FV.cptr()[v0*FV.shape(1)+a] * C0p1.cptr()[a0*C0p1.shape(1)*C0p1.shape(2)+a3*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+v0] * C0p1.cptr()[a1*C0p1.shape(1)*C0p1.shape(2)+a2*C0p1.shape(2)+P0];
                        }//a
                      }//v0
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +2 @D1(a0,a1) @cF(v0,a) @X(0')(1)(P,a2,a0) @T(0')(i,v0,P0) @X(0')(1)(P0,a2,a1)" << std::endl;
        // @S2(P,i,a) <<= +2 @D1(a0,a1) @cF(v0,a) @X(0')(1)(P,a2,a0) @T(0')(i,v0,P0) @X(0')(1)(P0,a2,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&i : irange(0L, _Ndocc)){
                    for(auto &&v0 : irange(0L, _Nvirt)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] +=  2.0 * d1.cptr()[a0*d1.shape(1)+a1] * FV.cptr()[v0*FV.shape(1)+a] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a0*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+v0] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a1*C0p1.shape(2)+P0];
                      }//a
                    }//v0
                  }//i
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +1 E0 @D2(a0,a1,a2,a3) @X(0')(1)(P,a0,a3) @T(0')(i,a,P0) @X(0')(1)(P0,a1,a2)" << std::endl;
        // @S2(P,i,a) <<= +1 E0 @D2(a0,a1,a2,a3) @X(0')(1)(P,a0,a3) @T(0')(i,a,P0) @X(0')(1)(P0,a1,a2)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] += h0_energy * d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * C0p1.cptr()[a0*C0p1.shape(1)*C0p1.shape(2)+a3*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p1.cptr()[a1*C0p1.shape(1)*C0p1.shape(2)+a2*C0p1.shape(2)+P0];
                      }//a
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -2 E0 @D1(a0,a1) @X(0')(1)(P,a2,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a2,a1)" << std::endl;
        // @S2(P,i,a) <<= -2 E0 @D1(a0,a1) @X(0')(1)(P,a2,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a2,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&i : irange(0L, _Ndocc)){
                    for(auto &&a : irange(0L, _Nvirt)){
                      S2.cptr()[i*S2.shape(1)+a] -=  2.0 * h0_energy * d1.cptr()[a0*d1.shape(1)+a1] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a0*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a1*C0p1.shape(2)+P0];
                    }//a
                  }//i
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +1 @D2(a0,a1,a2,a3) @cF(i,c0) @X(0')(1)(P,a0,a1) @T(0')(c0,a,P0) @X(0')(2)(P0,a3,a2)" << std::endl;
        // @S2(P,i,a) <<= +1 @D2(a0,a1,a2,a3) @cF(i,c0) @X(0')(1)(P,a0,a1) @T(0')(c0,a,P0) @X(0')(2)(P0,a3,a2)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&c0 : irange(0L, _Ndocc)){
                        for(auto &&a : irange(0L, _Nvirt)){
                          S2.cptr()[i*S2.shape(1)+a] += d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * Fcore.cptr()[i*Fcore.shape(1)+c0] * C0p1.cptr()[a0*C0p1.shape(1)*C0p1.shape(2)+a1*C0p1.shape(2)+P] * U2.cptr()[c0*U2.shape(1)+a] * C0p2.cptr()[a3*C0p2.shape(1)*C0p2.shape(2)+a2*C0p2.shape(2)+P0];
                        }//a
                      }//c0
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -2 Ecas @D2(a0,a1,a2,a3) @X(0')(1)(P,a0,a1) @T(0')(i,a,P0) @X(0')(2)(P0,a3,a2)" << std::endl;
        // @S2(P,i,a) <<= -2 Ecas @D2(a0,a1,a2,a3) @X(0')(1)(P,a0,a1) @T(0')(i,a,P0) @X(0')(2)(P0,a3,a2)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] -=  2.0 * Ecas * d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * C0p1.cptr()[a0*C0p1.shape(1)*C0p1.shape(2)+a1*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p2.cptr()[a3*C0p2.shape(1)*C0p2.shape(2)+a2*C0p2.shape(2)+P0];
                      }//a
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -1 @C3(a0,a1,a2,a3) @X(0')(1)(P,a3,a2) @T(0')(i,a,P0) @X(0')(2)(P0,a0,a1)" << std::endl;
        // @S2(P,i,a) <<= -1 @C3(a0,a1,a2,a3) @X(0')(1)(P,a3,a2) @T(0')(i,a,P0) @X(0')(2)(P0,a0,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] -= C3.cptr()[a0*C3.shape(1)*C3.shape(2)*C3.shape(3)+a1*C3.shape(2)*C3.shape(3)+a2*C3.shape(3)+a3] * C0p1.cptr()[a3*C0p1.shape(1)*C0p1.shape(2)+a2*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p2.cptr()[a0*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P0];
                      }//a
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(a4,a1) @X(0')(1)(P,a4,a0) @T(0')(i,a,P0) @X(0')(2)(P0,a2,a3)" << std::endl;
        // @S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(a4,a1) @X(0')(1)(P,a4,a0) @T(0')(i,a,P0) @X(0')(2)(P0,a2,a3)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&a4 : irange(0L, _Nact)){
                      for(auto &&i : irange(0L, _Ndocc)){
                        for(auto &&a : irange(0L, _Nvirt)){
                          S2.cptr()[i*S2.shape(1)+a] -= d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * Fact.cptr()[a4*Fact.shape(1)+a1] * C0p1.cptr()[a4*C0p1.shape(1)*C0p1.shape(2)+a0*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p2.cptr()[a2*C0p2.shape(1)*C0p2.shape(2)+a3*C0p2.shape(2)+P0];
                        }//a
                      }//i
                    }//a4
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -2 Ecas @D1(a0,a1) @X(0')(1)(P,a2,a0) @T(0')(i,a,P0) @X(0')(2)(P0,a2,a1)" << std::endl;
        // @S2(P,i,a) <<= -2 Ecas @D1(a0,a1) @X(0')(1)(P,a2,a0) @T(0')(i,a,P0) @X(0')(2)(P0,a2,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&i : irange(0L, _Ndocc)){
                    for(auto &&a : irange(0L, _Nvirt)){
                      S2.cptr()[i*S2.shape(1)+a] -=  2.0 * Ecas * d1.cptr()[a0*d1.shape(1)+a1] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a0*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p2.cptr()[a2*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P0];
                    }//a
                  }//i
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +1 @D1(a0,a1) @cF(i,c0) @X(0')(1)(P,a2,a0) @T(0')(c0,a,P0) @X(0')(2)(P0,a2,a1)" << std::endl;
        // @S2(P,i,a) <<= +1 @D1(a0,a1) @cF(i,c0) @X(0')(1)(P,a2,a0) @T(0')(c0,a,P0) @X(0')(2)(P0,a2,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&i : irange(0L, _Ndocc)){
                  for(auto &&c0 : irange(0L, _Ndocc)){
                    for(auto &&a2 : irange(0L, _Nact)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] += d1.cptr()[a0*d1.shape(1)+a1] * Fcore.cptr()[i*Fcore.shape(1)+c0] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a0*C0p1.shape(2)+P] * U2.cptr()[c0*U2.shape(1)+a] * C0p2.cptr()[a2*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P0];
                      }//a
                    }//a2
                  }//c0
                }//i
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(a2,a3) @X(0')(1)(P,a4,a0) @T(0')(i,a,P0) @X(0')(2)(P0,a4,a1)" << std::endl;
        // @S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(a2,a3) @X(0')(1)(P,a4,a0) @T(0')(i,a,P0) @X(0')(2)(P0,a4,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&a4 : irange(0L, _Nact)){
                      for(auto &&i : irange(0L, _Ndocc)){
                        for(auto &&a : irange(0L, _Nvirt)){
                          S2.cptr()[i*S2.shape(1)+a] -= d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * Fact.cptr()[a2*Fact.shape(1)+a3] * C0p1.cptr()[a4*C0p1.shape(1)*C0p1.shape(2)+a0*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p2.cptr()[a4*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P0];
                        }//a
                      }//i
                    }//a4
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(v0,a) @X(0')(1)(P,a0,a1) @T(0')(i,v0,P0) @X(0')(2)(P0,a3,a2)" << std::endl;
        // @S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(v0,a) @X(0')(1)(P,a0,a1) @T(0')(i,v0,P0) @X(0')(2)(P0,a3,a2)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&v0 : irange(0L, _Nvirt)){
                        for(auto &&a : irange(0L, _Nvirt)){
                          S2.cptr()[i*S2.shape(1)+a] -= d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * FV.cptr()[v0*FV.shape(1)+a] * C0p1.cptr()[a0*C0p1.shape(1)*C0p1.shape(2)+a1*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+v0] * C0p2.cptr()[a3*C0p2.shape(1)*C0p2.shape(2)+a2*C0p2.shape(2)+P0];
                        }//a
                      }//v0
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -1 @D1(a0,a1) @cF(v0,a) @X(0')(1)(P,a2,a0) @T(0')(i,v0,P0) @X(0')(2)(P0,a2,a1)" << std::endl;
        // @S2(P,i,a) <<= -1 @D1(a0,a1) @cF(v0,a) @X(0')(1)(P,a2,a0) @T(0')(i,v0,P0) @X(0')(2)(P0,a2,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&i : irange(0L, _Ndocc)){
                    for(auto &&v0 : irange(0L, _Nvirt)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] -= d1.cptr()[a0*d1.shape(1)+a1] * FV.cptr()[v0*FV.shape(1)+a] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a0*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+v0] * C0p2.cptr()[a2*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P0];
                      }//a
                    }//v0
                  }//i
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(a3,a4) @X(0')(1)(P,a0,a1) @T(0')(i,a,P0) @X(0')(2)(P0,a4,a2)" << std::endl;
        // @S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(a3,a4) @X(0')(1)(P,a0,a1) @T(0')(i,a,P0) @X(0')(2)(P0,a4,a2)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&a4 : irange(0L, _Nact)){
                      for(auto &&i : irange(0L, _Ndocc)){
                        for(auto &&a : irange(0L, _Nvirt)){
                          S2.cptr()[i*S2.shape(1)+a] -= d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * Fact.cptr()[a3*Fact.shape(1)+a4] * C0p1.cptr()[a0*C0p1.shape(1)*C0p1.shape(2)+a1*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p2.cptr()[a4*C0p2.shape(1)*C0p2.shape(2)+a2*C0p2.shape(2)+P0];
                        }//a
                      }//i
                    }//a4
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -1 @D1(a0,a1) @cF(a2,a3) @X(0')(1)(P,a2,a0) @T(0')(i,a,P0) @X(0')(2)(P0,a3,a1)" << std::endl;
        // @S2(P,i,a) <<= -1 @D1(a0,a1) @cF(a2,a3) @X(0')(1)(P,a2,a0) @T(0')(i,a,P0) @X(0')(2)(P0,a3,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] -= d1.cptr()[a0*d1.shape(1)+a1] * Fact.cptr()[a2*Fact.shape(1)+a3] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a0*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p2.cptr()[a3*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P0];
                      }//a
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +1 E0 @D2(a0,a1,a2,a3) @X(0')(1)(P,a0,a1) @T(0')(i,a,P0) @X(0')(2)(P0,a3,a2)" << std::endl;
        // @S2(P,i,a) <<= +1 E0 @D2(a0,a1,a2,a3) @X(0')(1)(P,a0,a1) @T(0')(i,a,P0) @X(0')(2)(P0,a3,a2)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] += h0_energy * d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * C0p1.cptr()[a0*C0p1.shape(1)*C0p1.shape(2)+a1*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p2.cptr()[a3*C0p2.shape(1)*C0p2.shape(2)+a2*C0p2.shape(2)+P0];
                      }//a
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +1 E0 @D1(a0,a1) @X(0')(1)(P,a2,a0) @T(0')(i,a,P0) @X(0')(2)(P0,a2,a1)" << std::endl;
        // @S2(P,i,a) <<= +1 E0 @D1(a0,a1) @X(0')(1)(P,a2,a0) @T(0')(i,a,P0) @X(0')(2)(P0,a2,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&i : irange(0L, _Ndocc)){
                    for(auto &&a : irange(0L, _Nvirt)){
                      S2.cptr()[i*S2.shape(1)+a] += h0_energy * d1.cptr()[a0*d1.shape(1)+a1] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a0*C0p1.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p2.cptr()[a2*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P0];
                    }//a
                  }//i
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -2 Ecas @D2(a0,a1,a2,a3) @X(0')(2)(P,a0,a1) @T(0')(i,a,P0) @X(0')(1)(P0,a3,a2)" << std::endl;
        // @S2(P,i,a) <<= -2 Ecas @D2(a0,a1,a2,a3) @X(0')(2)(P,a0,a1) @T(0')(i,a,P0) @X(0')(1)(P0,a3,a2)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] -=  2.0 * Ecas * d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * C0p2.cptr()[a0*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p1.cptr()[a3*C0p1.shape(1)*C0p1.shape(2)+a2*C0p1.shape(2)+P0];
                      }//a
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +1 @D2(a0,a1,a2,a3) @cF(i,c0) @X(0')(2)(P,a0,a1) @T(0')(c0,a,P0) @X(0')(1)(P0,a3,a2)" << std::endl;
        // @S2(P,i,a) <<= +1 @D2(a0,a1,a2,a3) @cF(i,c0) @X(0')(2)(P,a0,a1) @T(0')(c0,a,P0) @X(0')(1)(P0,a3,a2)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&c0 : irange(0L, _Ndocc)){
                        for(auto &&a : irange(0L, _Nvirt)){
                          S2.cptr()[i*S2.shape(1)+a] += d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * Fcore.cptr()[i*Fcore.shape(1)+c0] * C0p2.cptr()[a0*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P] * U2.cptr()[c0*U2.shape(1)+a] * C0p1.cptr()[a3*C0p1.shape(1)*C0p1.shape(2)+a2*C0p1.shape(2)+P0];
                        }//a
                      }//c0
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -1 @C3(a0,a1,a2,a3) @X(0')(2)(P,a3,a2) @T(0')(i,a,P0) @X(0')(1)(P0,a0,a1)" << std::endl;
        // @S2(P,i,a) <<= -1 @C3(a0,a1,a2,a3) @X(0')(2)(P,a3,a2) @T(0')(i,a,P0) @X(0')(1)(P0,a0,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] -= C3.cptr()[a0*C3.shape(1)*C3.shape(2)*C3.shape(3)+a1*C3.shape(2)*C3.shape(3)+a2*C3.shape(3)+a3] * C0p2.cptr()[a3*C0p2.shape(1)*C0p2.shape(2)+a2*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p1.cptr()[a0*C0p1.shape(1)*C0p1.shape(2)+a1*C0p1.shape(2)+P0];
                      }//a
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(a4,a1) @X(0')(2)(P,a4,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a2,a3)" << std::endl;
        // @S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(a4,a1) @X(0')(2)(P,a4,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a2,a3)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&a4 : irange(0L, _Nact)){
                      for(auto &&i : irange(0L, _Ndocc)){
                        for(auto &&a : irange(0L, _Nvirt)){
                          S2.cptr()[i*S2.shape(1)+a] -= d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * Fact.cptr()[a4*Fact.shape(1)+a1] * C0p2.cptr()[a4*C0p2.shape(1)*C0p2.shape(2)+a0*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a3*C0p1.shape(2)+P0];
                        }//a
                      }//i
                    }//a4
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -2 Ecas @D1(a0,a1) @X(0')(2)(P,a2,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a2,a1)" << std::endl;
        // @S2(P,i,a) <<= -2 Ecas @D1(a0,a1) @X(0')(2)(P,a2,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a2,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&i : irange(0L, _Ndocc)){
                    for(auto &&a : irange(0L, _Nvirt)){
                      S2.cptr()[i*S2.shape(1)+a] -=  2.0 * Ecas * d1.cptr()[a0*d1.shape(1)+a1] * C0p2.cptr()[a2*C0p2.shape(1)*C0p2.shape(2)+a0*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a1*C0p1.shape(2)+P0];
                    }//a
                  }//i
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +1 @D1(a0,a1) @cF(i,c0) @X(0')(2)(P,a2,a0) @T(0')(c0,a,P0) @X(0')(1)(P0,a2,a1)" << std::endl;
        // @S2(P,i,a) <<= +1 @D1(a0,a1) @cF(i,c0) @X(0')(2)(P,a2,a0) @T(0')(c0,a,P0) @X(0')(1)(P0,a2,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&i : irange(0L, _Ndocc)){
                  for(auto &&c0 : irange(0L, _Ndocc)){
                    for(auto &&a2 : irange(0L, _Nact)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] += d1.cptr()[a0*d1.shape(1)+a1] * Fcore.cptr()[i*Fcore.shape(1)+c0] * C0p2.cptr()[a2*C0p2.shape(1)*C0p2.shape(2)+a0*C0p2.shape(2)+P] * U2.cptr()[c0*U2.shape(1)+a] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a1*C0p1.shape(2)+P0];
                      }//a
                    }//a2
                  }//c0
                }//i
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(a2,a3) @X(0')(2)(P,a4,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a4,a1)" << std::endl;
        // @S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(a2,a3) @X(0')(2)(P,a4,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a4,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&a4 : irange(0L, _Nact)){
                      for(auto &&i : irange(0L, _Ndocc)){
                        for(auto &&a : irange(0L, _Nvirt)){
                          S2.cptr()[i*S2.shape(1)+a] -= d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * Fact.cptr()[a2*Fact.shape(1)+a3] * C0p2.cptr()[a4*C0p2.shape(1)*C0p2.shape(2)+a0*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p1.cptr()[a4*C0p1.shape(1)*C0p1.shape(2)+a1*C0p1.shape(2)+P0];
                        }//a
                      }//i
                    }//a4
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(a3,a4) @X(0')(2)(P,a0,a1) @T(0')(i,a,P0) @X(0')(1)(P0,a4,a2)" << std::endl;
        // @S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(a3,a4) @X(0')(2)(P,a0,a1) @T(0')(i,a,P0) @X(0')(1)(P0,a4,a2)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&a4 : irange(0L, _Nact)){
                      for(auto &&i : irange(0L, _Ndocc)){
                        for(auto &&a : irange(0L, _Nvirt)){
                          S2.cptr()[i*S2.shape(1)+a] -= d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * Fact.cptr()[a3*Fact.shape(1)+a4] * C0p2.cptr()[a0*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p1.cptr()[a4*C0p1.shape(1)*C0p1.shape(2)+a2*C0p1.shape(2)+P0];
                        }//a
                      }//i
                    }//a4
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -1 @D1(a0,a1) @cF(a2,a3) @X(0')(2)(P,a2,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a3,a1)" << std::endl;
        // @S2(P,i,a) <<= -1 @D1(a0,a1) @cF(a2,a3) @X(0')(2)(P,a2,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a3,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] -= d1.cptr()[a0*d1.shape(1)+a1] * Fact.cptr()[a2*Fact.shape(1)+a3] * C0p2.cptr()[a2*C0p2.shape(1)*C0p2.shape(2)+a0*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p1.cptr()[a3*C0p1.shape(1)*C0p1.shape(2)+a1*C0p1.shape(2)+P0];
                      }//a
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(v0,a) @X(0')(2)(P,a0,a1) @T(0')(i,v0,P0) @X(0')(1)(P0,a3,a2)" << std::endl;
        // @S2(P,i,a) <<= -1 @D2(a0,a1,a2,a3) @cF(v0,a) @X(0')(2)(P,a0,a1) @T(0')(i,v0,P0) @X(0')(1)(P0,a3,a2)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&v0 : irange(0L, _Nvirt)){
                        for(auto &&a : irange(0L, _Nvirt)){
                          S2.cptr()[i*S2.shape(1)+a] -= d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * FV.cptr()[v0*FV.shape(1)+a] * C0p2.cptr()[a0*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+v0] * C0p1.cptr()[a3*C0p1.shape(1)*C0p1.shape(2)+a2*C0p1.shape(2)+P0];
                        }//a
                      }//v0
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -1 @D1(a0,a1) @cF(v0,a) @X(0')(2)(P,a2,a0) @T(0')(i,v0,P0) @X(0')(1)(P0,a2,a1)" << std::endl;
        // @S2(P,i,a) <<= -1 @D1(a0,a1) @cF(v0,a) @X(0')(2)(P,a2,a0) @T(0')(i,v0,P0) @X(0')(1)(P0,a2,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&i : irange(0L, _Ndocc)){
                    for(auto &&v0 : irange(0L, _Nvirt)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] -= d1.cptr()[a0*d1.shape(1)+a1] * FV.cptr()[v0*FV.shape(1)+a] * C0p2.cptr()[a2*C0p2.shape(1)*C0p2.shape(2)+a0*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+v0] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a1*C0p1.shape(2)+P0];
                      }//a
                    }//v0
                  }//i
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +1 E0 @D2(a0,a1,a2,a3) @X(0')(2)(P,a0,a1) @T(0')(i,a,P0) @X(0')(1)(P0,a3,a2)" << std::endl;
        // @S2(P,i,a) <<= +1 E0 @D2(a0,a1,a2,a3) @X(0')(2)(P,a0,a1) @T(0')(i,a,P0) @X(0')(1)(P0,a3,a2)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] += h0_energy * d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * C0p2.cptr()[a0*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p1.cptr()[a3*C0p1.shape(1)*C0p1.shape(2)+a2*C0p1.shape(2)+P0];
                      }//a
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +1 E0 @D1(a0,a1) @X(0')(2)(P,a2,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a2,a1)" << std::endl;
        // @S2(P,i,a) <<= +1 E0 @D1(a0,a1) @X(0')(2)(P,a2,a0) @T(0')(i,a,P0) @X(0')(1)(P0,a2,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&i : irange(0L, _Ndocc)){
                    for(auto &&a : irange(0L, _Nvirt)){
                      S2.cptr()[i*S2.shape(1)+a] += h0_energy * d1.cptr()[a0*d1.shape(1)+a1] * C0p2.cptr()[a2*C0p2.shape(1)*C0p2.shape(2)+a0*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p1.cptr()[a2*C0p1.shape(1)*C0p1.shape(2)+a1*C0p1.shape(2)+P0];
                    }//a
                  }//i
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +4 Ecas @D2(a0,a1,a2,a3) @X(0')(2)(P,a0,a1) @T(0')(i,a,P0) @X(0')(2)(P0,a3,a2)" << std::endl;
        // @S2(P,i,a) <<= +4 Ecas @D2(a0,a1,a2,a3) @X(0')(2)(P,a0,a1) @T(0')(i,a,P0) @X(0')(2)(P0,a3,a2)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] +=  4.0 * Ecas * d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * C0p2.cptr()[a0*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p2.cptr()[a3*C0p2.shape(1)*C0p2.shape(2)+a2*C0p2.shape(2)+P0];
                      }//a
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -2 @D2(a0,a1,a2,a3) @cF(i,c0) @X(0')(2)(P,a0,a1) @T(0')(c0,a,P0) @X(0')(2)(P0,a3,a2)" << std::endl;
        // @S2(P,i,a) <<= -2 @D2(a0,a1,a2,a3) @cF(i,c0) @X(0')(2)(P,a0,a1) @T(0')(c0,a,P0) @X(0')(2)(P0,a3,a2)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&c0 : irange(0L, _Ndocc)){
                        for(auto &&a : irange(0L, _Nvirt)){
                          S2.cptr()[i*S2.shape(1)+a] -=  2.0 * d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * Fcore.cptr()[i*Fcore.shape(1)+c0] * C0p2.cptr()[a0*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P] * U2.cptr()[c0*U2.shape(1)+a] * C0p2.cptr()[a3*C0p2.shape(1)*C0p2.shape(2)+a2*C0p2.shape(2)+P0];
                        }//a
                      }//c0
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +2 @C3(a0,a1,a2,a3) @X(0')(2)(P,a3,a2) @T(0')(i,a,P0) @X(0')(2)(P0,a0,a1)" << std::endl;
        // @S2(P,i,a) <<= +2 @C3(a0,a1,a2,a3) @X(0')(2)(P,a3,a2) @T(0')(i,a,P0) @X(0')(2)(P0,a0,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] +=  2.0 * C3.cptr()[a0*C3.shape(1)*C3.shape(2)*C3.shape(3)+a1*C3.shape(2)*C3.shape(3)+a2*C3.shape(3)+a3] * C0p2.cptr()[a3*C0p2.shape(1)*C0p2.shape(2)+a2*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p2.cptr()[a0*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P0];
                      }//a
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +2 @D2(a0,a1,a2,a3) @cF(a4,a1) @X(0')(2)(P,a4,a0) @T(0')(i,a,P0) @X(0')(2)(P0,a2,a3)" << std::endl;
        // @S2(P,i,a) <<= +2 @D2(a0,a1,a2,a3) @cF(a4,a1) @X(0')(2)(P,a4,a0) @T(0')(i,a,P0) @X(0')(2)(P0,a2,a3)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&a4 : irange(0L, _Nact)){
                      for(auto &&i : irange(0L, _Ndocc)){
                        for(auto &&a : irange(0L, _Nvirt)){
                          S2.cptr()[i*S2.shape(1)+a] +=  2.0 * d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * Fact.cptr()[a4*Fact.shape(1)+a1] * C0p2.cptr()[a4*C0p2.shape(1)*C0p2.shape(2)+a0*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p2.cptr()[a2*C0p2.shape(1)*C0p2.shape(2)+a3*C0p2.shape(2)+P0];
                        }//a
                      }//i
                    }//a4
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +4 Ecas @D1(a0,a1) @X(0')(2)(P,a2,a0) @T(0')(i,a,P0) @X(0')(2)(P0,a2,a1)" << std::endl;
        // @S2(P,i,a) <<= +4 Ecas @D1(a0,a1) @X(0')(2)(P,a2,a0) @T(0')(i,a,P0) @X(0')(2)(P0,a2,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&i : irange(0L, _Ndocc)){
                    for(auto &&a : irange(0L, _Nvirt)){
                      S2.cptr()[i*S2.shape(1)+a] +=  4.0 * Ecas * d1.cptr()[a0*d1.shape(1)+a1] * C0p2.cptr()[a2*C0p2.shape(1)*C0p2.shape(2)+a0*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p2.cptr()[a2*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P0];
                    }//a
                  }//i
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -2 @D1(a0,a1) @cF(i,c0) @X(0')(2)(P,a2,a0) @T(0')(c0,a,P0) @X(0')(2)(P0,a2,a1)" << std::endl;
        // @S2(P,i,a) <<= -2 @D1(a0,a1) @cF(i,c0) @X(0')(2)(P,a2,a0) @T(0')(c0,a,P0) @X(0')(2)(P0,a2,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&i : irange(0L, _Ndocc)){
                  for(auto &&c0 : irange(0L, _Ndocc)){
                    for(auto &&a2 : irange(0L, _Nact)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] -=  2.0 * d1.cptr()[a0*d1.shape(1)+a1] * Fcore.cptr()[i*Fcore.shape(1)+c0] * C0p2.cptr()[a2*C0p2.shape(1)*C0p2.shape(2)+a0*C0p2.shape(2)+P] * U2.cptr()[c0*U2.shape(1)+a] * C0p2.cptr()[a2*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P0];
                      }//a
                    }//a2
                  }//c0
                }//i
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +2 @D2(a0,a1,a2,a3) @cF(a2,a3) @X(0')(2)(P,a4,a0) @T(0')(i,a,P0) @X(0')(2)(P0,a4,a1)" << std::endl;
        // @S2(P,i,a) <<= +2 @D2(a0,a1,a2,a3) @cF(a2,a3) @X(0')(2)(P,a4,a0) @T(0')(i,a,P0) @X(0')(2)(P0,a4,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&a4 : irange(0L, _Nact)){
                      for(auto &&i : irange(0L, _Ndocc)){
                        for(auto &&a : irange(0L, _Nvirt)){
                          S2.cptr()[i*S2.shape(1)+a] +=  2.0 * d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * Fact.cptr()[a2*Fact.shape(1)+a3] * C0p2.cptr()[a4*C0p2.shape(1)*C0p2.shape(2)+a0*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p2.cptr()[a4*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P0];
                        }//a
                      }//i
                    }//a4
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +2 @D2(a0,a1,a2,a3) @cF(v0,a) @X(0')(2)(P,a0,a1) @T(0')(i,v0,P0) @X(0')(2)(P0,a3,a2)" << std::endl;
        // @S2(P,i,a) <<= +2 @D2(a0,a1,a2,a3) @cF(v0,a) @X(0')(2)(P,a0,a1) @T(0')(i,v0,P0) @X(0')(2)(P0,a3,a2)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&v0 : irange(0L, _Nvirt)){
                        for(auto &&a : irange(0L, _Nvirt)){
                          S2.cptr()[i*S2.shape(1)+a] +=  2.0 * d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * FV.cptr()[v0*FV.shape(1)+a] * C0p2.cptr()[a0*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+v0] * C0p2.cptr()[a3*C0p2.shape(1)*C0p2.shape(2)+a2*C0p2.shape(2)+P0];
                        }//a
                      }//v0
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +2 @D1(a0,a1) @cF(v0,a) @X(0')(2)(P,a2,a0) @T(0')(i,v0,P0) @X(0')(2)(P0,a2,a1)" << std::endl;
        // @S2(P,i,a) <<= +2 @D1(a0,a1) @cF(v0,a) @X(0')(2)(P,a2,a0) @T(0')(i,v0,P0) @X(0')(2)(P0,a2,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&i : irange(0L, _Ndocc)){
                    for(auto &&v0 : irange(0L, _Nvirt)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] +=  2.0 * d1.cptr()[a0*d1.shape(1)+a1] * FV.cptr()[v0*FV.shape(1)+a] * C0p2.cptr()[a2*C0p2.shape(1)*C0p2.shape(2)+a0*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+v0] * C0p2.cptr()[a2*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P0];
                      }//a
                    }//v0
                  }//i
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +2 @D2(a0,a1,a2,a3) @cF(a3,a4) @X(0')(2)(P,a0,a1) @T(0')(i,a,P0) @X(0')(2)(P0,a4,a2)" << std::endl;
        // @S2(P,i,a) <<= +2 @D2(a0,a1,a2,a3) @cF(a3,a4) @X(0')(2)(P,a0,a1) @T(0')(i,a,P0) @X(0')(2)(P0,a4,a2)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&a4 : irange(0L, _Nact)){
                      for(auto &&i : irange(0L, _Ndocc)){
                        for(auto &&a : irange(0L, _Nvirt)){
                          S2.cptr()[i*S2.shape(1)+a] +=  2.0 * d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * Fact.cptr()[a3*Fact.shape(1)+a4] * C0p2.cptr()[a0*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p2.cptr()[a4*C0p2.shape(1)*C0p2.shape(2)+a2*C0p2.shape(2)+P0];
                        }//a
                      }//i
                    }//a4
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= +2 @D1(a0,a1) @cF(a2,a3) @X(0')(2)(P,a2,a0) @T(0')(i,a,P0) @X(0')(2)(P0,a3,a1)" << std::endl;
        // @S2(P,i,a) <<= +2 @D1(a0,a1) @cF(a2,a3) @X(0')(2)(P,a2,a0) @T(0')(i,a,P0) @X(0')(2)(P0,a3,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] +=  2.0 * d1.cptr()[a0*d1.shape(1)+a1] * Fact.cptr()[a2*Fact.shape(1)+a3] * C0p2.cptr()[a2*C0p2.shape(1)*C0p2.shape(2)+a0*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p2.cptr()[a3*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P0];
                      }//a
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -2 E0 @D2(a0,a1,a2,a3) @X(0')(2)(P,a0,a1) @T(0')(i,a,P0) @X(0')(2)(P0,a3,a2)" << std::endl;
        // @S2(P,i,a) <<= -2 E0 @D2(a0,a1,a2,a3) @X(0')(2)(P,a0,a1) @T(0')(i,a,P0) @X(0')(2)(P0,a3,a2)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&a3 : irange(0L, _Nact)){
                    for(auto &&i : irange(0L, _Ndocc)){
                      for(auto &&a : irange(0L, _Nvirt)){
                        S2.cptr()[i*S2.shape(1)+a] -=  2.0 * h0_energy * d2.cptr()[a0*d2.shape(1)*d2.shape(2)*d2.shape(3)+a1*d2.shape(2)*d2.shape(3)+a2*d2.shape(3)+a3] * C0p2.cptr()[a0*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p2.cptr()[a3*C0p2.shape(1)*C0p2.shape(2)+a2*C0p2.shape(2)+P0];
                      }//a
                    }//i
                  }//a3
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P
        std::cout << "@S2(P,i,a) <<= -2 E0 @D1(a0,a1) @X(0')(2)(P,a2,a0) @T(0')(i,a,P0) @X(0')(2)(P0,a2,a1)" << std::endl;
        // @S2(P,i,a) <<= -2 E0 @D1(a0,a1) @X(0')(2)(P,a2,a0) @T(0')(i,a,P0) @X(0')(2)(P0,a2,a1)
        for(auto &&P : irange(0L, Nrho0p)){
          const long I_P = RActList.GetIndex(P);
          if (I_P < 0) continue;
          orz::DTensor S2 = R2.GetTensor(P);
          for(auto &&P0 : irange(0L, Nrho0p)){
            const long I_P0 = TActList.GetIndex(P0);
            if (I_P0 < 0) continue;
            orz::DTensor U2 = T2.GetTensor(P0);
            for(auto &&a0 : irange(0L, _Nact)){
              for(auto &&a1 : irange(0L, _Nact)){
                for(auto &&a2 : irange(0L, _Nact)){
                  for(auto &&i : irange(0L, _Ndocc)){
                    for(auto &&a : irange(0L, _Nvirt)){
                      S2.cptr()[i*S2.shape(1)+a] -=  2.0 * h0_energy * d1.cptr()[a0*d1.shape(1)+a1] * C0p2.cptr()[a2*C0p2.shape(1)*C0p2.shape(2)+a0*C0p2.shape(2)+P] * U2.cptr()[i*U2.shape(1)+a] * C0p2.cptr()[a2*C0p2.shape(1)*C0p2.shape(2)+a1*C0p2.shape(2)+P0];
                    }//a
                  }//i
                }//a2
              }//a1
            }//a0
          }//P0
          R2.PutTensor(P, S2);
        }//P        
      }
      //-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN-GEN
      
      return;
    }

  } }
