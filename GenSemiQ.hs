
import System.Directory
import Data.List as List
import Data.Maybe as Maybe

import Data.Time
import System.Process
import System.Environment

_header :: String -> String
_header tcutpno =
  "\
   \! UHF DLPNO-CCSD def2-TZVPP def2-TZVPP/C ExtremeSCF NoAutoStart Angs \n\
   \ \n\
   \%maxcore 32000 \n\
   \ \n\
   \%mdci \n\
   \  UseQROs      true   \n\
   \  CIType       CCSD   \n\
   \  LevelShift   0.0    \n\
   \  MaxIter      50     \n\
   \  DoSingles    true   \n\
   \  Randomize    false  \n\
   \  TCutPNO      " ++ tcutpno   ++ " \n\
   \  TCutPairs    0.0    \n\
   \  TCutMKN      0.0    \n\
   \  TCutDO       0.0    \n\
   \  TCutCMO      0.0    \n\
   \  TCutCPAO     0.0    \n\
   \end\n\n"
      
_geom =
  "\
   \* xyz 0 2                                 \n\
   \  O      1.087086   -0.265839   -0.081996 \n\
   \  O     -1.407346    1.153453   -4.751394 \n\
   \  C      0.480715    0.127097   -1.226734 \n\
   \  C     -0.101543   -0.887934   -2.039620 \n\
   \  C      0.414306    1.492582   -1.623906 \n\
   \  C     -0.734296   -0.549432   -3.220844 \n\
   \  C     -0.217306    1.844122   -2.804881 \n\
   \  C     -0.832448    0.840535   -3.679032 \n\
   \  H     -0.029213   -1.931594   -1.696777 \n\
   \  H      0.871936    2.263549   -0.979470 \n\
   \  H     -1.190903   -1.315665   -3.865452 \n\
   \  H     -0.280357    2.893280   -3.131828 \n\
   \  H      1.440470    0.519476    0.380576 \n\
   \*\n\n"


names =
  [
    "1e-3"    ,
    "5e-3"    ,    
    "1e-4"    ,
    "5e-4"    ,    
    "1e-5"    ,
    "5e-5"    ,    
    "1e-6"    ,
    "5e-6"    ,    
    "1e-7"    ,
    "3.33e-7" ,    
    "1e-8"    ,
    "5e-8"    ,    
    "1e-9"    ,
    "5e-9"    ,    
    "1e-10"   ,
    "5e-10"   ,    
    "1e-11"   ,
    "5e-11"   ,    
    "1e-12"   ,
    "5e-12"   
  ]

_genInput :: String -> String
_genInput tcutpno = (_header tcutpno) ++ _geom

fileName :: String -> String
fileName name = "./semiquinone." ++ name ++ ".inp"

main = do
  let
    fileNames = map fileName names
    contents  = map (\x -> _genInput x) names
    nameCont  = zip fileNames contents
  res <- mapM (\x -> writeFile (fst x) (snd x)) nameCont
  putStrLn " Done.\n"

  
