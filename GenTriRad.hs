
import System.Directory
import Data.List as List
import Data.Maybe as Maybe

import Data.Time
import System.Process
import System.Environment

_header :: String -> String -> String -> String
_header tcutmkn tcutpairs tcutpno =
  "\
   \! UHF DLPNO-CCSD def2-TZVPP def2-TZVPP/C TightSCF Angs \n\
   \ \n\
   \%maxcore 32000 \n\
   \ \n\
   \%mdci \n\
   \  UseQROs      true    \n\
   \  CIType       CCSD    \n\
   \  LevelShift   0.0     \n\
   \  MaxIter      25      \n\
   \  DoSingles    true    \n\
   \  Randomize    false   \n\
   \  TCutPNO      " ++ tcutpno   ++ " \n\
   \  TCutPairs    " ++ tcutpairs ++ " \n\
   \  TCutMKN      " ++ tcutmkn   ++ " \n\
   \  TCutDO       1e-2    \n\
   \  TCutCMO      1e-3    \n\
   \  TCutCPAO     1e-3    \n\
   \end\n\n"
      
_geom :: String -> String
_geom name = "*xyzfile 0 5 " ++ name ++ ".xyz\n"

names =
  [
    "C10.tri",
    "C20.tri", 
    "C30.tri", 
    "C40.tri", 
    "C50.tri", 
    "C60.tri",
    "C70.tri",
    "C80.tri",
    "C90.tri",
    "C100.tri",
    "C120.tri",
    "C150.tri"
  ]

fileName :: String -> String
fileName name = "./" ++ name ++ ".inp"

_headerLoose :: String
_headerLoose = _header "1e-3" "1e-3" "1e-6"

_headerNormal :: String
_headerNormal = _header "1e-3" "1e-4" "3.33e-7"

main = do
  let
    fileNames = map fileName names
    geoms     = map _geom    names
    --contents  = map (\x -> _headerLoose ++ x) geoms -- Loose
    contents  = map (\x -> _headerNormal ++ x) geoms
    nameCont  = zip fileNames contents
  res <- mapM (\x -> writeFile (fst x) (snd x)) nameCont
  putStrLn " Done.\n"

  
