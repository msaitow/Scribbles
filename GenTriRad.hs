
import System.Directory
import Data.List as List
import Data.Maybe as Maybe

import Data.Time
import System.Process
import System.Environment

_headerPNO :: String -> String
_headerPNO val =
  "\
   \! UHF DLPNO-CCSD def2-TZVPP def2-TZVPP/C ExtremeSCF NoAutoStart Angs Pal4 \n\
   \ \n\
   \%maxcore 10000 \n\
   \ \n\
   \%mdci \n\
   \  UseQROs      true   \n\
   \  CIType       CCSD   \n\
   \  KCOpt       -20    \n\
   \  MaxIter      50     \n\
   \  Randomize    false  \n\
   \  TCutPairs    0.0    \m\
   \  TCutPNO      "  ++ val ++ "\n\"
   \end\n\n"

_headerMKN :: String -> String
_headerMKN val =
  "\
   \! UHF DLPNO-CCSD def2-TZVPP def2-TZVPP/C ExtremeSCF NoAutoStart Angs Pal4 \n\
   \ \n\
   \%maxcore 10000 \n\
   \ \n\
   \%mdci \n\
   \  UseQROs      true   \n\
   \  CIType       CCSD   \n\
   \  KCOpt       -20    \n\
   \  MaxIter      50     \n\
   \  Randomize    false  \n\
   \  TCutMKN      "  ++ val ++ "\n\"
   \end\n\n"

_headerPairs :: String -> String
_headerPairs val =
  "\
   \! UHF DLPNO-CCSD def2-TZVPP def2-TZVPP/C ExtremeSCF NoAutoStart Angs Pal4 \n\
   \ \n\
   \%maxcore 10000 \n\
   \ \n\
   \%mdci \n\
   \  UseQROs      true   \n\
   \  CIType       CCSD   \n\
   \  KCOpt       -20    \n\
   \  MaxIter      50     \n\
   \  Randomize    false  \n\
   \  TCutPairs  "  ++ val ++ "\n\"
   \end\n\n"

_headerDO :: String -> String
_headerDO val =
  "\
   \! UHF DLPNO-CCSD def2-TZVPP def2-TZVPP/C ExtremeSCF NoAutoStart Angs Pal4 \n\
   \ \n\
   \%maxcore 10000 \n\
   \ \n\
   \%mdci \n\
   \  UseQROs      true   \n\
   \  CIType       CCSD   \n\
   \  KCOpt       -20    \n\
   \  MaxIter      50     \n\
   \  Randomize    false  \n\
   \  TCutDO      "  ++ val ++ "\n\"
   \end\n\n"
                                  
_geom =
  "\
   \*xyzfile 0 2 trityl.opt.new.xyz \n\
   \\n\n"


mkns =
  [
    "1e-1" ,
    "1e-2" ,
    "1e-3" ,
    "1e-4" ,
    "1e-5" ,
    "1e-6" ,
    "1e-7" ,
    "1e-8"
    ]

pairs =
  [
    "1e-1" ,
    "1e-2" ,
    "1e-3" ,
    "1e-4" ,
    "1e-5" ,
    "1e-6" ,
    "1e-7" ,
    "1e-8"
    ]

pnos =
  [
    "1e-1" ,
    "1e-2" ,
    "1e-3" ,
    "1e-4" ,
    "1e-5" ,
    "1e-6" ,
    "1e-7" ,
    "1e-8" ,
    "1e-9"
    ]

dos =
  [
    "1e-1" ,
    "1e-2" ,
    "1e-3" ,
    "1e-4" ,
    "1e-5" ,
    "1e-6" ,
    "1e-7" ,
    "1e-8"
    ]


_genInput :: String -> String
_genInput tcut = (_headerPNO tcut) ++ _geom
--_genInput tcut = (_headerPairs tcut) ++ _geom
--_genInput tcut = (_headerDO tcut) ++ _geom
--_genInput tcut = (_headerMKN tcut) ++ _geom

fileName :: String -> String
fileName name = "./tri." ++ name ++ ".inp"

main = do
  let
    names     = pnos -- change!
    fileNames = map fileName names
    contents  = map (\x -> _genInput x) names
    nameCont  = zip fileNames contents
  res <- mapM (\x -> writeFile (fst x) (snd x)) nameCont
  putStrLn " Done.\n"

  
