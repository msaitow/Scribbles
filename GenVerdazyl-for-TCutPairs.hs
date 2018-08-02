
import System.Directory
import Data.List as List
import Data.Maybe as Maybe

import Data.Time
import System.Process
import System.Environment

_header :: String -> String
_header tcutpno =
  "\
   \! ROHF DLPNO-CCSD EPR-II autoaux TightSCF NoFrozenCore \n\
   \ \n\
   \%maxcore 16000 \n\
   \ \n\
   \%mdci \n\
   \  Randomize    false  \n\
   \  Density      Unrelaxed \n\
   \  TCutPairs    " ++ tcutpno   ++ " \n\
   \  TCutPNO      1e-10  \n\
   \  TCutMKN      0.0    \n\
   \  TCutDO       0.0    \n\
   \  TCutCMO      0.0    \n\
   \  TCutCPAO     0.0    \n\
   \end\n\n"
      
_geom =
  "\
   \* xyz 0 2                                 \n\
   \  C     -3.896933   -6.090437   -8.433495 \n\
   \  N     -3.764910   -4.995498   -7.655036 \n\
   \  N     -4.190231   -5.128386   -6.383405 \n\
   \  C     -4.283800   -6.486511   -5.854181 \n\
   \  N     -5.024894   -7.272058   -6.837544 \n\
   \  N     -4.620967   -7.185773   -8.120054 \n\
   \  H     -4.819969   -6.482341   -4.893003 \n\
   \  H     -3.260574   -6.912451   -5.718494 \n\
   \  C     -3.306537   -6.032618   -9.798314 \n\
   \  C     -3.468264   -7.111895  -10.696832 \n\
   \  C     -2.570086   -4.900682  -10.215943 \n\
   \  H     -4.032017   -7.995660  -10.365910 \n\
   \  H     -2.435216   -4.066746   -9.511857 \n\
   \  C     -2.916726   -7.053295  -11.985270 \n\
   \  C     -2.019894   -4.846553  -11.505168 \n\
   \  H     -3.053846   -7.901366  -12.674808 \n\
   \  H     -1.449925   -3.957038  -11.816607 \n\
   \  C     -2.191168   -5.920791  -12.396848 \n\
   \  H     -1.756877   -5.877388  -13.407672 \n\
   \  C     -4.450638   -3.974300   -5.620578 \n\
   \  C     -4.583839   -2.724146   -6.277270 \n\
   \  C     -4.571093   -4.026517   -4.208749 \n\
   \  H     -4.481101   -2.693068   -7.369995 \n\
   \  H     -4.430547   -4.969384   -3.661175 \n\
   \  C     -4.845044   -1.567403   -5.534501 \n\
   \  C     -4.846076   -2.857958   -3.482926 \n\
   \  H     -4.951444   -0.606105   -6.061487 \n\
   \  H     -4.940748   -2.919843   -2.387384 \n\
   \  C     -4.985856   -1.621341   -4.134935 \n\
   \  H     -5.199724   -0.708428   -3.559349 \n\
   \  C     -6.064880   -8.159428   -6.502321 \n\
   \  C     -6.902684   -8.666070   -7.528392 \n\
   \  C     -6.287200   -8.571946   -5.164306 \n\
   \  H     -6.719675   -8.342501   -8.561073 \n\
   \  H     -5.629868   -8.235250   -4.350003 \n\
   \  C     -7.939157   -9.551808   -7.212489 \n\
   \  C     -7.339749   -9.449962   -4.865543 \n\
   \  H     -8.582731   -9.932542   -8.020706 \n\
   \  H     -7.497313   -9.759669   -3.820800 \n\
   \  C     -8.172934   -9.947100   -5.881656 \n\
   \  H     -8.994475  -10.638065   -5.639362 \n\
   \*\n\n"

_footer =
  "\
   \%eprnmr nuclei = all N { aiso, adip, fgrad} \n\
   \end \n\
   \"
 
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
    "5e-7"    ,        
    "1e-8"    ,
    "5e-8"    ,    
    "1e-9"    ,
    "5e-9"    ,    
    "1e-10"   ,
    "5e-10"   ,    
    "1e-11"   ,
    "5e-11"   ,    
    "1e-12"   ,
    "5e-12"   ,
    "0"
  ]

_genInput :: String -> String
_genInput tcutpno = (_header tcutpno) ++ _geom ++ _footer

fileName :: String -> String
fileName name = "./verdazyl-" ++ name ++ ".inp"

main = do
  let
    fileNames = map fileName names
    contents  = map (\x -> _genInput x) names
    nameCont  = zip fileNames contents
  res <- mapM (\x -> writeFile (fst x) (snd x)) nameCont
  putStrLn " Done.\n"

 
