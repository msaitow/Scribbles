
import System.Directory
import Data.List as List
import Data.Maybe as Maybe

import Data.Time
import System.Process
import System.Environment

myBasis = ["def2-TZVPP", "def2-QZVPP"]
myPNOs  = ["NormalPNO", "TightPNO"]

genHeader :: String -> String -> String
genHeader basis pno =
  "\
   \! UHF DLPNO-CCSD " ++ basis ++ " " ++ basis ++ "/C VeryTightSCF NoAutoStart " ++ pno ++ "\n\
   \ \n\
   \%maxcore 4000 \n\
   \ \n\
   \%mdci \n\
   \      UseQROs      true \n\
   \      CIType       CCSD \n\
   \      MaxIter      35 \n\
   \      Randomize    false \n\
   \      end \n\n"  
                 
main = do
  let
    dirNames = ["martins." ++ x ++ "." ++ y | x <- myBasis, y <- myPNOs]
    headers  = [genHeader x y | x <- myBasis, y <- myPNOs]

  -- (1) make each directory
  mapM (\x -> callCommand $ "mkdir ./" ++ x) dirNames
  -- (2) generate header files
  mapM (\x -> writeFile (fst x) (snd x)) $ zip (map (\x -> "./" ++ x ++ "/header") dirNames) headers
  -- (3) copy gemetry
  mapM (\x -> callCommand $ "cp ./geom/*xyz " ++ x) dirNames
  -- (4) generate input
  mapM (\x -> callCommand $ "GenInput header " ++ x) dirNames
  -- (5) submit jobs
  mapM (\x -> callCommand $ "SubmitterPal " ++ x) dirNames
  print "End"
  
