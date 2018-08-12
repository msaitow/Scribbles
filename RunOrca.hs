
import System.Directory
import Data.List as List
import Data.Maybe as Maybe

import Data.Time
import System.Process
import System.Environment

import MSUtils
  
-- Variables
scratchDir = "/tmp/"
myOrca = "orca"

main = do
  -- Get the target directory
  args <- getArgs
  -- Get contents of the current repository
  if length args > 0 then do
    putStrLn $ " * Input file is " ++ (args !! 0)
    else do
    error "Give the name of your input file! Fowlloing are the instructions. \n 1) Input file name should not include any dot except the one at the end e.g., ch2-rhf.inp is OK but ch2.rhf.inp is NOT OK!\n 2) You should call RunOrca in the same directory as the input file e.g., RunOrca ./o2.inp (or simply o2.inp) is OK but RunOrca ../o2.inp is NOT OK!"

  -- Input label (label ++ ".inp" = input file name)      
  let myLabel = fst $ splitAtDotForward $ args !! 0

  -- Get current time to make foot note
  time <- getZonedTime
  let
    (LocalTime d t) = zonedTimeToLocalTime time
    myScratchDir = scratchDir ++ myLabel ++ "-" ++ (show t) ++ "-" ++ (show d) ++ "/"

  putStrLn $ " * Scratch directory : " ++ myScratchDir

  -- Create the scratch directory
  create <- callCommand $ "mkdir " ++ myScratchDir

  -- Copy the input file into the scratch
  if myLabel == "None"
    then do error " Input file name should end by \".inp\"!"
    else do callCommand $ "cp ./" ++ myLabel ++ ".inp " ++ myScratchDir

  -- Submit a calculation
  submit <- callCommand $ myOrca ++ " " ++ (myScratchDir ++ myLabel ++ ".inp") ++ " | tee ./" ++ myLabel ++ ".out"

  putStrLn "Calculation finished!"
  
  -- Convert the gbw file into molden file
  mv      <- callCommand $ "mv " ++ (myScratchDir ++ myLabel ++ ".gbw" ++ " ./")  
  convert <- callCommand $ myOrca ++ "_2mkl " ++ myLabel ++ " -molden"
  rename  <- callCommand $ "mv " ++ myLabel ++ ".molden.input" ++ " ./" ++ myLabel ++ ".molden"
  mv      <- callCommand $ "rm " ++ (myLabel ++ ".gbw")    

  putStrLn "Conversion of GBW file into molden file"
  
  -- Remove all the temp files
  --submit <- callCommand $ "rm -rf " ++ myScratchDir
  
