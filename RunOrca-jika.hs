
import System.Directory
import Data.List as List
import Data.Maybe as Maybe

import Data.Time
import System.Process
import System.Environment

-- Split the input into forward and afterward in terms of '.'
-- ex) splitAtDot "C90.inp" -> ("C90", ".inp")
-- NOTE: This splits, for example, "C90.loose.inp" into ("C90.", "loose.inp")
splitAtDotForward :: String -> (String, String)
splitAtDotForward x
  | myIndex /= Nothing = splitAt (Maybe.fromJust myIndex) x
  | otherwise          = ("None", "")
  where
    findDotIndex :: String -> Maybe Int
    findDotIndex kore = List.findIndex (== '.') kore

    myIndex = findDotIndex x

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
    error "Give the name of your input file\n"

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
  
