
import System.Directory
import Data.List as List
import Data.Maybe as Maybe

import Data.Time
import System.Process
import System.Environment

--import MSUtils

--------------------------------------------------------------------  
-- Variables
--------------------------------------------------------------------  
execName   = "orca"
scratchDir = "/work/$USER/"
mpiDir     = "/usr/local/orca/402/"
ldDir      = mpiDir ++ "orca_4_0_1_2_linux_x86-64_shared_openmpi202/"
--------------------------------------------------------------------

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

-- The return value should not contain any dot or slash
dotSlashRemover :: String -> String
dotSlashRemover s
  | b == "None" = a
  | b == ""     = if ( ('.' `elem` x) == False) && ( ('/' `elem` x) == False) then x else error $ "Input file should be in the same directory as you're now! >> (" ++ b ++ ", " ++ a ++ ", " ++ x ++ ")"
  | otherwise   = b               
  where
    (b, a) = splitAtDotForward s
    x      = fst $ splitAtDotForward (tail $ tail a)

main = do  
  
  putStrLn ""
  putStrLn " ==================="
  putStrLn " RunOrca            "
  putStrLn " ==================="
  putStrLn $ " * orca directory    : " ++ ldDir
  putStrLn $ " * mpi directory     : " ++ mpiDir
  
  -- Get the target directory
  args <- getArgs
  -- Get contents of the current repository
  if length args > 0 then do
    putStrLn $ " * Input file        : " ++ (args !! 0)
    else do
    error "Give the name of your input file! Fowlloing are the instructions. \n 1) Input file name should not include any dot except the one at the end e.g., ch2-rhf.inp is OK but ch2.rhf.inp is NOT OK!\n 2) You should call RunOrca in the same directory as the input file e.g., RunOrca ./o2.inp (or simply o2.inp) is OK but RunOrca ../o2.inp is NOT OK!"    

  -- Input label (label ++ ".inp" = input file name)
  let myLabel = if fst x == "None" then snd x else fst x
        where x = splitAtDotForward $ (dotSlashRemover $ args !! 0) ++ ".inp"     

  -- Get current time to make foot note
  time <- getZonedTime
  let
    (LocalTime d t) = zonedTimeToLocalTime time
    myScratchDir = scratchDir ++ myLabel ++ "-" ++ (show t) ++ "-" ++ (show d) ++ "/"

  putStrLn $ " * Scratch directory : " ++ myScratchDir

  -- Create the scratch directory
  create <- callCommand $ "mkdir -p " ++ myScratchDir
  
  -- Copy the input file into the scratch
  if myLabel == "None"
    then do error " Input file name should end by \".inp\"!"
    else do callCommand $ "cp ./" ++ myLabel ++ ".inp " ++ myScratchDir

  let myOrca = ldDir ++ execName

  -- Export proper exnvironmental variables
  let
    pathCommand  = "export PATH=" ++ mpiDir ++ "/bin:$PATH \n "
    ldCommand    = "export LD_LIBRARY_PATH=" ++ ldDir ++ ":$LD_LIBRARY_PATH \n "
    callORCA x = callCommand $ pathCommand ++ ldCommand ++ x

  -- Set Proper environmental variables
  setEnv "LD_LIBRARY_PATH" (mpiDir ++ "lib")
  
  -- Get contents of the current repository
  dirs <- getCurrentDirectory
  cons <- getDirectoryContents dirs
  -- Filter names of GBW files
  let gbwNames = map (\x -> (fst x) ++ (snd x)) $ filter (\x -> (snd x) == ".gbw") $ map (\x -> splitAtDotForward x) cons
  -- If there're GBW files in the current directory, we just copy them into scratch
  submitted <- mapM (\x -> callCommand $ "cp " ++ x ++ " " ++ myScratchDir) gbwNames
  
  -- Submit a calculation
  submit <- callORCA $ myOrca ++ " " ++ (myScratchDir ++ myLabel ++ ".inp") ++ " | tee ./" ++ myLabel ++ ".out"

  putStrLn "Calculation finished!"
  
  -- Convert the gbw file into molden file
  mv      <- callCommand $ "cp " ++ (myScratchDir ++ myLabel ++ ".gbw" ++ " ./")
  rm      <- callCommand $ "rm " ++ myScratchDir ++ "*.gbw"  
  convert <- callORCA    $ myOrca ++ "_2mkl " ++ myLabel ++ " -molden"
  rename  <- callCommand $ "mv " ++ myLabel ++ ".molden.input" ++ " ./" ++ myLabel ++ ".molden"
  --mv      <- callCommand $ "rm " ++ (myLabel ++ ".gbw")    

  putStrLn "Conversion of GBW file into molden file"
  
  -- Remove all the temp files
  --submit <- callCommand $ "rm -rf " ++ myScratchDir
  
