
import System.Directory
import Data.List as List
import Data.Maybe as Maybe

import Data.Time
import System.Process
import System.Environment

import MSUtils
  
-- Returns all the output and intermediate files whose input is there
-- For example:    
--   Input  [".", "..", "lct_entry.cpp", "lct_entry.h", "lct_pnosolver.cpp", "lct_LHS_V0_V0.h", "lct_LHS_V0_V0.cpp", "lct_LHS_V0_Vm1.cpp"]
--   Output ["lct_LHS_V0_V0.cpp", "lct_LHS_V0_Vm1.cpp"]     
returnFiles :: [String] -> [String]
returnFiles inFiles = lhsFiles
  where
    -- Screen "." and ".."
    actualFiles = filter (\x -> (x /= ".") && (x /= "..")) inFiles
    -- Files names splitted by "."
    splitted    = map (splitAtDotBackward) actualFiles
    -- Fetch input files    
    cppFiles  = map (\x -> (fst x)++(snd x)) $ List.nub $ filter (\x -> snd x == ".cpp") $ splitted
    -- C++ source files for LHS vectors
    lhsFiles  = [x | x<-cppFiles, List.isPrefixOf "lct_LHS_" x]

-- Do a sed command
doSed :: [String] -> [String]
doSed inFiles = map (\x -> sed ++ " " ++ x ++ cp ++ x) actualFiles
  where
    -- Screen "." and ".."
    actualFiles = filter (\x -> (x /= ".") && (x /= "..")) inFiles
    --
    sed = "sed -e \"s/ovl.GetTensor(/ovl.CopyTensor(/g\""
    cp  = " > ./sedded/"
  
main = do
  -- Get the target directory
  args <- getArgs
  -- Get contents of the current repository
  if length args > 0 then do
    repo <- setCurrentDirectory (args !! 0)
    putStrLn "Target directory set \n"
    else do
    repo <- setCurrentDirectory "./"
    putStrLn "Default directory is ./ \n"
  
  -- Get contents of the current repository
  dirs <- getCurrentDirectory
  cons <- getDirectoryContents dirs
  
  -- Submit jobs
  let jobs = doSed $ returnFiles cons

  print "Jobs .. "
  print jobs
      
  submitted <- mapM (\x -> callCommand x) jobs
  --submitted <- mapM (\x -> putStr x) jobs
  
  print "End"
  
