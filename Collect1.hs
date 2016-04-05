
import System.Directory
import Data.List as List
import Data.Maybe as Maybe

import Data.Time
import System.Process
import System.Environment

import MSUtils
import CollectUtils

main = do
  -- Get the target directory
  args <- getArgs
  -- Get contents of the current repository
  if length args > 0 then do
    repo <- setCurrentDirectory (args !! 0)
    putStrLn "Target directory set \n"
    else do
    putStrLn "This program reads all the input files in the gven target directory and prints the correlation energy\n"
    error "Aborting further operation"
    
  -- Get contents of the current repository
  dirs <- getCurrentDirectory
  cons <- getDirectoryContents dirs

  let outFiles = returnOutFiles cons
  results <- mapM getEnergy outFiles

  print "Outputfiles found .. "
  print outFiles
  print "Results"  
  print results
  print "End"

  
