
import System.Directory
import Data.List as List
import Data.Maybe as Maybe

import Data.Time
import System.Process
import System.Environment

import MSUtils
  
-- Returns all the output and intermediate files whose input is there
-- For example:    
--   Input  [".", "..", "C90.inp", "C90.tmp", "C90.h.tmp", "C90.out", "C120.out"]
--   Output ["C90.tmp", "C90.out"]     
returnFilesToBeRenamed :: [String] -> [String]
returnFilesToBeRenamed inFiles = withInput
  where
    -- Screen "." and ".."
    actualFiles = filter (\x -> (x /= ".") && (x /= "..")) inFiles
    -- Files names splitted by "."
    splitted    = map (splitAtDotBackward) actualFiles
    -- Fetch input files    
    inputFiles  = List.nub $ map (\x -> fst x) $ filter (\x -> snd x == ".inp") $ splitted
    -- Files which has the input
    withInput   = [(fst x)++(snd x) | x<-splitted, (fst x) `List.elem` inputFiles, (((snd x) == ".out") || ((snd x) == ".gbw"))]

-- Returns a list of pairs of names of files: [("original file name", "changed file name")]
-- where changed file name is "original file name" ++ footer
renameFiles :: String -> [String] -> [(String,String)]
renameFiles footer inFiles = zip inFiles changedNames
  where changedNames = map (++footer) inFiles

-- Move file
moveFile :: String -> (String, String) -> IO ()
moveFile path (iFile, oFile) = callCommand command
  where command = "cp " ++ path ++ iFile ++ " " ++ path ++ oFile  

-- Submit a job 
submitJobs :: [String] -> [String]
submitJobs inFiles = map (\x -> mySubOrca ++ " " ++ x) inputFiles 
  where
    -- Screen "." and ".."
    actualFiles = filter (\x -> (x /= ".") && (x /= "..")) inFiles
    -- Files names splitted by "."
    splitted    = map (splitAtDotBackward) actualFiles
    -- Fetch input files    
    inputFiles  = List.nub $ map (\x -> fst x) $ filter (\x -> snd x == ".inp") $ splitted
    -- Name of suborca
    mySubOrca = "~/bin/ms_suborca2"
  
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
  
--  -- Get current time to make foot note
--  time <- getZonedTime
--  let
--    (LocalTime d t) = zonedTimeToLocalTime time
--    myFooter = ".svd." ++ (show t) ++ "-" ++ (show d)
  
  -- Get contents of the current repository
  dirs <- getCurrentDirectory
  cons <- getDirectoryContents dirs
--  -- Get file names to be changed
--  let
--    fileNames = returnFilesToBeRenamed cons
--    newNames  = renameFiles myFooter fileNames
--
--  if False then do
--    putStr " Not moved"
--    else do
--    -- (1) Rename previous ouput files
--    moved <- mapM (\x -> moveFile ((show dirs) ++ "/") x) newNames 
--    print " Moved"

  -- Submit jobs
  let jobs = submitJobs cons

  print "~"
  print jobs
      
  submitted <- mapM (\x -> callCommand x) jobs
  
  print "End"
  
