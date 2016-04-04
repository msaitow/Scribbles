
import System.Directory
import Data.List as List
import Data.Maybe as Maybe

import Data.Time
import System.Process
import System.Environment

import MSUtils
  
-- Extract the correlation reference energy
-- Input lined contes of the orca output files
extractEnergies :: String -> (String, String)
extractEnergies contents = (zerothE, corrE)
  where
    tokenC  = fmap words $ lines contents
    zeroEcand = filter (\x ->  "E(0)" `elem` x) tokenC 
    zerothE = if length zeroEcand == 0 then "N.F." else last $ head zeroEcand
    corrEcand = filter (\x -> ("E(CORR)(corrected)" `elem` x) || ("E(CORR)" `elem` x) || ("E(CORR)(total)" `elem` x)) tokenC    
    corrE   = if length corrEcand == 0 then "N.F." else last $ head corrEcand

-- Generate the keyword for a given file
-- Input file name for a orca output              
genKey :: String -> String
genKey myName = (fst (splitAtDotForward myName)) ++ (snd (splitAtDotBackward myName))

-- Returns (keyword, (E0, Ecorr))
-- Input fileName 
getEnergy :: String -> IO (String, (String, String))
getEnergy myName = do
  cont <- readFile myName
  let keyword = genKey myName
      ene     =  extractEnergies cont
  return (keyword, ene)

-- Returns all the output
returnOutFiles :: [String] -> [String]
returnOutFiles inFiles = outFiles
  where
    -- Screen "." and ".."
    actualFiles = filter (\x -> (x /= ".") && (x /= "..")) inFiles
    -- Files names splitted by "."
    splitted    = map (splitAtDotBackward) actualFiles
    -- Fetch input files    
    outFiles  = List.nub $ map (\x -> (fst x) ++ (snd x)) $ filter (\x -> snd x == ".out") $ splitted

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

  
