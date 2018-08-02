
import System.Directory
import System.Environment

import AnalUtils

main = do
  -- Get the target directory
  args <- getArgs
  -- Get contents of the current repository
  if length args > 0 then do
    repo <- setCurrentDirectory $ (args !! 0) ++ "/canon"
    putStrLn "Target directory set \n"
    else do
    putStrLn "This program reads all the orca output files in the given target directory and its subdirectory and then it prints the correlation energy"
    putStrLn ""
    putStrLn "Argument: path to directory which should contain three subdirectories, canon/, cdlpno/ and odlpno/. All the output files located in each subdirectories"
    putStrLn "          are read the ODLPNO truncation errors are calculaed. Note that here errors in total pair correlation energies are analyzed."    
    putStrLn ""    
    error "Aborting further operation"

  -- (1) Canonical results
  -- Get contents of the current repository
  dirs <- getCurrentDirectory
  cons <- getDirectoryContents dirs
  let outFilesCanon = returnOutFiles cons
  results_canon <- mapM getEnergy outFilesCanon
  
  -- (3) ODLPNO results
  repo <- setCurrentDirectory $ (args !! 0) ++ "/odlpno"
  -- Get contents of the current repository
  dirs <- getCurrentDirectory
  cons <- getDirectoryContents dirs
  let outFilesDLPNO = returnOutFiles cons
  results_dlpno <- mapM getEnergy outFilesDLPNO

  -- Collect the results
  let cData = makeMaps results_canon results_dlpno results_dlpno
  
  print "Outputfiles found .. "
  --print outFiles
  print "ResultsCanon"  
  print results_canon
  print "ResultsODLPNO"  
  print results_dlpno
  print "Data"
  print cData

  printData  cData
  printData' cData
  analData   cData
  analData'  cData
  printRMSD  cData
  printRMSD' cData
  
  print "End"

  
