
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
    putStrLn "Argument: path to directory which should contain three subdirectories, canon/, lpno/ and dlpno/. All the output files located in each subdirectories"
    putStrLn "          are read and the LPNO and DLPNO truncation errors are calculaed. Note that here errors in strong pair correlation energies are analyzed."    
    putStrLn ""    
    error "Aborting further operation"

  -- (1) Canonical results
  -- Get contents of the current repository
  dirs <- getCurrentDirectory
  cons <- getDirectoryContents dirs
  let outFilesCanon = returnOutFiles cons
  results_canon <- mapM getEnergyStrong outFilesCanon
  
  -- (2) LPNO results
  repo <- setCurrentDirectory $ (args !! 0) ++ "/lpno"
  -- Get contents of the current repository
  dirs <- getCurrentDirectory
  cons <- getDirectoryContents dirs
  let outFilesLPNO = returnOutFiles cons
  results_lpno <- mapM getEnergyStrong outFilesLPNO

  -- (3) DLPNO results
  repo <- setCurrentDirectory $ (args !! 0) ++ "/dlpno"
  -- Get contents of the current repository
  dirs <- getCurrentDirectory
  cons <- getDirectoryContents dirs
  let outFilesDLPNO = returnOutFiles cons
  results_dlpno <- mapM getEnergyStrong outFilesDLPNO

  -- Collect the results
  let cData = makeMaps results_canon results_lpno results_dlpno
  
  print "Outputfiles found .. "
  --print outFiles
  print "ResultsCanon"  
  print results_canon
  print "ResultsLPNO"  
  print results_lpno
  print "ResultsDLPNO"  
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

  
