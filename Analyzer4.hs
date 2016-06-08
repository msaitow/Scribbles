
import System.Directory
import System.Environment
import Data.Maybe as Maybe
import Data.List as List
import Text.Printf as Printf
import MSUtils
import AnalUtils

makeMaps' :: [(String, (String, String))] -> [(String, (String, String))] -> [(String, String, String, String)]
makeMaps' cList oList = fmap (\x -> findPartner x oList) cList
  where
    -- the first argument has to be canonical results and the second and third have to be lpno and dlpno
    findPartner :: (String, (String, String)) -> [(String, (String, String))] -> (String, String, String, String)
    findPartner (tag, (eref, cEne)) odlpnos = (tag, eref, cEne, oEne)
      where
        oData = head $ filter (\x -> fst x == tag) odlpnos
        -- if the reference energies don't match, put smething else than the value .. Here, we trust the Eref for the canonical case
        -- R.N. stands for reference doesn't match .. 
        oEne = if (fst $ snd oData) == eref then (snd $ snd oData) else "R.N."

analDataMod :: [(String, String, String, String)] -> IO ()
analDataMod inp = do
  Printf.printf " Found %d entry .. \n" (length inp) 
  putStrLn ""
  putStrLn " ------------------------------------------"
  putStrLn " Tag ODLPNO(%) (Kcal/mol)"
  putStrLn " ------------------------------------------"  
  mapM_ analDatumMod inp
  putStrLn " ------------------------------------------"    
  putStrLn ""

-- Analyze and print data ..
analDatumMod :: (String, String, String, String) -> IO ()
analDatumMod (fName, eref, cdlpno, odlpno) = do
  let tag = fst $ splitAtDotBackward fName
  let odlpnoError = calcError2 (cdlpno,  odlpno)
  if      odlpnoError /= Nothing then Printf.printf " %10s : %9.2f (%5.2f) \n"  tag (fst $ Maybe.fromJust odlpnoError) (snd $ Maybe.fromJust odlpnoError)
  else                                Printf.printf " %10s : %9s           \n"  tag  "N.C."    

-- Print max deviation and RMSD in percent
printRMSDMod :: [(String, String, String, String)] -> IO ()
printRMSDMod cData = do
  let
    dDLPNOs = calcRMSD $ fmap (\x -> Maybe.fromJust x) $ filter(\x -> x /= Nothing) $ fmap calcError' $ fmap (\(_, _, x, c) -> (x, c)) $ cData 

  putStrLn   ""
  putStrLn   " -----------------------------------"
  putStrLn   "         Max (%Error)  RMSD (%Error)"
  putStrLn   " -----------------------------------"
  Printf.printf " DLPNO : %10.3f %10.3f \n" (fst dDLPNOs) (snd dDLPNOs)  
  putStrLn   " -----------------------------------"
  putStrLn   ""

-- Print max deviation and RMSD in Kcal/mol
printRMSDMod' :: [(String, String, String, String)] -> IO ()
printRMSDMod' cData = do
  let
    dDLPNOs = calcRMSD $ fmap (\x -> snd $ Maybe.fromJust x) $ filter(\x -> x /= Nothing) $ fmap calcError2 $ fmap (\(_, _, x, c) -> (x, c)) $ cData 

  putStrLn   ""
  putStrLn   " ---------------------------------------"
  putStrLn   "         Max (Kcal/mol)  RMSD (Kcal/mol)"
  putStrLn   " ---------------------------------------"
  Printf.printf " DLPNO : %10.3f %10.3f \n" (fst dDLPNOs) (snd dDLPNOs)    
  putStrLn   " ---------------------------------------"
  putStrLn   ""

main = do
  -- Get the target directory
  args <- getArgs
  -- Get contents of the current repository
  if length args > 0 then do
    repo <- setCurrentDirectory $ (args !! 0) ++ "/cdlpno"
    putStrLn "Target directory set \n"
    else do
    putStrLn "This program reads all the orca output files in the given target directory and its subdirectory and then it prints the correlation energy"
    putStrLn ""
    putStrLn "Argument: path to directory which should contain three subdirectories, cdlpno/ and odlpno/. All the output files located in each subdirectories"
    putStrLn "          are read and the closed and open-shell DLPNO truncation errors are calculaed for closed-shell molecules. Note that here errors in strong pair correlation energies are analyzed."    
    putStrLn ""    
    error "Aborting further operation"

  -- (1) Closed-shell DLPNO results
  -- Get contents of the current repository
  dirs <- getCurrentDirectory
  cons <- getDirectoryContents dirs
  let outFilesCDLPNO = returnOutFiles cons
  results_cdlpno <- mapM getEnergy outFilesCDLPNO
  
  -- (2) Open-shell DLPNO results
  repo <- setCurrentDirectory $ (args !! 0) ++ "/odlpno"
  -- Get contents of the current repository
  dirs <- getCurrentDirectory
  cons <- getDirectoryContents dirs
  let outFilesODLPNO = returnOutFiles cons
  results_odlpno <- mapM getEnergy outFilesODLPNO

  -- Collect the results
  let cData = makeMaps' results_cdlpno results_odlpno
  
  print "Outputfiles found .. "
  --print outFiles
  print "Results_CDLPNO"  
  print results_cdlpno
  print "Results_ODLPNO"  
  print results_odlpno
  print "Data"
  print cData

  analDataMod  cData
--  printData' cData
--  analData   cData
--  analData'  cData
  printRMSDMod  cData
  printRMSDMod' cData
  
  print "End"

  
