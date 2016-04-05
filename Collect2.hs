
import System.Directory
import Data.List as List
import Data.Maybe as Maybe
import System.Environment

import Text.Printf as Printf

import MSUtils
import CollectUtils

-- Just Print data
printData :: [(String, String, String, String, String)] -> IO ()
printData inp = do
  putStrLn ""
  putStrLn " -------------------------"
  putStrLn "                Ecorr     "
  putStrLn "          ----------------"
  putStrLn " Tag Eref Canon LPNO DLPNO"
  putStrLn " -------------------------"
  mapM_ printDatum inp
  putStrLn " -------------------------"
  putStrLn ""
  
-- Just print data ..
printDatum :: (String, String, String, String, String) -> IO ()
printDatum (fName, eref, ecorr, elpno, edlpno) = do
  let tag = fst $ splitAtDotBackward fName
  putStrLn $ " " ++ tag ++ " : " ++ eref ++ ", " ++ ecorr ++ ", " ++ elpno ++ ", " ++ edlpno

-- Just print data
printData' :: [(String, String, String, String, String)] -> IO ()
printData' inp = do
  putStrLn ""
  putStrLn " -------------------------"
  putStrLn "                Ecorr     "
  putStrLn "          ----------------"
  putStrLn " Tag Eref Canon LPNO DLPNO"
  putStrLn " -------------------------"
  mapM_ printDatum' inp
  putStrLn " -------------------------"
  putStrLn ""
  
-- Just print data ..
printDatum' :: (String, String, String, String, String) -> IO ()
printDatum' (fName, eref, ecorr, elpno, edlpno) = do
  let tag   = fst $ splitAtDotBackward fName
      
      erefD   :: Double
      erefD   = if eref   /= "N.F." then read eref   else 200.0
      ecorrD  :: Double
      ecorrD  = if ecorr  /= "N.F." then read ecorr  else 200.0
      elpnoD  :: Double
      elpnoD  = if elpno  /= "N.F." then read elpno  else 200.0
      edlpnoD :: Double
      edlpnoD = if edlpno /= "N.F." then read edlpno else 200.0
      
  --putStrLn $ " " ++ tag ++ " : " ++ (show erefD) ++ ", " ++ (show ecorrD) ++ ", " ++ (show elpnoD) ++ ", " ++ (show edlpnoD)
  Printf.printf " %10s : %20.10f %14.10f %14.10f %14.10f \n" tag erefD ecorrD elpnoD edlpnoD
  
-- Calculate the truncation error
calcError :: String -> String -> String
calcError eref epno
  | eref == "N.F." || epno == "N.F." = " N.C"
  | otherwise                       = show percent
  where
    percent = (erefD - (erefD - epnoD)) / erefD * 100.0
    erefD   = (read eref) :: Double    
    epnoD   = (read epno) :: Double

-- Calculate the truncation error
calcError' :: (String, String) -> Maybe Double
calcError' (eref, epno)
  | eref == "N.F." || epno == "N.F." = Nothing
  | otherwise                        = Just percent
  where
    percent =  - (erefD - epnoD) / erefD * 100.0
    erefD   = (read eref) :: Double    
    epnoD   = (read epno) :: Double

-- Calculate the truncation error
calcError2 :: (String, String) -> Maybe (Double, Double)
calcError2 (eref, epno)
  | eref == "N.F." || epno == "N.F." = Nothing
  | otherwise                        = Just (percent, error_kcal)
  where
    hartree2kcal = 627.5081068957
    percent      = (erefD - (erefD - epnoD)) / erefD * 100.0
    error_kcal   = hartree2kcal * (erefD - epnoD)
    erefD        = (read eref) :: Double    
    epnoD        = (read epno) :: Double

-- [Error in correlation energy (%)] -> (maximim error, RMSD)
calcRMSD :: [Double] -> (Double,Double)
calcRMSD inp = (max, rmsd)
  where
    rmsd = sqrt $ (sum $ fmap (\x -> x ^ 2) inp) / (fromIntegral $ length inp)
    max  = maximum $ fmap abs inp

-- Just Print data
analData :: [(String, String, String, String, String)] -> IO ()
analData inp = do
  putStrLn ""
  putStrLn " -------------------"
  putStrLn "        Ecorr (%)   "
  putStrLn "    ----------------"
  putStrLn " Tag LPNO DLPNO"
  putStrLn " -------------------"
  mapM_ analDatum inp
  putStrLn " -------------------"
  putStrLn ""

-- Analyze and print data ..
analDatum :: (String, String, String, String, String) -> IO ()
analDatum (fName, eref, ecorr, elpno, edlpno) = do
  let tag = fst $ splitAtDotBackward fName
  --putStrLn $ " " ++ tag ++ " : " ++ (calcError ecorr elpno) ++ ", " ++ (calcError ecorr edlpno)
  Printf.printf " %10s : %20s %20s\n" tag (calcError ecorr elpno) (calcError ecorr edlpno)

-- Just Print data
analData' :: [(String, String, String, String, String)] -> IO ()
analData' inp = do
  Printf.printf " Found %d entry .. \n" (length inp) 
  putStrLn ""
  putStrLn " ------------------------------------------"
  putStrLn " Tag LPNO(%) (Kcal/mol) DLPNO(%) (Kcal/mol)"
  putStrLn " ------------------------------------------"  
  mapM_ analDatum' inp
  putStrLn " ------------------------------------------"    
  putStrLn ""

-- Analyze and print data ..
analDatum' :: (String, String, String, String, String) -> IO ()
analDatum' (fName, eref, ecorr, elpno, edlpno) = do
  let tag = fst $ splitAtDotBackward fName
  let lpnoError  = calcError2 (ecorr,  elpno)
  let dlpnoError = calcError2 (ecorr, edlpno)  
  if      lpnoError /= Nothing && dlpnoError /= Nothing then Printf.printf " %10s : %9.2f (%5.2f) %9.2f (%5.2f) \n"  tag (fst $ Maybe.fromJust lpnoError) (snd $ Maybe.fromJust lpnoError)
                                                             (fst $ Maybe.fromJust dlpnoError) (snd $ Maybe.fromJust dlpnoError)   
  else if lpnoError == Nothing && dlpnoError /= Nothing then Printf.printf " %10s : %9s         %9.2f (%5.2f) \n"  tag "N.C." 
                                                             (fst $ Maybe.fromJust dlpnoError) (snd $ Maybe.fromJust dlpnoError)   
  else if lpnoError /= Nothing && dlpnoError == Nothing then Printf.printf " %10s : %9.2f %9s   \n"       tag (fst $ Maybe.fromJust lpnoError) "N.C."    
  else                                                       Printf.printf " %10s : %9s         %9s \n"           tag "N.C."                           "N.C."    

main = do
  -- Get the target directory
  args <- getArgs
  -- Get contents of the current repository
  if length args > 0 then do
    repo <- setCurrentDirectory $ (args !! 0) ++ "/canon"
    putStrLn "Target directory set \n"
    else do
    putStrLn "This program reads all the orca output files in the gven target directory and prints the correlation energy\n"
    error "Aborting further operation"

  -- (1) Canonical results
  -- Get contents of the current repository
  dirs <- getCurrentDirectory
  cons <- getDirectoryContents dirs
  let outFilesCanon = returnOutFiles cons
  results_canon <- mapM getEnergy outFilesCanon
  
  -- (2) LPNO results
  repo <- setCurrentDirectory $ (args !! 0) ++ "/lpno"
  -- Get contents of the current repository
  dirs <- getCurrentDirectory
  cons <- getDirectoryContents dirs
  let outFilesLPNO = returnOutFiles cons
  results_lpno <- mapM getEnergy outFilesLPNO

  -- (3) DLPNO results
  repo <- setCurrentDirectory $ (args !! 0) ++ "/dlpno"
  -- Get contents of the current repository
  dirs <- getCurrentDirectory
  cons <- getDirectoryContents dirs
  let outFilesDLPNO = returnOutFiles cons
  results_dlpno <- mapM getEnergy outFilesDLPNO

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

  let
    dLPNOs  = calcRMSD $ fmap (\x -> Maybe.fromJust x) $ filter(\x -> x /= Nothing) $ fmap calcError' $ fmap (\(_, _, x, b, _) -> (x, b)) $ cData 
    dDLPNOs = calcRMSD $ fmap (\x -> Maybe.fromJust x) $ filter(\x -> x /= Nothing) $ fmap calcError' $ fmap (\(_, _, x, _, c) -> (x, c)) $ cData 

  putStrLn   ""
  putStrLn   " -------------------"
  putStrLn   "         Max   RMSD "
  putStrLn   " -------------------"
  putStrLn $ " LPNO  : " ++ (show $ fst dLPNOs ) ++ ", " ++ (show $ snd dLPNOs )
  putStrLn $ " DLPNO : " ++ (show $ fst dDLPNOs) ++ ", " ++ (show $ snd dDLPNOs)  
  putStrLn   " -------------------"
  putStrLn   ""
  
  print "End"

  
