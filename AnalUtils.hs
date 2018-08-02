
module AnalUtils
(
  hartree2kcal,
  extractEnergiesStrong,
  extractEnergies,
  returnOutFiles,
  getEnergy,
  getEnergyStrong,
  makeMaps,
  printData,
  printData',
  calcError,
  calcError',
  calcError2,
  calcRMSD,
  printRMSD,
  printRMSD',
  analData,
  analData'
) where

import Data.Maybe as Maybe
import Data.List as List
import Text.Printf as Printf
import MSUtils

-- Conversion coefficient
hartree2kcal = 627.5081068957
                   
-- Extract the correlation energy for strong pairs
-- Input: lined contes of the orca output files
extractEnergiesStrong :: String -> (String, String)
extractEnergiesStrong contents = (zerothE, corrE)
  where
    tokenC  = fmap words $ lines contents
    zeroEcand = filter (\x ->  "E(0)" `elem` x) tokenC 
    zerothE = if length zeroEcand == 0 then "N.F." else last $ head zeroEcand
    corrEcand = filter (\x -> ("E(CORR)(strong-pairs)" `elem` x) || ("E(CORR)" `elem` x)) tokenC
    corrE   = if length corrEcand == 0 then "N.F." else last $ head corrEcand

-- Extract the correlation energy
-- Input: lined contes of the orca output files
extractEnergies :: String -> (String, String)
extractEnergies contents = (zerothE, corrE)
  where
    tokenC  = fmap words $ lines contents
    zeroEcand = filter (\x ->  "E(0)" `elem` x) tokenC 
    zerothE = if length zeroEcand == 0 then "N.F." else last $ head zeroEcand
    corrEcand = filter (\x -> ("E(CORR)(corrected)" `elem` x) || ("E(CORR)" `elem` x) || ("E(CORR)(total)" `elem` x)) tokenC
    corrE   = if length corrEcand == 0 then "N.F." else last $ head corrEcand

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

-- Generate the keyword for a given file
-- Input: file name for a orca output              
genKey :: String -> String
genKey myName = (fst (splitAtDotForward myName)) ++ (snd (splitAtDotBackward myName))

-- Returns (keyword, (E0, Ecorr))
-- Input: fileName 
getEnergy :: String -> IO (String, (String, String))
getEnergy myName = do
  cont <- readFile myName
  let keyword = genKey myName
      ene     =  extractEnergies cont
  return (keyword, ene)

-- Returns (keyword, (E0, Ecorr for strong pairs))
-- Input: fileName 
getEnergyStrong :: String -> IO (String, (String, String))
getEnergyStrong myName = do
  cont <- readFile myName
  let keyword = genKey myName
      ene     = extractEnergiesStrong cont
  return (keyword, ene)

-- make map [keyword, canonica, lpno. dlpno]
-- canonical -> lpno -> dlpno
makeMaps :: [(String, (String, String))] -> [(String, (String, String))] -> [(String, (String, String))] -> [(String, String, String, String, String)]
makeMaps cList lList dList = fmap (\x -> findPartners x lList dList) cList
  where
    -- the first argument has to be canonical results and the second and third have to be lpno and dlpno
    findPartners :: (String, (String, String)) -> [(String, (String, String))] -> [(String, (String, String))] -> (String, String, String, String, String)
    findPartners (tag, (eref, ecorr)) lpnos dlpnos = (tag, eref, ecorr, lEne, dEne)
      where
        lData = head $ filter (\x -> fst x == tag) lpnos
        dData = head $ filter (\x -> fst x == tag) dlpnos
        -- if the reference energies don't match, put smething else than the value .. Here, we trust the Eref for the canonical case
        -- R.N. stands for reference doesn't match ..
        lEne = (snd $ snd lData)
        dEne = (snd $ snd dData)
--        lEne = if (fst $ snd lData) == eref then (snd $ snd lData) else "R.N."
--        dEne = if (fst $ snd dData) == eref then (snd $ snd dData) else "R.N."

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

-- Print max deviation and RMSD in percent
printRMSD :: [(String, String, String, String, String)] -> IO ()
printRMSD cData = do
  let
    dLPNOs  = calcRMSD $ fmap (\x -> Maybe.fromJust x) $ filter(\x -> x /= Nothing) $ fmap calcError' $ fmap (\(_, _, x, b, _) -> (x, b)) $ cData 
    dDLPNOs = calcRMSD $ fmap (\x -> Maybe.fromJust x) $ filter(\x -> x /= Nothing) $ fmap calcError' $ fmap (\(_, _, x, _, c) -> (x, c)) $ cData 

  putStrLn   ""
  putStrLn   " -----------------------------------"
  putStrLn   "         Max (%Error)  RMSD (%Error)"
  putStrLn   " -----------------------------------"
  Printf.printf " LPNO  : %10.3f %10.3f \n" (fst dLPNOs ) (snd dLPNOs )
  Printf.printf " DLPNO : %10.3f %10.3f \n" (fst dDLPNOs) (snd dDLPNOs)  
--  putStrLn $ " LPNO  : " ++ (show $ fst dLPNOs ) ++ ", " ++ (show $ snd dLPNOs )
--  putStrLn $ " DLPNO : " ++ (show $ fst dDLPNOs) ++ ", " ++ (show $ snd dDLPNOs)  
  putStrLn   " -----------------------------------"
  putStrLn   ""

-- Print max deviation and RMSD in Kcal/mol
printRMSD' :: [(String, String, String, String, String)] -> IO ()
printRMSD' cData = do
  let
    dLPNOs  = calcRMSD $ fmap (\x -> snd $ Maybe.fromJust x) $ filter(\x -> x /= Nothing) $ fmap calcError2 $ fmap (\(_, _, x, b, _) -> (x, b)) $ cData 
    dDLPNOs = calcRMSD $ fmap (\x -> snd $ Maybe.fromJust x) $ filter(\x -> x /= Nothing) $ fmap calcError2 $ fmap (\(_, _, x, _, c) -> (x, c)) $ cData 

  putStrLn   ""
  putStrLn   " ---------------------------------------"
  putStrLn   "         Max (Kcal/mol)  RMSD (Kcal/mol)"
  putStrLn   " ---------------------------------------"
  Printf.printf " LPNO  : %10.3f %10.3f \n" (fst dLPNOs ) (snd dLPNOs )
  Printf.printf " DLPNO : %10.3f %10.3f \n" (fst dDLPNOs) (snd dDLPNOs)    
--  putStrLn $ " LPNO  : " ++ (show $ fst dLPNOs ) ++ ", " ++ (show $ snd dLPNOs )
--  putStrLn $ " DLPNO : " ++ (show $ fst dDLPNOs) ++ ", " ++ (show $ snd dDLPNOs)  
  putStrLn   " ---------------------------------------"
  putStrLn   ""

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
