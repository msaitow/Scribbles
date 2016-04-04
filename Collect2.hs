
import System.Directory
import Data.List as List
import Data.Maybe as Maybe
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
        lEne = if (fst $ snd lData) == eref then (snd $ snd lData) else "R.N"
        dEne = if (fst $ snd dData) == eref then (snd $ snd dData) else "R.N"

-- Just print data ..
printDatum :: (String, String, String, String, String) -> IO ()
printDatum (fName, eref, ecorr, elpno, edlpno) = do
  let tag = fst $ splitAtDotBackward fName
  putStrLn $ " " ++ tag ++ " : " ++ eref ++ ", " ++ ecorr ++ ", " ++ elpno ++ ", " ++ edlpno

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
      
  putStrLn $ " " ++ tag ++ " : " ++ (show erefD) ++ ", " ++ (show ecorrD) ++ ", " ++ (show elpnoD) ++ ", " ++ (show edlpnoD)

-- Calculate the truncation error
calcError :: String -> String -> String
calcError eref epno
  | eref == "N.F." || epno == "N.F." = " N.C"
  | otherwise                       = show percent
  where
    percent = (erefD - (erefD - epnoD)) / erefD * 100.0
    erefD   = (read eref) :: Double    
    epnoD   = (read epno) :: Double
    
-- Analyze and print data ..
analDatum :: (String, String, String, String, String) -> IO ()
analDatum (fName, eref, ecorr, elpno, edlpno) = do
  let tag = fst $ splitAtDotBackward fName
  putStrLn $ " " ++ tag ++ " : " ++ (calcError ecorr elpno) ++ ", " ++ (calcError ecorr edlpno)
  
main = do
  -- Get the target directory
  args <- getArgs
  -- Get contents of the current repository
  if length args > 0 then do
    repo <- setCurrentDirectory $ (args !! 0) ++ "/canon"
    putStrLn "Target directory set \n"
    else do
    putStrLn "This program reads all the input files in the gven target directory and prints the correlation energy\n"
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

  putStrLn ""
  putStrLn " -------------------------"
  putStrLn "                Ecorr     "
  putStrLn "          ----------------"
  putStrLn " Tag Eref Canon LPNO DLPNO"
  putStrLn " -------------------------"
  mapM_ printDatum cData
  putStrLn " -------------------------"
  putStrLn ""

--  putStrLn ""
--  putStrLn " -------------------------"
--  putStrLn "                Ecorr     "
--  putStrLn "          ----------------"
--  putStrLn " Tag Eref Canon LPNO DLPNO"
--  putStrLn " -------------------------"
--  mapM_ printDatum' cData
--  putStrLn " -------------------------"
--  putStrLn ""

  putStrLn ""
  putStrLn " -------------------"
  putStrLn "        Ecorr (%)   "
  putStrLn "    ----------------"
  putStrLn " Tag LPNO DLPNO"
  putStrLn " -------------------"
  mapM_ analDatum cData
  putStrLn " -------------------"
  putStrLn ""

  print "End"

  
