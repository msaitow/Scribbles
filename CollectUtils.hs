
module CollectUtils
(
  extractEnergiesStrong,
  extractEnergies,
  returnOutFiles,
  getEnergy,
  getEnergyStrong,
  makeMaps
) where

import Data.List as List  
import MSUtils

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
      ene     =  extractEnergiesStrong cont
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
        lEne = if (fst $ snd lData) == eref then (snd $ snd lData) else "R.N."
        dEne = if (fst $ snd dData) == eref then (snd $ snd dData) else "R.N."
