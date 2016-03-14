
import System.Directory
import Data.List as List
import Data.Maybe as Maybe

import Data.Time
import System.Process
import System.Environment

-- How to compile on clusters:
-- bash> ghc --make -L/home/saitow/lib ./ChangeHeader.hs

-- newline replaces the hotlne
replaceHeader :: String -> String -> FilePath -> IO ()
replaceHeader hotLine newLine fileName = do
  contents <- readFile fileName
  print (lines contents) -- debug
  let corrected     = unlines $ fmap (\x -> if x == hotLine then newLine else x) $ lines contents
  if (corrected == contents) then print " Nothing is found .. " 
  else writeFile (fileName ++ ".corrected.inp") corrected

-- This is supposed to work in the exactly the same way as the above one
-- but splits from the end so that "C90.loose.inp" is splitted into ("C90.loose.", "inp")
splitAtDotBackward :: String -> (String, String)
splitAtDotBackward x
  | myIndex /= Nothing = if (last $ fst splittedPair) == '.'
                         then ((init $ fst splittedPair), "." ++ (snd splittedPair))
                         else  splittedPair
  | otherwise          = ("None", "")
  where
    splittedPair = (reverse $ snd revsplitted, reverse $ fst revsplitted) 
    revsplitted  = splitAt (Maybe.fromJust myIndex) reversed    
    reversed     = reverse x

    findDotIndex :: String -> Maybe Int
    findDotIndex kore = List.findIndex (== '.') kore

    myIndex = findDotIndex reversed

-- Returns all the output and intermediate files whose input is there
-- For example:    
--   Input  [".", "..", "C90.inp", "C90.tmp", "C90.h.tmp", "C90.out", "C120.out"]
--   Output ["C90.tmp", "C90.out"]     
returnFilesToBeRenamed :: [String] -> [String]
returnFilesToBeRenamed inFiles = inputFiles
  where
    -- Screen "." and ".."
    actualFiles = filter (\x -> (x /= ".") && (x /= "..")) inFiles
    -- Files names splitted by "."
    splitted    = map (splitAtDotBackward) actualFiles
    -- Fetch input files    
    inputFiles  = List.nub $ map (\x -> (fst x) ++ (snd x)) $ filter (\x -> (snd x == ".inp") || (snd x == ".in")) $ splitted

main = do
  -- For angs
  let hotLine = "#! UHF  Angs LPNO-CCSD cc-pVQZ cc-pVQZ/C Decontract VeryTightSCF PModel NOFROZENCORE MiniPrint Conv"
      newLine = "! UHF Angs LPNO-CCSD cc-pVQZ cc-pVQZ/C Decontract DecontractAux RIJK VeryTightSCF NOFROZENCORE MiniPrint"
--  -- For bohr
--  let hotLine = "! UHF  Bohrs LPNO-CCSD cc-pVQZ cc-pVQZ/C Decontract VeryTightSCF PModel NOFROZENCORE MiniPrint Conv"
--      newLine = "! UHF  Bohrs LPNO-CCSD cc-pVQZ cc-pVQZ/C Decontract DecontractAux RIJK VeryTightSCF NOFROZENCORE MiniPrint"      
  -- Get the target directory
  args <- getArgs
  -- Get contents of the current repository
  if length args > 0 then do
    repo <- setCurrentDirectory (args !! 0)
    putStrLn "Target directory set \n"
  else do
    repo <- setCurrentDirectory "./"
    putStrLn "Default directory is ./ \n"
  
  -- Get contents of the current repository
  dirs <- getCurrentDirectory
  cons <- getDirectoryContents dirs
  let fileNames = returnFilesToBeRenamed cons

  -- Replace!
  ios <- mapM_ (\x -> replaceHeader hotLine newLine x) fileNames
   
  print fileNames 
  print "End"
