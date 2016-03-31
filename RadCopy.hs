
import System.Directory
import Data.List as List
import Data.Maybe as Maybe

import Control.Monad
import Data.Time
import System.Process
import System.Environment

-- How to compile on clusters:
-- bash> ghc --make -L/home/saitow/lib ./RadMove.hs

-- newline replaces the hotlne
radOrNot :: String -> FilePath -> IO (Maybe FilePath)
radOrNot hotLine fileName = do
  contents <- readFile fileName
  let isClosedLines = fmap (\x -> if x == hotLine then True else False) $ lines contents
      isClosed      = foldl (||) False isClosedLines
      name          = if isClosed == False then Just fileName else Nothing
  return name
  
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

-- Returns all the xyz files
returnFilesToBeRenamed :: [String] -> [String]
returnFilesToBeRenamed inFiles = inputFiles
  where
    -- Screen "." and ".."
    actualFiles = filter (\x -> (x /= ".") && (x /= "..")) inFiles
    -- Files names splitted by "."
    splitted    = map (splitAtDotBackward) actualFiles
    -- Fetch input files    
    inputFiles  = List.nub $ map (\x -> (fst x) ++ (snd x)) $ filter (\x -> (snd x == ".xyz")) $ splitted

-- Copy file
cpFile :: String -> String -> IO ()
cpFile path name = callCommand $ "cp " ++ path ++ name ++ " ./rad/"

main = do
  -- For closed-shells
  let hotLine = "0 1"
  -- Get the target directory
  args <- getArgs
  -- Get contents of the current repository
  if length args > 0 then do
    repo <- setCurrentDirectory (args !! 0)
    putStrLn "Target directory set \n"
  else do
    repo <- setCurrentDirectory "../"
    putStrLn "Default directory is ../ \n"
  
  -- Get contents of the current repository
  dirs <- getCurrentDirectory
  cons <- getDirectoryContents dirs
  let fileNames = returnFilesToBeRenamed cons

  -- Detect the radical geometries
  rads  <- mapM (\x -> radOrNot "0 1" x) fileNames
  let names = map (\x -> Maybe.fromJust x) $ filter (\x -> x /= Nothing) rads

  print "Target Files:"
  print fileNames 
  print "Radical XYZ Files:"
  print names 

  print "Copying Radicals into rad/"
  mapM (\x -> cpFile ((show dirs) ++ "/") x) names

  print "End"
