
import System.Directory
import Data.List as List
import Data.Maybe as Maybe
import Data.List.Split as S

import Data.Time
import System.Process
import System.Environment

import MSUtils

-- Returns all the xyz files
returnXYZFiles :: [String] -> [String]
returnXYZFiles inFiles = geomFiles
  where
    -- Screen "." and ".."
    actualFiles = filter (\x -> (x /= ".") && (x /= "..")) inFiles
    -- Files names splitted by "."
    splitted    = map (splitAtDotBackward) actualFiles
    -- Fetch geom files    
    geomFiles  = List.nub $ map (\x -> fst x) $ filter (\x -> snd x == ".xyz") $ splitted

-- Generate input file
genInputFile :: String -> String -> FilePath -> IO ()
genInputFile header label geom = do
  geoms  <- readFile $ geom ++ ".xyz"
  let
    --fileName  = geom ++ label ++ ".inp"
    fileName  = (foldl (\x y -> if x == "" then x++ y else x ++ "-" ++ y) "" $ S.splitOn "." $ geom ++ label) ++ ".inp"
    modifiedG = modGeom geoms
  putStr $ " fileName=" ++ fileName ++ " .. "
  writeFile fileName (header ++ "\n" ++ modifiedG) 
  putStrLn "done"

-- Modify the xyz file so thta it fits the orca input
-- ex) inGeo =
--             "1"
--             "0 2"
--             "Al         0.000000000      0.000000000      0.000000000"
--     outGeo =
--            "* xyz 0 2"
--            "    Al         0.000000000      0.000000000      0.000000000"
--            "*"  
modGeom :: String -> String
modGeom inGeo = unlines newDecomp
  where
    decomposed = tail $ lines inGeo
    myHead     = head decomposed -- 0 2
    myGeos     = tail decomposed -- Al 0 0 0
    myTail     = "*"
    newDecomp  = ("* xyz " ++ myHead) : (map (\x -> "    " ++ x) myGeos) ++ [myTail]
   
main = do
  -- Get the target directory
  args <- getArgs
  -- Get contents of the current repository
  if length args == 2 then do
    putStrLn $ "Hearder file     ... " ++ (show (args !! 0))     
    repo <- setCurrentDirectory (args !! 1)
    putStrLn $ "Target directory ... " ++ (show (args !! 1))
    else do
    putStrLn ""
    putStrLn " This program generates a set of orca input files with the same wavefunction / basis set settings for different geometries"
    putStrLn ""    
    putStrLn " Usage:"
    putStrLn "       GenInput HeaderFile GeomDirec"    
    putStrLn " "
    putStrLn " Arguments:"
    putStrLn "       HeaderFile - File which contains header part of the input"
    putStrLn "       GeomDirec  - Path to the directory which contains the xyz file"
    putStrLn ""
    error "Abort further execution .. "

  -- Names of the header file and path to the geometry directory
  let headerName = args !! 0
      geomDir    = args !! 1

  -- Get current time to make foot note
  time <- getZonedTime
  let
    (LocalTime d t) = zonedTimeToLocalTime time
    myFooter = ".gen." ++ (show t) ++ "-" ++ (show d)
    --myFooter = foldl (\x y -> if x == "" then x++ y else x ++ "_" ++ y) "" $ S.splitOn "." $ ".gen." ++ (show t) ++ "-" ++ (show d)
  
  -- Open the header file
  headers <- readFile headerName
  putStrLn " Hearder files .. \n"
  putStrLn " ==================================="
  putStr headers
  putStrLn " ==================================="
  putStrLn ""  
  -- Get contents of the current repository
  xyzDirs  <- getCurrentDirectory
  allFiles <- getDirectoryContents xyzDirs
  let xyzFiles = returnXYZFiles allFiles
  
  -- Generate the input files
  mapM (\x -> genInputFile headers myFooter x) xyzFiles

  putStrLn " End"
