
import System.Directory
import Data.List as List
import Data.Maybe as Maybe

import Data.Time
import System.Process
import System.Environment
    
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
    fileName  = geom ++ label ++ ".inp"
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
    repo <- setCurrentDirectory "./"
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
  
  -- Open the header file
  headers <- readFile headerName
  putStrLn " Hearder files .. \n"
  putStrLn " ==================================="
  putStr headers
  putStrLn " ==================================="
  
  -- Get contents of the current repository
  xyzDirs  <- getCurrentDirectory
  allFiles <- getDirectoryContents xyzDirs
  let xyzFiles = returnXYZFiles allFiles
  --print xyzFiles
  
  -- Generate the input files
  mapM (\x -> genInputFile headers myFooter x) xyzFiles
  
--  -- Get file names to be changed
--  let
--    fileNames = returnFilesToBeRenamed cons
--    newNames  = renameFiles myFooter fileNames
--
--  if False then do
--    putStr " Not moved"
--    else do
--    -- (1) Rename previous ouput files
--    moved <- mapM (\x -> moveFile ((show dirs) ++ "/") x) newNames 
--    print " Moved"
--    
  putStrLn " End"
--  print fileNames
--  print "~"
--  print newNames
  
