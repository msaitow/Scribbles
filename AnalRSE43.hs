
import System.Directory
import System.Environment

--
import Data.Maybe as Maybe
import Data.List as List
import Text.Printf as Printf
import MSUtils
--

import AnalUtils

-- Input : file name
-- Output: (Label, Etot)
-- where Label can be something like E1, E2, .., P1, P2, ..
getEtot :: String -> IO (String, String)
getEtot myName = do
  cont <- readFile myName
  let keyword = fst $ splitAtDotForward myName
      ene     = extractEtot cont
  return (keyword, ene)

-- Input:  unlined contents of the file
-- Output: E(TOT) ... XXX
---  XXX will be the output
extractEtot :: String -> String
extractEtot contents = etot
  where
    tokenC   = fmap words $ lines contents
    etotcand = filter (\x ->  "E(TOT)" `elem` x) tokenC 
    etot     = if length etotcand == 0 then "N.F." else last $ head etotcand

compileTableRSE43 :: [(String, String)] -> [(Maybe Double, Maybe Double, Maybe Double, Maybe Double, Maybe Double)]
compileTableRSE43 inData = fmap calcReactionEnergy arow
  where
    dData      = fmap (\x -> convert2Double x) inData
    numbers    = filter (/=14) [2..45] -- [E_x + P1 = E1 + P_x | x /= 14, x <- [2..45]]
                                       -- somehow there's no 14th .. 
    arow       = fmap constructOneRow numbers
    
    convert2Double :: (String, String) -> (String, Maybe Double)
    convert2Double (x, y) = (x, z)
      where z = if y == "N.F." then Nothing else Just $ read y

    constructOneRow :: Integer -> (Maybe Double, Maybe Double, Maybe Double, Maybe Double)
    constructOneRow num = (findEnergy keyEx, findEnergy "P1", findEnergy "E1", findEnergy keyPx) -- (Ex, P1, E, Px)
      where
        keyEx = "E" ++ (show num)
        keyPx = "P" ++ (show num)

        findEnergy :: String -> Maybe Double
        findEnergy key = if length candidates == 0 then Nothing else snd $ head candidates
          where candidates = filter (\x -> fst x == key) dData

    calcReactionEnergy :: (Maybe Double, Maybe Double, Maybe Double, Maybe Double) -> (Maybe Double, Maybe Double, Maybe Double, Maybe Double, Maybe Double)
    calcReactionEnergy (a, b, c, d)
      | a == Nothing = (a, b, c, d, Nothing)
      | b == Nothing = (a, b, c, d, Nothing)
      | c == Nothing = (a, b, c, d, Nothing)
      | d == Nothing = (a, b, c, d, Nothing)
      | otherwise    = (a, b, c, d, ereac  )
      where
        ereac        = Just $ hartree2kcal * (- (Maybe.fromJust a) - (Maybe.fromJust b) + (Maybe.fromJust c) + (Maybe.fromJust d))
        
printTable :: [(Maybe Double, Maybe Double, Maybe Double, Maybe Double, Maybe Double)] -> IO ()
printTable inp = do
  putStrLn ""
  putStrLn " --------------------------------------------------------------"
  putStrLn "                Etotal / Eh                                    "
  putStrLn " ------------------------------------------"
  putStrLn "    Ex(S)       P1(D)      E1(S)     Px(D)     Ereac / Kcal/mol"
  putStrLn " --------------------------------------------------------------"
  mapM_ printRow inp
  putStrLn " --------------------------------------------------------------"
  putStrLn ""


printRow :: (Maybe Double, Maybe Double, Maybe Double, Maybe Double, Maybe Double) -> IO ()
printRow (ex, p1, e1, px, ereac) = do
  printOne ex
  printOne p1
  printOne e1
  printOne px
  printOne ereac
  putStrLn ""

printOne :: Maybe Double -> IO ()
printOne elem = do
  if elem /= Nothing then Printf.printf " %9.2f " (Maybe.fromJust elem)
                     else Printf.printf " %9s "   ("N.F.")

printOne' :: Maybe Double -> IO ()
printOne' elem = do
  if elem /= Nothing then Printf.printf " (%9.2f)\n " (Maybe.fromJust elem)
                     else Printf.printf " (%9s)\n "   ("N.F.")

-- esref :: reference energies
-- es    :: dlpno energies for instance  
printErrors :: [Maybe Double] -> [Maybe Double] -> IO ()
printErrors erefs es = do
  --putStrLn " Errors .. "
  let errors = zipWith (subMaybe) es erefs
  mapM (printOne')  errors
  putStrLn ""

subMaybe :: Maybe Double -> Maybe Double -> Maybe Double
subMaybe a b
  | a == Nothing = Nothing
  | b == Nothing = Nothing
  | otherwise    = Just $ (Maybe.fromJust a) - (Maybe.fromJust b)

-- Print max deviation and RMSD in Kcal/mol
-- canonical, lpno, dlpno  
calc_printRMSD :: [Maybe Double] -> [Maybe Double] -> [Maybe Double] -> IO ()
calc_printRMSD cData lData dData = do
  let
    dLPNOs  = calcRMSD $ fmap (\x -> Maybe.fromJust x) $ filter(\x -> x /= Nothing) $ zipWith calcError3 lData cData 
    dDLPNOs = calcRMSD $ fmap (\x -> Maybe.fromJust x) $ filter(\x -> x /= Nothing) $ zipWith calcError3 dData cData 

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

calcError3 :: Maybe Double -> Maybe Double -> Maybe Double
calcError3 d1 d2
  | d1 == Nothing = Nothing
  | d2 == Nothing = Nothing
  | otherwise     = Just $ (Maybe.fromJust d2) - (Maybe.fromJust d1) -- d1 will be canonical number ..

main = do
  -- Get the target directory
  args <- getArgs
  -- Get contents of the current repository
  if length args > 0 then do
    repo <- setCurrentDirectory $ (args !! 0) ++ "/canon"
    putStrLn "Target directory set \n"
    else do
    putStrLn "This program reads all the orca output files for the RSE43 sets in the given target directory and its subdirectory and then it prints the reaction energies .. "
    putStrLn ""
    putStrLn "Argument: path to directory which should contain three subdirectories, canon/, lpno/ and dlpno/. All the output files located in each subdirectories"
    putStrLn "          are read and the LPNO and DLPNO truncation errors are calculaed. Note that here errors in total correlation energies are analyzed."    
    putStrLn ""    
    error "Aborting further operation"

  -- (1) Canonical results
  -- Get contents of the current repository
  dirs <- getCurrentDirectory
  cons <- getDirectoryContents dirs
  let outFilesCanon = returnOutFiles cons
  results_canon <- mapM getEtot outFilesCanon
  let tableCanon  = compileTableRSE43 results_canon
  let ereacsC = fmap (\(a, b, c, d, e) -> e) tableCanon
  
  -- (2) LPNO results
  repo <- setCurrentDirectory $ (args !! 0) ++ "/lpno"
  -- Get contents of the current repository
  dirs <- getCurrentDirectory
  cons <- getDirectoryContents dirs
  let outFilesLPNO = returnOutFiles cons
  results_lpno <- mapM getEtot outFilesLPNO
  let tableLPNO  = compileTableRSE43 results_lpno  
  let ereacsL = fmap (\(a, b, c, d, e) -> e) tableLPNO
      
  -- (3) DLPNO results
  repo <- setCurrentDirectory $ (args !! 0) ++ "/dlpno"
  -- Get contents of the current repository
  dirs <- getCurrentDirectory
  cons <- getDirectoryContents dirs
  let outFilesDLPNO = returnOutFiles cons
  results_dlpno <- mapM getEtot outFilesDLPNO
  let tableDLPNO  = compileTableRSE43 results_dlpno
  let ereacsD = fmap (\(a, b, c, d, e) -> e) tableDLPNO
      
  Printf.printf "Found %d entries in /canon ... \n" (length tableCanon)
  printTable tableCanon

  Printf.printf "Found %d entries in /lpno .. \n" (length tableLPNO)
  printTable tableLPNO

  Printf.printf "Found %d entries in /dlpno .. \n" (length tableDLPNO)
  printTable tableDLPNO

  putStrLn "LPNO errors .. \n"
  printErrors ereacsC ereacsL
    
  putStrLn "DLPNO errors .. \n"
  printErrors ereacsC ereacsD
  
  calc_printRMSD ereacsC ereacsL ereacsD
  
  print "End"

  
