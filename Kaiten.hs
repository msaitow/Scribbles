
import System.Directory
import Data.List as List
import Data.Maybe as Maybe
import Text.Printf as Printf

import Data.Time
import System.Process
import System.Environment

import MSUtils

-- Generate input file
genGeom :: String -> Int -> IO ()
genGeom geom n = do
  geoms  <- readFile geom
  let
    tag       = fst $ splitAtDotBackward geom
    decomp    = lines geoms
    numAtoms :: Double
    numAtoms  = (read (decomp !! 0))
    charge   :: Double    
    charge    = (read ((words (decomp !! 1)) !! 0))
    mult     :: Double
    mult      = (read ((words (decomp !! 1)) !! 1))
    body      = readMolStructure $ unlines $ tail $ tail $ lines geoms
    mol       = Molecule body charge mult numAtoms

    -- Rotation angles
    d_theta = 2 * pi / (fromIntegral n)
    angles  = take (n-1) $ map (\x -> d_theta * x) [1.0,2.0..]
    d_zaxis = 3.0
    dist    = take (n-1) $ map (\x -> d_zaxis * x) [1.0,2.0..]
    
    -- concatenated rotated structures
    multimer = foldl (highSpinCoupling) mol $ zipWith (\x y -> transMol (kaitenMol mol x) y) angles dist     
  
  --putStr $ show multimer
  writeFile (tag ++ ".rot" ++ (show n) ++ ".xyz") (show multimer)
  putStrLn "done"

readMolStructure :: String -> [AtomCoord]
readMolStructure inLine = map readline decomposed
  where
    decomposed = map words $ lines inLine    
    readline :: [String] -> AtomCoord
    readline aLine = AtomCoord (aLine !! 0) (read (aLine !! 1)) (read (aLine !! 2)) (read (aLine !! 3))
  
-- Coordination of an atom
-- Eg. "X" 0.0 0.1 0.2    
data AtomCoord = AtomCoord { aLabel :: String,
                             aX     :: Double,
                             aY     :: Double,
                             aZ     :: Double }

instance Show AtomCoord where
  show = showAtom

showAtom :: AtomCoord -> String
showAtom atom = Printf.printf " %4s %12.6f %12.6f %12.6f " label x y z
  where
    label = aLabel atom
    x     = aX atom
    y     = aY atom
    z     = aZ atom
  
data Molecule = Molecule { molCoord :: [AtomCoord],
                           charge   :: Double,
                           mult     :: Double,
                           nAtom    :: Double}

instance Show Molecule where
  show mol = (Printf.printf "  %3.0f\n  %3.0f %3.0f\n" (nAtom mol) (charge mol) (mult mol)) ++ (foldl (++) "" $ map (\x -> show x ++ "\n") $ molCoord mol)

highSpinCoupling :: Molecule -> Molecule -> Molecule
highSpinCoupling mol1 mol2 = Molecule atomList charge3 mult3 nAtom3
  where    
    atomList = (molCoord mol1) ++ (molCoord mol2)
    charge3  = (charge mol1) + (charge mol2)
    nAtom3   = (nAtom mol1) + (nAtom mol2)
    nsomo1   = ((mult mol1) - 1) / 2.0 
    nsomo2   = ((mult mol2) - 1) / 2.0
    mult3    = 2 * (nsomo1+nsomo2) + 1
    
-- KatenMol
kaitenMol :: Molecule -> Double -> Molecule
kaitenMol inMol theta = Molecule (map (\x -> kaiten_N x theta) body)  (charge inMol) (mult inMol) (nAtom inMol)
  where
    body = molCoord inMol
      
-- TransMol
transMol :: Molecule -> Double -> Molecule
transMol inMol theta = Molecule (map (\x -> trans_Z x theta) body)  (charge inMol) (mult inMol) (nAtom inMol)
  where
    body = molCoord inMol

-- Rotates the given coordination by theta angle around Z-axis
kaiten_N :: AtomCoord -> Double -> AtomCoord
kaiten_N inAtom theta = outAtom
  where
    outAtom = AtomCoord (aLabel inAtom) rX rY (aZ inAtom)
    rX      = (cos theta) * (aX inAtom) - (sin theta) * (aY inAtom)
    rY      = (sin theta) * (aX inAtom) + (cos theta) * (aY inAtom)

-- Translate the atom along with Z-axis
trans_Z :: AtomCoord -> Double -> AtomCoord
trans_Z inAtom n = outAtom
  where
    outAtom = AtomCoord (aLabel inAtom) (aX inAtom) (aY inAtom) ((aZ inAtom) + n)

main = do
  -- Get the target directory
  args <- getArgs
  -- Get contents of the current repository
  if length args == 2 then do
    putStrLn $ "Target geometry    ... " ++ (show (args !! 0))     
    --repo <- setCurrentDirectory (args !! 1)
    putStrLn $ "Number of monomers ... " ++ (show (args !! 1))
    else do
    putStrLn ""
    putStrLn " This program generates xyz files which contains multimer of given molecule"
    putStrLn ""    
    putStrLn " Usage:"
    putStrLn "       Kaiten xyzFile_of_monomer Number_of_monomers_in_multimer"    
    putStrLn " "
    putStrLn ""
    error "Abort further execution .. "

  -- Names of the header file and path to the geometry directory
  let xyzFile  = args !! 0
      numMonos = args !! 1

  -- List of integers
  let nums = [1..(read numMonos)]
      
  mapM (\x -> genGeom xyzFile x) nums
  putStrLn " End"
