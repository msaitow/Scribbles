
import System.Directory
import Data.List as List
import Data.Maybe as Maybe
import Text.Printf as Printf

import Data.Time
import System.Process
import System.Environment

import MSUtils
    
-- Whether the hydrogens are ligning on which plane
data Plane = XY | YZ | ZX deriving (Show, Eq)

-- A point in 3D
data Point = Point { x :: Float,
                     y :: Float,
                     z :: Float } deriving (Show)

-- A water molecule
data Water = Water { o  :: Point,
                     h1 :: Point,
                     h2 :: Point }

printWater :: Water -> String
printWater water =  (Printf.printf " O  %14.10f %14.10f %14.10f\n" (x atom_o ) (y atom_o ) (z atom_o ))
                 ++ (Printf.printf " H  %14.10f %14.10f %14.10f\n" (x atom_h1) (y atom_h1) (z atom_h1))
                 ++ (Printf.printf " H  %14.10f %14.10f %14.10f\n" (x atom_h2) (y atom_h2) (z atom_h2))
  where
    atom_o  = o  water
    atom_h1 = h1 water
    atom_h2 = h2 water     

instance Show Water where
  show = printWater
  
-- Generate a water
-- theta_h_o_h should be in angle
-- dist_h_ should be in angsatoms
genWaterOnAxis :: Point -> Float -> Float -> Plane -> Water
genWaterOnAxis origin theta_h_o_h dist_h_o myPlane
  | myPlane == XY && onX == True = Water (origin) (Point (x origin + dy) (y origin + dx) (z origin)) (Point (x origin + dy) (y origin - dx) (z origin))
  | myPlane == XY && onY == True = Water (origin) (Point (x origin + dx) (y origin + dy) (z origin)) (Point (x origin - dx) (y origin + dy) (z origin))
  | myPlane == YZ && onY == True = Water (origin) (Point (x origin) (y origin + dy) (z origin + dx)) (Point (x origin) (y origin + dy) (z origin - dx))
  | myPlane == YZ && onZ == True = Water (origin) (Point (x origin) (y origin + dx) (z origin + dy)) (Point (x origin) (y origin - dx) (z origin + dy))
  | myPlane == ZX && onZ == True = Water (origin) (Point (x origin + dx) (y origin) (z origin + dy)) (Point (x origin - dx) (y origin) (z origin + dy))
  | myPlane == ZX && onX == True = Water (origin) (Point (x origin + dy) (y origin) (z origin + dx)) (Point (x origin + dy) (y origin) (z origin - dx))                           
  | otherwise                    = error "Fuck .."
  where
    -- To idicate whether the origin in on which axis
    onX          = (y origin == 0.0) && (z origin == 0.0)
    onY          = (x origin == 0.0) && (z origin == 0.0)
    onZ          = (x origin == 0.0) && (y origin == 0.0)
    -- Whether the point is on the positive side or not
    isPositive
      | onX       = ((x origin) >= 0.0)
      | onY       = ((y origin) >= 0.0)
      | onZ       = ((z origin) >= 0.0)
      | otherwise = False
      
    -- Parameters for the geometry
    theta_in_rad = 2 * pi * theta_h_o_h / 360.0
    phi          = pi / 2.0 - theta_in_rad / 2.0
    dx           = dist_h_o * (cos phi)
    dy           = if isPositive then dist_h_o * (sin phi)
                   else - dist_h_o * (sin phi)

-- From F. Neese Coord. Chem. Rev. 251, 288 (2007).
-- Symbol, bond length, charge and multiplicity
divalent = [
  ("V" , 2.128, "2 4"),
  ("Cr", 2.052, "2 5"),
  ("Mn", 2.192, "2 6"),
  ("Fe", 2.114, "2 5"),
  ("Co", 2.106, "2 4"),
  ("Ni", 2.061, "2 3"),
  ("Cu", 1.964, "2 2")
  ]

trivalent = [
  ("Ti", 2.028, "3 2"),  
  ("V" , 1.992, "3 3"),
  ("Cr", 1.959, "3 4"),
  ("Mn", 1.991, "3 5"),
  ("Fe", 1.995, "3 6"),
  ("Co", 1.873, "3 5")
  ]

-- Angle and distance are from F. Neese Coord. Chem. Rev. 251, 288 (2007).
genMyWater :: Point -> Plane -> Water
genMyWater origin myPlane = genWaterOnAxis origin 107.7 0.9817 myPlane

-- Return structure 
resStructure :: (String, Float, String) -> String
resStructure (name, length, label) =
  " 19\n" ++
  label ++ "\n" ++
  (Printf.printf " %2s %14.10f %14.10f %14.10f\n" name myZero myZero myZero) ++
  waters
  where
    myZero = 0.0 :: Float
    origins = [(Point length 0.0 0.0), (Point (-length) 0.0 0.0),
               (Point 0.0 length 0.0), (Point 0.0 (-length) 0.0),
               (Point 0.0 0.0 length), (Point 0.0 0.0 (-length))]
    coordinations = [XY, XY,
                     YZ, YZ,
                     ZX, ZX]
    waters = concat $ map show $ zipWith genMyWater origins coordinations

genFile :: (String, String) -> IO ()
genFile (name, structure) = writeFile name structure
    
main = do
  -- Names
  let dinames  = map (\x -> x ++ ".div.xyz") $ map (\(x, _, _) ->  x) divalent
      trinames = map (\x -> x ++ ".tri.xyz") $ map (\(x, _, _) ->  x) trivalent
  -- Structures
  let diStruct  = map (\x -> resStructure x) divalent
      triStruct = map (\x -> resStructure x) trivalent
  -- Inputs
  let diInputs  = zipWith (\x y -> (x,y))  dinames  diStruct
      triInputs = zipWith (\x y -> (x,y)) trinames triStruct

  mapM genFile diInputs
  mapM genFile triInputs
    
  putStrLn " Done"
