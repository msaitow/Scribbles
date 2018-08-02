
module AlkUtils
( 
  printElement,
  printAlkane,
  printLiakosLinearTri,
  printLiakosLinearSing  
) where
  
import Text.Printf as Printf

-- Print triplet linear polyene in Liako's format
printLiakosLinearTri :: Integer -> IO ()
printLiakosLinearTri myNum = do
  let myName   = "C.tri." ++ (show myNum) ++ ".xyz"
      myLength = 3 * myNum
      geom     = (Printf.printf " %d \n" myLength) ++ " 0 3 \n" ++ (printAlkane True myNum) :: String
  putStr $ " fileName=" ++ myName ++ " .. "
  writeFile myName geom
  putStrLn " done"

-- Print singlet linear polyene in Liako's format
printLiakosLinearSing :: Integer -> IO ()
printLiakosLinearSing myNum = do
  let myName   = "C.sing." ++ (show myNum) ++ ".xyz"
      myLength = 3 * myNum + 2
      geom     = (Printf.printf " %d \n" myLength) ++ " 0 1 \n" ++ (printAlkane False myNum) :: String
  putStr $ " fileName=" ++ myName ++ " .. "
  writeFile myName geom
  putStrLn " done"

-- Print one element ..
printElement :: String -> Double -> Double -> Double -> String
printElement label x y z = Printf.printf " %4s %12.6f %12.6f %12.6f \n" label x y z

-- Retunrn the contents of the xyz files for linear alkane chain .. 
printAlkane :: Bool -> Integer -> String
printAlkane isRadical clength
  | clength == 0 = error "Chain length cannot be zero .. "
  | clength <  0 = error "Chain length cannot be smaller than 0 .. "
  | otherwise    = if isRadical then chain else initH ++ chain ++ lastH
  where
    -- angle between two legs of a tetrahedron
    angle_deg = 109.4712               :: Double
    angle     = angle_deg * pi / 180.0 :: Double
    -- bons lengths (Angs)
    rch       = 1.09                   :: Double
    rcc       = 1.55                   :: Double
    dxh       = cos(angle/2.0) * rch   :: Double
    
    -- the first H atom
    initH     = printElement "H" (-cos(angle/2.0)*rch)  (0.0) (-sin(angle/2.0)*rch)
    lastH     = printElement "H"   xh                    yh     zh 
      where
        xc = if even (clength-1) then 0.0 else -cos(angle/2.0) * rcc :: Double
        yc = 0.0                                                     :: Double
        zc = (fromIntegral $ clength-1) * sin(angle/2.0) * rcc       :: Double
        xh = if even (clength-1) then xc - dxh else xc + dxh         :: Double
        yh = yc                                                      :: Double
        zh = zc + sin(angle/2.0) * rch                               :: Double
                
    cNums     = [0..clength-1]
    ch2s   = fmap (makeCH2) cNums
    chain  = foldr (++) "" ch2s
    
    makeCH2 :: Integer -> String
    makeCH2 myNum = cartC ++ cartH1 ++ cartH2
      where
        -- Coordinate for C atom        
        xc    = if even myNum then 0.0 else -cos(angle/2.0) * rcc :: Double
        yc    = 0.0                                               :: Double
        zc    = (fromIntegral myNum) * sin(angle/2.0) * rcc       :: Double
        cartC = printElement "C" xc yc zc
        -- Coordinate for H atoms
        xh    = if even myNum then xc + dxh else xc - dxh         :: Double
        yh    = sin(angle/2.0) * rch                              :: Double
        zh    = zc                                                :: Double
        cartH1 = printElement "H"  xh   yh  zh
        cartH2 = printElement "H"  xh (-yh)  zh
