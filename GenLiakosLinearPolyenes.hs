import System.Directory
import Data.List as List

import Data.Time
import System.Process
import System.Environment

import AlkUtils

main = do
  let lengths = filter (\x -> mod x 5 == 0) [1..150]
  mapM (\x -> printLiakosLinearTri x) lengths

  putStrLn " End"
  
