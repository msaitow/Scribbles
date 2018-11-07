import System.Directory
import Data.List as List

import Data.Time
import System.Process
import System.Environment

import AlkUtils

main = do
  let lengths = filter (\x -> mod x 5 == 0) [5..100]
  mapM (\x -> printLiakosLinearSing x) lengths

  putStrLn " End"
  
