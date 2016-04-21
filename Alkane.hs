import System.Directory
import Data.List as List

import Data.Time
import System.Process
import System.Environment

import AlkUtils
        
main = do
  -- Get the chain length
  args <- getArgs
  let clength = read (args !! 0) :: Integer
      mypoly = printAlkane False clength
  putStrLn mypoly
  putStrLn "End .. "
      
