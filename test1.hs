
import System.Directory

main = do
  d <- getCurrentDirectory
  c <- getDirectoryContents d
  print c
  
