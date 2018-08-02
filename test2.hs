
import System.Directory
import System.Environment

main=do
  args <- getArgs
  -- Get contents of the current repository
  if length args > 0 then do
    repo <- setCurrentDirectory (args !! 0)
    putStrLn "Set .. \n"
    else do
    repo <- setCurrentDirectory "./"
    putStrLn "Default .. \n"
    
  dirs <- getCurrentDirectory
  cons <- getDirectoryContents dirs
  putStrLn "CurrentDirectory"
  print cons
  
  putStrLn "The arguments are:"
  mapM putStrLn args

  
