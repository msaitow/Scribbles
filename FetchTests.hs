
import System.Directory
import System.Process
import System.Environment

testDir    = "/home1/saitow/orca/testsuite/InputFiles"
numMyTests = [717..734]

copyFiles :: Integer -> IO ()
copyFiles fileNum = callCommand command
  where command = "cp " ++ testDir ++ "/ORCA_Test_" ++ (show fileNum) ++ ".inp" ++ " ./"
  
main = do
  mapM (\x -> copyFiles x) numMyTests
  putStrLn " End"

  
    
