
import Data.Time.Clock
import Data.Time.Clock.POSIX

main = do
 now <- getCurrentTime
 putStrLn $ show $ posixSecondsToUTCTime $ utcTimeToPOSIXSeconds now
 
