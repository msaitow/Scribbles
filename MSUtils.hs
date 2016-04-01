
module MSUtils
(
  splitAtDotForward,
  splitAtDotBackward
) where

import Data.List as List
import Data.Maybe as Maybe

-- Split the input into forward and afterward in terms of '.'
-- ex) splitAtDot "C90.inp" -> ("C90", ".inp")
-- NOTE: This splits, for example, "C90.loose.inp" into ("C90.", "loose.inp")
splitAtDotForward :: String -> (String, String)
splitAtDotForward x
  | myIndex /= Nothing = splitAt (Maybe.fromJust myIndex) x
  | otherwise          = ("None", "")
  where
    findDotIndex :: String -> Maybe Int
    findDotIndex kore = List.findIndex (== '.') kore

    myIndex = findDotIndex x

-- Splits from the end so that "C90.loose.inp" is splitted into ("C90.loose.", "inp")
splitAtDotBackward :: String -> (String, String)
splitAtDotBackward x
  | myIndex /= Nothing = if (last $ fst splittedPair) == '.'
                         then ((init $ fst splittedPair), "." ++ (snd splittedPair))
                         else  splittedPair
  | otherwise          = ("None", "")
  where
    splittedPair = (reverse $ snd revsplitted, reverse $ fst revsplitted) 
    revsplitted  = splitAt (Maybe.fromJust myIndex) reversed    
    reversed     = reverse x

    findDotIndex :: String -> Maybe Int
    findDotIndex kore = List.findIndex (== '.') kore

    myIndex = findDotIndex reversed
