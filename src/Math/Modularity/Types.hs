{- Math.Modularity.Types
Gregory W. Schwartz

Collects the types used in modularity.
-}

{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE DeriveGeneric #-}

module Math.Modularity.Types where

-- Remote
import GHC.Generics (Generic)

-- Local

newtype Q = Q { unQ :: Double } deriving (Eq, Ord, Read, Show, Num, Generic)
