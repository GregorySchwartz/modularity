{- Math.Modularity.Types
Gregory W. Schwartz

Collects the types used in modularity.
-}

{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Math.Modularity.Types where

-- Remote

-- Local

newtype Q = Q { unQ :: Double } deriving (Eq, Ord, Read, Show, Num)
