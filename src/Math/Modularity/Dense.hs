{- Math.Modularity
Gregory W. Schwartz

Collects the functions pertaining to finding the Newman-Girvan modularity of an
adjacency matrix.
-}

module Math.Modularity.Dense
    ( getModularity
    , Q (..)
    ) where

-- Remote
import Data.Bool (bool)
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Numeric.LinearAlgebra as H

-- Local
import Math.Modularity.Types

type LabelVector     = H.Vector Double
type AdjacencyMatrix = H.Matrix Double

-- | Find modularity from a vector of community labels (0 or 1) corresponding to
-- rows in the adjacency matrix. Needs 0s on the diagonal for the adjacency
-- matrix.
getModularity :: LabelVector -> AdjacencyMatrix -> Q
getModularity moduleVec mat = Q $ (1 / (2 * m)) * sumQ mat
  where
    sumQ = V.sum
         . V.imap (\ v
                  -> H.sumElements
                   . VS.imap (\w x -> inner v w x * delta v w)
                  )
         . V.fromList
         . H.toRows
    inner v w x = x - ((k v * k w) / (2 * m))
    delta v w = ((s v * s w) + 1) / 2
    m = (/ 2) . H.sumElements $ mat -- Symmetric matrix so divide by 2.
    d = H.vector . fmap H.sumElements . H.toRows $ mat
    s = bool (-1) 1 . (== 0) . H.atIndex moduleVec
    k = H.atIndex d
