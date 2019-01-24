{- Math.Modularity
Gregory W. Schwartz

Collects the functions pertaining to finding the Newman-Girvan modularity of an
adjacency matrix.
-}

module Math.Modularity.Dense
    ( getModularity
    , getBModularity
    , Q (..)
    ) where

-- Remote
import Data.Bool (bool)
import Math.Clustering.Spectral.Dense (B (..), getB)
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

-- | Find modularity from a vector of community labels (0 or 1) corresponding to
-- rows in the normalized matrix B. See Shu et al., "Efficient Spectral
-- Neighborhood Blocking for Entity Resolution", 2011.
-- L = sum_i^n sum_j^n A(i,j) - n = 1^TA1 - n = (B^T1)^T(B^T1) - n.
getBModularity :: LabelVector -> B -> Q
getBModularity moduleVec (B b) = Q . sum . fmap inner $ [first, second]
  where
    inner v = (a v v / l) - ((a v (ones n) / l) ** 2)
    first  = H.fromColumns [moduleVec]
    second = H.fromColumns [H.cmap (bool 1 0 . (== 1)) moduleVec]
    l    = a (ones n) (ones n)
    a :: H.Matrix Double -> H.Matrix Double -> Double
    a oneL oneR = ( flip H.atIndex (0, 0)
                  $ (H.tr (partA oneL)) H.<> (partA oneR)
                  )
                - (H.sumElements oneL)
    partA one = (H.tr b) <> one
    n    = H.rows b
    ones x = (x H.>< 1) [1,1..]
