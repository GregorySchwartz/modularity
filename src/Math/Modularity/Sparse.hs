{- Math.Sparse.Modularity
Gregory W. Schwartz

Collects the functions pertaining to finding the Newman-Girvan modularity of a
sparse adjacency matrix.
-}

{-# LANGUAGE BangPatterns #-}

module Math.Modularity.Sparse
    ( getModularity
    , getBModularity
    , Q (..)
    , testModularity
    ) where

-- Remote
import Data.Bool (bool)
import Math.Clustering.Spectral.Sparse (B (..), getB)
import qualified Data.Sparse.Common as S
import qualified Numeric.LinearAlgebra.Sparse as S
import qualified Data.Vector as V

-- Local
import Math.Modularity.Types

type LabelVector     = S.SpVector Double
type AdjacencyMatrix = S.SpMatrix Double

-- | Find modularity from a vector of community labels (0 or 1) corresponding to
-- rows in the adjacency matrix. Needs 0s on the diagonal for the adjacency
-- matrix.
getModularity :: LabelVector -> AdjacencyMatrix -> Q
getModularity moduleVec mat = Q $ (1 / (2 * m)) * sumQ mat
  where
    sumQ :: S.SpMatrix Double -> Double
    sumQ = sum
         . fmap (\ (!v, !xs)
                  -> sum
                   . fmap (\(!w, !x) -> inner v w x * delta v w)
                   . zip [0,1..]
                   . S.toDenseListSV
                   $ xs
                )
         . zip [0,1..]
         . S.toRowsL
    inner v w x = x - ((k v * k w) / (2 * m))
    delta v w = ((s v * s w) + 1) / 2
    m = (/ 2) . sum $ mat -- Symmetric matrix so divide by 2.
    d = S.sparsifySV . S.vr . fmap sum . S.toRowsL $ mat
    s = bool (-1) 1 . (== 0) . flip S.lookupDenseSV moduleVec
    k = flip S.lookupDenseSV d

-- | Find modularity from a vector of community labels (0 or 1) corresponding to
-- rows in the normalized matrix B. See Shu et al., "Efficient Spectral
-- Neighborhood Blocking for Entity Resolution", 2011.
-- L = sum_i^n sum_j^n A(i,j) - n = 1^TA1 - n = (B^T1)^T(B^T1) - n.
getBModularity :: LabelVector -> B -> Q
getBModularity moduleVec (B b) = Q . sum . fmap inner $ [first, second]
  where
    inner v = (a v v / l) - ((a v (ones n) / l) ** 2)
    first  = S.fromColsL [moduleVec]
    second = S.fromColsL
           . (:[])
           . S.sparsifySV
           . S.vr
           . fmap (bool 1 0 . (== 1))
           . S.toDenseListSV
           $ moduleVec
    l    = a (ones n) (ones n)
    a :: S.SpMatrix Double -> S.SpMatrix Double -> Double
    a oneL oneR = ( flip S.lookupWD_SM (0, 0)
                  $ (S.transposeSM (partA oneL)) S.#~# (partA oneR)
                  )
                - (sum oneL)
    partA one = (S.transposeSM b) S.#~# one
    n    = S.nrows b

-- | Get a column vector of ones.
ones :: Int -> S.SpMatrix Double
ones n = S.fromColsL [S.onesSV n]

-- | Set the diagonal of a sparse matrix to 0.
setDiag0 :: S.SpMatrix Double -> S.SpMatrix Double
setDiag0 mat = S.fromListSM (S.dimSM mat)
             . fmap (\(!x, !y, !z) -> if x == y then (x, y, 0) else (x, y, z))
             . S.toListSM
             $ mat

-- | Test whether getModularity BB^T is the same as getBModularity B.
testModularity :: (Bool, Q, Q)
testModularity = (modA == modB, modA, modB)
  where
    items = S.fromListDenseSV 4 ([1,1,0,0] :: [Double])
    b     = getB $ S.fromListDenseSM 4 ([1,1,0,0,0,0,1,1] :: [Double])
    a     = setDiag0 $ (unB b) S.#~# S.transposeSM (unB b)
    modA  = getModularity items a
    modB  = getBModularity items b
