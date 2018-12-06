{- Math.Sparse.Modularity.Eigen.FeatureMatrix
Gregory W. Schwartz

Collects the functions pertaining to finding the Newman-Girvan modularity of a
sparse adjacency matrix.
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE TupleSections #-}

module Math.Modularity.Eigen.Sparse
    ( getModularity
    , getBModularity
    , Q (..)
    , testModularity
    ) where

-- Remote
import Data.Bool (bool)
import Math.Clustering.Spectral.Eigen.FeatureMatrix (B (..), getB)
import qualified Data.Eigen.SparseMatrix as S
import qualified Data.Vector.Storable as VS

-- Local
import Math.Modularity.Types

type LabelVector     = S.SparseMatrixXd
type AdjacencyMatrix = S.SparseMatrixXd

-- | Find modularity from a vector of community labels (0 or 1) corresponding to
-- rows in the adjacency matrix. Needs 0s on the diagonal for the adjacency
-- matrix.
getModularity :: LabelVector -> AdjacencyMatrix -> Q
getModularity moduleVec mat = Q $ (1 / (2 * m)) * sumQ mat
  where
    sumQ :: S.SparseMatrixXd -> Double
    sumQ = S.getSum
         . S._imap (\ i j v -> inner i j v * delta i j)
    inner v w x = x - ((k v * k w) / (2 * m))
    delta v w = ((s v * s w) + 1) / 2
    m = (/ 2) . S.getSum $ mat -- Symmetric matrix so divide by 2.
    d = S.getColSums mat
    s = bool (-1) 1 . (== 0) . (S.!) moduleVec . (,0)
    k = (S.!) d . (0,)

-- | Find modularity from a vector of community labels (0 or 1) corresponding to
-- rows in the normalized matrix B. See Shu et al., "Efficient Spectral
-- Neighborhood Blocking for Entity Resolution", 2011.
-- L = sum_i^n sum_j^n A(i,j) - n = 1^TA1 - n = (B^T1)^T(B^T1) - n.
getBModularity :: LabelVector -> B -> Q
getBModularity moduleVec (B b) = Q . sum . fmap inner $ [first, second]
  where
    inner v = (a v v / l) - ((a v (S.ones n) / l) ** 2)
    first  = moduleVec
    second = S.fromDenseList
           . (fmap . fmap) (bool 1 0 . (== 1))
           . S.toDenseList
           $ moduleVec
    l    = a (S.ones n) (S.ones n)
    a :: S.SparseMatrixXd -> S.SparseMatrixXd -> Double
    a oneL oneR = ( flip (S.!) (0, 0)
                  $ (S.transpose (partA oneL)) * (partA oneR)
                  )
                - (S.getSum oneL)
    partA one = (S.transpose b) * one
    n    = S.rows b

-- | Set the diagonal of a sparse matrix to 0.
setDiag0 :: S.SparseMatrixXd -> S.SparseMatrixXd
setDiag0 = S._imap (\x y z -> if x == y then 0 else z)

-- | Test whether getModularity BB^T is the same as getBModularity B.
testModularity :: (Bool, Q, Q)
testModularity = (modA == modB, modA, modB)
  where
    items = S.fromDenseList (fmap (:[]) [1,1,0,0] :: [[Double]])
    b     = getB True $ S.fromDenseList ([[1,1],[0,0],[0,0],[1,1]] :: [[Double]])
    a     = setDiag0 $ (unB b) * S.transpose (unB b)
    modA  = getModularity items a
    modB  = getBModularity items b
