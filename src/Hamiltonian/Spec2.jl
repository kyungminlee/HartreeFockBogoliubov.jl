module Spec2

abstract HermitianTensorElement{S <: Number, K}

immutable DiagonalHermitianTensorElement{S <: Real, K, N1}
          <: HermitianTensorElement{S, K}
  amplitude ::S
  sub ::NTuple{N1, K}
  latticepoint ::Vector{Int64}
end

immutable OffdiagonalHermitianTensorElement{S <: Number, K, N1, N2}
          <: HermitianTensorElement{S, K}
  amplitude ::S
  rowsub ::NTuple{N1, K}
  rowlatticepoint ::Vector{Int64}
  colsub ::NTuple{N2, K}
  collatticepoint ::Vector{Int64}
end

type HermitianTensor{S <: Number, K}
  elements::Vector{HermitianTensorElement{S, K}}
end

end
