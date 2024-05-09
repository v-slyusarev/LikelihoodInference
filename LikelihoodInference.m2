-- -*- coding: utf-8 -*-
newPackage(
    "LikelihoodInference",
    Version => "0.0.1",
    Date => "May 1, 2024",
    Authors => {
      {Name => "Vladislav Sliusarev", Email => "vnvdvc@nmsu.edu", HomePage => "https://www.linkedin.com/in/vladislav-sliusarev/"}
    },
    Headline => "Algebraic methods for Maximum Likelihood Estimates and Likelihood Ratio Tests",
    Keywords => {"Statistics", "Likelihood", "Probability"},
    DebuggingMode => false
)

export {
  "monomialParametrization",
  "expRandomCensoringMLE", "lambdaS",
  "independenceModelMLE",
  "iterativeProportionalScalingMLE", "err",
  "MultivariateNormalDistribution", "multivariateNormalDistribution","meanVector", "covarianceMatrix",
  "saturatedGaussianModelMLE",
  "parametricCurveGaussianModelMLE",
  "likelihoodIdeal",
  "implicitModelMLE"
}

needsPackage "NumericalAlgebraicGeometry";
needsPackage "RealRoots";

---------------------
-- Private methods --
---------------------

zipWithOperator := method(TypicalValue => List);
zipWithOperator (List, List, Function) := List => (xs, ys, f) -> (
  if (#xs != #ys) then error("Zipped lists must have equal number of elements");
  return toList((0..#xs-1) / (i -> f(xs#i, ys#i)));
);

norm1 := method(TypicalValue => RR);
norm1 List := RR => xs -> sum((toList xs) / abs);
norm1 Matrix := RR => m -> norm1 flatten entries m;
norm1 Vector := RR => v -> norm1 entries v;

validateLogLinearModelData := method();
validateLogLinearModelData(Matrix, List) := (A,u) -> (
  d := numgens target A;
  k := numgens source A;
  
  if (#u != k) then error(
    "Dimension mismatch. Expected an (d*k)-matrix A and a k-vector u, but the dimensions of A are (" | d | "," | k | ") and the dimension of u is " | (toString(#u)) | "."
  );

  a := norm1 A_{1};

  if any(entries transpose A, col -> norm1(col) != a) then error(
    "All column sums of A must be equal."
  );
)

--------------------------------------------------------
-- Maximum Likelihood Estimate For Independence Model --
--------------------------------------------------------

independenceModelMLE = method(TypicalValue => Sequence);
independenceModelMLE(Matrix) := u -> (
  totalCounts := sum(flatten entries u);

  rows := entries u;
  rowSums := rows / sum;
  rowDistribution := rowSums / (x -> x / totalCounts);

  columns := entries transpose u;
  columnSums := columns / sum;
  columnDistribution := columnSums / (x -> x / totalCounts);

  return (rowDistribution, columnDistribution);
);

-------------------------------------------------------
-- Monomial Parametrization For Affine Toric Variety --
-------------------------------------------------------

monomialParametrization = method(TypicalValue => List);
monomialParametrization(Matrix, Ring) := (A, R) -> (
  variables := gens R;

  if(numgens target A != #variables) then error(
    "The number of variables must match the number of columns, but the given matrix has " | numgens target A | " columns and there are " | #variables | " variables."
  );

  columns := entries transpose A;

  return columns / (a -> product(zipWithOperator(variables, a, power)));
)


---------------------------------------------------------------------------------------
-- Maximum Likelihood Estimate For Exponential Distribution Random Censoring Problem --
---------------------------------------------------------------------------------------

expRandomCensoringMLE = method(TypicalValue => List, Options => true);
expRandomCensoringMLE(ZZ, ZZ, ZZ, ZZ) := { lambdaS => 1 } >> opt -> (u0,u1,u2,u12) -> (
  
  l := symbol l;
  R := QQ[l_1, l_2];

  f1 := (u1+u12)*(l_1+l_2+2)*(l_1+1)*(l_1+l_2+1) +
          (u12)*l_1*(l_1+1)*(l_1+l_2+1) -
          (u2+u12)*l_1*(l_1+l_2+2)*(l_1+l_2+1) -
          (u0+u1+u2+u12)*l_1*(l_1+l_2+2)*(l_1+1);

  f2 := (u2+u12)*(l_1+l_2+2)*(l_2+1)*(l_1+l_2+1) +
          (u12)*l_2*(l_2+1)*(l_1+l_2+1) -
          (u1+u12)*l_2*(l_1+l_2+2)*(l_1+l_2+1) -
          (u0+u1+u2+u12)*l_2*(l_1+l_2+2)*(l_2+1);

  g := l_1*l_2*(l_1+1)*(l_2+1)*(l_1+l_2+1)*(l_1+l_2+2);

  I := ideal{f1, f2};
  J := ideal{g};
  K := saturate(I,J);

  allSolutions := solveSystem flatten entries gens K;
  positiveSolutions := select(allSolutions, point -> all(coordinates point, value -> value > 0));
  scaledSolutions := positiveSolutions / (point -> (coordinates point) / (x -> x * opt.lambdaS));
  return scaledSolutions;
);

----------------------------------------------
-- Iterative Proportional Scaling Algorithm --
----------------------------------------------

iterativeProportionalScalingMLE = method(TypicalValue => List, Options => true);
iterativeProportionalScalingMLE(Matrix, List) := { err => .01 } >> opt -> (A, u) -> (
  validateLogLinearModelData(A, u);

  d := numgens target A;
  k := numgens source A;
  
  a := norm1 A_{1};

  v := new MutableList from (k:(norm1(u) / k));

  Au := A * vector(u);
  Av := A * vector(toList v);

  theta := symbol theta;
  R := QQ[theta_1..theta_d];

  phi := monomialParametrization(A, R);
  phiAu := phi / (monomial -> sub (monomial, transpose matrix Au));

  while norm1(Av - Au) > opt.err do(
    phiAv := phi / (monomial -> sub (monomial, transpose matrix Av));
    
    for i from 0 to k-1 do(
      ratio := ((phiAu#i)/(phiAv#i))^(1/a);
      v#i = v#i * ratio;
    );

    Av = A * vector(toList v);
  );


  return toList v;
);

--------------------------------------
-- Multivariate Normal Distribution --
--------------------------------------

MultivariateNormalDistribution = new Type of HashTable;

meanVectorFieldName := "meanVector";
covarianceMatrixFieldName := "covarianceMatrix";

multivariateNormalDistribution = method(TypicalValue => MultivariateNormalDistribution);
multivariateNormalDistribution(Vector,Matrix) := (mu,Sigma) -> (
  if (transpose Sigma != Sigma) then error(
    "Expected a symmetric matrix but received:\n"|Sigma
  );

  if (#(entries mu) != numgens source Sigma) then error(
    "The dimension of the vector and the matrix must be the same."
  );

  return new MultivariateNormalDistribution from (
    new HashTable from {
      meanVectorFieldName => mu,
      covarianceMatrixFieldName => Sigma
    }
  );
);

meanVector = method(TypicalValue => Vector);
meanVector MultivariateNormalDistribution := N -> (N#meanVectorFieldName);

covarianceMatrix = method(TypicalValue => Matrix);
covarianceMatrix MultivariateNormalDistribution := N -> (N#covarianceMatrixFieldName);


validateVectorSample := method();
validateVectorSample List := X -> (
  dim := #(X#0);

  if (any(X, x -> #x != dim)) then error(
    "All sample vectors must have the same dimension."
  );
);


----------------------------------------------------------
-- Saturated Gaussian Model Maximum Likelihood Estimate --
----------------------------------------------------------

saturatedGaussianModelMLE = method(TypicalValue => MultivariateNormalDistribution);
saturatedGaussianModelMLE(List) := X -> (
  validateVectorSample(X);

  n := #X;

  meanVector := (1/n) * sum( X / vector);
  deviations := X / (x -> vector(x) - meanVector);
  covarianceMatrix := (1/n) * sum(deviations / (x -> (matrix x) * (transpose matrix x)));
  return multivariateNormalDistribution(meanVector, covarianceMatrix);
)

--------------------------------------------------------------------------
-- Rational Parametric Curve Gaussian Model Maximum Likelihood Estimate --
--------------------------------------------------------------------------

solveForReals := (f, err) -> (
  intervals := realRootIsolation(f, err);
  return intervals / (interval -> numeric (interval#0 + interval#1)/2);
)

parametricCurveGaussianModelMLE = method(TypicalValue => MultivariateNormalDistribution, Options => true);
parametricCurveGaussianModelMLE(List, List) := { err => .01 } >> opt -> (g,X) -> (
  validateVectorSample(X);

  m := #(X#0);

  if (#g != m) then error(
    "The parametrization and the sample must have the same dimension."
  );

  ts := gens ring g#0;

  if (any(g, gi -> gens ring gi != ts)) then error(
    "The components of the parametrization must be polynomials in the same variable."
  );

  if (#ts != 1) then error(
    "The parametrization must be a list of polynomials in one variable."
  );

  d := degree g#0;

  if (any(g, gi -> degree gi != d)) then error(
    "All components of the parametrization must have the same degree."
  );

  n := #X;
  g = vector g;
  t := ts#0;

  meanX := (1/n) * sum(X / vector);
  h := matrix(g - meanX);
  distance := trace ((transpose h) * h);
  derivativeOfDistance := diff(t, distance);
  solutions := solveForReals(derivativeOfDistance,opt.err);
  solutionDistances := solutions / (t0 -> sub(distance, {t=>t0}));

  meanVector := sub(g, {t=>solutions#(maxPosition solutionDistances)} );

  covarianceMatrix := id_(QQ^m);

  return multivariateNormalDistribution(meanVector, covarianceMatrix);
)

--------------------------------------------
-- Likelihood Ideal For An Implicit Model --
--------------------------------------------

augmentedJacobian := method(TypicalValue => Matrix);
augmentedJacobian Ideal := P -> (
  R := ring P;
  return (matrix {gens R}) || ((transpose jacobian P) * (diagonalMatrix vars R));
);

likelihoodIdeal = method(TypicalValue => Ideal);
likelihoodIdeal (List, List) := (g,u) -> (
  P := ideal g;
  r := #g;
  R := ring P;

  n := #u;

  if(#(gens R) != n) then error(
    "The number of variables must match the dimension of the data vector."
  );

  -- Step 1
  c := codim P;

  -- Step 2
  Jaug := augmentedJacobian(P);

  G := fold(
    (m1,m2)->m1|m2,
    entries id_(R^(r+1)) / (column -> fold(
          (m1,m2)->m1|m2,
          #g:(transpose matrix {column}) * diagonalMatrix(g)
      )
    )
  );
  

  JaugG := Jaug | G;
  generatorsOfM := (entries transpose gens ker JaugG) / (column -> take(column, n));

  -- Step 3
  generatorsOfI'u := generatorsOfM / (phi -> sum zipWithOperator(u, phi, times));
  I'u := ideal generatorsOfI'u;

  -- Step 4
  Q := minors(c+1, Jaug);
  Iu := saturate(I'u, (sum gens R)*(product gens R)*Q);
  return Iu;
)

------------------------------------------------
-- Implicit Model Maximum Likelihood Estimate --
------------------------------------------------
implicitModelMLE = method(TypicalValue => List);
implicitModelMLE(List, List) := (g,u) -> (
  R := ring ideal g;

  Iu := likelihoodIdeal(g,u);
  cirticalPointEquations := (flatten entries transpose gens Iu) | {(sum gens R) -1 };
  allSolutions := (solveSystem cirticalPointEquations) / (point -> (coordinates point)); 
  positiveSolutions := select(allSolutions, coordinates -> all(coordinates, value -> value > 0));


  Jaug := augmentedJacobian ideal g;

  S := RR[gens R];
  Lp := (product zipWithOperator(gens S, u, power)) / ((sum gens S)^(sum u));
  -- Lp := sub(Lp, matrix({ (gens R) / (p -> p_S)}));

  logLp := log Lp;

  << Lp << endl;

  return positiveSolutions;

)