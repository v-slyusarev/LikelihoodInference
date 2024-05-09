loadPackage "LikelihoodInference";

-------------------------------------------------------------------------
-- Saturated Gaussian Model Maximum Likelihood Estimate
-------------------------------------------------------------------------
N = saturatedGaussianModelMLE({ {1,2},{-1,-2}, {0,0}});
<< meanVector N << "    " << covarianceMatrix N << endl;
-------------------------------------------------------------------------

<< endl << endl << endl;

-------------------------------------------------------------------------
-- Rational Parametric Curve Gaussian Model Maximum Likelihood Estimate
-------------------------------------------------------------------------
R = QQ[t];
N = parametricCurveGaussianModelMLE(
  {(t+1)^2, t^2},
  { {1,2}, {-1,-2}, {1,2}}
);
<< meanVector N << "    " << covarianceMatrix N << endl;
-------------------------------------------------------------------------