loadPackage "LikelihoodInference";

-------------------------------------------------------------------------
-- Iterative Proportional Scaling Maximum Likelihood Estimate
-------------------------------------------------------------------------
A = matrix({
  {2,1,0},
  {1,2,3}
});
u = {10,20,30};
<< iterativeProportionalScalingMLE(A, u, err=>0.001) << endl;
-------------------------------------------------------------------------