## Overview
The repo serves the following three development objects:
1. implementation of conic optimization algorithms at principle level;
2. test open-source convex solvers' interface and performance in Julia and C/C++;
3. customized solvers based on open-source projects, for embeded system (CPU + Parallel acceleration).

Based on solvers' classification in three dimensions, {First-Order, Second-Order}x{IPM, Exterior Point Method(Penalty)}x{Primal, Primal-dual},
only the three types of solvers and responding projects are focused, for real-time trajectory generation's OCP on embeded system:
1. Homogeneous Self-dual Embedding(HSDE) IPM for Linear SOCP :   ECOS
2. Homogeneous Embedding(HE) IPM for Quadratic SOCP:             Clarabel 
3. HE Douglas-Rachford Splitting(DRS) for Quadratic SOCP:             SCS \
   and a similar variant, Proportional-integral projected gradient method (PIPG)

HE DRS is preferred for real-time platform, due to the advantages:
1. Infeasibility Certification concurrent with convergence, by Homogeneous Embedding;
2. Converge rapidly to local modest precision solution, by first-order splitting iteration;
3. Low computational load of single iteration, by conical projection and fixed LDLT decomposition, \
   compare to updated decomposition of each iteration in IPM's newton direction;
4. Warm-start, beneficial for MPC.
   
