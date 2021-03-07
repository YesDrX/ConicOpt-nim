# ConicOpt-nim
Conic Optimization Library for Nim

# Requirements
* blas

# Install
```
git clone https://github.com/YesDrX/ConicOpt-nim
cd ConicOpt-nim
nimble install
```

# Example
* Code
```nim
import conicOpt

when isMainModule:
  import sequtils
  import strformat
  import strutils
  var
    env = Env()
  env.appendvars(3)
  env.putvarboundslice(0, 3, ra.repeat(3), 0.0.repeat(3), 0.6.repeat(3))
  env.putclist(@[0, 1, 2], @[0.1, 0.2, 0.3])
  env.putobjsense(maximize)
  env.appendcons(1)
  env.putaijlist(@[0, 0, 0], @[0, 1, 2], @[1.0, 1.0, 1.0])
  env.putconbound(0, fx, 1.0, 1.0)
  env.setup()
  env.update_settings("eps", 1e-15)
  env.update_settings("verbose", 1)
  env.optimize()
  var sol = env.get_solution()
  echo sol
  env.clean()
```
* Compile
```
nim c --run -d:debug example.nim
```
* Output
```
=======================================================================================
Details of the conic problem:
---------------------------------------------------------------------------------------
A = 
[
[1.0, 1.0, 1.0],
[-1.0, 0.0, 0.0],
[1.0, 0.0, 0.0],
[0.0, -1.0, 0.0],
[0.0, 1.0, 0.0],
[0.0, 0.0, -1.0],
[0.0, 0.0, 1.0],
]
b     = @[1.0, -0.0, 0.6, -0.0, 0.6, -0.0, 0.6]
c     = @[0.1, 0.2, 0.3]
f     = 1     # rows for equality constraints
l     = 6     # rows for Inequality linear constraints
qsize = 0     # number of Quradic Cones (SOC)
q     = @[]   # size of each Quradic Cone (SOC)
ep    = 0     # number of exponential cones (EPC)
psize = 0     # number of power cones (POWC)
p     = @[]   # power parameter of power cones
m     = 7     # rows for A
---------------------------------------------------------------------------------------
Equality Constraints:
x0 + x1 + x2 = 1.0
---------------------------------------------------------------------------------------
Inequality Constraints:
 - x0 <= -0.0
x0 <= 0.6
 - x1 <= -0.0
x1 <= 0.6
 - x2 <= -0.0
x2 <= 0.6
---------------------------------------------------------------------------------------
Obj:
max : x0 * 0.1 + x1 * 0.2 + x2 * 0.3
=======================================================================================
---------------------------------------------------------------------------------------
        SCS v1.3.3 - Superlinear Splitting Conic Solver (SuperSCS)
        Web: https://kul-forbes.github.io/scs
        (c) P. Sopasakis, K. Menounou, P. Patrinos, KU Leuven, 2017-8
        (c) Brendan O'Donoghue, Stanford University, 2012-2016
---------------------------------------------------------------------------------------
Lin-sys: sparse-direct, nnz in A = 9
eps = 1.00e-15, alpha = 1.50, max_iters = 10000, normalize = 1, scale = 1.00
do_super_scs = 1, direction = 150, memory = 3
Variables n = 3, constraints m = 7
Cones:  primal zero / dual free vars: 1
        linear vars: 6
Setup time: 2.20e-04s

Running SuperSCS...
Allocated Memory: 36.86kB
---------------------------------------------------------------------------------------
 Iter | pri res | dua res | rel gap | pri obj | dua obj | kap/tau |   FPR   | time (s)
---------------------------------------------------------------------------------------
     0| 5.46e+00  4.14e-01  4.92e-01 -2.67e+00 -5.78e-01  0.00e+00  2.16e+00  9.28e-05 
    20| 1.77e-03  2.89e-03  9.00e-04 -2.60e-01 -2.59e-01  0.00e+00  1.09e-02  5.04e-04 
    40| 2.66e-04  2.71e-04  4.53e-05 -2.60e-01 -2.60e-01  0.00e+00  1.04e-03  7.06e-04 
    60| 4.15e-05  1.54e-05  9.85e-08 -2.60e-01 -2.60e-01  0.00e+00  1.07e-04  9.08e-04 
    80| 5.07e-06  7.93e-07  1.13e-06 -2.60e-01 -2.60e-01  0.00e+00  1.19e-05  1.09e-03 
   100| 3.18e-07  2.76e-07  1.53e-07 -2.60e-01 -2.60e-01  0.00e+00  1.25e-06  1.27e-03 
   120| 7.19e-09  4.00e-08  1.03e-08 -2.60e-01 -2.60e-01  0.00e+00  1.40e-07  1.46e-03 
   140| 4.53e-09  2.87e-09  3.74e-10 -2.60e-01 -2.60e-01  0.00e+00  1.37e-08  1.62e-03 
   160| 6.80e-10  5.68e-11  1.33e-10 -2.60e-01 -2.60e-01  0.00e+00  1.62e-09  1.79e-03 
   180| 5.64e-11  2.56e-11  1.89e-11 -2.60e-01 -2.60e-01  0.00e+00  1.66e-10  1.96e-03 
   200| 3.00e-12  4.64e-12  2.09e-12 -2.60e-01 -2.60e-01  0.00e+00  1.69e-11  2.13e-03 
   220| 3.47e-13  4.97e-13  6.76e-14 -2.60e-01 -2.60e-01  0.00e+00  1.89e-12  2.30e-03 
   240| 7.53e-14  2.15e-14  1.28e-14 -2.60e-01 -2.60e-01  0.00e+00  1.99e-13  2.47e-03 
   260| 9.33e-15  2.07e-15  2.34e-15 -2.60e-01 -2.60e-01  0.00e+00  2.14e-14  2.64e-03 
   280| 5.65e-16  2.23e-16  1.05e-16 -2.60e-01 -2.60e-01  0.00e+00  7.80e-16  2.80e-03 
---------------------------------------------------------------------------------------
Status: Solved
Timing: Solve time: 00:00:00.2
        Lin-sys: nnz in L factor: 19, avg solve time: 6.73e-07s
        Cones: avg projection time: 6.73e-07s
---------------------------------------------------------------------------------------
Error metrics:
dist(s, K) = -0.0000e+00, dist(y, K*) = 0.0000e+00, s'y/|s||y| = 0.0000e+00
|Ax + s - b|_2 / (1 + |b|_2) = 5.6550e-16
|A'y + c|_2 / (1 + |c|_2) = 2.2263e-16
|c'x + b'y| / (1 + |c'x| + |b'y|) = 1.0499e-16
---------------------------------------------------------------------------------------
c'x = -0.2600, -b'y = -0.2600
=======================================================================================
(primal_obj: 0.2599999999999997, dual_obj: 0.2599999999999999, primal_x: @[-4.520440012763929e-16, 0.4, 0.5999999999999993], dual_y: @[0.2000000000000001, 0.1000000000000001, 0.0], slack_s: @[0.0, 0.0, 0.6000000000000002, 0.4, 0.1999999999999998, 0.5999999999999994, 0.0])
```
