import superscs
import algorithm
import sequtils
import strformat
import strutils

const
  sqrt2 = 1.4142135623730951
  inv_sqrt2 = 0.7071067811865475

type
  BoundKey* = enum
    fx # uj =...= lj
    fr # free
    lo # lj<=...
    ra # lj<=...<=uj
    up # ...<=uj
  
  ConeType* = enum
    quad  # quadratic cone
    rquad # rotated quadratic cone
    pexp  # primal exponential cone
    # dexp  # dual exponential cone
    ppow  # primal power cone
    # dpow  # dual power cone
    # zero  # zero cone

  ObjSense* = enum
    minimize
    maximize
    
  ConeData* = object
    ct* : ConeType
    conepar* : float
    submem* : seq[int]
    
  Env* = object
    m* : cint
    n* : cint
    data* : ptr ScsData
    cone* : ScsCone
    sol* : ptr ScsSolution
    info* : ptr ScsInfo
          
    A* : ScsAMatrix #CSC format: https://kul-forbes.github.io/scs/page_sparse_matrices.html
    b* : seq[float]
    c* : seq[float]

    # lower and upper bounds of variables
    var_idx* : seq[cint]
    var_boundKeys* : seq[BoundKey]
    var_lowerBounds* : seq[float]
    var_upperBounds* : seq[float]

    # linear constraints
    numConstraints*: int      
    A_rowIdx* : seq[cint]
    A_colIdx* : seq[cint]
    A_values* : seq[float]
    A_boundkeys* : seq[BoundKey]
    A_lowerbounds* : seq[float]
    A_upperbounds* : seq[float]

    # cones
    cones* : seq[ConeData]
    q*: seq[cint]  # for quadratic cones
    p*: seq[float] # for power cones

    # maximize or minimize, if maximize flip sign of c
    objsense* : ObjSense
  
  OptimizationSolution* = object
    primal_obj* : float
    dual_obj*   : float
    primal_x*   : seq[float]
    dual_y*     : seq[float]
    slack_s*    : seq[float]

proc appendvars*(env: var Env, n : SomeInteger) =
  when defined(debug):
    assert env.n == 0
  env.n = n.cint
  env.c = newSeq[float](n)

proc putvarbound*(env: var Env, idx: SomeInteger, boundKey: BoundKey, lowerBound: float, upperBound: float) =
  when defined(debug):
    case boundKey:
      of fx:
        assert lowerBound == upperBound
      of fr:
        discard
      of lo:
        discard
      of ra:
        assert lowerBound <= upperBound
      of up:
        discard
  
  env.var_idx.add(idx.cint)
  env.var_boundKeys.add(boundKey)
  env.var_lowerBounds.add(lowerBound)
  env.var_upperBounds.add(upperBound)

proc putvarboundslice*(env: var Env, start_idx: SomeInteger, end_idx: SomeInteger, boundKeys: openArray[BoundKey], lowerBounds: openArray[float], upperBounds: openArray[float]) =
  when defined(debug):
    assert (end_idx.int - start_idx.int).int == boundKeys.len, fmt"start_idx={start_idx}, end_idx={end_idx}, boundKeys={boundKeys}"
    assert (end_idx.int - start_idx.int).int == lowerBounds.len, fmt"start_idx={start_idx}, end_idx={end_idx}, boundKeys={boundKeys}"
    assert (end_idx.int - start_idx.int).int == upperBounds.len, fmt"start_idx={start_idx}, end_idx={end_idx}, boundKeys={boundKeys}"

  for idx in start_idx .. end_idx-1:
    env.putvarbound(idx, boundKeys[idx], lowerBounds[idx], upperBounds[idx])

proc appendcons*(env: var Env, n : SomeInteger) =
  when defined(debug):
    assert env.numConstraints == 0
  env.numConstraints += n.cint
  env.A_boundkeys = newSeq[BoundKey](n)
  env.A_lowerbounds = newSeq[float](n)
  env.A_upperbounds = newSeq[float](n)

proc putaij*(env: var Env, i : SomeInteger, j : SomeInteger, aij: float) =
  when defined(debug):
    assert i.cint < env.numConstraints
    assert j.cint < env.n
  if aij != 0.0:
    env.A_rowIdx.add(i.cint)
    env.A_colIdx.add(j.cint)
    env.A_values.add(aij)

proc putaijlist*(env: var Env, subi : openArray[SomeInteger], subj : openArray[SomeInteger], valij : openArray[float]) =
  when defined(debug):
    assert subi.len == subj.len
    assert subi.len == valij.len
  for idx in 0 .. subi.len-1:
    env.putaij(subi[idx], subj[idx], valij[idx])

proc putarow*(env: var Env, i: SomeInteger, subj : openArray[SomeInteger], valij : openArray[float]) =
  when defined(debug):
    assert subj.len == valij.len
  for idx in 0 .. subj.len-1:
    env.putaij(i, subj[idx], valij[idx])

proc putconbound*(env: var Env, i: SomeInteger, boundKey : BoundKey, lowerBound : float, upperBound : float) =
  when defined(debug):
    assert lowerBound == upperBound
    assert i < env.numConstraints
  env.A_boundkeys[i] = boundKey
  env.A_lowerbounds[i] = lowerBound
  env.A_upperbounds[i] = upperBound      

proc putconboundlist*(env: var Env, subi: openArray[SomeInteger], boundKeys: openArray[BoundKey], lowerBounds: openArray[float], upperBounds: openArray[float]) =
  when defined(debug):
    assert subi.len == boundKeys.len
    assert subi.len == lowerBounds.len
    assert subi.len == upperBounds.len
  for idx in subi:
    env.putconbound(idx, boundKeys[idx], lowerBounds[idx], upperBounds[idx])

proc appendcone*(env : var Env, coneType: ConeType, conepar: float = 0.0, submem : seq[int]) =
  when defined(debug):
    case coneType:
      of quad:
        assert submem.len >= 1
      of rquad:
        assert submem.len >= 2
      of pexp:
        assert submem.len == 3
      of ppow:
        assert submem.len >= 2
  env.cones.add(ConeData(ct: coneType, conepar: conepar, submem: submem))

proc putcj*(env : var Env, j : SomeInteger, cj: float) =
  when defined(debug):
    assert j < env.n
  env.c[j] = cj

proc putclist*(env : var Env, subj: openArray[SomeInteger], val : openArray[float]) = 
  for j in subj:
    env.putcj(j, val[j])

proc putobjsense*(env : var Env, objsense : ObjSense) =
  env.objsense = objsense

proc tmp_cmp[T](arr1, arr2 : seq[T], i, j: int): int =
  if i == j:
    return 0
  elif arr1[i] < arr1[j] or (arr1[i] == arr1[j] and arr2[i] < arr2[j]):
    return -1
  elif arr1[i] > arr1[j] or (arr1[i] == arr1[j] and arr2[i] > arr2[j]):
    return 1
  else:
    var
      err : ref OSError
    new(err)
    err.msg = fmt"""
  
  Duplicate index, i={i}, j={j}, arr1[i]={arr1[i]}, arr1[j]={arr1[j]}, arr2[i]={arr2[i]}, arr2[j]={arr2[j]}
  Why does this happen? Remember to create a copy for variables in cones.
  For example, to create the following two cones for variables (x1,x2,x3)
    1) (x1, x2, x3) in SOC(3)
    2) (x1, x2, x3) in PowerC
  Remember to create a copy for x1, x2, x3
  """
    raise err

proc argsort[T](arr1, arr2 : seq[T]): seq[int]=
  result = (0 .. arr1.len-1).toSeq
  sort(result, proc (i, j : int): int = tmp_cmp[T](arr1, arr2, i, j))

proc coo_to_csc(
  n_row, n_col : SomeInteger,
  Ai : openArray[SomeInteger], Aj : openArray[SomeInteger], Ax : openArray[float],
  ) : (seq[cint], seq[cint], seq[float]) = 
  #https://github.com/scipy/scipy/blob/3b36a574dc657d1ca116f6e230be694f3de31afc/scipy/sparse/sparsetools/coo.h
  #https://kul-forbes.github.io/scs/page_sparse_matrices.html
  var
    nnz = Ax.len
    Bp = newSeq[cint](n_col + 1)
    Bj = newSeq[cint](nnz)
    Bx = newSeq[float](nnz)
    sorted_idx = argsort(Aj.toSeq, Ai.toSeq)

  for new_idx, idx in sorted_idx:
  # for idx in 0 .. nnz-1:
    Bp[Aj[idx] + 1] += 1.cint
    Bj[new_idx] = Ai[idx].cint
    Bx[new_idx] = Ax[idx]
  
  for idx in 1 .. n_col:
    Bp[idx] += Bp[idx-1]
  
  return (Bp, Bj, Bx)

proc print_coo(n_row, n_col: SomeInteger, Ai, Aj : openArray[SomeInteger], Ax: openArray[float]):seq[seq[float]] =
  var
    A : seq[seq[float]]
  for idx in 0 .. n_row-1:
    A.add(newSeq[float](n_col))

  for idx in 0 .. Ai.len-1:
    A[Ai[idx]][Aj[idx]] = Ax[idx]
  
  echo "["
  for row in A:
    var rowstr = ($row).string
    echo rowstr[1 .. rowstr.len-1], ","
  echo "]"

  return A

proc get_row_expr(Arow: seq[float], flip : bool = false): string =
  var
    adj = if flip: -1.0 else: 1.0
    tmp : float
    
  for j, aij in Arow:
    tmp = aij * adj
    if tmp != 0.0:
      if tmp > 0:
        if result.len > 0:
          result &= " + "
        if tmp == sqrt2:
          result &= fmt"x{j} * sqrt(2)"
        elif tmp == inv_sqrt2:
          result &= fmt"x{j} * 1/sqrt(2)"
        elif tmp != 1.0:
          result &= fmt"x{j} * {tmp}"
        else:
          result &= fmt"x{j}"
      else:
        result &= " - "
        if tmp == -sqrt2:
          result &= fmt"x{j} * sqrt(2)"
        elif tmp == -inv_sqrt2:
          result &= fmt"x{j} * 1/sqrt(2)"
        elif tmp != -1.0:
          result &= fmt"x{j} * {-tmp}"
        else:
          result &= fmt"x{j}"

proc setup*(env: var Env) =
  var
    f : cint = 0.cint
    l : cint = 0.cint
    m : cint = 0.cint
    qsize : cint = 0.cint
    sum_qsize : cint = 0.cint
    ssize : cint = 0.cint
    sum_ssize : cint = 0.cint
    ep : cint = 0.cint
    ed : cint = 0.cint
    psize : cint = 0.cint

    nq : cint = 0.cint
    ns : cint = 0.cint

    rows : seq[cint]
    cols : seq[cint]
    vals : seq[float]
  
  ####################################################################################
  # Fixed Constraints (Zero Cones)
  ####################################################################################

  # fx var bounds
  for i in 0 .. env.var_boundKeys.len-1:
    if (env.var_boundKeys[i] == fx) or (env.var_boundKeys[i] == ra and env.var_lowerBounds[i] == env.var_upperBounds[i]):
      rows.add(f)
      cols.add(env.var_idx[i])
      vals.add(1.0)
      env.b.add(env.var_lowerBounds[i])
      f += 1.cint
  
  # fx linear constraints bounds
  var
    lp_cols : seq[cint]
    lp_vals : seq[float]
  for i in 0 .. env.A_boundkeys.len-1:
    if (env.A_boundkeys[i] == fx) or (env.A_boundkeys[i] == ra and env.A_lowerbounds[i] == env.A_upperbounds[i]):
      lp_cols = @[]
      lp_vals = @[]
      for idx in 0 .. env.A_rowIdx.len-1:
        if env.A_rowIdx[idx] == i:
          lp_cols.add(env.A_colIdx[idx])
          lp_vals.add(env.A_values[idx])
      if lp_cols.len > 0:
        for idx in 0 .. lp_cols.len-1:
          rows.add(f)
          cols.add(lp_cols[idx])
          vals.add(lp_vals[idx])
        env.b.add(env.A_lowerbounds[i])
        f += 1.cint

  m = f
  ####################################################################################
  # LP Cones
  ####################################################################################
  # var bounds
  for i in 0 .. env.var_boundKeys.len-1:
    if (env.var_boundKeys[i] == lo):
      rows.add(f + l)
      cols.add(env.var_idx[i])
      vals.add(-1.0)
      env.b.add(-env.var_lowerBounds[i])
      l += 1.cint
    elif (env.var_boundKeys[i] == up):
      rows.add(f + l)
      cols.add(env.var_idx[i])
      vals.add(1.0)
      env.b.add(env.var_upperBounds[i])
      l += 1.cint
    elif env.var_boundKeys[i] == ra and env.var_lowerBounds[i] != env.var_upperBounds[i]:
      rows.add(f + l)
      cols.add(env.var_idx[i])
      vals.add(-1.0)
      env.b.add(-env.var_lowerBounds[i])
      l += 1.cint

      rows.add(f + l)
      cols.add(env.var_idx[i])
      vals.add(1.0)
      env.b.add(env.var_upperBounds[i])
      l += 1.cint

  # linear constraints bounds
  for i in 0 .. env.A_boundkeys.len-1:
    if env.A_boundkeys[i] == lo:
      lp_cols = @[]
      lp_vals = @[]
      for idx in 0 .. env.A_rowIdx.len-1:
        if env.A_rowIdx[idx] == i:
          lp_cols.add(env.A_colIdx[idx])
          lp_vals.add(env.A_values[idx])
      if lp_cols.len > 0:
        for idx in 0 .. lp_cols.len-1:
          rows.add(f + l)
          cols.add(lp_cols[idx])
          vals.add(-lp_vals[idx])
        env.b.add(-env.A_lowerbounds[i])
        l += 1.cint          
    elif env.A_boundkeys[i] == up:
      lp_cols = @[]
      lp_vals = @[]
      for idx in 0 .. env.A_rowIdx.len-1:
        if env.A_rowIdx[idx] == i:
          lp_cols.add(env.A_colIdx[idx])
          lp_vals.add(env.A_values[idx])
      if lp_cols.len > 0:
        for idx in 0 .. lp_cols.len-1:
          rows.add(f + l)
          cols.add(lp_cols[idx])
          vals.add(lp_vals[idx])
        env.b.add(env.A_upperbounds[i])
        l += 1.cint
    elif env.A_boundkeys[i] == ra and env.A_lowerbounds[i] != env.A_upperbounds[i]:
      lp_cols = @[]
      lp_vals = @[]
      for idx in 0 .. env.A_rowIdx.len-1:
        if env.A_rowIdx[idx] == i:
          lp_cols.add(env.A_colIdx[idx])
          lp_vals.add(env.A_values[idx])
      if lp_cols.len > 0:
        for idx in 0 .. lp_cols.len-1:
          rows.add(f)
          cols.add(lp_cols[idx])
          vals.add(-lp_vals[idx])
        env.b.add(-env.A_lowerbounds[i])
        l += 1.cint

        for idx in 0 .. lp_cols.len-1:
          rows.add(f)
          cols.add(lp_cols[idx])
          vals.add(lp_vals[idx])
        env.b.add(env.A_upperbounds[i])
        l += 1.cint

  m = f + l
  ####################################################################################
  # SO Cones
  ####################################################################################
  for cone in env.cones:
    if cone.ct == quad:
      nq = 0.cint
      for idx in cone.submem:
        rows.add(f + l + sum_qsize)
        cols.add(idx.cint)
        vals.add(-1.0)
        env.b.add(0.0)
        nq += 1.cint
      qsize += 1.cint
      sum_qsize += nq
      env.q.add(nq)
    elif cone.ct == rquad: # convert rquad cone to quad cone
      nq = 0.cint

      rows.add(f + l + sum_qsize)
      cols.add(cone.submem[0].cint)
      vals.add(-inv_sqrt2)
      rows.add(f + l + sum_qsize)
      cols.add(cone.submem[1].cint)
      vals.add(-inv_sqrt2)
      env.b.add(0.0)        
      nq += 1.cint

      rows.add(f + l + sum_qsize + 1.cint)
      cols.add(cone.submem[0].cint)
      vals.add(-inv_sqrt2)
      rows.add(f + l + sum_qsize + 1.cint)
      cols.add(cone.submem[1].cint)
      vals.add(inv_sqrt2)
      env.b.add(0.0)        
      nq += 1.cint

      for idx in 2 .. cone.submem.len-1:
        rows.add(f + l + sum_qsize + idx.cint)
        cols.add(cone.submem[idx].cint)
        vals.add(-1.0)
        nq += 1.cint
      
      qsize += 1.cint
      sum_qsize += nq
      env.q.add(nq)
  
  m = f + l + sum_qsize
  ####################################################################################
  # PSD Cones : Not Implemented
  ####################################################################################
  ssize = 0.cint
  sum_ssize = 0.cint

  m = f + l + sum_qsize + sum_ssize
  ####################################################################################
  # Exponential Cones
  ####################################################################################
  for cone in env.cones:
    if cone.ct == pexp:
      for idx in 0 .. 2:
        rows.add(f + l + sum_qsize + sum_ssize + ep * 3.cint + idx.cint)
        cols.add(cone.submem[idx].cint)
        vals.add(-1.0)
      ep += 1.cint
  
  m = f + l + sum_qsize + sum_ssize + ep * 3.cint
  ####################################################################################
  # Dual Exponential Cones : Not Implemented
  ####################################################################################
  ed = 0.cint

  m = f + l + sum_qsize + sum_ssize + ep * 3.cint + ed * 3.cint
  ####################################################################################
  # Power Cones
  ####################################################################################
  for cone in env.cones:
    if cone.ct == ppow:
      for idx in 0 .. 2:
        rows.add(f + l + sum_qsize + sum_ssize + ep * 3.cint + ed * 3.cint + psize * 3.cint + idx.cint)
        cols.add(cone.submem[idx].cint)
        vals.add(-1.0)
      env.p.add(cone.conepar)
      psize += 1.cint
  m = f + l + sum_qsize + sum_ssize + ep * 3.cint + ed * 3.cint + psize * 3.cint

  ####################################################################################
  # Convert CSR to CSC for A
  ####################################################################################
  when defined(debug):
    echo '='.repeat(87)
    echo "Details of the conic problem:"

    echo '-'.repeat(87)
    echo "A = "
    var A_dense = print_coo(m, env.n, rows, cols, vals)
    echo fmt"b     = {env.b}"
    echo fmt"c     = {env.c}"
    echo fmt"f     = {f:<5} # rows for equality constraints"
    echo fmt"l     = {l:<5} # rows for Inequality linear constraints"
    echo fmt"qsize = {qsize:<5} # number of Quradic Cones (SOC)"
    echo fmt"q     = {env.q:<5} # size of each Quradic Cone (SOC)"
    echo fmt"ep    = {ep:<5} # number of exponential cones (EPC)"
    echo fmt"psize = {psize:<5} # number of power cones (POWC)"
    echo fmt"p     = {env.p:<5} # power parameter of power cones"
    echo fmt"m     = {m:<5} # rows for A"

    var
      expr : string

    if f > 0:
      echo '-'.repeat(87)
      echo "Equality Constraints:"

      for idx in 0 .. f-1:
        expr = get_row_expr(A_dense[idx])
        expr &= fmt" = {env.b[idx]}"
        echo expr

    if l > 0:
      echo '-'.repeat(87)
      echo "Inequality Constraints:"
      for idx in f .. f+l-1:
        expr =  get_row_expr(A_dense[idx])
        expr &= fmt" <= {env.b[idx]}"
        echo expr
    
    if sum_qsize > 0:
      echo '-'.repeat(87)
      echo "SO Cones:"
      var start_idx = f + l
      var end_idx : cint
      var item_expr : string
      for nq in env.q:
        end_idx = start_idx + nq
        expr = ""
        for idx in start_idx .. end_idx-1:
          item_expr = get_row_expr(A_dense[idx], flip=true)
          if idx == start_idx:
            expr = item_expr & " >= norm2("
          else:
            expr &= item_expr
            if idx < end_idx-1: expr &= ","
        expr &= fmt")"
        echo expr
        start_idx = end_idx

    if ep > 0:  
      echo '-'.repeat(87)
      echo "Exponential Cones:"
      var
        start_idx = f + l + sum_qsize + sum_ssize
      for idx in 0 .. ep-1:
        var
          ep_start = start_idx + idx * 3
          ep_end = start_idx + idx * 3 + 3
          item_exprs = newSeq[string](3)
          tmp_item_expr : string
        for row_idx in ep_start .. ep_end-1:
          item_exprs[row_idx - ep_start] = get_row_expr(A_dense[row_idx], flip=true)
        expr = fmt"({item_exprs[1]}) * exp(({item_exprs[0]})/({item_exprs[1]})) <= {item_exprs[2]}"
        echo expr

    if psize > 0:
      echo '-'.repeat(87)
      echo "Power Cones:"
      var
        start_idx = f + l + sum_qsize + sum_ssize + ep * 3.cint
      for idx in 0 .. psize-1:
        var
          p_start = start_idx + idx * 3
          p_end = start_idx + idx * 3 + 3
          item_exprs = newSeq[string](3)
          tmp_item_expr : string
          alpha = env.p[idx]
        for row_idx in p_start .. p_end-1:
          item_exprs[row_idx - p_start] = get_row_expr(A_dense[row_idx], flip=true)
        expr = fmt"({item_exprs[0]})^({alpha}) * ({item_exprs[1]})^({1-alpha}) >= |{item_exprs[2]}|"
        echo expr

    echo '-'.repeat(87)
    echo "Obj:"
    expr = get_row_expr(env.c)
    if env.objsense == ObjSense.maximize:
      expr = "max : " & expr
    else:
      expr = "min : " & expr
    echo expr
    echo '='.repeat(87)

  ####################################################################################
  ####################################################################################
  # Set Up ScSData
  ####################################################################################
  env.data = scs_init_data()
  env.m = m
  env.data.m = env.m
  env.data.n = env.n

  var
    (Ap, Ai, Ax) = coo_to_csc(env.m, env.n, rows, cols, vals)
  env.A.x = Ax[0].addr
  env.A.p = Ap[0].addr
  env.A.i = Ai[0].addr
  env.A.m = env.m
  env.A.n = env.n

  env.data.A = env.A.addr
  env.data.b = env.b[0].addr

  if env.objsense == maximize:
    for idx in 0 .. env.c.len-1:
      env.c[idx] = -env.c[idx]
  
  env.data.c = env.c[0].addr
  
  scs_set_default_settings(env.data)

  env.data.stgs.eps = 1e-9
  env.data.stgs.rho_x = 1.0
  env.data.stgs.verbose = 0
  env.data.stgs.sse = 0.7
  env.data.stgs.direction = restarted_broyden

  scs_set_restarted_broyden_settings(env.data, 50.cint)
  scs_set_anderson_settings(env.data, 3.cint)
  scs_set_tolerance(env.data, 1e-10)

  env.cone.f = f
  env.cone.l = l
  env.cone.qsize = qsize
  if env.q.len > 0:
    env.cone.q = env.q[0].addr
  else:
    env.cone.q = nil
  env.cone.ssize = 0.cint
  env.cone.s = nil
  env.cone.ep = ep
  env.cone.ed = ed
  env.cone.psize = psize
  if env.p.len > 0:
    env.cone.p = env.p[0].addr
  else:
    env.cone.p = nil

  env.sol = scs_init_sol()
  env.info = scs_init_info()

proc update_settings*(env: var Env, key: string, value: float) =
  case key:
    of "normalize":
      env.data.stgs.normalize = value.cint
    of "scale":
      env.data.stgs.scale = value
    of "rho_x":
      env.data.stgs.rho_x = value
    of "max_time_milliseconds":
      env.data.stgs.max_time_milliseconds = value
    of "max_iters":
      env.data.stgs.max_iters = value.cint
    of "previous_max_iters":
      env.data.stgs.previous_max_iters = value.cint
    of "eps":
      env.data.stgs.eps = value
    of "alpha":
      env.data.stgs.alpha = value
    of "cg_rate":
      env.data.stgs.cg_rate = value
    of "verbose":
      env.data.stgs.verbose = value.cint
    of "warm_start":
      env.data.stgs.warm_start = value.cint
    of "do_super_scs":
      env.data.stgs.do_super_scs = value.cint
    of "k0":
      env.data.stgs.k0 = value.cint
    of "c_bl":
      env.data.stgs.c_bl = value
    of "k1":
      env.data.stgs.k1 = value.cint
    of "k2":
      env.data.stgs.k2 = value.cint
    of "c1":
      env.data.stgs.c1 = value
    of "sse":
      env.data.stgs.sse = value
    of "ls":
      env.data.stgs.ls = value.cint
    of "beta":
      env.data.stgs.beta = value
    of "sigma":
      env.data.stgs.sigma = value
    of "direction":
      env.data.stgs.direction = value.ScsDirectionType
    of "thetabar":
      env.data.stgs.thetabar = value
    of "memory":
      env.data.stgs.memory = value.cint
    of "tRule":
      env.data.stgs.tRule = value.cint
    of "broyden_init_scaling":
      env.data.stgs.broyden_init_scaling = value.cint
    of "do_record_progress":
      env.data.stgs.do_record_progress = value.cint
    of "do_override_streams":
      env.data.stgs.do_override_streams = value.cint
    else:
      raise newException(OSError, fmt"""

    Unknown optimization setting key : {key}
    Check available settings at https://github.com/twvacek/superscs/blob/bfe3e9f968ff6bbb14b5fb3b8e1ebb4f5f233a9a/include/scs.h#L290
    """)

proc update_settings*(env: var Env, key: string, value: SomeInteger) =
  env.update_settings(key, value.float)

proc optimize*(env : var Env) =
  discard scs(env.data, env.cone.addr, env.sol, env.info)
  
proc get_solution*(env : var Env): OptimizationSolution =
  var
    ptr_sol_x = cast[ptr UncheckedArray[float]](env.sol.x)
    ptr_sol_y = cast[ptr UncheckedArray[float]](env.sol.y)
    ptr_sol_s = cast[ptr UncheckedArray[float]](env.sol.s)
    sol_x = newSeq[float](env.n)
    sol_y = newSeq[float](env.n)
    sol_s = newSeq[float](env.m)
    primal_obj = env.info.pobj
    dual_obj = env.info.dobj
  
  if env.objsense == maximize:
    primal_obj = -primal_obj
    dual_obj = -dual_obj
  
  for idx in 0 .. env.n-1:
    sol_x[idx] = ptr_sol_x[idx]
    sol_y[idx] = ptr_sol_y[idx]
  for idx in 0 .. env.m-1:
    sol_s[idx] = ptr_sol_s[idx]
  
  return OptimizationSolution(primal_obj: primal_obj, dual_obj: dual_obj, primal_x: sol_x, dual_y: sol_y, slack_s: sol_s)

proc clean*(env : var Env) =
  deallocShared(env.data.stgs)
  deallocShared(env.data)
  scs_free_sol(env.sol)
  scs_free_info(env.info)
