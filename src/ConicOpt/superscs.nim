const 
  libPath = currentSourcePath()[0 .. ^22] & "lib"

when defined(windows):
  const
    libscsdir* = libPath & "\\scsdir.dll"
elif defined(macosx):
  const
    libscsdir* = libPath & "/libscsdir.dylib"
else:
  const
    libscsdir* = libPath & "/libscsdir.so"

{.emit:"""/*INCLUDESECTION*/
#include <time.h>
""".}

const
  SCS_INFEASIBLE_INACCURATE* = -7.cint
  SCS_UNBOUNDED_INACCURATE* = -6.cint
  SCS_SIGINT* = -5.cint
  SCS_FAILED* = -4.cint
  SCS_INDETERMINATE* = -3.cint
  SCS_INFEASIBLE* = -2.cint
  SCS_UNBOUNDED* = -1.cint
  SCS_UNFINISHED* = 0.cint
  SCS_SOLVED* = 1.cint
  SCS_SOLVED_INACCURATE* = 2.cint

type
  ScsDirectionType* {.pure.} = enum
    restarted_broyden = 100.cint
    anderson_acceleration = 150.cint
    fixed_point_residual = 200.cint
    full_broyden = 300.cint

  scs_direction_cache* {.bycopy.} = object
    S*: ptr float
    U*: ptr float
    S_minus_Y*: ptr float
    t*: ptr float
    ls_wspace_length*: cint
    ls_wspace*: ptr float
    mem_cursor*: cint
    mem*: cint
    current_mem*: cint

  scs_work* {.bycopy.} = object
    m*: cint
    n*: cint
    l*: cint
    u*: ptr float
    v*: ptr float
    u_t*: ptr float
    u_prev*: ptr float
    u_b*: ptr float
    h*: ptr float
    g*: ptr float
    pr*: ptr float
    dr*: ptr float
    b*: ptr float
    c*: ptr float
    R*: ptr float
    R_prev*: ptr float
    dir*: ptr float
    H*: ptr float
    dut*: ptr float
    wu*: ptr float
    wu_t*: ptr float
    wu_b*: ptr float
    Rwu*: ptr float
    nrmR_con*: float
    Sk*: ptr float
    Yk*: ptr float
    stepsize*: float
    s_b*: ptr float
    gTh*: float
    sc_b*: float
    sc_c*: float
    nm_b*: float
    nm_c*: float
    kap_b*: float
    A*: ptr ScsAMatrix
    p*: ptr ScsPrivWorkspace
    stgs*: ptr ScsSettings
    scal*: ptr ScsScaling
    coneWork*: ptr ScsConeWork
    direction_cache*: ptr ScsDirectionCache

  scs_data* {.bycopy.} = object
    m*: cint
    n*: cint
    A*: ptr ScsAMatrix
    b*: ptr float
    c*: ptr float
    stgs*: ptr ScsSettings

  scs_settings* {.bycopy.} = object
    normalize*: cint
    scale*: float
    rho_x*: float
    max_time_milliseconds*: float
    max_iters*: cint
    previous_max_iters*: cint
    eps*: float
    alpha*: float
    cg_rate*: float
    verbose*: cint
    warm_start*: cint
    do_super_scs*: cint
    k0*: cint
    c_bl*: float
    k1*: cint
    k2*: cint
    c1*: float
    sse*: float
    ls*: cint
    beta*: float
    sigma*: float
    direction*: ScsDirectionType
    thetabar*: float
    memory*: cint
    tRule*: cint
    broyden_init_scaling*: cint
    do_record_progress*: cint
    do_override_streams*: cint
    output_stream*: ptr FILE

  scs_solution* {.bycopy.} = object
    x*: ptr float
    y*: ptr float
    s*: ptr float

  scs_info* {.bycopy.} = object
    status*: array[32, char]
    iter*: cint
    statusVal*: cint
    history_length*: cint
    cg_total_iters*: cint
    pobj*: float
    dobj*: float
    resPri*: float
    resDual*: float
    resInfeas*: float
    resUnbdd*: float
    relGap*: float
    setupTime*: float
    solveTime*: float
    linsys_total_solve_time_ms*: float
    progress_relgap*: ptr float
    progress_respri*: ptr float
    progress_resdual*: ptr float
    progress_pcost*: ptr float
    progress_dcost*: ptr float
    progress_norm_fpr*: ptr float
    progress_time*: ptr float
    progress_iter*: ptr cint
    progress_mode*: ptr cint
    progress_ls*: ptr cint
    allocated_memory*: culonglong

  scs_scaling* {.bycopy.} = object
    D*: ptr float
    E*: ptr float
    meanNormRowA*: float
    meanNormColA*: float

  scs_cone* {.bycopy.} = object
    f*: cint
    l*: cint
    q*: ptr cint
    qsize*: cint
    s*: ptr cint
    ssize*: cint
    ep*: cint
    ed*: cint
    p*: ptr float
    psize*: cint

  ScsConeWork* {.bycopy.} = object
    total_cone_time*: float

  scs_a_data_matrix* {.bycopy.} = object # https://kul-forbes.github.io/scs/page_sparse_matrices.html
    x* : ptr float
    i* : ptr cint
    p* : ptr cint
    m* : cint
    n* : cint

  ScsAMatrix* = scs_a_data_matrix

  scs_cs* {.bycopy.} = object
    nzmax* : cint
    m* : cint
    n* : cint
    p* : ptr cint
    i* : ptr cint
    x* : ptr float
    nz* : cint

  scs_private_data {.bycopy.} = object
    L : ptr scs_cs
    D : ptr float
    P : ptr cint
    bp: ptr float
    totalSolveTime : float

  ScsPrivWorkspace* = scs_private_data
  
  timespec* {.importc: "timespec", nodecl.} = object

  scs_timer* {.bycopy.} = object
    tic*: timespec
    toc*: timespec

  ScsTimer* = scs_timer
  ScsData* = scs_data
  ScsSettings* = scs_settings
  ScsSolution* = scs_solution
  ScsInfo* = scs_info
  ScsScaling* = scs_scaling
  ScsWork* = scs_work
  ScsCone* = scs_cone
  ScsDirectionCache* = scs_direction_cache

proc scs_get_cone_boundaries*(k: ptr ScsCone; boundaries: ptr ptr cint): cint {.cdecl, importc: "scs_get_cone_boundaries", dynlib: libscsdir.}
proc scs_init_conework*(k: ptr ScsCone): ptr ScsConeWork {.cdecl, importc: "scs_init_conework", dynlib: libscsdir.}
proc scs_get_cone_header*(k: ptr ScsCone): cstring {.cdecl, importc: "scs_get_cone_header", dynlib: libscsdir.}
proc scs_validate_cones*(d: ptr ScsData; k: ptr ScsCone): cint {.cdecl, importc: "scs_validate_cones", dynlib: libscsdir.}
proc scs_project_dual_cone*(x: ptr float; k: ptr ScsCone; c: ptr ScsConeWork; warm_start: ptr float; iter: cint): cint {.cdecl, importc: "scs_project_dual_cone", dynlib: libscsdir.}
proc scs_finish_cone*(coneWork: ptr ScsConeWork) {.cdecl, importc: "scs_finish_cone", dynlib: libscsdir.}
proc scs_get_cone_summary*(info: ptr ScsInfo; c: ptr ScsConeWork): cstring {.cdecl, importc: "scs_get_cone_summary", dynlib: libscsdir.}
proc scs_svd_workspace_size*(m: cint; n: cint): cint {.cdecl, importc: "scs_svd_workspace_size", dynlib: libscsdir.}
proc scs_qr_workspace_size*(m: cint; n: cint): cint {.cdecl, importc: "scs_qr_workspace_size", dynlib: libscsdir.}
proc scs_qrls*(m: cint; n: cint; A: ptr float; b: ptr float; wspace: ptr float; wsize: cint): cint {.cdecl, importc: "scs_qrls", dynlib: libscsdir.}
proc scs_svdls*(m: cint; n: cint; A: ptr float; b: ptr float;  wspace: ptr float; wsize: cint; rcond: float;  singular_values: ptr float; rank: ptr cint): cint {.cdecl, importc: "scs_svdls", dynlib: libscsdir.}
proc scs_set_as_scaled_array*(x: ptr float; a: ptr float; b: float;   len: cint) {.cdecl, importc: "scs_set_as_scaled_array", dynlib: libscsdir.}
proc scs_scale_array*(a: ptr float; b: float; len: cint) {.cdecl, importc: "scs_scale_array", dynlib: libscsdir.}
proc scs_inner_product*(x: ptr float; y: ptr float; len: cint): float {.cdecl, importc: "scs_inner_product", dynlib: libscsdir.}
proc scs_norm_squared*(v: ptr float; len: cint): float {.cdecl, importc: "scs_norm_squared", dynlib: libscsdir.}
proc scs_norm*(v: ptr float; len: cint): float {.cdecl, importc: "scs_norm", dynlib: libscsdir.}
proc scs_norm_infinity*(a: ptr float; l: cint): float {.cdecl, importc: "scs_norm_infinity", dynlib: libscsdir.}
proc scs_add_scaled_array*(a: ptr float; b: ptr float; n: cint; sc: float) {.cdecl, importc: "scs_add_scaled_array", dynlib: libscsdir.}
proc scs_add_array*(a: ptr float; b: ptr float; n: cint) {.cdecl, importc: "scs_add_array", dynlib: libscsdir.}
proc scs_axpy*(x: ptr float; u: ptr float; v: ptr float; a: float; b: float; n: cint) {.cdecl, importc: "scs_axpy", dynlib: libscsdir.}
proc scs_subtract_array*(a: ptr float; b: ptr float; n: cint) {.cdecl, importc: "scs_subtract_array", dynlib: libscsdir.}
proc scs_norm_difference*(a: ptr float; b: ptr float; l: cint): float {.cdecl, importc: "scs_norm_difference", dynlib: libscsdir.}
proc scs_norm_infinity_difference*(a: ptr float; b: ptr float; l: cint): float {.cdecl, importc: "scs_norm_infinity_difference", dynlib: libscsdir.}
proc scs_matrix_multiply*(rows_A: cint; cols_B: cint; cols_A: cint; alpha: float; A: ptr float; beta: float; B: ptr float; C: ptr float) {.cdecl, importc: "scs_matrix_multiply", dynlib: libscsdir.}
proc scs_matrix_transpose_multiply*(rows_A: cint; cols_B: cint; cols_A: cint; alpha: float; A: ptr float; beta: float; B: ptr float; C: ptr float) {.cdecl, importc: "scs_matrix_transpose_multiply", dynlib: libscsdir.}
proc scs_cgls_malloc_workspace*(m: cint; n: cint): ptr float {.cdecl, importc: "scs_cgls_malloc_workspace", dynlib: libscsdir.}
proc scs_cgls*(m: cint; n: cint; A: ptr float; b: ptr float; x: ptr float; tol: float; maxiter: ptr cint; wspace: ptr float): cint {.cdecl, importc: "scs_cgls", dynlib: libscsdir.}
proc scs_init_priv*(A: ptr ScsAMatrix; stgs: ptr ScsSettings): ptr ScsPrivWorkspace {.cdecl, importc: "scs_init_priv", dynlib: libscsdir.}
proc scs_solve_lin_sys*(A: ptr ScsAMatrix; stgs: ptr ScsSettings;  p: ptr ScsPrivWorkspace; b: ptr float; s: ptr float;  iter: cint): cint {.cdecl, importc: "scs_solve_lin_sys", dynlib: libscsdir.}
proc scs_free_priv*(p: ptr ScsPrivWorkspace) {.cdecl, importc: "scs_free_priv", dynlib: libscsdir.}
proc scs_accum_by_a_trans*(A: ptr ScsAMatrix; p: ptr ScsPrivWorkspace;  x: ptr float; y: ptr float) {.cdecl, importc: "scs_accum_by_a_trans", dynlib: libscsdir.}
proc scs_accum_by_a*(A: ptr ScsAMatrix; p: ptr ScsPrivWorkspace; x: ptr float; y: ptr float) {.cdecl, importc: "scs_accum_by_a", dynlib: libscsdir.}
proc scs_validate_linsys*(A: ptr ScsAMatrix): cint {.cdecl, importc: "scs_validate_linsys", dynlib: libscsdir.}
proc scs_get_linsys_method*(A: ptr ScsAMatrix; stgs: ptr ScsSettings): cstring {.cdecl, importc: "scs_get_linsys_method", dynlib: libscsdir.}
proc scs_get_linsys_summary*(p: ptr ScsPrivWorkspace; info: ptr ScsInfo): cstring {.cdecl, importc: "scs_get_linsys_summary", dynlib: libscsdir.}
proc scs_normalize_a*(A: ptr ScsAMatrix; stgs: ptr ScsSettings; k: ptr ScsCone;  scal: ptr ScsScaling) {.cdecl, importc: "scs_normalize_a", dynlib: libscsdir.}
proc scs_unnormalize_a*(A: ptr ScsAMatrix; stgs: ptr ScsSettings; scal: ptr ScsScaling) {.cdecl, importc: "scs_unnormalize_a", dynlib: libscsdir.}
proc scs_free_a_matrix*(A: ptr ScsAMatrix) {.cdecl, importc: "scs_free_a_matrix", dynlib: libscsdir.}
proc scs_linsys_is_indirect*(): cint {.cdecl, importc: "scs_linsys_is_indirect", dynlib: libscsdir.}
proc scs_linsys_total_cg_iters*(priv: ptr ScsPrivWorkspace): cint {.cdecl, importc: "scs_linsys_total_cg_iters", dynlib: libscsdir.}
proc scs_linsys_total_solve_time_ms*(priv: ptr ScsPrivWorkspace): float {.cdecl, importc: "scs_linsys_total_solve_time_ms", dynlib: libscsdir.}
proc scs_tic*(timer: ptr ScsTimer) {.cdecl, importc: "scs_tic", dynlib: libscsdir.}
proc scs_toc*(timer: ptr ScsTimer): float {.cdecl, importc: "scs_toc", dynlib: libscsdir.}
proc scs_strtoc*(str: cstring; timer: ptr ScsTimer): float {.cdecl, importc: "scs_strtoc", dynlib: libscsdir.}
proc scs_toc_quiet*(timer: ptr ScsTimer): float {.cdecl, importc: "scs_toc_quiet", dynlib: libscsdir.}
proc scs_print_cone_data*(cone: ptr ScsCone) {.cdecl, importc: "scs_print_cone_data", dynlib: libscsdir.}
proc scs_print_data*(data: ptr ScsData) {.cdecl, importc: "scs_print_data", dynlib: libscsdir.}
proc scs_print_work*(work: ptr ScsWork) {.cdecl, importc: "scs_print_work", dynlib: libscsdir.}
proc scs_print_array*(arr: ptr float; n: cint; name: cstring) {.cdecl, importc: "scs_print_array", dynlib: libscsdir.}
proc scs_set_default_settings*(data: ptr ScsData) {.cdecl, importc: "scs_set_default_settings", dynlib: libscsdir.}
proc scs_set_restarted_broyden_settings*(data: ptr ScsData; broyden_memory: cint) {.cdecl, importc: "scs_set_restarted_broyden_settings", dynlib: libscsdir.}
proc scs_set_anderson_settings*(data: ptr ScsData; anderson_memory: cint) {.cdecl, importc: "scs_set_anderson_settings", dynlib: libscsdir.}
proc scs_set_tolerance*(data: ptr ScsData; tolerance: float) {.cdecl, importc: "scs_set_tolerance", dynlib: libscsdir.}
proc scs_set_memory*(data: ptr ScsData; memory: cint) {.cdecl, importc: "scs_set_memory", dynlib: libscsdir.}
proc scs_free_sol*(sol: ptr ScsSolution) {.cdecl, importc: "scs_free_sol", dynlib: libscsdir.}
proc scs_free_data*(data: ptr ScsData) {.cdecl, importc: "scs_free_data", dynlib: libscsdir.}
proc scs_free_cone*(cone: ptr ScsCone) {.cdecl, importc: "scs_free_cone", dynlib: libscsdir.}
proc scs_free_data_cone*(d: ptr ScsData; k: ptr ScsCone) {.cdecl, importc: "scs_free_data_cone", dynlib: libscsdir.}
proc scs_free_info*(info: ptr ScsInfo) {.cdecl, importc: "scs_free_info", dynlib: libscsdir.}
proc scs_special_print*(print_mode: cint; stream: ptr FILE; format: cstring): cint {.varargs, cdecl, importc: "scs_special_print", dynlib: libscsdir.}
proc scs_init_sol*(): ptr ScsSolution {.cdecl, importc: "scs_init_sol", dynlib: libscsdir.}
proc scs_init_info*(): ptr ScsInfo {.cdecl, importc: "scs_init_info", dynlib: libscsdir.}
proc scs_init_data*(): ptr ScsData {.cdecl, importc: "scs_init_data", dynlib: libscsdir.}
proc scs*(d: ptr ScsData; k: ptr ScsCone; sol: ptr ScsSolution; info: ptr ScsInfo): cint {.cdecl, importc: "scs", dynlib: libscsdir.}
proc scs_version*(): cstring {.cdecl, importc: "scs_version", dynlib: libscsdir.}
proc scs_millis_to_time*(time: float; hours: ptr cint; minutes: ptr cint;   secs: ptr cint; sec_rest: ptr float) {.cdecl, importc: "scs_millis_to_time", dynlib: libscsdir.}
