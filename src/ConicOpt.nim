import ConicOpt/api
export api

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
  # env.appendcone(rquad, 0.0, @[0, 1, 2])
  # env.appendcone(pexp, 0.0, @[0, 1, 2])
  # env.appendcone(ppow, 0.3, @[0, 1, 2])
  env.setup()
  env.update_settings("eps", 1e-15)
  env.update_settings("verbose", 0)
  env.optimize()
  var sol = env.get_solution()
  echo sol
  env.clean()