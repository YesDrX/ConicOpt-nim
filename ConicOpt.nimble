version       = "0.1.0"
author        = "YesDrX"
description   = "A nim library for conic (convex cone) optimization"
license       = "MIT"
backend       = "c"
srcDir        = "src"
installExt    = @["nim","so","dylib","dll","h"]

import os, strutils

when defined(window):
    const
        libname = "scsdir.dll"
elif defined(macos):
    const
        libname = "libscsdir.dynlib"
else:
    const
        libname = "libscsdir.so"

before install:
    var cmd : string
    echo "Building SuperScS Library : https://github.com/twvacek/superscs"
    cd("superscs")
    cmd = "make"
    echo cmd
    exec cmd
    cd("..")
    cmd = "cp -vrf superscs/out/" & libname & " src/lib/"
    echo cmd
    exec cmd