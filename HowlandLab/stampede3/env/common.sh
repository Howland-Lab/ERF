#!/bin/bash
# Sourced by the three environment files.  ERF_HOME = this ERF clone (derived from this file's location if unset).
HL_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)          # HowlandLab/stampede3
export ERF_HOME=${ERF_HOME:-$(cd "$HL_DIR/../.." && pwd)}
export HL_DIR
export HL_BUILD=$HL_DIR/_build                                    # build trees (git-ignored)
