#!/usr/bin/env bash
set -euo pipefail

if [ ! -x /opt/salome/salome ]; then
  echo "SALOME launcher not found at /opt/salome/salome" >&2
  exit 1
fi

exec /opt/salome/salome "$@"
