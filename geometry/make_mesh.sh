#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_salome.sh gui
#   ./run_salome.sh
#   ./run_salome.sh tui
#   ./run_salome.sh /work/geometry/build_geometry.py
#
# Behavior:
#   - First arg "gui"  -> run SALOME GUI
#   - Otherwise        -> run headless/TUI (salome -t <script>)
#     Default script:  /work/geometry/build_geometry.py

IMAGE="trophime/salome:9.9.0-focal"
DEFAULT_TUI_SCRIPT="/work/geometry/build_geometry.py"

MODE="${1:-tui}"

# Resolve project root: parent directory of this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Detect NVIDIA support
if command -v nvidia-smi >/dev/null 2>&1; then
  USE_NVIDIA=1
else
  USE_NVIDIA=0
fi

DISPLAY_ENV="${DISPLAY:-}"
HOST_XAUTHORITY="${XAUTHORITY:-$HOME/.Xauthority}"

# If a .py file is passed, treat it as TUI script
TUI_SCRIPT="$DEFAULT_TUI_SCRIPT"
if [[ "${1:-}" == *.py ]]; then
  TUI_SCRIPT="$1"
  MODE="tui"
fi

DOCKER_ARGS=(
  --rm
  --user "$(id -u):$(id -g)"
  --shm-size=2g
  --ipc=host
  --net=host
  -v "$PROJECT_ROOT:/work"
)

# GPU support
if [[ "$USE_NVIDIA" -eq 1 ]]; then
  DOCKER_ARGS+=(
    --device nvidia.com/gpu=all
    -e NVIDIA_VISIBLE_DEVICES=all
    -e NVIDIA_DRIVER_CAPABILITIES=all
  )
fi

GUI_ARGS=(
  -e "DISPLAY=$DISPLAY_ENV"
  -e QT_X11_NO_MITSHM=1
  -e "XAUTHORITY=$HOST_XAUTHORITY"
  -v /tmp/.X11-unix:/tmp/.X11-unix:ro
  -v "$HOST_XAUTHORITY:$HOST_XAUTHORITY:ro"
)

if [[ "$MODE" == "gui" ]]; then
  if [[ -z "$DISPLAY_ENV" ]]; then
    echo "ERROR: DISPLAY is not set. Cannot run GUI mode."
    exit 1
  fi

  docker run "${DOCKER_ARGS[@]}" "${GUI_ARGS[@]}" \
    "$IMAGE" salome
else
  docker run "${DOCKER_ARGS[@]}" \
    "$IMAGE" salome -t "$TUI_SCRIPT"
fi
