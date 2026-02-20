#!/usr/bin/env bash
set -euo pipefail


# -----------------------------
# Argument parsing
# -----------------------------

MODE="tui"
DEFAULT_TUI_SCRIPT="/work/geometry/build_geometry.py"
BASE_PATH_ENV=""   # will map to BASE_PATH inside container

while [[ $# -gt 0 ]]; do
  case "$1" in
    gui)
      MODE="gui"
      shift
      ;;
    tui)
      MODE="tui"
      shift
      ;;
    --config_path)
      if [[ -z "${2:-}" ]]; then
        echo "ERROR: --config_path requires a value"
        exit 1
      fi
      BASE_PATH_ENV="$2"
      shift 2
      ;;
    *.py)
      MODE="tui"
      TUI_SCRIPT="$1"
      shift
      ;;
    *)
      echo "ERROR: Unknown argument '$1'"
      echo "Usage:"
      echo "  ./run_salome.sh [gui|tui] [--config_path PATH] [script.py]"
      exit 1
      ;;
  esac
done

IMAGE="trophime/salome:9.9.0-focal"
DEFAULT_TUI_SCRIPT="/work/geometry/build_geometry.py"

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

# -----------------------------
# Argument parsing
# -----------------------------

MODE="tui"
TUI_SCRIPT="$DEFAULT_TUI_SCRIPT"

if [[ $# -ge 1 ]]; then
  case "$1" in
    gui)
      MODE="gui"
      ;;
    tui)
      MODE="tui"
      ;;
    *.py)
      MODE="tui"
      TUI_SCRIPT="$1"
      ;;
    *)
      echo "ERROR: Unknown argument '$1'"
      echo "Usage:"
      echo "  ./run_salome.sh"
      echo "  ./run_salome.sh tui"
      echo "  ./run_salome.sh gui"
      echo "  ./run_salome.sh /path/to/script.py"
      exit 1
      ;;
  esac
fi

# -----------------------------
# Docker arguments
# -----------------------------

DOCKER_ARGS=(
  --rm
  --user "$(id -u):$(id -g)"
  --shm-size=2g
  --ipc=host
  --net=host
  -v "$PROJECT_ROOT:/work"
  -w /work/geometry 
)

if [[ -n "$BASE_PATH_ENV" ]]; then
  DOCKER_ARGS+=(
    -e "BASE_PATH=$BASE_PATH_ENV"
  )
fi

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

DOCKER_ARGS+=(
  -it
)
# -----------------------------
# Run SALOME
# -----------------------------

if [[ "$MODE" == "gui" ]]; then
  if [[ -z "$DISPLAY_ENV" ]]; then
    echo "ERROR: DISPLAY is not set. Cannot run GUI mode."
    exit 1
  fi

  docker run "${DOCKER_ARGS[@]}" "${GUI_ARGS[@]}" \
    "$IMAGE" salome
else
  echo "Running SALOME in TUI mode"
  echo "TUI script: $TUI_SCRIPT"
  echo "Working dir: /work/geometry"

  docker run "${DOCKER_ARGS[@]}" \
    "$IMAGE" salome -t "$TUI_SCRIPT"
fi
