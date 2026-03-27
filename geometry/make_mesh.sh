#!/usr/bin/env bash
set -euo pipefail

SALOME_VERSION="${SALOME_VERSION:-9.15.0}"
SALOME_URL="${SALOME_URL:-https://www.salome-platform.org/?download_id=2885&sdm_process_download=1}"
IMAGE="${SALOME_IMAGE:-mfem-salome:${SALOME_VERSION}}"
DEFAULT_TUI_SCRIPT="/work/geometry/build_geometry.py"

# Resolve project root: parent directory of this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
DOCKERFILE_PATH="$SCRIPT_DIR/Dockerfile.salome"
ENTRYPOINT_PATH="$SCRIPT_DIR/salome-container-entrypoint.sh"
BUILD_LABEL="mfem.salome.build_signature"

compute_build_signature() {
  local sources_hash
  sources_hash="$(sha256sum "$DOCKERFILE_PATH" "$ENTRYPOINT_PATH" | sha256sum | awk '{print $1}')"
  printf '%s\n%s\n%s\n' "$SALOME_VERSION" "$SALOME_URL" "$sources_hash" | sha256sum | awk '{print $1}'
}

BUILD_SIGNATURE="$(compute_build_signature)"

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
BASE_PATH_ENV=""

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
      [[ -n "${2:-}" ]] || { echo "ERROR: --config_path requires a value"; exit 1; }
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
      echo "  ./make_mesh.sh [gui|tui] [--config_path PATH] [script.py]"
      exit 1
      ;;
  esac
done

ensure_image() {
  local current_signature=""

  if docker image inspect "$IMAGE" >/dev/null 2>&1; then
    current_signature="$(docker image inspect --format "{{ index .Config.Labels \"$BUILD_LABEL\" }}" "$IMAGE" 2>/dev/null || true)"
    if [[ "$current_signature" == "$BUILD_SIGNATURE" ]]; then
      return
    fi

    echo "Rebuilding SALOME image: $IMAGE"
    echo "Build inputs changed since the existing image was created"
  else
    echo "Building SALOME image: $IMAGE"
  fi

  echo "SALOME version: $SALOME_VERSION"
  echo "SALOME URL: $SALOME_URL"

  docker build \
    -f "$DOCKERFILE_PATH" \
    -t "$IMAGE" \
    --label "$BUILD_LABEL=$BUILD_SIGNATURE" \
    --build-arg "SALOME_VERSION=$SALOME_VERSION" \
    --build-arg "SALOME_URL=$SALOME_URL" \
    --build-arg "USER_ID=$(id -u)" \
    --build-arg "GROUP_ID=$(id -g)" \
    "$PROJECT_ROOT"
}
# -----------------------------
# Docker arguments
# -----------------------------

DOCKER_ARGS=(
  --rm
  --shm-size=2g
  --ipc=host
  --net=host
  -e "HOST_UID=$(id -u)"
  -e "HOST_GID=$(id -g)"
  -v "$PROJECT_ROOT:/work"
  -w /work/geometry
)

if [[ -n "${SALOME_CONTAINMENT_TRACE:-}" ]]; then
  DOCKER_ARGS+=(-e "SALOME_CONTAINMENT_TRACE=$SALOME_CONTAINMENT_TRACE")
fi

if [[ -n "${SALOME_CONTAINMENT_TRACE_INNER:-}" ]]; then
  DOCKER_ARGS+=(-e "SALOME_CONTAINMENT_TRACE_INNER=$SALOME_CONTAINMENT_TRACE_INNER")
fi

if [[ -n "${SALOME_CONTAINMENT_ALLOW_FEATURE_POINT_CLOUD:-}" ]]; then
  DOCKER_ARGS+=(-e "SALOME_CONTAINMENT_ALLOW_FEATURE_POINT_CLOUD=$SALOME_CONTAINMENT_ALLOW_FEATURE_POINT_CLOUD")
fi

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

# -----------------------------
# Run SALOME
# -----------------------------

if [[ "$MODE" == "gui" ]]; then
  ensure_image
  if [[ -z "$DISPLAY_ENV" ]]; then
    echo "ERROR: DISPLAY is not set. Cannot run GUI mode."
    exit 1
  fi

  docker run "${DOCKER_ARGS[@]}" -it "${GUI_ARGS[@]}" \
    "$IMAGE" --gui
else
  ensure_image
  echo "Running SALOME in TUI mode"
  echo "TUI script: $TUI_SCRIPT"
  echo "Working dir: /work/geometry"

  docker run "${DOCKER_ARGS[@]}" \
    "$IMAGE" -t "$TUI_SCRIPT"
fi
