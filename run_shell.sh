#!/usr/bin/env bash
set -euo pipefail

# ---------------- User data ----------------
export DOCKER_UID="$(id -u)"
export GID="$(id -g)"
export XEMFEM_PATH="$(pwd)"

# ---------------- X detection & stubbing ----------------
if [ -n "${DISPLAY:-}" ] && [ -S /tmp/.X11-unix/X0 ]; then
  echo "[mfem] X detected"
  export ENABLE_X=1
  export XAUTHORITY="${XAUTHORITY:-$HOME/.Xauthority}"
else
  echo "[mfem] X NOT detected, using stubs"
  export ENABLE_X=0
  export DISPLAY=
  export XAUTHORITY=/tmp/.Xauthority.stub
  touch "$XAUTHORITY"
fi

# ---------------- GPU group IDs ----------------
export VIDEO_GID="$(getent group video | cut -d: -f3 || true)"
export RENDER_GID="$(getent group render | cut -d: -f3 || true)"

# ---------------- Environment detection ----------------
# We will choose EXACTLY ONE compose service to run.
SERVICE="mfem-shell"

# ---- Detect NixOS ----
IS_NIXOS=0
if [ -f /etc/NIXOS ] || grep -qi nixos /etc/os-release 2>/dev/null; then
  IS_NIXOS=1
  echo "[mfem] Host OS: NixOS"
else
  echo "[mfem] Host OS: non-Nix"
fi

# ---- Detect GPU availability ----
HAS_GPU=0
HAS_NVIDIA=0

if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1; then
  HAS_GPU=1
  HAS_NVIDIA=1
  echo "[mfem] NVIDIA GPU detected"
elif [ -e /dev/dri/renderD128 ]; then
  HAS_GPU=1
  HAS_NVIDIA=0
  echo "[mfem] DRM render node detected (non-NVIDIA)"
else
  HAS_GPU=0
  HAS_NVIDIA=0
  echo "[mfem] No GPU detected"
fi

# ---------------- Dispatch (service selection) ----------------
if [ "${HAS_GPU:-0}" -eq 1 ]; then
  if [ "${IS_NIXOS:-0}" -eq 1 ] && [ "${HAS_NVIDIA:-0}" -eq 1 ]; then
    SERVICE="mfem-shell-nixos"     # NixOS + NVIDIA CDI
  elif [ "${HAS_NVIDIA:-0}" -eq 1 ]; then
    SERVICE="mfem-shell-gpu"       # non-Nix + NVIDIA (gpus: all)
  else
    SERVICE="mfem-shell-drm"       # non-NVIDIA GPU via /dev/dri
  fi
else
  SERVICE="mfem-shell"             # no GPU
fi

# ---------------- Compose invocation ----------------
exec docker compose \
  -f ./docker/docker-compose.yml \
  run --rm --build --service-ports "$SERVICE"
