#!/usr/bin/env bash
set -euo pipefail

# ---------------- User data ----------------
export DOCKER_UID=$(id -u)
export GID=$(id -g)
export XEMFEM_PATH=$(pwd)

# ---------------- X detection & stubbing ----------------
if [ -n "${DISPLAY:-}" ] && [ -S /tmp/.X11-unix/X0 ]; then
    echo "[mfem] X detected"
    export ENABLE_X=1
    export XAUTHORITY=${XAUTHORITY:-$HOME/.Xauthority}
else
    echo "[mfem] X NOT detected, using stubs"
    export ENABLE_X=0
    export DISPLAY=
    export XAUTHORITY=/tmp/.Xauthority.stub
    touch "$XAUTHORITY"
fi

# ---------------- GPU group IDs ----------------
export VIDEO_GID=$(getent group video | cut -d: -f3 || echo "")
export RENDER_GID=$(getent group render | cut -d: -f3 || echo "")

# ---------------- Environment detection ----------------

COMPOSE_PROFILES=("shell")

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

# NVIDIA GPU present?
if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1; then
    HAS_GPU=1
    echo "[mfem] NVIDIA GPU detected"
# Fallback: DRM render node (covers Mesa/iGPU cases)
elif [ -e /dev/dri/renderD128 ]; then
    HAS_GPU=1
    echo "[mfem] DRM render node detected"
else
    echo "[mfem] No GPU detected"
fi

# ---------------- Profile selection ----------------

if [ "$HAS_GPU" -eq 1 ]; then
    COMPOSE_PROFILES+=("gpu")
    if [ "$IS_NIXOS" -eq 1 ]; then
        COMPOSE_PROFILES+=("nixos")
        echo "[mfem] Dispatching: GPU + NixOS (CDI)"
    else
        COMPOSE_PROFILES+=("classic")
        echo "[mfem] Dispatching: GPU + classic (gpus=all)"
    fi
else
    COMPOSE_PROFILES+=("no-gpu")
    echo "[mfem] Dispatching: no-GPU (software rendering)"
fi

# ---------------- Compose invocation ----------------

PROFILE_ARGS=()
for p in "${COMPOSE_PROFILES[@]}"; do
    PROFILE_ARGS+=(--profile "$p")
done

exec docker compose \
    -f ./docker/docker-compose.yml \
    "${PROFILE_ARGS[@]}" \
    run --rm --build --service-ports mfem-shell
