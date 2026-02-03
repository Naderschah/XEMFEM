#!/usr/bin/env bash
set -euo pipefail

# ---- user-configurable (no floating tags) ----
FREECAD_VERSION="${FREECAD_VERSION:-1.0.2}"   # pick an explicit version directory from the repo
IMAGE_NAME="local/freecad-cli:${FREECAD_VERSION}"
REPO_URL="https://github.com/amrit3701/docker-freecad-cli.git"

# Prefer Ubuntu when both exist (Dockerfile.ubuntu), else fall back to Dockerfile / Dockerfile.debian
PREFER_UBUNTU="${PREFER_UBUNTU:-1}"

ARCH="$(uname -m)"
case "$ARCH" in
  x86_64|amd64) ARCH_DIR="amd64" ;;
  aarch64|arm64) ARCH_DIR="arm64" ;;
  *)
    echo "Unsupported arch: $ARCH" >&2
    exit 1
    ;;
esac

WORKDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CACHE_DIR="${WORKDIR}/.freecad-docker-src"
OUT_DIR="${WORKDIR}/.freecad-docker"

mkdir -p "$OUT_DIR"

if [[ ! -d "$CACHE_DIR/.git" ]]; then
  git clone --depth 1 "$REPO_URL" "$CACHE_DIR"
else
  git -C "$CACHE_DIR" fetch --depth 1 origin
  git -C "$CACHE_DIR" reset --hard origin/master
fi

BASE_DIR="${CACHE_DIR}/${FREECAD_VERSION}/${ARCH_DIR}"
if [[ ! -d "$BASE_DIR" ]]; then
  echo "Could not find version/arch directory: $BASE_DIR" >&2
  echo "Available versions in repo:" >&2
  ls -1 "$CACHE_DIR" | sed 's/^/  - /' >&2
  exit 1
fi

# ---- Choose Dockerfile: prefer Ubuntu if present ----
CANDIDATES=()
if [[ "${PREFER_UBUNTU}" == "1" ]]; then
  CANDIDATES+=("${BASE_DIR}/Dockerfile.ubuntu")
fi
CANDIDATES+=(
  "${BASE_DIR}/Dockerfile"
  "${BASE_DIR}/Dockerfile.debian"
)

SELECTED=""
for f in "${CANDIDATES[@]}"; do
  if [[ -f "$f" ]]; then
    SELECTED="$f"
    break
  fi
done

if [[ -z "$SELECTED" ]]; then
  echo "No Dockerfile found in: $BASE_DIR" >&2
  echo "Expected one of:" >&2
  for f in "${CANDIDATES[@]}"; do echo "  - $f" >&2; done
  echo "Files present:" >&2
  ls -la "$BASE_DIR" >&2
  exit 1
fi

echo "Using Dockerfile: $SELECTED"
cp -f "$SELECTED" "${OUT_DIR}/Dockerfile"

echo "Building image: ${IMAGE_NAME}"
docker build -t "${IMAGE_NAME}" -f "${OUT_DIR}/Dockerfile" "${OUT_DIR}"

echo "Verifying FreeCAD command exists in the image..."
docker run --rm -t "${IMAGE_NAME}" bash -lc 'command -v FreeCADCmd || command -v freecadcmd || (echo "No FreeCADCmd/freecadcmd found" >&2; exit 1)'

echo "OK."
echo "Next:"
echo "  docker compose run --rm occ"
