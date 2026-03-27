#!/usr/bin/env bash
set -euo pipefail

TARGET_UID="${HOST_UID:-}"
TARGET_GID="${HOST_GID:-}"
SALOME_LAUNCHER="${SALOME_LAUNCHER:-}"

resolve_salome_launcher() {
  local candidate

  if [ -n "$SALOME_LAUNCHER" ] && [ -x "$SALOME_LAUNCHER" ]; then
    printf '%s\n' "$SALOME_LAUNCHER"
    return 0
  fi

  for candidate in \
    /opt/salome/salome \
    "${SALOME_ROOT:-}/salome" \
    /opt/SALOME-*/salome \
    /opt/salome/run_salome.sh \
    "${SALOME_ROOT:-}/run_salome.sh" \
    /opt/SALOME-*/run_salome.sh; do
    if [ -n "$candidate" ] && [ -x "$candidate" ]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done

  return 1
}

ensure_machine_id() {
  if [ -s /etc/machine-id ]; then
    return 0
  fi

  if command -v dbus-uuidgen >/dev/null 2>&1; then
    rm -f /etc/machine-id
    dbus-uuidgen --ensure=/etc/machine-id
  fi
}

if ! SALOME_LAUNCHER="$(resolve_salome_launcher)"; then
  echo "SALOME launcher not found. Checked /opt/salome and /opt/SALOME-*/{salome,run_salome.sh}" >&2
  exit 1
fi

if [ "${1:-}" = "--gui" ]; then
  shift
fi

ensure_machine_id

if [ "$(id -u)" -eq 0 ]; then
  if [ -n "$TARGET_UID" ] || [ -n "$TARGET_GID" ]; then
    if [ -z "$TARGET_UID" ] || [ -z "$TARGET_GID" ]; then
      echo "HOST_UID and HOST_GID must be set together" >&2
      exit 1
    fi

    case "$TARGET_UID:$TARGET_GID" in
      (*[!0-9:]*)
        echo "HOST_UID and HOST_GID must be numeric, got ${TARGET_UID}:${TARGET_GID}" >&2
        exit 1
        ;;
    esac

    current_uid="$(id -u salome)"
    current_gid="$(id -g salome)"
    current_group="$(id -g -n salome)"

    if [ "$current_gid" != "$TARGET_GID" ]; then
      existing_group="$(getent group "$TARGET_GID" | cut -d: -f1 || true)"
      if [ -n "$existing_group" ]; then
        usermod -g "$existing_group" salome
      else
        groupmod -o -g "$TARGET_GID" "$current_group"
      fi
    fi

    if [ "$current_uid" != "$TARGET_UID" ]; then
      existing_user="$(getent passwd "$TARGET_UID" | cut -d: -f1 || true)"
      if [ -n "$existing_user" ] && [ "$existing_user" != "salome" ]; then
        echo "HOST_UID $TARGET_UID already belongs to user $existing_user" >&2
        exit 1
      fi
      usermod -o -u "$TARGET_UID" salome
    fi
  fi

  chown -R salome:"$(id -g -n salome)" /home/salome
  export HOME=/home/salome
  export USER=salome
  export LOGNAME=salome
  exec runuser -u salome -- "$SALOME_LAUNCHER" "$@"
fi

exec "$SALOME_LAUNCHER" "$@"
