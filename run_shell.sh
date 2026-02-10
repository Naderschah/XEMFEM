# Get user data 
export DOCKER_UID=$(id -u)
export GID=$(id -g)
export XEMFEM_PATH=$(pwd)

# ---- X detection & stubbing ----
if [ -n "${DISPLAY:-}" ] && [ -S /tmp/.X11-unix/X0 ]; then
    echo "X detected"
    export ENABLE_X=1
    export XAUTHORITY=${XAUTHORITY:-$HOME/.Xauthority}
else
    echo "X NOT detected, using stubs"
    export ENABLE_X=0
    export DISPLAY=
    export XAUTHORITY=/tmp/.Xauthority.stub
    touch "$XAUTHORITY"
fi

docker compose -f ./docker/docker-compose.yml run --rm --build --service-ports mfem-shell

