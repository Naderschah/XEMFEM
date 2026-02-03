# Get user data 
#export UID=$(id -u) 
export GID=$(id -g)
export XEMFEM_PATH=$(pwd)
docker compose -f ./docker/docker-compose.yml run --rm --build --service-ports mfem-shell

