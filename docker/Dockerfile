FROM ubuntu:22.04


ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PYTHONUNBUFFERED=1


RUN apt-get update && apt-get install -y --no-install-recommends \
    # --- Core build / system tools ---
    build-essential cmake git wget ca-certificates curl \
    pkg-config unzip locales gfortran gdb vim-common libyaml-cpp-dev \
    \
    # --- Python toolchain ---
    python3 python3-dev python3-pip python3-venv python3-wheel \
    swig \
    \
    # --- MPI + parallel math stacks ---
    openmpi-bin libopenmpi-dev libmetis-dev \
    libparmetis-dev libsuperlu-dist-dev \
    \
    # --- BLAS/LAPACK (keep one optimized stack; OpenBLAS is fine) ---
    libopenblas-dev libblas-dev liblapack-dev \
    \
    # --- Linear algebra / scientific libs ---
    libblas-dev liblapack-dev libopenblas-dev libatlas-base-dev \
    libmetis-dev libsuperlu-dev zlib1g-dev \
    \
    # --- Graphics & visualization stack ---
    gmsh libgmsh-dev paraview netgen \
    \
    # --- X11 + OpenGL runtime libraries (needed by gmsh, paraview, netgen) ---
    libx11-dev libxrender1 libxext6 libsm6 libxft2 \
    libgl1-mesa-dev libglu1-mesa-dev libglew-dev libfreetype-dev \
    libfontconfig1-dev libpng-dev libsdl2-dev libglm-dev \
    x11-apps \
    && ln -sf $(which python3) /usr/local/bin/python \
    && locale-gen en_US.UTF-8 \
    && rm -rf /var/lib/apt/lists/*



    
# Create user
ARG USER_ID=1000
ARG GROUP_ID=100


RUN group_name=$(getent group $GROUP_ID | cut -d: -f1) \
    && if [ -z "$group_name" ]; then \
         groupadd -g $GROUP_ID felix; \
       else \
         echo "Using existing group: $group_name"; \
       fi \
    && useradd -m -s /bin/bash -u $USER_ID -g $GROUP_ID felix


# Switch to felix
USER felix
WORKDIR /home/felix

# Python packages
COPY --chown=felix:users requirements.txt /home/felix/requirements.txt
RUN python3 -m pip install --upgrade pip setuptools wheel && \
    pip install --user -r requirements.txt
ENV PATH=/home/felix/.local/bin:${PATH}

# ===== hypre build with OpenMP =====
ENV HYPRE_TAG=v2.33.0
ENV HYPRE_DIR=/home/felix/hypre_src
RUN git clone --depth 1 https://github.com/hypre-space/hypre.git --branch ${HYPRE_TAG} ${HYPRE_DIR} && \
    mkdir -p /home/felix/hypre && cd /home/felix/hypre && \
    cmake ${HYPRE_DIR}/src \
      -DCMAKE_BUILD_TYPE=Debug \
      -DHYPRE_ENABLE_PRINT_ERRORS=ON \
      -DHYPRE_ENABLE_MPI=ON \
      -DHYPRE_ENABLE_OPENMP=ON \
      -DBUILD_SHARED_LIBS=ON \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_CXX_COMPILER=mpicxx \
      -DCMAKE_INSTALL_PREFIX=${HYPRE_DIR}/install && \
    cmake --build . -j"$(nproc)" && \
    cmake --install .

ENV HYPRE_ROOT=${HYPRE_DIR}/install
ENV LD_LIBRARY_PATH=${HYPRE_ROOT}/lib:${LD_LIBRARY_PATH}
ENV LIBRARY_PATH=${HYPRE_ROOT}/lib:${LIBRARY_PATH}
ENV C_INCLUDE_PATH=${HYPRE_ROOT}/include:${C_INCLUDE_PATH}
ENV CPLUS_INCLUDE_PATH=${HYPRE_ROOT}/include:${CPLUS_INCLUDE_PATH}


# ===== MFEM C++ build (now with MPI + OpenMP + ParMETIS + SuperLU_DIST) =====
# FInd SuperLU libararies and turn on 
ENV MFEM_TAG=v4.8
ENV CC=mpicc CXX=mpicxx
RUN git clone --depth 1 --branch ${MFEM_TAG} https://github.com/mfem/mfem.git /home/felix/mfem && \
    mkdir /home/felix/mfem/build && cd /home/felix/mfem/build && \
    cmake -S .. -B . \
      -DCMAKE_BUILD_TYPE=Debug \
      -DMFEM_USE_MPI=ON \
      -DMFEM_USE_OPENMP=ON \
      -DMFEM_USE_METIS=ON \
      -DMFEM_USE_SUPERLU=OFF \
      -DMFEM_USE_PETSC=OFF \
      -DMFEM_USE_GZSTREAM=ON \
      -DMFEM_USE_SUNDIALS=OFF \
      -DCMAKE_INSTALL_PREFIX=/home/felix/mfem-install && \
    cmake --build . --target install -j$(nproc)

ENV PATH=/home/felix/mfem-install/bin:${PATH}
ENV LD_LIBRARY_PATH=/home/felix/mfem-install/lib:${LD_LIBRARY_PATH}
ENV CPLUS_INCLUDE_PATH=/home/felix/mfem-install/include:${CPLUS_INCLUDE_PATH}
ENV LIBRARY_PATH=/home/felix/mfem-install/lib:${LIBRARY_PATH}

# PyMFEM
RUN pip install --user mfem

RUN git clone --depth 1 https://github.com/glvis/glvis.git /home/felix/glvis && \
    cd /home/felix/glvis && \
    make MFEM_DIR=/home/felix/mfem/build && \
    mkdir -p /home/felix/.local/bin && \
    cp glvis /home/felix/.local/bin/


WORKDIR /home/felix/work
EXPOSE 8888

