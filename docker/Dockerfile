FROM ubuntu:22.04


ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PYTHONUNBUFFERED=1


RUN apt-get update && apt-get install -y --no-install-recommends \
    # --- Core build / system tools ---
    build-essential cmake git wget ca-certificates curl ninja-build \
    pkg-config unzip locales gfortran gdb vim-common libyaml-cpp-dev \
    \
    # ----- QT buttons for paraview ------
    adwaita-icon-theme libqt5svg5 qt5-image-formats-plugins \
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
    # ---- HDF5 used for interpolator output -----
    libhdf5-dev \
    \
    # --- X11 + OpenGL runtime libraries (needed by gmsh, paraview, netgen) ---
    libx11-dev libxrender1 libxext6 libsm6 libxft2 \
    libgl1-mesa-dev libglu1-mesa-dev libglew-dev libfreetype-dev \
    mesa-utils mesa-common-dev freeglut3-dev \
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



# ====== VTK Toolchain for plotting and paraview -- very heavy to compile ==========
ARG VTK_VERSION=v9.5.0
RUN mkdir -p $HOME/vtk \
 && git clone --branch ${VTK_VERSION} --depth 1 --shallow-submodules \
        https://gitlab.kitware.com/vtk/vtk.git $HOME/vtk/source \
 && mkdir -p $HOME/vtk/build $HOME/vtk/install \
 && cd $HOME/vtk/build \
 && cmake ../source -GNinja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$HOME/vtk/install \
    -DVTK_BUILD_TESTING=OFF \
    -DVTK_BUILD_EXAMPLES=OFF \
    -DVTK_BUILD_ALL_MODULES=OFF \
    \
    # Core / data model
    -DVTK_MODULE_ENABLE_VTK_CommonCore=YES \
    -DVTK_MODULE_ENABLE_VTK_CommonDataModel=YES \
    \
    # IO for .vtu / .pvtu / PNG
    -DVTK_MODULE_ENABLE_VTK_IOXML=YES \
    -DVTK_MODULE_ENABLE_VTK_IOImage=YES \
    \
    # Filters we use (streamlines, geometry, etc.)
    -DVTK_MODULE_ENABLE_VTK_FiltersCore=YES \
    -DVTK_MODULE_ENABLE_VTK_FiltersGeometry=YES \
    -DVTK_MODULE_ENABLE_VTK_FiltersGeneral=YES \
    -DVTK_MODULE_ENABLE_VTK_FiltersSources=YES \
    -DVTK_MODULE_ENABLE_VTK_FiltersFlowPaths=YES \
    \
    # Rendering + annotation + charts for scalar bars, cube axes
    -DVTK_MODULE_ENABLE_VTK_RenderingCore=YES \
    -DVTK_MODULE_ENABLE_VTK_RenderingOpenGL2=YES \
    -DVTK_MODULE_ENABLE_VTK_RenderingFreeType=YES \
    -DVTK_MODULE_ENABLE_VTK_RenderingContext2D=YES \
    -DVTK_MODULE_ENABLE_VTK_RenderingContextOpenGL2=YES \
    -DVTK_MODULE_ENABLE_VTK_RenderingAnnotation=YES \
    -DVTK_MODULE_ENABLE_VTK_ChartsCore=YES \
    -DVTK_MODULE_ENABLE_VTK_InteractionStyle=YES \
    \
    # Headless / offscreen setup
    -DVTK_RENDERING_BACKEND=OpenGL2 \
    -DVTK_OPENGL_HAS_OSMESA=ON \
    -DVTK_USE_X=OFF \
    -DVTK_DEFAULT_RENDER_WINDOW_OFFSCREEN=ON \
 && ninja \
 && ninja install

ENV VTK_DIR="/home/felix/vtk/install/lib/cmake/vtk-9.5"
ENV CMAKE_PREFIX_PATH="/home/felix/vtk/install:${CMAKE_PREFIX_PATH}"
ENV LD_LIBRARY_PATH="/home/felix/vtk/install/lib:${LD_LIBRARY_PATH}"

# I do not understand why but the cmake files are lower AND upper case
RUN ln -s $VTK_DIR/VTK-targets.cmake $VTK_DIR/vtk-targets.cmake
RUN ln -s $VTK_DIR/VTK-vtk-module-properties.cmake $VTK_DIR/vtk-vtk-module-properties.cmake
RUN ln -s $VTK_DIR/VTK-vtk-module-find-packages.cmake $VTK_DIR/vtk-vtk-module-find-packages.cmake

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

# ===== GSLIB for MFEM 4.8 (correct layout) =====
ENV GSLIB_TAG=v1.0.9
ENV GSLIB_SRC=/home/felix/gslib_src
ENV GSLIB_PREFIX=/home/felix/gslib

# ===== GSLIB for MFEM 4.8 (copy everything relevant) =====
ENV GSLIB_TAG=v1.0.9
ENV GSLIB_SRC=/home/felix/gslib_src
ENV GSLIB_PREFIX=/home/felix/gslib

RUN curl -L https://github.com/gslib/gslib/archive/${GSLIB_TAG}.tar.gz \
    | tar -xz -C /home/felix && \
    mv /home/felix/gslib-1.0.9 ${GSLIB_SRC} && \
    cd ${GSLIB_SRC} && \
    make clean && \
    make CC=mpicc && \
    \
    mkdir -p ${GSLIB_PREFIX} && \
    \
    # Copy the full build outputs (generated headers include config.h)
    cp -a ${GSLIB_SRC}/build/include ${GSLIB_PREFIX}/ && \
    cp -a ${GSLIB_SRC}/build/lib ${GSLIB_PREFIX}/ && \
    \
    # Also copy the wrapper header used by MFEM's FindGSLIB logic / mfem glue
    cp -v ${GSLIB_SRC}/src/gslib.h ${GSLIB_PREFIX}/include/ && \
    cp -v ${GSLIB_SRC}/build/include/gslib/* ${GSLIB_PREFIX}/include/  && \
    \
    # MFEM expects libgslib.a; gslib builds libgs.a
    ln -sf ${GSLIB_PREFIX}/lib/libgs.a ${GSLIB_PREFIX}/lib/libgslib.a && \
    \
    # Sanity
    test -f ${GSLIB_PREFIX}/include/gslib.h && \
    test -f ${GSLIB_PREFIX}/include/gslib/gslib.h && \
    test -f ${GSLIB_PREFIX}/include/gslib/config.h && \
    test -f ${GSLIB_PREFIX}/lib/libgs.a



# ===== MFEM C++ build (now with MPI + OpenMP + ParMETIS + SuperLU_DIST) =====
# FInd SuperLU libararies and turn on 
ENV MFEM_TAG=v4.8
ENV CC=mpicc CXX=mpicxx

RUN git clone --depth 1 --branch ${MFEM_TAG} https://github.com/mfem/mfem.git /home/felix/mfem && \
    mkdir -p /home/felix/mfem/build && cd /home/felix/mfem/build && \
    cmake -S .. -B . \
      -DCMAKE_BUILD_TYPE=Debug \
      -DMFEM_USE_MPI=ON \
      -DMFEM_USE_OPENMP=ON \
      -DMFEM_USE_METIS=ON \
      -DMFEM_USE_GSLIB=ON \
      -DGSLIB_DIR=/home/felix/gslib \
      -DMFEM_USE_SUPERLU=OFF \
      -DMFEM_USE_PETSC=OFF \
      -DMFEM_USE_GZSTREAM=ON \
      -DMFEM_USE_SUNDIALS=OFF \
      -DCMAKE_INSTALL_PREFIX=/home/felix/mfem-install && \
    cmake --build . --target install -j"$(nproc)"

ENV PATH=/home/felix/mfem-install/bin:${PATH}
ENV LD_LIBRARY_PATH=/home/felix/mfem-install/lib:${LD_LIBRARY_PATH}
ENV CPLUS_INCLUDE_PATH=/home/felix/mfem-install/include:${CPLUS_INCLUDE_PATH}
ENV LIBRARY_PATH=/home/felix/mfem-install/lib:${LIBRARY_PATH}


# PyMFEM
RUN pip install --user mfem

# GLvis
RUN git clone --depth 1 https://github.com/glvis/glvis.git /home/felix/glvis && \
    cd /home/felix/glvis && \
    make MFEM_DIR=/home/felix/mfem/build && \
    mkdir -p /home/felix/.local/bin && \
    cp glvis /home/felix/.local/bin/


# Nelder-Mead header keeps being annoying
RUN git clone --depth 1 https://github.com/YibaiMeng/nelder-mead.git /home/felix/third_party/nelder-mead


WORKDIR /home/felix/work
EXPOSE 8888

