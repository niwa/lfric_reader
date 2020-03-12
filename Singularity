Bootstrap: docker
From: ubuntu:18.04

%environment

    # Required for auto-loading the netCDFLFRicReader plugin
    export LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib/netCDFLFRicReader
    export PV_PLUGIN_PATH=/usr/local/lib/netCDFLFRicReader

%post -c /bin/bash

    apt-get update && apt-get upgrade -y

    # Set timezone to UTC non-interactively
    ln -fs /usr/share/zoneinfo/UTC /etc/localtime
    apt-get install -y tzdata
    dpkg-reconfigure --frontend noninteractive tzdata

    apt-get install -y xz-utils pkg-config python3-dev python3-numpy python3-matplotlib \
                       wget make cmake libmpich-dev libosmesa6-dev git

    # Select ParaView major.minor version and patch to form "major.minor.patch"
    pv_version=5.8
    pv_patch=0
    pv_sha256=219e4107abf40317ce054408e9c3b22fb935d464238c1c00c0161f1c8697a3f9

    # Download and verity ParaView sources
    pv_url="https://www.paraview.org/paraview-downloads/download.php"
    pv_url_params="submit=Download&version=v${pv_version}&type=source&os=Sources&downloadFile="
    pv_filename=ParaView-v${pv_version}.${pv_patch}.tar.xz
    wget "${pv_url}?${pv_url_params}${pv_filename}" -O ${pv_filename}
    echo "${pv_sha256} ${pv_filename}" | sha256sum --check

    # Build ParaView without GUI and X ("headless") - this will provide ParaView Server
    # and ParaView Python
    tar -xf ${pv_filename}
    mkdir -p ParaView-v${pv_version}.${pv_patch}/build
    cd ParaView-v${pv_version}.${pv_patch}/build
    cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=/usr/local \
    -DPARAVIEW_BUILD_EDITION=CANONICAL \
    -DPARAVIEW_BUILD_SHARED_LIBS=ON \
    -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON \
    -DPARAVIEW_USE_QT=OFF \
    -DPARAVIEW_USE_PYTHON=ON \
    -DPARAVIEW_USE_MPI=ON \
    -DPARAVIEW_USE_VTKM=OFF \
    -DPARAVIEW_ENABLE_RAYTRACING=OFF \
    -DPARAVIEW_ENABLE_VISITBRIDGE=OFF \
    -DPARAVIEW_PLUGIN_ENABLE_NetCDFTimeAnnotationPlugin=ON \
    -DPARAVIEW_PLUGIN_ENABLE_StreamLinesRepresentation=ON \
    -DPARAVIEW_PLUGIN_ENABLE_StreamingParticles=ON \
    -DPARAVIEW_PLUGIN_ENABLE_SurfaceLIC=ON \
    -DVTK_USE_X=OFF \
    -DVTK_OPENGL_HAS_OSMESA=ON \
    -DVTK_SMP_IMPLEMENTATION_TYPE=OpenMP \
    -DVTK_MODULE_ENABLE_VTK_GeovisCore=YES \
    -DVTK_MODULE_ENABLE_VTK_IONetCDF=YES
    make -j 4
    make install
    cd ../..

    # Build netCDFLFRicReader plugin
    git clone https://github.com/tinyendian/lfric_reader.git
    mkdir -p lfric_reader/src/cxx/build
    cd lfric_reader/src/cxx/build
    cmake .. -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=/usr/local
    make -j 4
    make install
