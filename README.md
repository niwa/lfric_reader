# LFRic netCDF file reader

This project implements a file reader for LFRic output files in netCDF UGRID format.

## Building the reader

To build the reader, you need to have an installation of ParaView on your system. The directories with the ParaView installation need to allow write access as a plugin description file will need to be written at build time. You can either install the reader directly in your ParaView installation, or install it in any other place that can be accessed by ParaView (or ParaView Server, if you run a client-server configuration).

Use the following commands to set up the build system and build the reader:
```
PV_ROOT_PATH=/path/to/paraview/installation
export PATH=$PV_ROOT_PATH/bin:$PATH
export CPATH=$PV_ROOT_PATH/include${CPATH:+:$CPATH}
export LD_LIBRARY_PATH=$PV_ROOT_PATH/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export LIBRARY_PATH=$PV_ROOT_PATH/lib${LIBRARY_PATH:+:$LIBRARY_PATH}

cd src/cxx
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/install/dir
make
```

## Testing the reader

To run the built-in tests, add the flag ```-DBUILD_TESTING=ON``` in CMake, rebuild the code with ```make```, and type command
```
ctest
```

## Using the reader

The reader is a "server side" plugin, so it needs to be accessed by the ParaView Server (or by the built-in server if you simply run the ParaView GUI). Using the ParaView GUI, navigate to "Tools > Manage Plugins...". If you only run the GUI, look for the LFRic reader plugin in the "Local Plugins" section, or in the "Remote Plugins" section if you run a separate ParaView server. If the reader is not shown in the list, click on "Load New...", navigate to the reader installation directory, and select the shared library ```.so``` file. The reader plugin should now appear in the list.

Click on "Load Selected" after selecting the reader plugin, if it is not already loaded. You can also tick the "Auto Load" checkbox in the plugin details section, so that the reader will be automatically available next time you launch ParaView.

ParaView should now allow you to choose the LFRic netCDF file reader when you open a netCDF file.

## Python scripting

The LFRic reader can easily by used in Python scripts like any other ParaView reader. An example is included in this repository.

## Debugging the reader

The reader will send error messages to the error console in the ParaView GUI, and to the command line (if available). If you want to see more detailed output to track down a problem, you can also rebuild the reader with `-DCMAKE_BUILD_TYPE=Debug`.

## Features and limitations

The reader currently includes the following features:
* Visualisation of LFRic fields that are stored in the "full-level" and "half-level" grids
* Fields that are stored on "edge-half-level" grids can be visualised as points and glyphs (e.g., arrows)
* Visualisation using Cartesian coordinates and longitude-latitude-radius coordinates
* Automatic handling of periodic grids (vertices are automatically mirrored)
* Support for very large data files using parallel operation

Limitations:
* Fields on "full-level" grids are interpolated onto "half-level" grids using nearest-neighbour averaging, resulting in cell-based VTK data

## Parallelisation

The LFRic reader supports MPI-based parallelisation with VTK's parallel pipelines, allowing users to visualise very large datafiles on high-performance systems. This requires using a separate ParaView Server. When the server is launched on more than one MPI process, the reader will automatically only load a subset of the data from the netCDF file on each process. This enables parallel read from disk, reduces memory consumption per process, and parallel pipeline operation for faster computation.

Grids are automatically partitioned along the vertical axis, which is the only structured axis in the data. Users are free to choose any number of MPI processes up to the number of available vertical layers.

To operate the LFRic reader in parallel, you will need to build ParaView with MPI support (```-DPARAVIEW_USE_MPI=ON```) and launch a separate ParaView Server using
```
mpiexec -np <number of processes> pvserver
```
The server will report its port number (the default is 11111) and wait for a connection. Launch the GUI and connect it to the server by clicking on "Connect" in the toolbar. Create a new connection, entering the server hostname, e.g., "localhost" if you run the server on the same machine as the GUI or if you use an ssh tunnel, as well as the server port number. If the server connection is successful, the server name in the "Pipeline Browser" should change from "builtin:" to the server name, and you will see the file space accessible to the server when you try to open a file.

## Singularity Container

The repository includes a Singularity container definition file ```Singularity``` that provides a recent version of ParaView together with the LFRic reader plugin. Note that only a "headless" version of ParaView will currently be built providing ParaView server and ParaView Python. Please refer to the instructions below on how to use the container.

Build the container first using
```
singularity build lfric_reader.sif Singularity
```
Note that this usually requires root privileges. Once container build has finished, you can run a ParaView Python script using
```
singularity exec lfric_reader.sif pvpython myscript.py
```
If you want to use the container with the ParaView GUI, you'll need to download a ParaView executable from paraview.org (make sure that you download the exact same version that is used in the container!) and launch the container using
```
singularity exec lfric_reader.sif pvserver
```
Lauch the GUI and connect to the server from the GUI by clicking on "Connect" in the toolbar. If the container runs on the same host, create a connection to ```localhost``` with port number ```11111``` and connect to the server. If it runs on a different host, you'll need to create an ssh tunnel first and forward port ```11111``` to localhost, then create a connection to ```localhost``` - `11111`.

It is also possible to run the container with MPI parallelisation. If a single host is sufficient, launch the container using
```
singularity exec lfric_reader.sif mpiexec -np <number of ranks> pvserver
```
If more than one host is needed, launch the container with the same MPI distribution that was used to build ParaView inside (currently MPICH - it should be straightforward to switch to, e.g., OpenMPI in the Singularity definition file and rebuild the container),
```
mpiexec -np <number of ranks> singularity exec lfric_reader.sif pvserver
```
This should also work with a job submission system such as Slurm, provided that the underlying MPI distribution matches the one used inside the container.
