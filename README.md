# LFRic netCDF file reader

This project implements a file reader for LFRic output files in netCDF UGRID format.

## Building the reader

To build the reader, you need to have an installation of ParaView on your system. The directories with the ParaView installation need to allow write access as a plugin description file will be written at build time.

Use the following commands to set up the build system:
```
PV_ROOT_PATH=/path/to/paraview/installation
export PATH=$PV_ROOT_PATH/bin:$PATH
export CPATH=$PV_ROOT_PATH/include${CPATH:+:$CPATH}
export LD_LIBRARY_PATH=$PV_ROOT_PATH/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export LIBRARY_PATH=$PV_ROOT_PATH/lib${LIBRARY_PATH:+:$LIBRARY_PATH}

cd lfric_reader
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

## Using the reader

You can either install the reader in your ParaView installation, or install it in a place that can be accessed by a ParaView. The reader is a "server side" plugin. Using the ParaView GUI, navigate to "Tools > Manage Plugins...". If you only run the GUI, look for the LFRic reader plugin in the "Local Plugins" section, or in the "Remote Plugins" section if you run a separate ParaView server. If the reader is not shown in the list, click on "Load New...", navigate to the reader installation directory, and select the shared library ```.so``` file. The reader plugin should now appear in the list.

Click on "Load Selected" after selecting the reader plugin, if it is not already loaded. You can also tick the "Auto Load" checkbox in the plugin details section, so that the reader becomes automatically available at runtime.

ParaView should now allow you to choose the LFRic netCDF file reader when you open a netCDF file.

## Debugging the reader

The reader will send error messages to the error console in the ParaView GUI, and to the command line if available. You can get additional debugging output by setting the `DEBUG` macro at build time.

## Features and limitations

The reader currently includes the following features:
* LFRic fields that are stored in the "full-level" and "half-level" grids
* Visualisation using Cartesian coordinates and longitude-latitude-radius coordinates
* Automatic handling of periodic grids (vertices are automatically mirrored)

Limitations:
* Fileds on "full-level" grids are interpolated onto "half-level" grids using nearest-neighbour averaging, resulting in cell-based VTK data
* Fields that are stored on "edge" grids are not yet fully supported
