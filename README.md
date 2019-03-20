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
make test
```

Note: Some of the tests require sample netCDF files that are not currently included and will thus fail.

## Using the reader

The reader is a "server side" plugin, so it needs to be accessed by the ParaView Server (or by the built-in server if you simply run the ParaView GUI). Using the ParaView GUI, navigate to "Tools > Manage Plugins...". If you only run the GUI, look for the LFRic reader plugin in the "Local Plugins" section, or in the "Remote Plugins" section if you run a separate ParaView server. If the reader is not shown in the list, click on "Load New...", navigate to the reader installation directory, and select the shared library ```.so``` file. The reader plugin should now appear in the list.

Click on "Load Selected" after selecting the reader plugin, if it is not already loaded. You can also tick the "Auto Load" checkbox in the plugin details section, so that the reader will be automatically available next time you launch ParaView.

ParaView should now allow you to choose the LFRic netCDF file reader when you open a netCDF file.

## Debugging the reader

The reader will send error messages to the error console in the ParaView GUI, and to the command line (if available). If you want to see more detailed output to track down a problem, you can also rebuild the reader with `-DCMAKE_BUILD_TYPE=Debug`.

## Features and limitations

The reader currently includes the following features:
* Visualisation of LFRic fields that are stored in the "full-level" and "half-level" grids
* Visualisation using Cartesian coordinates and longitude-latitude-radius coordinates
* Automatic handling of periodic grids (vertices are automatically mirrored)

Limitations:
* Fields on "full-level" grids are interpolated onto "half-level" grids using nearest-neighbour averaging, resulting in cell-based VTK data
* Fields that are stored on "edge" grids are not yet fully supported
* Vector fields are not yet recognised and appear as scalar component fields
* Vector components are not yet transformed to Cartesian coordinates