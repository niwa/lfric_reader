import os

import paraview.simple as pvs


# Load the LFRic reader plugin - it may be necessary to provide the full path to
# libnetCDFLFRicReader.so, depending on the installation. Set remote=True if you
# connect ParaView with a separate server.
pluginPath = os.environ.get('PV_NC_PLUGIN_PATH', 'libnetCDFLFRicReader.so')
pvs.LoadPlugin(pluginPath, remote=False)

# Create a new data source with the LFRic output data, load the pressure field,
# and use Cartesian coordinates rather than lon-lat-rad
filePath = ''
data = pvs.NetCDFLFRicReader(FileName=filePath+'lfric_output.nc')
data.CellArrayStatus = ['pressure']
data.UseCartesiancoordinates = 1

# Add a clip filter to the pipeline and set plane orientation
clip = pvs.Clip(Input=data)
clip.ClipType = 'Plane'
clip.ClipType.Normal = [0.0, 0.0, 1.0]

# Set up a new render view for displaying the grid and data
renderView = pvs.CreateView('RenderView')
renderView.ViewSize = [1600, 800]
renderView.CameraPosition = [-43.0, 39.0, 32.0]
renderView.CameraViewUp = [0.76, 0.44, 0.49]

# Define a colour look-up table using (data value, r, g, b) tuples
# Use value range from first timestep
lut = pvs.GetColorTransferFunction('pressure')
valueRange = data.CellData['pressure'].GetRange()
lut.RGBPoints = [valueRange[0], 0.23, 0.30, 0.75,
                 valueRange[1], 0.71, 0.02, 0.15]

# Show clip filter in the render view, render grid surfaces,
# and colour surfaces using the pressure field and our LUT
clipDisplay = pvs.Show(clip, renderView)
clipDisplay.Representation = 'Surface With Edges'
clipDisplay.ColorArrayName = ['CELLS', 'pressure']
clipDisplay.LookupTable = lut

# Loop over time steps and render an image for each one
timeSteps = data.TimestepValues
for istep in range(0,len(timeSteps)):

  print("Rendering time step %i" % istep)

  # Set time step and render image
  renderView.ViewTime = timeSteps[istep]
  pvs.Render()

  # Write image file
  istepstr = "%i" % istep
  pvs.SaveScreenshot("lfric_clip_" + istepstr.zfill(2) + ".png")
