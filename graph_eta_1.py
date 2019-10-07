#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
eta_1_Aelpvd = PVDReader(FileName='/home/mochong/MduanNjw/ActiveEL/Output/eta_1_Ael.pvd')
#eta_1_Aelpvd = PVDReader(FileName='/home/mochong/MduanNjw/ActiveEL/Output/eta_-1_Ael.pvd')
# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [1148, 802]

# show data in view
eta_1_AelpvdDisplay = Show(eta_1_Aelpvd, renderView1)
# trace defaults for the display properties.
eta_1_AelpvdDisplay.Representation = 'Surface'

# reset view to fit data
renderView1.ResetCamera()

# show color bar/color legend
eta_1_AelpvdDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'VelocityPressureDirector1'
velocityPressureDirector1LUT = GetColorTransferFunction('VelocityPressureDirector1')

# create a new 'Glyph'
glyph1 = Glyph(Input=eta_1_Aelpvd,
    GlyphType='Arrow')

# Properties modified on glyph1
glyph1.Scalars = ['POINTS', 'None']
glyph1.Vectors = ['POINTS', "['Velocity', 'Pressure', 'Director'][2]"]
glyph1.ScaleFactor = 0.19600000000000004

# show data in view
glyph1Display = Show(glyph1, renderView1)
# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

animationScene1.GoToNext()

animationScene1.GoToNext()

animationScene1.GoToNext()

animationScene1.GoToNext()

animationScene1.GoToNext()

animationScene1.GoToNext()

animationScene1.GoToNext()

animationScene1.GoToNext()

animationScene1.GoToNext()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [0.5, 3.5, 13.724888285812703]
renderView1.CameraFocalPoint = [0.5, 3.5, 0.0625]
renderView1.CameraParallelScale = 3.536086289953909

#### uncomment the following to render all views
RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
