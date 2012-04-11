from subprocess import call
from paraview.simple import *
from paraview import servermanager
from os import remove

call("tutorials/tutorial1")
connection = servermanager.Connect()
grid = servermanager.sources.XMLUnstructuredGridReader(FileName="tutorial1.vtu")
grid.UpdatePipeline()
Show(grid)

dp = GetDisplayProperties(grid)
dp.Representation = 'Wireframe'
dp.LineWidth = 5
dp.AmbientColor = [1, 0, 0]

view = GetActiveView()
view.Background = [1, 1, 1]
camera = GetActiveCamera()
camera.SetViewUp(0.5,0.3,0.7)
camera.SetViewAngle(30)
camera.SetFocalPoint(1,1,0.5)

Render()

WriteImage("Figure/tutorial1.png")

remove("tutorial1.vtu")
