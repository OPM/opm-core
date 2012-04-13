from subprocess import call
from paraview.simple import *
from paraview import servermanager
from os import remove

connection = servermanager.Connect()

# tutorial 1
call("tutorials/tutorial1")
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
camera.SetPosition(20, 20, 110)
camera.SetViewUp(0.5,0.3,0.7)
camera.SetViewAngle(30)
camera.SetFocalPoint(1,1,0.5)
Render()
WriteImage("Figure/tutorial1.png")
Hide(grid)
# remove("tutorial1.vtu")

# tutorial 2
call("tutorials/tutorial2")
grid = servermanager.sources.XMLUnstructuredGridReader(FileName="tutorial2.vtu")
grid.UpdatePipeline()
Show(grid)
dp = GetDisplayProperties(grid)
dp.Representation = 'Surface'
data = grid.GetCellDataInformation()
pressure = data.GetArray('pressure')
pmin = pressure.GetRange()[0]
pmax = pressure.GetRange()[1]
dp.LookupTable = MakeBlueToRedLT(0.5*(pmin+pmax)-0.2*(pmax-pmin), 0.5*(pmin+pmax)-0.2*(pmax+pmin))
dp.ColorArrayName = 'pressure'
view = GetActiveView()
view.Background = [1, 1, 1]
camera = GetActiveCamera()
camera.SetPosition(20, 20, 110)
camera.SetViewUp(0, 1, 0)
camera.SetViewAngle(30)
camera.SetFocalPoint(20, 20, 0.5)
Render()
WriteImage("Figure/tutorial2.png")
# remove("tutorial2.vtu")
