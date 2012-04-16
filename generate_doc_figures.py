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
camera.SetPosition(4, -6, 5)
camera.SetViewUp(-0.19, 0.4, 0.9)
camera.SetViewAngle(30)
camera.SetFocalPoint(1.5, 1.5, 1)
Render()
WriteImage("Figure/tutorial1.png")
Hide(grid)
remove("tutorial1.vtu")

# tutorial 2
call("tutorials/tutorial2")
grid = servermanager.sources.XMLUnstructuredGridReader(FileName="tutorial2.vtu")
grid.UpdatePipeline()
Show(grid)
dp = GetDisplayProperties(grid)
dp.Representation = 'Surface'
dp.ColorArrayName = 'pressure'
pres = grid.CellData.GetArray(0)
pres_lookuptable = GetLookupTableForArray( "pressure", 1, RGBPoints=[pres.GetRange()[0], 1, 0, 0, pres.GetRange()[1], 0, 0, 1] )
dp.LookupTable = pres_lookuptable
view = GetActiveView()
view.Background = [1, 1, 1]
camera = GetActiveCamera()
camera.SetPosition(20, 20, 110)
camera.SetViewUp(0, 1, 0)
camera.SetViewAngle(30)
camera.SetFocalPoint(20, 20, 0.5)
Render()
WriteImage("Figure/tutorial2.png")
Hide(grid)
# remove("tutorial2.vtu")

# # tutorial 3
call("tutorials/tutorial3")
cases = ["000", "005", "010", "015", "019"]
for case in cases:
    grid = servermanager.sources.XMLUnstructuredGridReader(FileName="tutorial3-"+case+".vtu")
    grid.UpdatePipeline()
    Show(grid)
    dp = GetDisplayProperties(grid)
    dp.Representation = 'Surface'
    dp.ColorArrayName = 'saturation'
    sat = grid.CellData.GetArray(1)
    sat_lookuptable = GetLookupTableForArray( "saturation", 1, RGBPoints=[0, 1, 0, 0, 1, 0, 0, 1])
    dp.LookupTable = sat_lookuptable
    view.Background = [1, 1, 1]
    camera = GetActiveCamera()
    camera.SetPosition(100, 100, 550)
    camera.SetViewUp(0, 1, 0)
    camera.SetViewAngle(30)
    camera.SetFocalPoint(100, 100, 5)
    Render()
    WriteImage("Figure/tutorial3-"+case+".png")
Hide(grid)
#  remove("tutorial3.vtu")

