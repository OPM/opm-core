from subprocess import call
from paraview.simple import *
from paraview import servermanager
from os import remove, mkdir, curdir
from os.path import join, isdir

figure_path = curdir
tutorial_data_path = curdir
tutorial_path = "tutorials"

if not isdir(figure_path):
    mkdir(figure_path)
    
connection = servermanager.Connect()

# tutorial 1
call(join(tutorial_path,"tutorial1"))
grid = servermanager.sources.XMLUnstructuredGridReader(FileName=join(tutorial_data_path,"tutorial1.vtu"))
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
WriteImage(join(figure_path,"tutorial1.png"))
Hide(grid)
remove(join(tutorial_data_path,"tutorial1.vtu"))

# tutorial 2
call(join(tutorial_path,"tutorial2"))
grid = servermanager.sources.XMLUnstructuredGridReader(FileName=join(tutorial_data_path,"tutorial2.vtu"))
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
WriteImage(join(figure_path,"tutorial2.png"))
Hide(grid)
remove(join(tutorial_data_path,"tutorial2.vtu"))

# tutorial 3
call(join(tutorial_path,"tutorial3"))
cases = ["000", "005", "010", "015", "019"]
for case in cases:
    grid = servermanager.sources.XMLUnstructuredGridReader(FileName=join(tutorial_data_path,"tutorial3-"+case+".vtu"))
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
    WriteImage(join(figure_path,"tutorial3-"+case+".png"))
Hide(grid)
# remove(join(tutorial_data_path,"tutorial3.vtu"))

