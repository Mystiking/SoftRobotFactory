import PySimpleGUI as sg
import numpy as np
import pyvista as pv
import tetgen
import trimesh
import os, sys, time, ast
# Mesh utility functionality
from utils import ComputeSurfaceMesh, ComputeCablePoints, ComputeEndEffectorPoint

supported_files = ['stl', 'obj']

def TextLabel(text):
    return sg.Text(text+':', justification='r', size=(15,1))

sg.theme('DarkAmber')

layout0 = [[sg.Text('Soft Robot Configuration Creation', font='Any 18')],
           [sg.Text('File paths', font='Any 14')],
           [TextLabel('Mesh filename'), sg.Input(key='meshFileName'), sg.FileBrowse(target='meshFileName')],
           [TextLabel('Out folder'), sg.Input(key='outfolder'), sg.FolderBrowse(target='outfolder')],
           [TextLabel('Cfg file'), sg.Input(key='cfgfile'), sg.FileBrowse(target='cfgfile')],
           [sg.Text('TetGen & Misc', font='Any 14')],
           [sg.Frame(layout=[ [TextLabel('Show Robot'), sg.Checkbox('', key='showRobot', size=(10,1), default=True)],
                              [TextLabel('nobisect'), sg.Input('False', key='nobisect')],
                              [TextLabel('mindihedral'), sg.Input(15, key='mindihedral')],
                              [TextLabel('optmaxdihedral'), sg.Input(140, key='optmaxdihedral')],
                              [TextLabel('minratio'), sg.Input(1.5, key='minratio')],
                              [TextLabel('# cables'), sg.Input(0, key='numCables')],
                              [TextLabel('# air chambers'), sg.InputText('0', key='numAirChambers')],
                              [TextLabel('Direchlet Conditions'), sg.Checkbox('', key='direchlet', size=(10,1), default=False)],
                              [TextLabel('Has end-effector'), sg.Checkbox('', key='endEffector', size=(10,1), default=False)]
                            ], title='Options', relief=sg.RELIEF_SUNKEN,)],
           [sg.Button('Create Soft Robot'), sg.Exit()]]

window0 = sg.Window('Soft Robot Creator v.0.1', layout0)

ns = None
elems = None
tet = None
showRobot = False
surfacemesh = None
numCables = 0
numAirChambers = 0
direchlet_conditions = False
endEffector = False
meshname = ''
outfolder = ''
cfgfile = ''
while True:
    event, values = window0.read()
    if event in ('Create Soft Robot'):
        # Sanitize
        if values['meshFileName'] == '' or\
           not os.path.exists(values['meshFileName']) or\
           not (values['meshFileName'][-3:] in supported_files):
            print('Unsupported or non-existant surface mesh file.')
            window0.close()
            exit(0)
        if values['outfolder'] == '':
            print('No output folder specified.')
            window0.close()
            exit(0)
        outfolder = values['outfolder']
        cfgfile = values['cfgfile']
        try:
            values['nobisect'] = bool(values['nobisect'])
        except:
            print('nobisect should be a boolean.')
            window0.close()
            exit(0)
        try:
            values['mindihedral'] = float(values['mindihedral'])
        except:
            print('mindihedral should be a number.')
            window0.close()
            exit(0)
        try:
            values['optmaxdihedral'] = float(values['optmaxdihedral'])
        except:
            print('optmaxdihedral should be a number.')
            window0.close()
            exit(0)
        try:
            values['minratio'] = float(values['minratio'])
        except:
            print('minratio should be a number.')
            window0.close()
            exit(0)
        try:
            numCables = int(values['numCables'])
        except:
            print('numCables should be an integer.')
            window0.close()
            exit(0)
        try:
            numAirChambers = int(values['numAirChambers'])
        except:
            print('numAirChambers should be an integer.')
            window0.close()
            exit(0)
        direchlet_conditions = bool(values['direchlet'])
        endEffector = bool(values['endEffector'])

        # Perform the tetgen tetrahedralise step
        mesh = trimesh.load(values['meshFileName'])
        tet = tetgen.TetGen(mesh.vertices, mesh.faces)
        tet.make_manifold()
        ns, elems = tet.tetrahedralize(order=1,
                                       nobisect=values['nobisect'],
                                       mindihedral=values['mindihedral'],
                                       optmaxdihedral=values['optmaxdihedral'],
                                       minratio=values['minratio'])
        showRobot = values['showRobot']
        meshname = values['meshFileName'].split('/')[-1][:-4]
        break
    if event in (None, 'Exit'):
        window0.close()
        exit(0)

window0.close()

if showRobot:
    surfacemesh = ComputeSurfaceMesh(tet, ns, elems)
    grid = tet.grid
    grid.plot()
    surfacemesh.show()

layout1 = [[sg.Text('Do you want to continue?', font='Any 18')],
           [sg.Button('Continue'), sg.Exit()]]

window1 = sg.Window('Soft Robot Creator v.0.1', layout1)
while True:
    event, values = window1.read()
    if event in ('Continue'):
        break
    if event in (None, 'Exit'):
        window1.close()
        exit(0)

window1.close()

# Cable, air chamber and dirichlet conditions set-up

layout2 = [[sg.Text('Soft Robot Configuration Creation', font='Any 18')]]
if numCables > 0:
    if cfgfile != '':
        cablePointsFromCfg = open(cfgfile).readlines()[0]
        layout2.append([sg.Frame(layout=[
                              [TextLabel('stiffness ([float])'), sg.Input('[' + ('1e3, ' * numCables)[:-2] + ']', key='stiffness')],
                              [TextLabel('damping ([float])'), sg.Input('[' + ('1e-2, ' * numCables)[:-2] + ']', key='damping')],
                              [TextLabel('# cable points ([int])'), sg.Input('[' + ('{0}, '.format(len(ast.literal_eval(cablePointsFromCfg)[0])) * numCables)[:-2] + ']', key='numCablePoints')],
                              [TextLabel('cable points ([[float]])'), sg.Input(cablePointsFromCfg, key='cablePoints')],
                            ], title='Cable Options', relief=sg.RELIEF_SUNKEN,)])
    else:
        layout2.append([sg.Frame(layout=[
                              [TextLabel('stiffness ([float])'), sg.Input('[' + ('1e3, ' * numCables)[:-2] + ']', key='stiffness')],
                              [TextLabel('damping ([float])'), sg.Input('[' + ('1e-2, ' * numCables)[:-2] + ']', key='damping')],
                              [TextLabel('# cable points ([int])'), sg.Input('[' + ('0, ' * numCables)[:-2] + ']', key='numCablePoints')],
                              [TextLabel('cable points ([[float]])'), sg.Input('[' + ('[], ' * numCables)[:-2] + ']', key='cablePoints')],
                            ], title='Cable Options', relief=sg.RELIEF_SUNKEN,)])
if direchlet_conditions:
    if cfgfile != '':
        dcFromCfg = ast.literal_eval(open(cfgfile).readlines()[1])
        layout2.append([sg.Frame(layout=[
                              [TextLabel('xmin ([float])'), sg.Input('[' + str(dcFromCfg[0]) + ']', key='xmin')],
                              [TextLabel('ymin ([float])'), sg.Input('[' + str(dcFromCfg[1]) + ']', key='ymin')],
                              [TextLabel('zmin ([float])'), sg.Input('[' + str(dcFromCfg[2]) + ']', key='zmin')],
                              [TextLabel('xmax ([float])'), sg.Input('[' + str(dcFromCfg[3]) + ']', key='xmax')],
                              [TextLabel('ymax ([float])'), sg.Input('[' + str(dcFromCfg[4]) + ']', key='ymax')],
                              [TextLabel('zmax ([float])'), sg.Input('[' + str(dcFromCfg[5]) + ']', key='zmax')]
                            ], title='Direchlet Conditions', relief=sg.RELIEF_SUNKEN,)])
    else:
        layout2.append([sg.Frame(layout=[
                          [TextLabel('xmin ([float])'), sg.Input('[]', key='xmin')],
                          [TextLabel('ymin ([float])'), sg.Input('[]', key='ymin')],
                          [TextLabel('zmin ([float])'), sg.Input('[]', key='zmin')],
                          [TextLabel('xmax ([float])'), sg.Input('[]', key='xmax')],
                          [TextLabel('ymax ([float])'), sg.Input('[]', key='ymax')],
                          [TextLabel('zmax ([float])'), sg.Input('[]', key='zmax')]
                        ], title='Direchlet Conditions', relief=sg.RELIEF_SUNKEN,)])
if endEffector:
    if cfgfile != '':
        endEffectorFromCfg = open(cfgfile).readlines()[2]
        layout2.append([sg.Frame(layout=[
                              [TextLabel('End effector ([float])'), sg.Input(endEffectorFromCfg, key='endEffectorPoint')],
                            ], title='End Effector', relief=sg.RELIEF_SUNKEN,)])
    else:
        layout2.append([sg.Frame(layout=[
                          [TextLabel('End effector ([float])'), sg.Input('[]', key='endEffectorPoint')],
                        ], title='End Effector', relief=sg.RELIEF_SUNKEN,)])
layout2.append([sg.Button('Create'), sg.Exit()])

# Cables
stiffness = None
damping = None
numCablePoints = None
cablePoints = None
# Air Chambers TBD
# Dirichlet Conditions
xmin = None
ymin = None
zmin = None
xmax = None
ymax = None
zmax = None

window2 = sg.Window('Soft Robot Creator v.0.1', layout2)
while True:
    event, values = window2.read()
    if event in ('Create'):
        # Sanitize
        if numCables > 0:
            try:
                stiffness = ast.literal_eval(values['stiffness'])
            except:
                print('stiffness should be a list of stiffness coefficients.')
                window2.close()
                exit(0)
            try:
                damping = ast.literal_eval(values['damping'])
            except:
                print('damping should be a list of damping coefficients.')
                window2.close()
                exit(0)
            try:
                numCablePoints = ast.literal_eval(values['numCablePoints'])
            except:
                print('# cable points should be a list of the amount of points in each cable.')
                window2.close()
                exit(0)
            try:
                cablePoints = ast.literal_eval(values['cablePoints'])
            except:
                print('cablePoints should be a list of cable point lists.')
                window2.close()
                exit(0)
        if numAirChambers > 0:
            # TODO: Add air-chamber support
            pass
        if direchlet_conditions:
            try:
                xmin = ast.literal_eval(values['xmin'])
            except:
                print('xmin should be a list of min x-values for dirichlet bounding boxes.')
                window2.close()
                exit(0)
            try:
                ymin = ast.literal_eval(values['ymin'])
            except:
                print('ymin should be a list of min y-values for dirichlet bounding boxes.')
                window2.close()
                exit(0)
            try:
                zmin = ast.literal_eval(values['zmin'])
            except:
                print('zmin should be a list of min z-values for dirichlet bounding boxes.')
                window2.close()
                exit(0)
            try:
                xmax = ast.literal_eval(values['xmax'])
            except:
                print('xmax should be a list of max x-values for dirichlet bounding boxes.')
                window2.close()
                exit(0)
            try:
                ymax = ast.literal_eval(values['ymax'])
            except:
                print('ymax should be a list of max y-values for dirichlet bounding boxes.')
                window2.close()
                exit(0)
            try:
                zmax = ast.literal_eval(values['zmax'])
            except:
                print('zmax should be a list of max z-values for dirichlet bounding boxes.')
                window2.close()
                exit(0)
        if endEffector:
            try:
                endEffectorPoint = ast.literal_eval(values['endEffectorPoint'])
            except:
                print('End Effector should be a point on the mesh.')
                window2.close()
                exit(0)
        # Compute remaining info
        if numCables > 0:
            barycentricCoordinates, cableParticles, cableLengths = ComputeCablePoints(cablePoints, surfacemesh)
        if numAirChambers > 0:
            # TODO: Add air-chamber support
            pass
        if endEffector:
            endEffectorBC, endEffectorPs = ComputeEndEffectorPoint(endEffectorPoint, surfacemesh)
        # Save all the correct info
        if not os.path.exists(outfolder):
            os.mkdir(outfolder)
        meshname = outfolder + meshname if outfolder[-1] == '/' else outfolder + '/' + meshname
        ## .tet file (volume mesh)
        with open(meshname + '.tet', 'w') as tetfile:
            for v in ns:
                tetfile.write('v ' + ' '.join([str(val) for val in v]) + '\r\n')
            for t in elems:
                tetfile.write('t ' + ' '.join([str(val) for val in t[:4]]) + '\r\n')
        ## .obj file (Render mesh)
        with open(meshname + '.obj', 'w') as objfile:
            for vi in range(len(surfacemesh.vertices)):
                objfile.write("vn " + " ".join([str(val) for val in surfacemesh.vertex_normals[vi]]) + "\r\n")
                objfile.write("v " + " ".join([str(val) for val in surfacemesh.vertices[vi]]) + "\r\n")
            for fi in range(len(surfacemesh.faces)):
                objfile.write("f " + " ".join([str(val+1) + "//" + str(val+1) for val in surfacemesh.faces[fi]]) + "\r\n")
        ## .objmax file (Used to compute surface -> particle mapping)
        with open(meshname + '.objmax', 'w') as objmaxfile:
            for vi in range(len(surfacemesh.vertices)):
                objmaxfile.write("v " + " ".join([str(val) for val in surfacemesh.vertices[vi]]) + "\r\n")
            for fi in range(len(surfacemesh.faces)):
                objmaxfile.write("f " + " ".join([str(val+1) for val in surfacemesh.faces[fi]]) + "\r\n")
        ## .cfg file (Robot configuration file consisting of cables, air chambers, direchlet conditions and end effectors)
        with open(meshname + '.cfg', 'w') as cfgfile:
            # Writing cable data
            if numCables > 0:
                for c in range(numCables):
                    cfgfile.write("c {0} {1} {2} {3}\n".format(stiffness[c], damping[c], numCablePoints[c], cableLengths[c]))
                    for i in range(numCablePoints[c]):
                        b = barycentricCoordinates[c][i]
                        p = cableParticles[c][i]
                        cfgfile.write("p {0} {1} {2} {3} {4} {5}\n".format(b[0], b[1], b[2], p[0], p[1], p[2]))
            # Writing air-chamber data
            if numAirChambers > 0:
                # TODO: Add air-chamber support
                pass
            # Writing dirichlet conditions
            if direchlet_conditions:
                for vi in range(len(surfacemesh.vertices)):
                    v = surfacemesh.vertices[vi]
                    is_affected = False
                    for j in range(len(xmin)):
                        if xmin[j] <= v[0] and v[0] <= xmax[j]\
                            and ymin[j] <= v[1] and v[1] <= ymax[j]\
                            and zmin[j] <= v[2] and v[2] <=zmax[j]:
                            is_affected = True
                            break
                    if is_affected:
                        cfgfile.write('d {0} {1} {2}\n'.format(vi, 0.0, 0.0)) # Fix that we write bogus values
            # Writing end-effector information
            if endEffector:
                bc = endEffectorBC
                cp = endEffectorPs
                cfgfile.write('e {0} {1} {2} {3} {4} {5}'.format(bc[0], bc[1], bc[2], cp[0], cp[1], cp[2]))
        break
    if event in (None, 'Exit'):
        window1.close()
        exit(0)

window2.close()