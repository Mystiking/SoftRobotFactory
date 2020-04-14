import numpy as np
import pyvista as pv
import tetgen
import trimesh

def ComputeSurfaceMesh(tet, ns, elems):
    # Compute surface of the Tet mesh
    objmesh = trimesh.Trimesh(vertices=tet.node, faces=tet.f)
    facetBounds = []
    for facet, facetn in zip(objmesh.facets, objmesh.facets_normal):
        bounds = np.array([np.inf, np.inf, np.inf, -np.inf, -np.inf, -np.inf])
        for f in facet:
            for v in objmesh.vertices[objmesh.faces[f]]:
                bounds[0] = min(bounds[0], v[0])
                bounds[1] = min(bounds[1], v[1])
                bounds[2] = min(bounds[2], v[2])
                bounds[3] = max(bounds[3], v[0])
                bounds[4] = max(bounds[4], v[1])
                bounds[5] = max(bounds[5], v[2])
            facetBounds.append(bounds)
        triPoints = [objmesh.vertices[objmesh.faces[f][0]] for f in facet] 
        triPoint = objmesh.vertices[objmesh.faces[facet[0]][0]]

    isWithinBounds = lambda v, bounds: bounds[0] <= v[0] and v[0] <= bounds[3] and bounds[1] <= v[1] and v[1] <= bounds[4] and bounds[2] <= v[2] and v[2] <= bounds[5]
    isTriangleWithinBounds = lambda t, bounds: isWithinBounds(t[0], bounds) and isWithinBounds(t[1], bounds) and isWithinBounds(t[2], bounds)

    surfaceVerts = []
    surfaceTris = []
    for t in elems:
        tris = [[t[0], t[2], t[1]],
                        [t[1], t[2], t[3]],
                        [t[0], t[1], t[3]],
                        [t[0], t[3], t[2]]]
        for tri in tris:
            vs = np.array([[ns[tri[0]][0], ns[tri[0]][1], ns[tri[0]][2]], [ns[tri[1]][0], ns[tri[1]][1], ns[tri[1]][2]], [ns[tri[2]][0], ns[tri[2]][1], ns[tri[2]][2]]])
            withinBounds = False
            for bounds in facetBounds:
                withinBounds = withinBounds or isTriangleWithinBounds(vs, bounds)
            if withinBounds:
                index = []
                 #surfaceVerts.append(ns[tri[0]].tolist())
                 #print(ns[tri[0]].tolist() in surfaceVerts)
                if not (ns[tri[0]].tolist() in    surfaceVerts):
                    surfaceVerts.append(ns[tri[0]].tolist())
                    index.append(len(surfaceVerts)-1)
                else:
                    index.append(surfaceVerts.index(ns[tri[0]].tolist()))
                if not (ns[tri[1]].tolist() in    surfaceVerts):
                    surfaceVerts.append(ns[tri[1]].tolist())
                    index.append(len(surfaceVerts)-1)
                else:
                    index.append(surfaceVerts.index(ns[tri[1]].tolist()))
                if not (ns[tri[2]].tolist() in    surfaceVerts):
                    surfaceVerts.append(ns[tri[2]].tolist())
                    index.append(len(surfaceVerts)-1)
                else:
                    index.append(surfaceVerts.index(ns[tri[2]].tolist()))
                surfaceTris.append(index)
    surfaceMesh = trimesh.Trimesh(vertices=np.array(surfaceVerts), faces=np.array(surfaceTris))
    return surfaceMesh

def triArea(v0, v1, v2):
    l0 = np.linalg.norm(v1 - v0)
    l1 = np.linalg.norm(v2 - v1)
    l2 = np.linalg.norm(v0 - v2)
    p = (l0 + l1 + l2) / 2
    return np.sqrt(p * (p - l0) * (p - l1) * (p - l2))

def isInside(triangle, point):
    Ainv = 1.0 / triArea(triangle[0], triangle[1], triangle[2])
    A1 = triArea(point, triangle[1], triangle[2]) * Ainv
    A2 = triArea(triangle[0], point, triangle[2]) * Ainv
    A3 = triArea(triangle[0], triangle[1], point) * Ainv
    return 0 <= A1 and A1 <= 1 and 0 <= A2 and A2 <= 1 and 0 <= A3 and A3 <= 1 and abs((A1 + A2 + A3) - 1) < 1e-6

def computeBarycentric(triangle, point):
    Ainv = 1.0 / triArea(triangle[0], triangle[1], triangle[2])
    A1 = triArea(point, triangle[1], triangle[2]) * Ainv
    A2 = triArea(triangle[0], point, triangle[2]) * Ainv
    A3 = triArea(triangle[0], triangle[1], point) * Ainv
    return A1, A2, A3

def ComputeCablePoints(cables, surfaceMesh):
    cableLengths = []
    for c in cables:
        cableLengths.append(sum([np.linalg.norm(np.array(c[i]) - np.array(c[i+1])) for i in range(len(c) - 1)]))

    barycentricCoordinates = []
    cablepointParticles = []

    minDist = np.inf
    minTri = -1
    minTriIndex = -1
    inside = False
    for cablePoints in cables:
        bC = []
        cP = []
        for p in cablePoints:
            for triIndex in surfaceMesh.faces:
                tri = [np.array(surfaceMesh.vertices[triIndex[0]]), np.array(surfaceMesh.vertices[triIndex[1]]), np.array(surfaceMesh.vertices[triIndex[2]])]
                dist = min(np.linalg.norm(np.array(surfaceMesh.vertices[triIndex[0]]) - np.array(p)), min(np.linalg.norm(np.array(surfaceMesh.vertices[triIndex[1]]) - np.array(p)), np.linalg.norm(np.array(surfaceMesh.vertices[triIndex[2]]) - np.array(p))))
                if dist < minDist:
                    minDist = dist
                    minTri = tri
                    minTriIndex = triIndex
                if isInside(tri, np.array(p)):
                    inside = True
                    bC.append(computeBarycentric(tri, np.array(p)))
                    cP.append(triIndex)
                    break
            if not inside:
                bC.append(computeBarycentric(minTri, np.array(p)))
                cP.append(minTriIndex)
            inside = False
        barycentricCoordinates.append(bC)
        cablepointParticles.append(cP)
    return barycentricCoordinates, cablepointParticles, cableLengths

def ComputeEndEffectorPoint(p, surfaceMesh):
    minDist = np.inf
    minTri = -1
    minTriIndex = -1
    inside = False
    bc = -1
    cp = -1
    for triIndex in surfaceMesh.faces:
        tri = [np.array(surfaceMesh.vertices[triIndex[0]]), np.array(surfaceMesh.vertices[triIndex[1]]), np.array(surfaceMesh.vertices[triIndex[2]])]
        dist = min(np.linalg.norm(np.array(surfaceMesh.vertices[triIndex[0]]) - np.array(p)), min(np.linalg.norm(np.array(surfaceMesh.vertices[triIndex[1]]) - np.array(p)), np.linalg.norm(np.array(surfaceMesh.vertices[triIndex[2]]) - np.array(p))))
        if dist < minDist:
            minDist = dist
            minTri = tri
            minTriIndex = triIndex
        if isInside(tri, np.array(p)):
            inside = True
            bc = computeBarycentric(tri, np.array(p))
            cp = triIndex
            break
    if not inside:
        bc = computeBarycentric(minTri, np.array(p))
        cp = minTriIndex
    return bc, cp