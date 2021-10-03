import pymem3dg as dg

# icoFace, icoVertex = dg.getIcosphere(1,5)
# icoFace, icoVertex = dg.getCylinder(1,20, 10)
mesh = "./frame7.ply"

g = dg.System(mesh)

g.saveRichData("frame7_geometry.obj", True)