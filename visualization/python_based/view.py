import pymem3dg as dg
import numpy as np
import polyscope as ps

outputMesh = "stress_frame0.ply"
vmin = -5
vmax = 10

ps.init()
face, vertex = dg.readMesh(outputMesh)
mesh = ps.register_surface_mesh("frame0", vertex, face)
mesh.add_scalar_quantity("mean_curvature", dg.readData(
    outputMesh, 'vertex', 'mean_curvature'), enabled=True, vminmax=(vmin, vmax))

ps.show()