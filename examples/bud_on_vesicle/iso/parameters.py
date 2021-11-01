import pymem3dg as dg
def parameters():
    p = dg.Parameters()

    p.proteinMobility = 3
    p.temperature = 0

    p.point.pt = [0, 0, 1]
    p.point.isFloatVertex = True
    
    p.proteinDistribution.protein0 = [0.1]
    
    p.boundary.shapeBoundaryCondition = "none"
    p.boundary.proteinBoundaryCondition = "none"
    
    p.variation.isProteinVariation = True
    p.variation.isShapeVariation = True
    p.variation.radius = -1
    
    p.bending.Kb = 8.22e-5
    p.bending.Kbc = 0  # 8.22e-4 #DEFINITION OF LARGE AND SMALL VALUE
    p.bending.H0c = 10
    
    p.tension.isConstantSurfaceTension = False
    p.tension.Ksg = 1
    p.tension.A_res = 0
    p.tension.At = 12.5025
    p.tension.lambdaSG = 0
    
    p.adsorption.epsilon = -1e-3
    
    p.osmotic.isPreferredVolume = True
    p.osmotic.isConstantOsmoticPressure = False
    p.osmotic.Kv = 0.5
    p.osmotic.V_res = 0
    p.osmotic.n = 1
    p.osmotic.Vt = 0.95 * 4.15889  # 1 * 4 * 3.1416 / 3
    p.osmotic.cam = -1
    p.osmotic.lambdaV = 0
    
    p.dirichlet.eta = 0.1
    
    p.dpd.gamma = 0
    
    p.external.Kf = 0
    return p;


if __name__ == '__main__':
    p = parameters()
