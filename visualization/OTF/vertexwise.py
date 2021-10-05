import matplotlib.pyplot as plt
import pymem3dg as dg
import polyscope as ps
import numpy as np
import netCDF4 as nc
import imp


def matplotlibStyle():
    """ Formatting style of matplotlib """
    plt.rcParams['font.sans-serif'] = "Arial"
    plt.rcParams['font.family'] = "sans-serif"
    # mpl.rcParams.update({'font.size': 8})
    SMALL_SIZE = 6
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 9
    plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    plt.rc('pdf', fonttype=42)


def getColorMap(fig, range, orient="horizontal", num_points=100):
    x = np.linspace(0, 1, num=num_points)
    y = np.linspace(0, 1, num=num_points)
    fig.gca().set_visible(False)
    Temp = np.linspace(range[0], range[1], num=num_points)
    plt.scatter(x, y, c=Temp)
    plt.colorbar(orientation=orient)
    plt.tight_layout()
    return fig


def readMeshByPly(ply):
    face, vertex = dg.readMesh(ply)
    return face, vertex


def readMeshDataByPly(ply, dataType, dataName):
    """ example: proteinDensity = readMeshDataByPly(meshList[i], "vertex", "protein_density") """
    data = dg.readData(ply + ".ply",
                       dataType, dataName)
    return data


def readMeshByNc(ncFrame):
    trajNc, frame = ncFrame
    with nc.Dataset(trajNc + ".nc") as ds:
        coordinates = np.array(
            ds.groups['Trajectory'].variables['coordinates'][frame])
        topology = np.array(
            ds.groups['Trajectory'].variables['topology'][frame])
        # -1 to dynamically allocate size
        coordinates = np.reshape(coordinates, (-1, 3))
        topology = np.reshape(topology, (-1, 3))
    return topology, coordinates


def readMeshDataByNc(ncFrame, group, variable, num_col):
    """ example: proteinDensity = readMeshDataByPly(meshList[i], "Trajectory", "proteindensity") """
    trajNc, frame = ncFrame
    with nc.Dataset(trajNc + ".nc") as ds:
        data = np.array(
            ds.groups[group].variables[variable][frame])
        # -1 to dynamically allocate size
        data = np.reshape(data, (-1, num_col))
    return np.squeeze(data)


def constructSystemByNc(parameters, ncFrame):
    trajNc, frame = ncFrame
    return dg.System(trajNc + ".nc", frame, parameters, 0, True)


def rowwiseScaling(scaling, matrix):
    if np.shape(matrix)[0] == np.size(matrix):
        return matrix * scaling
    else:
        return matrix * scaling[:, None]


def rowwiseNorm(matrix):
    if np.shape(matrix)[0] == np.size(matrix):
        return np.abs(matrix)
    else:
        return np.sum(matrix**2, axis=1)**(0.5)


def rowwiseNormalize(matrix):
    # print("norm:", rowwiseNorm(matrix))
    return rowwiseScaling(rowwiseNorm(matrix)**(-1), matrix)


def rowwiseDotProduct(a, b):
    return np.sum(a*b, axis=1)


if __name__ == "__main__":
    """ working directory """
    wd = "../../results/bud_on_vesicle/"

    """ initialize polyscope """
    ps.init()

    """ list all parameters """
    parametersList = [wd + '0p7_neg/parameters.py',
                      wd + '0p7_neg/parameters.py', wd + '0p7/parameters.py']

    """ list all meshes """
    """ based on .ply file """
    # meshList = ['0p6/frame3', '0p6/frame7', '0p8/frame10',
    #                  '0p8/frame20', '1/frame70', '1/frame140', '1p2/frame166', '1p2/frame332']

    """ based on .nc file """
    meshList = [(wd + '0p7_neg/traj', 14),
                (wd + '0p7_neg/traj', 27), (wd + '0p7/traj', 200)]

    """ loop over meshes """
    for i in range(len(meshList)):

        """ construct a tag """
        tag = meshList[i][0] + "{}".format(meshList[i][1])

        """ read mesh """
        face, vertex = readMeshByNc(meshList[i])
        mesh = ps.register_surface_mesh(tag, vertex, face)

        """ read parameters """
        parameterFile = imp.load_source("module.name", parametersList[i])

        """ construct system """
        system = constructSystemByNc(parameterFile.parameters(), meshList[i])

        """ read protein density """
        proteinDensity = system.getProteinDensity()
        limit = (np.min(proteinDensity), np.max(proteinDensity))
        mesh.add_scalar_quantity(
            "proteinDensity", proteinDensity, vminmax=limit)

        """ read spontaneous curvature """
        spontaneousCurvature = system.getSpontaneousCurvature()
        limit = (np.min(spontaneousCurvature), np.max(spontaneousCurvature))
        mesh.add_scalar_quantity(
            "spontaneousCurvature", spontaneousCurvature, vminmax=limit)

        """ compute forces """
        system.computePhysicalForces()

        """ read line capillary force """
        lineCapillaryForce = system.forces.getLineCapillaryForce()
        mesh.add_vector_quantity("lineCapillaryForce", lineCapillaryForce)
        mesh.add_scalar_quantity(
            "lineCapillaryForce_norm", rowwiseNorm(lineCapillaryForce))

        """ read bending force """
        bendingForce = system.forces.getBendingForce()
        mesh.add_vector_quantity("bendingForce", bendingForce)
        mesh.add_scalar_quantity(
            "bendingForce_norm", rowwiseNorm(bendingForce))

        """ color map """
        fig, axs = plt.subplots(1, 1)
        fig.set_size_inches(1, 1)
        fig = getColorMap(fig, limit)
        fig.savefig(tag + ".pdf", transparent=True)
        # plt.subplots_adjust(left=0.164, bottom=0.07, right=0.988, top=0.988)
        # fig.savefig(meshFile +"_horizontal.pdf", transparent=True)

    """ configure Polyscope """
    ps.set_up_dir("z_up")
    ps.show()
