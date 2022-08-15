import numpy as np
import math
import itk
import PythonVersorRigid3DPerspectiveTransform as T
import matplotlib.pyplot as plt


def GetProjectedPointsTRTP(points, transformClass, size):
    '''
    Project 3D TomoRecon points to 2D points using PythonVersorRigid3DPerspectiveTransform class
    https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection 
    Parameters:
        points (array): nx3 array where n is the number of points to be projected
        transformClass (class): constructed PythonVersorRigid3DPerspectiveTransform class
        size (array): 1x3 array of pixel dimensions of tomosynthesis recon image
    Returns:
        array: nx2 array of projected points where n is the number of points
        array: nx3 array of the point of intersections in 3D space between the emitter and the detection plane
        array: nx1 array of the t value in the parametric equation defines in the link
    '''
    projectedPoints = []
    point_of_intersection = []
    t_vals = []
    for point in points:
        projected, num, t = transformClass.TransformPoint(point)
        projectedPoint = [((projected[0]+(size[0]/2))/0.194),
                         ((projected[1]+(size[1]/2))/0.194)
                         ]
        point_of_intersection = np.append(point_of_intersection, num)
        t_vals = np.append(t_vals, t)
        projectedPoints.append(projectedPoint)
    projectedPoints = np.array(projectedPoints)
    projectedPoints = projectedPoints.reshape((len(points), 2))
    point_of_intersection = point_of_intersection.reshape((len(points), 3))
    return projectedPoints, point_of_intersection, t_vals


def TransformAllPoints3D(translationMatrix, points):
    '''
    This function uses the dot product to transform 3D points given a translation matrix
    Parameters:
        translationMatrix (array): 4x4 transformation matrix
        points (array): nx3 array of 3D points to be transformed
    Returns:
        array: nx3 array of transformed 3D points
    '''
    new_source_points = []
    for point in points:
        point = np.append(point, 1)
        new_point = np.dot(translationMatrix, point.T)
        new_source_points = np.append(new_source_points, new_point[0:3])
    new_source_points = new_source_points.reshape((len(points), 3))
    return new_source_points
 
    
def TransformAllPointsAllParameters(x, points):
    '''
    Function to transform 3D points (preprocess source points)
    transforms based on 6 parameters defined in a vector x
    Parameters: 
        x (array): 1x6 array where the inputs define the x translation, y translation
                   z translation, z rotation, y rotation and x rotation respectively
        points (array): nx3 array of 3D points to be transformed 
    Returns:
        array: nx3 array of transformed 3D points 
    '''
    # rotation matrix in z direction
    rotationMatrixZ = np.array([[math.cos(x[3]), -math.sin(x[3]), 0, 0],
                                [math.sin(x[3]), math.cos(x[3]), 0, 0],
                                [0, 0, 1, 0],
                                [0, 0, 0, 1]
    ])
    # rotation matrix in y direction
    rotationMatrixY = np.array([[math.cos(x[4]), 0, math.sin(x[4]), 0],
                                [0, 1, 0, 0],
                                [-math.sin(x[4]), 0, math.cos(x[4]), 0],
                                [0, 0, 0, 1]
    ])
    # rotation matrix in x direction
    rotationMatrixX = np.array([[1, 0, 0, 0],
                                [0, math.cos(x[5]), -math.sin(x[5]), 0],
                                [0, math.sin(x[5]), math.cos(x[5]), 0],
                                [0, 0, 0, 1]
    ])
    # translation matrix in x, y, z direction
    translationMatrixD = np.array([[1, 0, 0, x[0]],
                                   [0, 1, 0, x[1]],
                                   [0, 0, 1, x[2]],
                                   [0, 0, 0, 1]
    ])
    # use transformPoints3D to transform each point, returns array of all transformed points
    points = TransformAllPoints3D(rotationMatrixZ, points)
    points = TransformAllPoints3D(rotationMatrixY, points)
    points = TransformAllPoints3D(rotationMatrixX, points)
    points = TransformAllPoints3D(translationMatrixD, points)
    return points


def TransformOnePointAllParameters(x, points):
    '''
    Function to transform 3D points (preprocess source points)
    transforms based on 6 parameters defined in a vector x
    Parameters: 
        x (array): 1x6 array where the inputs define the x translation, y translation
                   z translation, z rotation, y rotation and x rotation respectively
        points (array): nx3 array of 3D points to be transformed 
    Returns:
        array: nx3 array of transformed 3D points 
    '''
    # rotation matrix in z direction
    rotationMatrixZ = np.array([[math.cos(x[3]), -math.sin(x[3]), 0, 0],
                                [math.sin(x[3]), math.cos(x[3]), 0, 0],
                                [0, 0, 1, 0],
                                [0, 0, 0, 1]
    ])
    # rotation matrix in y direction
    rotationMatrixY = np.array([[math.cos(x[4]), 0, math.sin(x[4]), 0],
                                [0, 1, 0, 0],
                                [-math.sin(x[4]), 0, math.cos(x[4]), 0],
                                [0, 0, 0, 1]
    ])
    # rotation matrix in x direction
    rotationMatrixX = np.array([[1, 0, 0, 0],
                                [0, math.cos(x[5]), -math.sin(x[5]), 0],
                                [0, math.sin(x[5]), math.cos(x[5]), 0],
                                [0, 0, 0, 1]
    ])
    # translation matrix in x, y, z direction
    translationMatrixD = np.array([[1, 0, 0, x[0]],
                                   [0, 1, 0, x[1]],
                                   [0, 0, 1, x[2]],
                                   [0, 0, 0, 1]
    ])
    # use transformPoints3D to transform each point, returns array of all transformed points
    points = TransformOnePoint3D(rotationMatrixZ, points)
    points = TransformOnePoint3D(rotationMatrixY, points)
    points = TransformOnePoint3D(rotationMatrixX, points)
    points = TransformOnePoint3D(translationMatrixD, points)
    return points


def model(x, u, transformClass, imageSize):
    '''
    Optimization model used to project one point from 3D tomosynthesis reconstruction points
    to 2D coordinates
    Paramters:
        x (array): 1x6 array where the inputs define the x translation, y translation
                   z translation, z rotation, y rotation and x rotation respectively
        u (array): 1x3 array that defines one 3D point
        transformClass (class): constructed PythonVersorRigid3DPerspectiveTransform class
        imageSize (array): 1x3 array of pixel dimensions of tomosynthesis recon image
    Returns:
        array: 1x2 array that defines the projected 2D point
    '''
    # rotation matrix in z direction
    rotationMatrixZ = np.array([[math.cos(x[3]), -math.sin(x[3]), 0, 0],
                                [math.sin(x[3]), math.cos(x[3]), 0, 0],
                                [0, 0, 1, 0],
                                [0, 0, 0, 1]
    ])
    # rotation matrix in y direction
    rotationMatrixY = np.array([[math.cos(x[4]), 0, math.sin(x[4]), 0],
                                [0, 1, 0, 0],
                                [-math.sin(x[4]), 0, math.cos(x[4]), 0],
                                [0, 0, 0, 1]
    ])
    # rotation matrix in x direction
    rotationMatrixX = np.array([[1, 0, 0, 0],
                                [0, math.cos(x[5]), -math.sin(x[5]), 0],
                                [0, math.sin(x[5]), math.cos(x[5]), 0],
                                [0, 0, 0, 1]
    ])
    # translation matrix in x, y, z direction
    translationMatrixD = np.array([[1, 0, 0, x[0]],
                                   [0, 1, 0, x[1]],
                                   [0, 0, 1, x[2]],
                                   [0, 0, 0, 1]
    ])
    # use transformPoints3D function to transform each point in a pointset using a tranformation matrix
    point = TransformOnePoint3D(rotationMatrixZ, u)
    point = TransformOnePoint3D(rotationMatrixY, point)
    point = TransformOnePoint3D(rotationMatrixX, point)
    point = TransformOnePoint3D(translationMatrixD, point)
    projected, _, _ = transformClass.TransformPoint(point)
    point = [((projected[0]+(imageSize[0]/2))/0.194),
             ((projected[1]+(imageSize[1]/2))/0.194)
            ]
    return point


def TransformOnePoint3D(translationMatrix, point):
    '''
    This function uses the dot product to transform one 3D point given a translation matrix
    Parameters:
        translationMatrix (array): 4x4 transformation matrix
        points (array): 1x3 array of 3D point to be transformed
    Returns
        array: 1x3 array of 3D transformed point
    '''
    point = np.append(point, 1)
    new_point = np.dot(translationMatrix, point.T)
    new_source_point = new_point[0:3]
    return new_source_point


def modelCT(x, u, transformClass, size, spacing):
    '''
    Optimization model used to project one point from 3D CT points to 2D coordinates
    Paramters:
        x (array): 1x6 array where the inputs define the x translation, y translation
                    z translation, z rotation, y rotation and x rotation respectively
        u (array): 1x3 array that defines one 3D point
        transformClass (class): constructed PythonVersorRigid3DPerspectiveTransform class
        imageSize (array): 1x3 array of pixel dimensions of tomo image
        spacing (float): value that decribes the x, y spacing of the CT scan for the purpose of scaling the source points
    Returns:
        array: 1x2 array that defines the projected 2D point
    '''
    # rotation matrix in z direction
    rotationMatrixZ = np.array([[math.cos(x[3]), -math.sin(x[3]), 0, 0],
                                [math.sin(x[3]), math.cos(x[3]), 0, 0],
                                [0, 0, 1, 0],
                                [0, 0, 0, 1]
    ])
    # rotation matrix in y direction
    rotationMatrixY = np.array([[math.cos(x[4]), 0, math.sin(x[4]), 0],
                                [0, 1, 0, 0],
                                [-math.sin(x[4]), 0, math.cos(x[4]), 0],
                                [0, 0, 0, 1]
    ])
    # rotation matrix in x direction
    rotationMatrixX = np.array([[1, 0, 0, 0],
                                [0, math.cos(x[5]), -math.sin(x[5]), 0],
                                [0, math.sin(x[5]), math.cos(x[5]), 0],
                                [0, 0, 0, 1]
    ])
    # translation matrix in x, y, z direction
    translationMatrixD = np.array([[1, 0, 0, x[0]],
                                   [0, 1, 0, x[1]],
                                   [0, 0, 1, x[2]],
                                   [0, 0, 0, 1]
    ])
    # use transformPoints3D function to transform each point in a pointset using a tranformation matrix
    point = TransformOnePoint3D(rotationMatrixZ, u)
    point = TransformOnePoint3D(rotationMatrixY, point)
    point = TransformOnePoint3D(rotationMatrixX, point)
    point = TransformOnePoint3D(translationMatrixD, point)
    # project transformed points
    projected, _, _ = transformClass.TransformPoint(point)
    projected_point =  [projected[0]/(spacing*0.194)+size[0]/2, projected[1]/(spacing*0.194)+size[1]/2]
    # return projected points
    return projected_point


def GetProjectedPointsCTTP(points, transformClass, size, spacing):
    '''
    Project 3D CT points to 2D points using PythonVersorRigid3DPerspectiveTransform
    https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection 
    Parameters:
        points (array): nx3 array where n is the number of points to be projected
        transformClass (class): constructed PythonVersorRigid3DPerspectiveTransform class
        size (array): 1x3 array of pixel dimensions of tomosynthesis projection image
        spacing (float): value that decribes the x, y spacing of the CT scan for the purpose of scaling the source points
    Returns:
        array: nx2 array of projected points where n is the number of points
        array: nx3 array of the point of intersections in 3D space between the emitter and the detection plane
        array: nx1 array of the t value in the parametric equation defines in the link
    '''
    projectedPoints = []
    point_of_intersection = []
    t_vals = []
    for point in points:
        # tranform each point and store output values
        projected, num, t = transformClass.TransformPoint(point)
        projectedPoints = np.append(projectedPoints, 
                                    [projected[0]/(spacing*0.194)+size[0]/2, 
                                     projected[1]/(spacing*0.194)+size[1]/2])
        point_of_intersection = np.append(point_of_intersection, num)
        t_vals = np.append(t_vals, t)
    # reshape (np.append flattens array)
    projectedPoints = projectedPoints.reshape((len(points), 2))
    point_of_intersection = point_of_intersection.reshape((len(points), 3))
    # return list of each output value of each point
    return projectedPoints, point_of_intersection[:,:2], t_vals


def MakeVesselMask(y, source, size, destDir):
    '''
    Function to make x number of toy vessel images with different geometries
    using transformed tomosynthesis reconstruction source points, where x is the number of emitter positions
    Parameters:
        y (array): parsed geo.txt file that contains emitter positions
        source (array): 3D tomoRecon source points
        size (array): 1x3 array of pixel dimensions of tomosynthesis projection image
        destDir (string): string that descibes directory of where the mask .dcm files should be written 
    Returns:
        None
    '''
    CenterOfRotation = np.array([0, 0, 0])                     # not used in the transform point (using center of field)
    PlaneCenter = np.array([size[0]/2, size[1]/2, 0])          # center of detector
    PlaneNormal = np.array([0, 0, 1])                          # vector normal
    xDirection = np.array([size[0]/2+1, size[1]/2, 0])         # vector x direction
    yDirection = np.array([size[0]/2, size[1]/2+1, 0])         # vector y direction
    count = 1
    for line, i in zip(y, range(len(y))):
        new_points = []
        EmitterPosition = [(size[0]/2)*0.194-line[0], (size[1]/2)*0.194+line[1], line[2]]
        # update projection class
        ProjP = T.VersorRigid3DPerspectiveTransform(CenterOfRotation, EmitterPosition,
                                                    PlaneCenter, PlaneNormal, 
                                                    xDirection, yDirection) 
        # project point
        # make blank array with Tomo Projection size
        array = np.zeros(shape=[size[1], size[0]], dtype=np.uint8)
        # for each projected point, draw 10x12 point in image
        for point3D in source:
            projected,_ ,_ = ProjP.TransformPoint(np.array([point3D[0], point3D[1], point3D[2]]))
            point = [((projected[0]+(size[0]/2))/0.194),
                     ((projected[1]+(size[1]/2))/0.194)
                    ]
            new_points.append(point)
        for point in new_points:
            if point[0]>=1536:
                None
            elif point[0]<0:
                None
            elif point[1]>=2048:
                None
            elif point[1]<0:
                None
            else:
                ind_1lb = int(point[1])-10
                ind_1ub = int(point[1])+10
                ind_2lb = int(point[0])-12
                ind_2ub = int(point[0])+12
                array[ind_1lb:ind_1ub, ind_2lb:ind_2ub] = 255

        print("line {} of 29 complete".format(count))
        count = count + 1
        img = itk.GetImageFromArray(array)
        new_suffix = f'{i+1:02}.dcm'
        itk.imwrite(img, destDir+"/mask_ct_"+new_suffix)


def ProjectPositionCT(x, line, pnts, imgPath, size):
    '''
    Function that plots CT points over an image given the tranformation matrix and given emitter position
    and image path. Displays plot over iamge
    Parameters:
        x (array): 1x6 array where the inputs define the x translation, y translation
                   z translation, z rotation, y rotation and x rotation respectively
        line (array): 1x3 array that contains an emitter position (most likely one line of the geo.txt file)
        pnts (array): 3D CT source points
        imgPath (string): path to the tomosynthesis projection image to be plotted over
        size (array): 1x3 array of pixel dimensions of tomosynthesis projection image
    Returns:
        None
    '''
    CenterOfRotation = np.array([0, 0, 0])                     # not used in the transform point (using center of field)
    PlaneCenter = np.array([0, 0, 0])                          # center of detector
    PlaneNormal = np.array([0, 0, 1])                          # vector normal
    xDirection = np.array([1, 0, 0])                           # vector x direction
    yDirection = np.array([0, 1, 0])                           # vector y direction
    EmitterPosition = [-line[0], line[1], line[2]]
    ProjP = T.VersorRigid3DPerspectiveTransform(CenterOfRotation, EmitterPosition, 
                                                PlaneCenter, PlaneNormal, 
                                                xDirection, yDirection) 
    transformed = TransformAllPointsAllParameters(x, pnts) 
    projected, _, _ = GetProjectedPointsCTTP(transformed, ProjP, size)
    plt.close()
    plt.rcParams["figure.figsize"] = (12,8.5)
    img = itk.imread(imgPath) # read image
    img = np.squeeze(img) 
    fig, ax = plt.subplots()
    ax.imshow(img)
    for projectedpoint in projected:
        line1 = plt.scatter(projectedpoint[0], projectedpoint[1], 10, c='y')
    plt.show()


def ProjectPositionTR(x, line, pnts, imgPath, size):
    '''
    Function that plots tomoRecon points over an image given the tranformation matrix and given emitter position
    and image path. Displays plot over iamge
    Parameters:
        x (array): 1x6 array where the inputs define the x translation, y translation
                   z translation, z rotation, y rotation and x rotation respectively
        line (array): 1x3 array that contains an emitter position (most likely one line of the geo.txt file)
        pnts (array): 3D tomorecon source points
        imgPath (string): path to the tomosynthesis projection image to be plotted over
        size (array): 1x3 array of pixel dimensions of tomosynthesis projection image
    Returns:
        None
    '''
    CenterOfRotation = np.array([0, 0, 0])                     # not used in the transform point (using center of field)
    PlaneCenter = np.array([size[0]/2, size[1]/2, 0])          # center of detector
    PlaneNormal = np.array([0, 0, 1])                          # vector normal
    xDirection = np.array([size[0]/2+1, size[1]/2, 0])         # vector x direction
    yDirection = np.array([size[0]/2, size[1]/2+1, 0])         # vector y direction
    EmitterPosition = [(size[0]/2)*0.194-line[0], (size[1]/2)*0.194+line[1], 
                       line[2]]
    ProjP = T.VersorRigid3DPerspectiveTransform(CenterOfRotation, EmitterPosition, 
                                                PlaneCenter, PlaneNormal,
                                                xDirection, yDirection) 
    transformed = TransformAllPointsAllParameters(x, pnts) 
    projected, _, _ = GetProjectedPointsCTTP(transformed, ProjP, size)
    plt.close()
    plt.rcParams["figure.figsize"] = (12,8.5)
    img = itk.imread(imgPath) # read image
    img = np.squeeze(img) 
    fig, ax = plt.subplots()
    ax.imshow(img)
    for projectedpoint in projected:
        line1 = plt.scatter(projectedpoint[0], projectedpoint[1], 10, c='y')
    plt.show()


def MakeVesselOverlay(y, source, size, destDir, spacing):
    '''
    Function to write x number of vessel overlay images with different geometries
    using source vessel points
    Parameters:
        y (array): nx3 array that contains an emitter position per line (most likely one line of the geo.txt file)
        source (array): 3D CT source points
        size (array): 1x3 array of pixel dimensions of tomosynthesis projection image
        destDir (string): path to directory where user wants the overlay to be written
        spacing (float): value that decribes the x, y spacing of the CT scan for the purpose of scaling the source points
    Returns:
        None
    '''
    # Set up VersorRigid3DPerspectiveTransform for projecting CT Data
    CenterOfRotation = np.array([0, 0, 0])                     # not used in the transform point (using center of field)
    PlaneCenter = np.array([0, 0, 0])                          # center of detector
    PlaneNormal = np.array([0, 0, 1])                          # vector normal
    xDirection = np.array([1, 0, 0])                           # vector x direction
    yDirection = np.array([0, 1, 0])                           # vector y direction
    count = 1
    for line, i in zip(y, range(len(y))):
        new_points = []
        EmitterPosition = [-line[0], line[1], line[2]]
        # update projection class
        ProjP = T.VersorRigid3DPerspectiveTransform(CenterOfRotation, EmitterPosition, 
                                                    PlaneCenter, PlaneNormal, 
                                                    xDirection, yDirection) 
        # project point
        # make blank array with Tomo Projection size
        array = np.zeros(shape=[size[1], size[0]], dtype=np.uint8)
        for point3D in source:
            projected,_ ,_ = ProjP.TransformPoint(np.array([point3D[0], point3D[1], point3D[2]]))
            projected_point = [projected[0]/(spacing*0.194)+size[0]/2, projected[1]/(spacing*0.194)+size[1]/2]
            new_points.append(projected_point)
        for point in new_points:
            if point[0]>=1536:
                None
            elif point[0]<0:
                None
            elif point[1]>=2048:
                None
            elif point[1]<0:
                None
            else:
                ind_1lb = int(point[1])-3
                ind_1ub = int(point[1])+3
                ind_2lb = int(point[0])-3
                ind_2ub = int(point[0])+3
                array[ind_1lb:ind_1ub, ind_2lb:ind_2ub] = 255

        print("line {} of 29 complete".format(count))
        count = count + 1
        img = itk.GetImageFromArray(array)
        new_suffix = f'{i+1:02}.dcm'
        itk.imwrite(img, destDir+"/vessOverlay_"+new_suffix)


def CT_TomoProjectionRegistration(x, u, y, mask, x_scale, size, spacing):
    '''
    For each point and emitter position, evaluate projected point based on voxel intensity.
    Optimization wants to align points with greatest average voxel vals
    Parameters:
        x (array): 1x6 array where the inputs define the x translation, y translation
                   z translation, z rotation, y rotation and x rotation respectively
        u (array): 3D CT source points
        y (array): nx3 array that contains an emitter position per line (most likely one line of the geo.txt file)
        mask (list): list of strings of the location of the images to register the source points with
        x_scale (array): 1x6 array where the inputs define the x translation, y translation
                         z translation, z rotation, y rotation and x rotation scales respectively
        size (array): 1x3 array of pixel dimensions of tomosynthesis projection image
        spacing (float): value that decribes the x, y spacing of the CT scan for the purpose of scaling the source points
    Returns:
        int: 6000 minus the average voxel value of the projected points average voxel values at each emitter position
    '''
    CenterOfRotation = np.array([0, 0, 0])                     # not used in the transform point (using center of field)
    PlaneCenter = np.array([0, 0, 0])                          # center of detector
    PlaneNormal = np.array([0, 0, 1])                          # vector normal
    xDirection = np.array([1, 0, 0])                           # vector x direction
    yDirection = np.array([0, 1, 0])                           # vector y direction
    totVal = 0
    count = 1
    for point in u:
        voxVal = 0
        for line, path in zip(y, mask):
            EmitterPosition = [-line[0], line[1], line[2]]
            # update projection class
            ProjP = T.VersorRigid3DPerspectiveTransform(CenterOfRotation, EmitterPosition, 
                                                        PlaneCenter, PlaneNormal,
                                                        xDirection, yDirection) 
            # project point
            proj = modelCT(x/x_scale, point, ProjP, size, spacing)
            # get val at projected voxel
            im = itk.imread(path)
            im.SetSpacing((0.194, 0.194))
            im.Update()
            # adjust for projections outside of index of tomo image
            if int(proj[0]) >= size[0]:
                print("out of bounds")
            elif int(proj[1]) >= size[1]:
                print("out of bounds")
            elif int(proj[0]) <= 0:
                print("out of bounds")
            elif int(proj[1]) <= 0:
                print("out of bounds")
            elif im[int(proj[1]), int(proj[0])] == 0:
                None
            else:
                voxVal = voxVal + im[int(proj[1]), int(proj[0])]
        print("Optimized point {} of {}".format(count, len(u)))
        count = count+1
        totVal = totVal + voxVal/len(y)
    returnVal = (6000 - totVal/len(u))
    print (returnVal, x/x_scale)
    return returnVal