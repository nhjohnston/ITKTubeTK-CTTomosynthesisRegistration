import itk
import math 
import numpy as np

class VersorRigid3DPerspectiveTransform:
    def __init__(self, CenterOfRotation, EmitterPosition, PlaneCenter, PlaneNormal, xDirection, yDirection):
        self.CenterOfRotation = CenterOfRotation
        self.EmitterPosition = EmitterPosition
        self.PlaneCenter = PlaneCenter
        self.PlaneNormal = PlaneNormal
        self.xDirection = xDirection
        self.yDirection = yDirection

    def TransformPoint(self, point):
        # https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection

        # t = ((p01 x p02) * (Ia - p0)) / (-Iab * (p01 x p02))
        # u = ((p02 x -Iab) * (Ia - p0)) / (-Iab * (p01 x p02))
        # v = ((-Iab x p01) * (Ia - p0)) / (-Iab * (p01 x p02))

        # where 

        # * = dot product
        # x = cross product

        # plane_origin = p0
        # plane_x_axis direction = p01 = p1 - p0 
        # plane_y_axis direction = p02 = p2 - p0 
        # emitter_center = Ia
        # point = Ib
        # line = Iab = Ib - Ia

        # u = plane x-axis intersect of point
        # v = plane y-axis intersect of point
        # t = ratio of intersect along line Iab

        # Intersection exists if 
        # (p01 x p02) * Iab != 0

        # Intersection is a line if
        # (p01 x p02) * Iab == 1 
    
        Ia = point
        Ib = self.EmitterPosition
        p0 = self.PlaneCenter
        p01 = self.xDirection - p0
        p02 = self.yDirection - p0
        Iab = -(Ib - Ia)  ## actually storing -Iab, for convenience
        p0Ia = Ia - p0
        transformed_point = np.array([0, 0])
        denom = np.dot(Iab, np.cross(p01, p02))
        transformed_point[0] = np.dot(np.cross(p02, Iab), p0Ia) /  denom
        transformed_point[1] = np.dot(np.cross(Iab, p01), p0Ia) / denom

        m_LastTransformedPointDistance = np.dot(np.cross(p01, p02), p0Ia) / denom
        pointOfIntersection = Ia + Iab*m_LastTransformedPointDistance

        return transformed_point, pointOfIntersection, m_LastTransformedPointDistance
    
    def ComputeJacobianWithRespectToParameters(self, angles, point):
        # compute derivatives with respect to rotation
        x = angles[0]
        y = angles[1]
        z = angles[2]

        cz = math.cos(z*0.5)
        sz = math.sin(z*0.5)
        cy = math.cos(y*0.5)
        sy = math.sin(z*0.5)
        cx = math.cos(x*0.5)
        sx = math.sin(z*0.5)

        vx = sx * cy * cz + cx * sy * sz
        vy = cx * sy * cz + sx * cy * sz
        vz = cx * cy * sz + sx * sy * cz
        vw = cx * cy * cz + sx * sy * sz

        jacobian = np.empty([3, 6])

        px = point[0] - self.CenterOfRotation[0]
        py = point[1] - self.CenterOfRotation[1]
        pz = point[2] - self.CenterOfRotation[2]

        vxx = vx * vx
        vyy = vy * vy
        vzz = vz * vz
        vww = vw * vw

        vxy = vx * vy
        vxz = vx * vz
        vxw = vx * vw

        vyz = vy * vz
        vyw = vy * vw

        vzw = vz * vw

        # compute Jacobian with respect to quaternion parameters
        jacobian[0][0] = 2.0 * ((vyw + vxz) * py + (vzw - vxy) * pz) / vw
        jacobian[1][0] = 2.0 * ((vyw - vxz) * px - 2 * vxw * py + (vxx - vww) * pz) / vw
        jacobian[2][0] = 2.0 * ((vzw + vxy) * px + (vww - vxx) * py - 2 * vxw * pz) / vw

        jacobian[0][1] = 2.0 * (-2 * vyw * px + (vxw + vyz) * py + (vww - vyy) * pz) / vw
        jacobian[1][1] = 2.0 * ((vxw - vyz) * px + (vzw + vxy) * pz) / vw
        jacobian[2][1] = 2.0 * ((vyy - vww) * px + (vzw - vxy) * py - 2 * vyw * pz) / vw

        jacobian[0][2] = 2.0 * (-2 * vzw * px + (vzz - vww) * py + (vxw - vyz) * pz) / vw
        jacobian[1][2] = 2.0 * ((vww - vzz) * px - 2 * vzw * py + (vyw + vxz) * pz) / vw
        jacobian[2][2] = 2.0 * ((vxw + vyz) * px + (vyw - vxz) * py) / vw

        jacobian[0][3] = 1.0
        jacobian[1][4] = 1.0
        jacobian[2][5] = 1.0
        return jacobian
