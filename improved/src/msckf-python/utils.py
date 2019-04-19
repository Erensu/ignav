import numpy as np


# quaternion representation: [x, y, z, w]
# JPL convention
def skew(vec):
    """
    Create a skew-symmetric matrix from a 3-element vector.
    """
    x, y, z = vec
    return np.array([
        [0, -z, y],
        [z, 0, -x],
        [-y, x, 0]])


def to_rotation(q):
    """
    Convert a quaternion to the corresponding rotation matrix.
    Pay attention to the convention used. The function follows the
    conversion in "Indirect Kalman Filter for 3D Attitude Estimation:
    A Tutorial for Quaternion Algebra", Equation (78).
    The input quaternion should be in the form [q1, q2, q3, q4(scalar)]
    """
    q = q / np.linalg.norm(q)
    vec = q[:3]
    w = q[3]

    R = (2*w*w-1)*np.identity(3) - 2*w*skew(vec) + 2*vec[:, None]*vec
    return R


def to_quaternion(R):
    """
    Convert a rotation matrix to a quaternion.
    Pay attention to the convention used. The function follows the
    conversion in "Indirect Kalman Filter for 3D Attitude Estimation:
    A Tutorial for Quaternion Algebra", Equation (78).
    The input quaternion should be in the form [q1, q2, q3, q4(scalar)]
    """
    if R[2, 2] < 0:
        if R[0, 0] > R[1, 1]:
            t = 1 + R[0,0] - R[1,1] - R[2,2]
            q = [t, R[0, 1]+R[1, 0], R[2, 0]+R[0, 2], R[1, 2]-R[2, 1]]
        else:
            t = 1 - R[0,0] + R[1,1] - R[2,2]
            q = [R[0, 1]+R[1, 0], t, R[2, 1]+R[1, 2], R[2, 0]-R[0, 2]]
    else:
        if R[0, 0] < -R[1, 1]:
            t = 1 - R[0,0] - R[1,1] + R[2,2]
            q = [R[0, 2]+R[2, 0], R[2, 1]+R[1, 2], t, R[0, 1]-R[1, 0]]
        else:
            t = 1 + R[0,0] + R[1,1] + R[2,2]
            q = [R[1, 2]-R[2, 1], R[2, 0]-R[0, 2], R[0, 1]-R[1, 0], t]

    q = np.array(q) # * 0.5 / np.sqrt(t)
    return q / np.linalg.norm(q)


def quaternion_normalize(q):
    """
    Normalize the given quaternion to unit quaternion.
    """
    return q / np.linalg.norm(q)


def quaternion_conjugate(q):
    """
    Conjugate of a quaternion.
    """
    return np.array([*-q[:3], q[3]])


def quaternion_multiplication(q1, q2):
    """
    Perform q1 * q2
    """
    q1 = q1 / np.linalg.norm(q1)
    q2 = q2 / np.linalg.norm(q2)

    L = np.array([
        [ q1[3],  q1[2], -q1[1], q1[0]],
        [-q1[2],  q1[3],  q1[0], q1[1]],
        [ q1[1], -q1[0],  q1[3], q1[2]],
        [-q1[0], -q1[1], -q1[2], q1[3]]
    ])

    q = L @ q2
    return q / np.linalg.norm(q)


def small_angle_quaternion(dtheta):
    """
    Convert the vector part of a quaternion to a full quaternion.
    This function is useful to convert delta quaternion which is  
    usually a 3x1 vector to a full quaternion.
    For more details, check Equation (238) and (239) in "Indirect Kalman 
    Filter for 3D Attitude Estimation: A Tutorial for quaternion Algebra".
    """
    dq = dtheta / 2.
    dq_square_norm = dq @ dq

    if dq_square_norm <= 1:
        q = np.array([*dq, np.sqrt(1-dq_square_norm)])
    else:
        q = np.array([*dq, 1.])
        q /= np.sqrt(1+dq_square_norm)
    return q


def from_two_vectors(v0, v1):
    """
    Rotation quaternion from v0 to v1.
    """
    v0 = v0 / np.linalg.norm(v0)
    v1 = v1 / np.linalg.norm(v1)
    d = v0 @ v1

    # if dot == -1, vectors are nearly opposite
    if d < -0.999999:
        axis = np.cross([1,0,0], v0)
        if np.linalg.norm(axis) < 0.000001:
            axis = np.cross([0,1,0], v0)
        q = np.array([*axis, 0.])
    elif d > 0.999999:
        q = np.array([0., 0., 0., 1.])
    else:
        s = np.sqrt((1+d)*2)
        axis = np.cross(v0, v1)
        vec = axis / s
        w = 0.5 * s
        q = np.array([*vec, w])
        
    q = q / np.linalg.norm(q)
    return quaternion_conjugate(q)   # hamilton -> JPL


class Isometry3d(object):
    """
    3d rigid transform.
    """
    def __init__(self, R, t):
        self.R = R
        self.t = t

    def matrix(self):
        m = np.identity(4)
        m[:3, :3] = self.R
        m[:3, 3] = self.t
        return m

    def inverse(self):
        return Isometry3d(self.R.T, -self.R.T @ self.t)

    def __mul__(self, T1):
        R = self.R @ T1.R
        t = self.R @ T1.t + self.t
        return Isometry3d(R, t)