import unittest
import numpy as np
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from utils import *


class TestUtils(unittest.TestCase):
    def test_skew(self):
        """
        Test skew-symmetric matrix.
        """
        w = np.array([1.0, 2.0, 3.0])
        w_hat = skew(w)
        zero_vec = w_hat @ w

        self.assertEqual(np.linalg.matrix_rank(w_hat), 2)
        self.assertTrue(np.all(zero_vec == 0))

    def test_normalize(self):
        """
        Test quaternion normalize.
        """
        q = np.array([1.0, 2.0, 3.0, 4.0])
        q = quaternion_normalize(q)

        self.assertAlmostEqual(np.linalg.norm(q), 1.0)

    def test_to_rotation(self):
        """
        Test converting quaternion to rotation matrix.
        """
        q = np.array([-1, 1, 3, 2])
        q = q / np.linalg.norm(q)
        R_gt = np.array([
            [-1/3.0, -14/15., -2/15.0],
            [2/3.0, -1/3.0, 2/3.0],
            [-2/3.0, 2/15.0, 11/15.0]]).T
        R = to_rotation(q)

        zero_matrix = R - R_gt
        self.assertAlmostEqual(np.linalg.norm(zero_matrix), 0.0)

        for _ in range(20):
            q = np.random.randn(4)
            q /= np.linalg.norm(q)
            q_inv = quaternion_conjugate(q)

            R = to_rotation(q)
            R_inv = to_rotation(q_inv)

            zero_matrix = R @ R_inv - np.identity(3)
            self.assertAlmostEqual(np.linalg.norm(zero_matrix), 0.0)

            # orthogonal matrix
            zero_matrix = R @ R.T - np.identity(3)
            self.assertAlmostEqual(np.linalg.norm(zero_matrix), 0.0)

    def test_to_quaternion(self):
        """
        Test converting rotation matrix quaternion.
        """
        R = np.identity(3)
        q = to_quaternion(R)
        zero_vec = q - np.array([0., 0., 0., 1.])
        self.assertAlmostEqual(np.linalg.norm(zero_vec), 0.0)

        for _ in range(20):
            q = np.random.randn(4)
            q /= np.linalg.norm(q)

            R = to_rotation(q)
            R2 = to_rotation(to_quaternion(R))
            zero_matrix = R - R2
            self.assertAlmostEqual(np.linalg.norm(zero_matrix), 0.0)

    def test_multiplication(self):
        """
        Test quaternion multiplication.
        """
        for _ in range(20):
            q1 = np.random.randn(4)
            q2 = np.random.randn(4)
            q1 /= np.linalg.norm(q1)
            q2 /= np.linalg.norm(q2)
            q_prod = quaternion_multiplication(q1, q2)

            R1 = to_rotation(q1)
            R2 = to_rotation(q2)
            R_prod = R1 @ R2
            R_prod_cp = to_rotation(q_prod)

            zero_matrix = R_prod - R_prod_cp
            self.assertAlmostEqual(np.linalg.norm(zero_matrix), 0.0)

    def test_from_two_vectors(self):
        """
        Test rotation quaternion from two vectors.
        """
        for _ in range(20):
            v0 = np.random.randn(3)
            v1 = np.random.randn(3)
            v0 /= np.linalg.norm(v0)
            v1 /= np.linalg.norm(v1)

            q = from_two_vectors(v0, v1)
            R = to_rotation(q)

            zero_vec = R @ v0 - v1
            self.assertAlmostEqual(np.linalg.norm(zero_vec), 0.0)

            q_inv = from_two_vectors(v1, v0)
            R_inv = to_rotation(q_inv)
            zero_matrix = R @ R_inv - np.identity(3)
            self.assertAlmostEqual(np.linalg.norm(zero_matrix), 0.0)

    def test_small_angle_quaternion(self):
        pass

    def test_isometry3d(self):
        """
        Test Isometry3d.
        """
        q = np.random.randn(4)
        R = to_rotation(q)
        t = np.random.randn(3)

        T = Isometry3d(R, t)
        T_inv = T.inverse()
        T_identity = T * T_inv
        zero_matrix = T_identity.R - np.identity(3)
        zero_vec = T_identity.t

        self.assertAlmostEqual(np.linalg.norm(zero_matrix), 0.0)
        self.assertAlmostEqual(np.linalg.norm(zero_vec), 0.0)

        T1 = Isometry3d(
            to_rotation(np.random.randn(4)), 
            np.random.randn(3))
        T2 = Isometry3d(
            to_rotation(np.random.randn(4)), 
            np.random.randn(3))

        T3 = T1 * T2
        zero_matrix = T3.matrix() - T1.matrix() @ T2.matrix()
        self.assertAlmostEqual(np.linalg.norm(zero_matrix), 0.0)


if __name__ == '__main__':
    unittest.main()