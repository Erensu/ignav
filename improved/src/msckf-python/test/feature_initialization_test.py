import unittest
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from config import OptimizationConfigEuRoC
from utils import to_quaternion, to_rotation, Isometry3d
from feature import Feature
from msckf import CAMState

class TestFeature(unittest.TestCase):
    def test_feature_initialization(self):
        """
        Test feature initialization.
        """
        optimization_config = OptimizationConfigEuRoC()
        Feature.R_cam0_cam1 = np.identity(3)
        Feature.t_cam0_cam1 = np.zeros(3)

        # feature = np.array([0.5, 0., 0.])
        feature = np.random.random(3) * 0.5

        # Add 6 camera poses, all of which are able to see the
        # feature at the origin. For simplicity, the six camera
        # view are located at the six intersections between a
        # unit sphere and the coordinate system. And the z axes
        # of the camera frames are facing the origin.
        cam_poses = [
            Isometry3d(np.array([
                [0., 0., -1.],
                [1., 0., 0.],
                [0., -1., 0.]]
            ), np.array([1., 0., 0.])),
            Isometry3d(np.array([
                [-1., 0., 0.],
                [0., 0., -1.],
                [0., -1., 0.]]
            ), np.array([0., 1., 0.])),
            Isometry3d(np.array([
                [0., 0., 1.],
                [-1., 0., 0.],
                [0., -1., 0.]]
            ), np.array([-1., 0., 0.])),
            Isometry3d(np.array([
                [1., 0., 0.],
                [0., 0., 1.],
                [0., -1., 0.]]
            ), np.array([0., -1., 0.])),
            Isometry3d(np.array([
                [0., -1., 0.],
                [-1., 0., 0.],
                [0., 0., -1.]]
            ), np.array([0., 0., 1.])),
            Isometry3d(np.array([
                [1., 0., 0.],
                [0., 1., 0.],
                [0., 0., 1.]]
            ), np.array([0., 0., -1.])),
        ]

        # Set the camera states
        cam_states = dict()
        for i in range(6):
            cam_state = CAMState(i)
            cam_state.timestamp = i
            cam_state.orientation = to_quaternion(cam_poses[i].R.T)
            cam_state.position = cam_poses[i].t
            cam_states[i] = cam_state

        # Compute measurements.
        measurements = []
        for i in range(6):
            cam_pose_inv = cam_poses[i].inverse()
            p = cam_pose_inv.R @ feature + cam_pose_inv.t
            u, v = p[:2] / p[2] + np.random.randn(2) * 0.01
            measurements.append(np.array([u, v, u, v]))

        # Initialize a feature object.
        feature_object = Feature(0, optimization_config=optimization_config)
        for i in range(6):
            feature_object.observations[i] = measurements[i]

        # Compute the 3d position of the feature.
        status = feature_object.initialize_position(cam_states)

        # Check the difference between the computed 3d
        # feature position and the groud truth.
        print('status:', status)
        print('ground truth position:\n', feature)
        print('estimated position:\n', feature_object.position)
        e = np.linalg.norm(feature - feature_object.position)
        print('error norm:', e)
        self.assertTrue(e < 0.05)

if __name__ == '__main__':
    unittest.main()