import numpy as np
from utils import Isometry3d, to_rotation


class Feature(object):
    # id for next feature
    next_id = 0

    # Takes a vector from the cam0 frame to the cam1 frame.
    R_cam0_cam1 = None
    t_cam0_cam1 = None

    def __init__(self, new_id=0, optimization_config=None):
        # An unique identifier for the feature.
        self.id = new_id

        # Store the observations of the features in the
        # state_id(key)-image_coordinates(value) manner.
        self.observations = dict()   # <StateID, vector4d>

        # 3d postion of the feature in the world frame.
        self.position = np.zeros(3)

        # A indicator to show if the 3d postion of the feature
        # has been initialized or not.
        self.is_initialized = False

        # Optimization configuration for solving the 3d position.
        self.optimization_config = optimization_config

    def cost(self, T_c0_ci, x, z):
        """
        Compute the cost of the camera observations

        Arguments:
            T_c0_ci: A rigid body transformation takes a vector in c0 frame
                to ci frame. (Isometry3d)
            x: The current estimation. (vec3)
            z: The ith measurement of the feature j in ci frame. (vec2)

        Returns:
            e: The cost of this observation. (double)
        """
        # Compute hi1, hi2, and hi3 as Equation (37).
        alpha, beta, rho = x
        h = T_c0_ci.R @ np.array([alpha, beta, 1.0]) + rho * T_c0_ci.t

        # Predict the feature observation in ci frame.
        z_hat = h[:2] / h[2]

        # Compute the residual.
        e = ((z_hat - z)**2).sum()
        return e

    def jacobian(self, T_c0_ci, x, z):
        """
        Compute the Jacobian of the camera observation

        Arguments:
            T_c0_ci: A rigid body transformation takes a vector in c0 frame
                to ci frame. (Isometry3d)
            x: The current estimation. (vec3)
            z: The ith measurement of the feature j in ci frame. (vec2)

        Returns:
            J: The computed Jacobian. (Matrix23)
            r: The computed residual. (vec2)
            w: Weight induced by huber kernel. (double)
        """
        # Compute hi1, hi2, and hi3 as Equation (37).
        alpha, beta, rho = x
        h = T_c0_ci.R @ np.array([alpha, beta, 1.0]) + rho * T_c0_ci.t
        h1, h2, h3 = h

        # Compute the Jacobian.
        W = np.zeros((3, 3))
        W[:, :2] = T_c0_ci.R[:, :2]
        W[:, 2] = T_c0_ci.t

        J = np.zeros((2, 3))
        J[0] = W[0]/h3 - W[2]*h1/(h3*h3)
        J[1] = W[1]/h3 - W[2]*h2/(h3*h3)

        # Compute the residual.
        z_hat = np.array([h1/h3, h2/h3])
        r = z_hat - z

        # Compute the weight based on the residual.
        e = np.linalg.norm(r)
        if e <= self.optimization_config.huber_epsilon:
            w = 1.0
        else:
            w = self.optimization_config.huber_epsilon / (2*e)

        return J, r, w

    def generate_initial_guess(self, T_c1_c2, z1, z2):
        """
        Compute the initial guess of the feature's 3d position using 
        only two views.

        Arguments:
            T_c1_c2: A rigid body transformation taking a vector from c2 frame 
                to c1 frame. (Isometry3d)
            z1: feature observation in c1 frame. (vec2)
            z2: feature observation in c2 frame. (vec2)

        Returns:
            p: Computed feature position in c1 frame. (vec3)
        """
        # Construct a least square problem to solve the depth.
        m = T_c1_c2.R @ np.array([*z1, 1.0])
        a = m[:2] - z2*m[2]                   # vec2
        b = z2*T_c1_c2.t[2] - T_c1_c2.t[:2]   # vec2

        # Solve for the depth.
        depth = a @ b / (a @ a)
        
        p = np.array([*z1, 1.0]) * depth
        return p

    def check_motion(self, cam_states):
        """
        Check the input camera poses to ensure there is enough translation 
        to triangulate the feature

        Arguments:
            cam_states: input camera poses. (dict of <CAMStateID, CAMState>)

        Returns:
            True if the translation between the input camera poses 
                is sufficient. (bool)
        """
        if self.optimization_config.translation_threshold < 0:
            return True

        observation_ids = list(self.observations.keys())
        first_id = observation_ids[0]
        last_id = observation_ids[-1]

        first_cam_pose = Isometry3d(
            to_rotation(cam_states[first_id].orientation).T,
            cam_states[first_id].position)

        last_cam_pose = Isometry3d(
            to_rotation(cam_states[last_id].orientation).T,
            cam_states[last_id].position)

        # Get the direction of the feature when it is first observed.
        # This direction is represented in the world frame.
        feature_direction = np.array([*self.observations[first_id][:2], 1.0])
        feature_direction = feature_direction / np.linalg.norm(feature_direction)
        feature_direction = first_cam_pose.R @ feature_direction

        # Compute the translation between the first frame and the last frame. 
        # We assume the first frame and the last frame will provide the 
        # largest motion to speed up the checking process.
        translation = last_cam_pose.t - first_cam_pose.t
        parallel = translation @ feature_direction
        orthogonal_translation = translation - parallel * feature_direction

        return (np.linalg.norm(orthogonal_translation) > 
            self.optimization_config.translation_threshold)

    def initialize_position(self, cam_states):
        """
        Intialize the feature position based on all current available 
        measurements.

        The computed 3d position is used to set the position member variable. 
        Note the resulted position is in world frame.

        Arguments:
            cam_states: A dict containing the camera poses with its ID as the 
                associated key value. (dict of <CAMStateID, CAMState>)

        Returns:
            True if the estimated 3d position of the feature is valid. (bool)
        """
        cam_poses = []     # [Isometry3d]
        measurements = []  # [vec2]

        T_cam1_cam0 = Isometry3d(
            Feature.R_cam0_cam1, Feature.t_cam0_cam1).inverse()

        for cam_id, m in self.observations.items():
            try:
                cam_state = cam_states[cam_id]
            except KeyError:
                continue
            
            # Add measurements.
            measurements.append(m[:2])
            measurements.append(m[2:])

            # This camera pose will take a vector from this camera frame
            # to the world frame.
            cam0_pose = Isometry3d(
                to_rotation(cam_state.orientation).T, cam_state.position)
            cam1_pose = cam0_pose * T_cam1_cam0

            cam_poses.append(cam0_pose)
            cam_poses.append(cam1_pose)

        # All camera poses should be modified such that it takes a vector 
        # from the first camera frame in the buffer to this camera frame.
        T_c0_w = cam_poses[0]
        cam_poses_tmp = []
        for pose in cam_poses:
            cam_poses_tmp.append(pose.inverse() * T_c0_w)
        cam_poses = cam_poses_tmp

        # Generate initial guess
        initial_position = self.generate_initial_guess(
            cam_poses[-2], measurements[0], measurements[-2])
        solution = np.array([*initial_position[:2], 1.0]) / initial_position[2]

        # Apply Levenberg-Marquart method to solve for the 3d position.
        lambd = self.optimization_config.initial_damping
        inner_loop_count = 0
        outer_loop_count = 0
        is_cost_reduced = False
        delta_norm = float('inf')

        # Compute the initial cost.
        total_cost = 0.0
        # for i, cam_pose in enumerate(cam_poses):
        for cam_pose, measurement in zip(cam_poses, measurements):
            total_cost += self.cost(cam_pose, solution, measurement)

        # Outer loop.
        while (outer_loop_count < 
            self.optimization_config.outer_loop_max_iteration
            and delta_norm > 
            self.optimization_config.estimation_precision):

            A = np.zeros((3, 3))
            b = np.zeros(3)
            for cam_pose, measurement in zip(cam_poses, measurements):
                J, r, w = self.jacobian(cam_pose, solution, measurement)
                if w == 1.0:
                    A += J.T @ J
                    b += J.T @ r
                else:
                    A += w * w * J.T @ J
                    b += w * w * J.T @ r

            # Inner loop.
            # Solve for the delta that can reduce the total cost.
            while (inner_loop_count < 
                self.optimization_config.inner_loop_max_iteration
                and not is_cost_reduced):

                delta = np.linalg.solve(A + lambd * np.identity(3), b)   # vec3
                new_solution = solution - delta
                delta_norm = np.linalg.norm(delta)

                new_cost = 0.0
                for cam_pose, measurement in zip(cam_poses, measurements):
                    new_cost += self.cost(
                        cam_pose, new_solution, measurement)

                if new_cost < total_cost:
                    is_cost_reduced = True
                    solution = new_solution
                    total_cost = new_cost
                    lambd = max(lambd/10., 1e-10)
                else:
                    is_cost_reduced = False
                    lambd = min(lambd*10., 1e12)
                
                inner_loop_count += 1
            inner_loop_count = 0
            outer_loop_count += 1

        # Covert the feature position from inverse depth
        # representation to its 3d coordinate.
        final_position = np.array([*solution[:2], 1.0]) / solution[2]

        # Check if the solution is valid. Make sure the feature
        # is in front of every camera frame observing it.
        is_valid_solution = True
        for pose in cam_poses:
            position = pose.R @ final_position + pose.t
            if position[2] <= 0:
                is_valid_solution = False
                break

        # Convert the feature position to the world frame.
        self.position = T_c0_w.R @ final_position + T_c0_w.t

        self.is_initialized = is_valid_solution
        return is_valid_solution

    def initialize_position_mono(self, cam_states):
        """
        Intialize the feature position based on all current available
        measurements.

        The computed 3d position is used to set the position member variable.
        Note the resulted position is in world frame.

        Arguments:
            cam_states: A dict containing the camera poses with its ID as the
                associated key value. (dict of <CAMStateID, CAMState>)

        Returns:
            True if the estimated 3d position of the feature is valid. (bool)
        """
        cam_poses = []  # [Isometry3d]
        measurements = []  # [vec2]

        for cam_id, m in self.observations.items():
            try:
                cam_state = cam_states[cam_id]
            except KeyError:
                continue

            # Add measurements.
            measurements.append(m[:2])

            # This camera pose will take a vector from this camera frame
            # to the world frame.
            cam0_pose = Isometry3d(to_rotation(cam_state.orientation).T, cam_state.position)
            cam_poses.append(cam0_pose)

        # All camera poses should be modified such that it takes a vector
        # from the first camera frame in the buffer to this camera frame.
        T_c0_w = cam_poses[0]
        cam_poses_tmp = []
        for pose in cam_poses:
            cam_poses_tmp.append(pose.inverse() * T_c0_w)
        cam_poses = cam_poses_tmp

        # Generate initial guess
        initial_position = self.generate_initial_guess(cam_poses[-2], measurements[0], measurements[-2])
        solution = np.array([*initial_position[:2], 1.0]) / initial_position[2]

        # Apply Levenberg-Marquart method to solve for the 3d position.
        lambd = self.optimization_config.initial_damping
        inner_loop_count = 0
        outer_loop_count = 0
        is_cost_reduced = False
        delta_norm = float('inf')

        # Compute the initial cost.
        total_cost = 0.0
        # for i, cam_pose in enumerate(cam_poses):
        for cam_pose, measurement in zip(cam_poses, measurements):
            total_cost += self.cost(cam_pose, solution, measurement)

        # Outer loop.
        while (outer_loop_count <
               self.optimization_config.outer_loop_max_iteration
               and delta_norm >
               self.optimization_config.estimation_precision):

            A = np.zeros((3, 3))
            b = np.zeros(3)
            for cam_pose, measurement in zip(cam_poses, measurements):
                J, r, w = self.jacobian(cam_pose, solution, measurement)
                if w == 1.0:
                    A += J.T @ J
                    b += J.T @ r
                else:
                    A += w * w * J.T @ J
                    b += w * w * J.T @ r

            # Inner loop.
            # Solve for the delta that can reduce the total cost.
            while (inner_loop_count <
                   self.optimization_config.inner_loop_max_iteration
                   and not is_cost_reduced):

                delta = np.linalg.solve(A + lambd * np.identity(3), b)  # vec3
                new_solution = solution - delta
                delta_norm = np.linalg.norm(delta)

                new_cost = 0.0
                for cam_pose, measurement in zip(cam_poses, measurements):
                    new_cost += self.cost(
                        cam_pose, new_solution, measurement)

                if new_cost < total_cost:
                    is_cost_reduced = True
                    solution = new_solution
                    total_cost = new_cost
                    lambd = max(lambd / 10., 1e-10)
                else:
                    is_cost_reduced = False
                    lambd = min(lambd * 10., 1e12)

                inner_loop_count += 1
            inner_loop_count = 0
            outer_loop_count += 1

        # Covert the feature position from inverse depth
        # representation to its 3d coordinate.
        final_position = np.array([*solution[:2], 1.0]) / solution[2]

        # Check if the solution is valid. Make sure the feature
        # is in front of every camera frame observing it.
        is_valid_solution = True
        for pose in cam_poses:
            position = pose.R @ final_position + pose.t
            if position[2] <= 0:
                is_valid_solution = False
                break

        # Convert the feature position to the world frame.
        self.position = T_c0_w.R @ final_position + T_c0_w.t

        self.is_initialized = is_valid_solution
        return is_valid_solution