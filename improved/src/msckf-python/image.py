import numpy as np
import cv2
import time
from itertools import chain, compress
from collections import defaultdict, namedtuple


class FeatureMetaData(object):
    """
    Contain necessary information of a feature for easy access.
    """
    def __init__(self):
        self.id = None           # int
        self.response = None     # float
        self.lifetime = None     # int
        self.cam0_point = None   # vec2
        self.cam1_point = None   # vec2


class FeatureMeasurement(object):
    """
    Stereo measurement of a feature.
    """
    def __init__(self):
        self.id = None
        self.u0 = None
        self.v0 = None
        self.u1 = None
        self.v1 = None


class ImageProcessor(object):
    """
    Detect and track features in image sequences.
    """
    def __init__(self, config):
        self.config = config

        # Indicate if this is the first image message.
        self.is_first_img = True

        # ID for the next new feature.
        self.next_feature_id = 0

        # Feature detector
        self.detector = cv2.FastFeatureDetector_create(self.config.fast_threshold)

        # IMU message buffer.
        self.imu_msg_buffer = []

        # Previous and current images
        self.cam0_prev_img_msg = None
        self.cam0_curr_img_msg = None
        self.cam1_curr_img_msg = None

        # Pyramids for previous and current image
        self.prev_cam0_pyramid = None
        self.curr_cam0_pyramid = None
        self.curr_cam1_pyramid = None

        # Features in the previous and current image.
        # list of lists of FeatureMetaData
        self.prev_features = [[] for _ in range(self.config.grid_num)]  # Don't use [[]] * N
        self.curr_features = [[] for _ in range(self.config.grid_num)]

        # Number of features after each outlier removal step.
        # keys: before_tracking, after_tracking, after_matching, after_ransac
        self.num_features = defaultdict(int)

        # load config
        # Camera calibration parameters
        self.cam0_resolution = config.cam0_resolution   # vec2
        self.cam0_intrinsics = config.cam0_intrinsics   # vec4
        self.cam0_distortion_model = config.cam0_distortion_model     # string
        self.cam0_distortion_coeffs = config.cam0_distortion_coeffs   # vec4

        self.cam1_resolution = config.cam1_resolution   # vec2
        self.cam1_intrinsics = config.cam1_intrinsics   # vec4
        self.cam1_distortion_model = config.cam1_distortion_model     # string
        self.cam1_distortion_coeffs = config.cam1_distortion_coeffs   # vec4

        # Take a vector from cam0 frame to the IMU frame.
        self.T_cam0_imu = np.linalg.inv(config.T_imu_cam0)
        self.R_cam0_imu = self.T_cam0_imu[:3, :3]
        self.t_cam0_imu = self.T_cam0_imu[:3, 3]

        # Take a vector from cam1 frame to the IMU frame.
        self.T_cam1_imu = np.linalg.inv(config.T_imu_cam1)
        self.R_cam1_imu = self.T_cam1_imu[:3, :3]
        self.t_cam1_imu = self.T_cam1_imu[:3, 3]

    def stereo_callback(self, stereo_msg):
        """
        Callback function for the stereo images.
        """
        self.cam0_curr_img_msg = stereo_msg.cam0_msg
        self.cam1_curr_img_msg = stereo_msg.cam1_msg

        # Build the image pyramids once since they're used at multiple places.
        self.create_image_pyramids()

        # Detect features in the first frame.
        if self.is_first_img:
            self.initialize_first_frame()
            self.is_first_img = False

            # Draw results.
            # self.draw_features_stereo()
        else:
            # Track the feature in the previous image.
            self.track_features()

            # Add new features into the current image.
            self.add_new_features()
            self.prune_features()

            # Draw results.
            # self.draw_features_stereo()

        # print('===image process elapsed:', time.time() - start, f'({stereo_msg.timestamp})')

        try:
            return self.publish()
        finally:
            self.cam0_prev_img_msg = self.cam0_curr_img_msg
            self.prev_features = self.curr_features
            self.prev_cam0_pyramid = self.curr_cam0_pyramid

            # Initialize the current features to empty vectors.
            self.curr_features = [[] for _ in range(self.config.grid_num)]

    def imu_callback(self, msg):
        """
        Callback function for the imu message.
        """
        self.imu_msg_buffer.append(msg)

    def create_image_pyramids(self):
        """
        Create image pyramids used for KLT tracking.
        (Seems doesn't work in python)
        """
        curr_cam0_img = self.cam0_curr_img_msg.image
        self.curr_cam0_pyramid = curr_cam0_img

        curr_cam1_img = self.cam1_curr_img_msg.image
        self.curr_cam1_pyramid = curr_cam1_img

    def initialize_first_frame(self):
        """
        Initialize the image processing sequence, which is basically detect 
        new features on the first set of stereo images.
        """
        img = self.cam0_curr_img_msg.image
        grid_height, grid_width = self.get_grid_size(img)

        # Detect new features on the frist image.
        new_features = self.detector.detect(img)

        # Find the stereo matched points for the newly detected features.
        cam0_points = [kp.pt for kp in new_features]
        cam1_points, inlier_markers = self.stereo_match(cam0_points)

        cam0_inliers, cam1_inliers = [], []
        response_inliers = []
        for i, inlier in enumerate(inlier_markers):
            if not inlier:
                continue
            cam0_inliers.append(cam0_points[i])
            cam1_inliers.append(cam1_points[i])
            response_inliers.append(new_features[i].response)

        # Group the features into grids
        grid_new_features = [[] for _ in range(self.config.grid_num)]

        for i in range(len(cam0_inliers)):
            cam0_point = cam0_inliers[i]
            cam1_point = cam1_inliers[i]
            response = response_inliers[i]

            row = int(cam0_point[1] / grid_height)
            col = int(cam0_point[0] / grid_width)
            code = row*self.config.grid_col + col

            new_feature = FeatureMetaData()
            new_feature.response = response
            new_feature.cam0_point = cam0_point
            new_feature.cam1_point = cam1_point
            grid_new_features[code].append(new_feature)

        # Sort the new features in each grid based on its response.
        # And collect new features within each grid with high response.
        for i, new_features in enumerate(grid_new_features):
            for feature in sorted(new_features, key=lambda x:x.response, reverse=True)[:self.config.grid_min_feature_num]:

                self.curr_features[i].append(feature)
                self.curr_features[i][-1].id = self.next_feature_id
                self.curr_features[i][-1].lifetime = 1
                self.next_feature_id += 1

    def track_features(self):
        """
        Tracker features on the newly received stereo images.
        """
        img = self.cam0_curr_img_msg.image
        grid_height, grid_width = self.get_grid_size(img)

        # Compute a rough relative rotation which takes a vector 
        # from the previous frame to the current frame.
        cam0_R_p_c, cam1_R_p_c = self.integrate_imu_data()

        # Organize the features in the previous image.
        prev_ids = []
        prev_lifetime = []
        prev_cam0_points = []
        prev_cam1_points = []

        for feature in chain.from_iterable(self.prev_features):
            prev_ids.append(feature.id)
            prev_lifetime.append(feature.lifetime)
            prev_cam0_points.append(feature.cam0_point)
            prev_cam1_points.append(feature.cam1_point)
        prev_cam0_points = np.array(prev_cam0_points, dtype=np.float32)

        # Number of the features before tracking.
        self.num_features['before_tracking'] = len(prev_cam0_points)

        # Abort tracking if there is no features in the previous frame.
        if len(prev_cam0_points) == 0:
            return

        # Track features using LK optical flow method.
        curr_cam0_points = self.predict_feature_tracking(
            prev_cam0_points,
            cam0_R_p_c,
            self.cam0_intrinsics)

        curr_cam0_points, track_inliers, _ = cv2.calcOpticalFlowPyrLK(
            self.prev_cam0_pyramid,
            self.curr_cam0_pyramid,
            prev_cam0_points.astype(np.float32), 
            curr_cam0_points.astype(np.float32), 
            **self.config.lk_params)
            
        # Mark those tracked points out of the image region as untracked.
        for i, point in enumerate(curr_cam0_points):
            if not track_inliers[i]:
                continue
            if (point[0] < 0 or point[0] > img.shape[1]-1 or 
                point[1] < 0 or point[1] > img.shape[0]-1):
                track_inliers[i] = 0

        # Collect the tracked points.
        prev_tracked_ids = select(prev_ids, track_inliers)
        prev_tracked_lifetime = select(prev_lifetime, track_inliers)
        prev_tracked_cam0_points = select(prev_cam0_points, track_inliers)
        prev_tracked_cam1_points = select(prev_cam1_points, track_inliers)
        curr_tracked_cam0_points = select(curr_cam0_points, track_inliers)

        # Number of features left after tracking.
        self.num_features['after_tracking'] = len(curr_tracked_cam0_points)

        # Outlier removal involves three steps, which forms a close
        # loop between the previous and current frames of cam0 (left)
        # and cam1 (right). Assuming the stereo matching between the
        # previous cam0 and cam1 images are correct, the three steps are:
        #
        # prev frames cam0 ----------> cam1
        #              |                |
        #              |ransac          |ransac
        #              |   stereo match |
        # curr frames cam0 ----------> cam1
        #
        # 1) Stereo matching between current images of cam0 and cam1.
        # 2) RANSAC between previous and current images of cam0.
        # 3) RANSAC between previous and current images of cam1.
        #
        # For Step 3, tracking between the images is no longer needed.
        # The stereo matching results are directly used in the RANSAC.

        # Step 1: stereo matching.
        curr_cam1_points, match_inliers = self.stereo_match(curr_tracked_cam0_points)

        prev_matched_ids = select(prev_tracked_ids, match_inliers)
        prev_matched_lifetime = select(prev_tracked_lifetime, match_inliers)
        prev_matched_cam0_points = select(prev_tracked_cam0_points, match_inliers)
        prev_matched_cam1_points = select(prev_tracked_cam1_points, match_inliers)
        curr_matched_cam0_points = select(curr_tracked_cam0_points, match_inliers)
        curr_matched_cam1_points = select(curr_cam1_points, match_inliers)

        # Number of features left after stereo matching.
        self.num_features['after_matching'] = len(curr_matched_cam0_points)

        cam0_ransac_inliers = [1] * len(prev_matched_cam0_points)
        cam1_ransac_inliers = [1] * len(prev_matched_cam1_points)

        # Number of features after ransac.
        after_ransac = 0
        for i in range(len(cam0_ransac_inliers)):
            if not (cam0_ransac_inliers[i] and cam1_ransac_inliers[i]):
                continue 
            row = int(curr_matched_cam0_points[i][1] / grid_height)
            col = int(curr_matched_cam0_points[i][0] / grid_width)
            code = row * self.config.grid_col + col

            grid_new_feature = FeatureMetaData()
            grid_new_feature.id = prev_matched_ids[i]
            grid_new_feature.lifetime = prev_matched_lifetime[i] + 1
            grid_new_feature.cam0_point = curr_matched_cam0_points[i]
            grid_new_feature.cam1_point = curr_matched_cam1_points[i]
            prev_matched_lifetime[i] += 1

            self.curr_features[code].append(grid_new_feature)
            after_ransac += 1
        self.num_features['after_ransac'] = after_ransac

    def track_features_mono(self):
        """
        Tracker features on the newly received mono images.
        """

    def add_new_features(self):
        """
        Detect new features on the image to ensure that the features are 
        uniformly distributed on the image.
        """
        curr_img = self.cam0_curr_img_msg.image
        grid_height, grid_width = self.get_grid_size(curr_img)

        # Create a mask to avoid redetecting existing features.
        mask = np.ones(curr_img.shape[:2], dtype='uint8')

        for feature in chain.from_iterable(self.curr_features):
            x, y = map(int, feature.cam0_point)
            mask[y-3:y+4, x-3:x+4] = 0

        # Detect new features.
        new_features = self.detector.detect(curr_img, mask=mask)

        # Collect the new detected features based on the grid.
        # Select the ones with top response within each grid afterwards.
        new_feature_sieve = [[] for _ in range(self.config.grid_num)]
        for feature in new_features:
            row = int(feature.pt[1] / grid_height)
            col = int(feature.pt[0] / grid_width)
            code = row * self.config.grid_col + col
            new_feature_sieve[code].append(feature)

        new_features = []
        for features in new_feature_sieve:
            if len(features) > self.config.grid_max_feature_num:
                features = sorted(features, key=lambda x:x.response, reverse=True)[:self.config.grid_max_feature_num]
            new_features.append(features)

        new_features = list(chain.from_iterable(new_features))

        # Find the stereo matched points for the newly detected features.
        cam0_points = [kp.pt for kp in new_features]
        cam1_points, inlier_markers = self.stereo_match(cam0_points)

        cam0_inliers, cam1_inliers, response_inliers = [], [], []
        for i, inlier in enumerate(inlier_markers):
            if not inlier:
                continue
            cam0_inliers.append(cam0_points[i])
            cam1_inliers.append(cam1_points[i])
            response_inliers.append(new_features[i].response)

        # Group the features into grids
        grid_new_features = [[] for _ in range(self.config.grid_num)]
        for i in range(len(cam0_inliers)):
            cam0_point = cam0_inliers[i]
            cam1_point = cam1_inliers[i]
            response = response_inliers[i]

            row = int(cam0_point[1] / grid_height)
            col = int(cam0_point[0] / grid_width)
            code = row*self.config.grid_col + col

            new_feature = FeatureMetaData()
            new_feature.response = response
            new_feature.cam0_point = cam0_point
            new_feature.cam1_point = cam1_point
            grid_new_features[code].append(new_feature)

        # Sort the new features in each grid based on its response.
        # And collect new features within each grid with high response.
        for i, new_features in enumerate(grid_new_features):
            for feature in sorted(new_features, key=lambda x:x.response, reverse=True)[:self.config.grid_min_feature_num]:
                self.curr_features[i].append(feature)
                self.curr_features[i][-1].id = self.next_feature_id
                self.curr_features[i][-1].lifetime = 1
                self.next_feature_id += 1

    def prune_features(self):
        """
        Remove some of the features of a grid in case there are too many 
        features inside of that grid, which ensures the number of features 
        within each grid is bounded.
        """
        for i, features in enumerate(self.curr_features):
            # Continue if the number of features in this grid does
            # not exceed the upper bound.
            if len(features) <= self.config.grid_max_feature_num:
                continue
            self.curr_features[i] = sorted(features, key=lambda x:x.lifetime, 
                reverse=True)[:self.config.grid_max_feature_num]

    def publish(self):
        """
        Publish the features on the current image including both the 
        tracked and newly detected ones.
        """
        curr_ids = []
        curr_cam0_points = []
        curr_cam1_points = []
        for feature in chain.from_iterable(self.curr_features):
            curr_ids.append(feature.id)
            curr_cam0_points.append(feature.cam0_point)
            curr_cam1_points.append(feature.cam1_point)

        curr_cam0_points_undistorted = self.undistort_points(
            curr_cam0_points,
            self.cam0_intrinsics,
            self.cam0_distortion_model,
            self.cam0_distortion_coeffs)

        curr_cam1_points_undistorted = self.undistort_points(
            curr_cam1_points,
            self.cam1_intrinsics,
            self.cam1_distortion_model,
            self.cam1_distortion_coeffs)

        features = []
        for i in range(len(curr_ids)):
            fm = FeatureMeasurement()
            fm.id = curr_ids[i]
            fm.u0 = curr_cam0_points_undistorted[i][0]
            fm.v0 = curr_cam0_points_undistorted[i][1]
            fm.u1 = curr_cam1_points_undistorted[i][0]
            fm.v1 = curr_cam1_points_undistorted[i][1]
            features.append(fm)

        feature_msg = namedtuple('feature_msg', ['timestamp', 'features'])(
            self.cam0_curr_img_msg.timestamp, features)
        return feature_msg

    def integrate_imu_data(self):
        """
        Integrates the IMU gyro readings between the two consecutive images, 
        which is used for both tracking prediction and 2-point RANSAC.

        Returns:
            cam0_R_p_c: a rotation matrix which takes a vector from previous 
                cam0 frame to current cam0 frame.
            cam1_R_p_c: a rotation matrix which takes a vector from previous 
                cam1 frame to current cam1 frame.
        """
        # Find the start and the end limit within the imu msg buffer.
        idx_begin = None
        for i, msg in enumerate(self.imu_msg_buffer):
            if msg.timestamp >= self.cam0_prev_img_msg.timestamp - 0.01:
                idx_begin = i
                break

        idx_end = None
        for i, msg in enumerate(self.imu_msg_buffer):
            if msg.timestamp >= self.cam0_curr_img_msg.timestamp - 0.004:
                idx_end = i
                break

        if idx_begin is None or idx_end is None:
            return np.identity(3), np.identity(3)

        # Compute the mean angular velocity in the IMU frame.
        mean_ang_vel = np.zeros(3)
        for i in range(idx_begin, idx_end):
            mean_ang_vel += self.imu_msg_buffer[i].angular_velocity

        if idx_end > idx_begin:
            mean_ang_vel /= (idx_end - idx_begin)

        # Transform the mean angular velocity from the IMU frame to the 
        # cam0 and cam1 frames.
        cam0_mean_ang_vel = self.R_cam0_imu.T @ mean_ang_vel
        cam1_mean_ang_vel = self.R_cam1_imu.T @ mean_ang_vel

        # Compute the relative rotation.
        dt = self.cam0_curr_img_msg.timestamp - self.cam0_prev_img_msg.timestamp
        cam0_R_p_c = cv2.Rodrigues(cam0_mean_ang_vel * dt)[0].T
        cam1_R_p_c = cv2.Rodrigues(cam1_mean_ang_vel * dt)[0].T

        # Delete the useless and used imu messages.
        self.imu_msg_buffer = self.imu_msg_buffer[idx_end:]
        return cam0_R_p_c, cam1_R_p_c

    def rescale_points(self, pts1, pts2):
        """
        Arguments:
            pts1: first set of points.
            pts2: second set of points.

        Returns:
            pts1: scaled first set of points.
            pts2: scaled second set of points.
            scaling_factor: scaling factor
        """
        scaling_factor = 0
        for pt1, pt2 in zip(pts1, pts2):
            scaling_factor += np.linalg.norm(pt1)
            scaling_factor += np.linalg.norm(pt2)

        scaling_factor = (len(pts1) + len(pts2)) / scaling_factor * np.sqrt(2)

        for i in range(len(pts1)):
            pts1[i] *= scaling_factor
            pts2[i] *= scaling_factor

        return pts1, pts2, scaling_factor

    def get_grid_size(self, img):
        """
        # Size of each grid.
        """
        grid_height = int(np.ceil(img.shape[0] / self.config.grid_row))
        grid_width  = int(np.ceil(img.shape[1] / self.config.grid_col))
        return grid_height, grid_width

    def predict_feature_tracking(self, input_pts, R_p_c, intrinsics):
        """
        predictFeatureTracking Compensates the rotation between consecutive 
        camera frames so that feature tracking would be more robust and fast.

        Arguments:
            input_pts: features in the previous image to be tracked.
            R_p_c: a rotation matrix takes a vector in the previous camera 
                frame to the current camera frame. (matrix33)
            intrinsics: intrinsic matrix of the camera. (vec3)

        Returns:
            compensated_pts: predicted locations of the features in the 
                current image based on the provided rotation.
        """
        # Return directly if there are no input features.
        if len(input_pts) == 0:
            return []

        # Intrinsic matrix.
        K = np.array([
            [intrinsics[0], 0.0, intrinsics[2]],
            [0.0, intrinsics[1], intrinsics[3]],
            [0.0, 0.0, 1.0]])
        H = K @ R_p_c @ np.linalg.inv(K)

        compensated_pts = []
        for i in range(len(input_pts)):
            p1 = np.array([*input_pts[i], 1.0])
            p2 = H @ p1
            compensated_pts.append(p2[:2] / p2[2])
        return np.array(compensated_pts, dtype=np.float32)

    def stereo_match(self, cam0_points):
        """
        Matches features with stereo image pairs.

        Arguments:
            cam0_points: points in the primary image.

        Returns:
            cam1_points: points in the secondary image.
            inlier_markers: 1 if the match is valid, 0 otherwise.
        """
        cam0_points = np.array(cam0_points)
        if len(cam0_points) == 0:
            return []

        R_cam0_cam1 = self.R_cam1_imu.T @ self.R_cam0_imu
        cam0_points_undistorted = self.undistort_points(
            cam0_points,
            self.cam0_intrinsics,
            self.cam0_distortion_model,
            self.cam0_distortion_coeffs,
            R_cam0_cam1)

        cam1_points = self.distort_points(
            cam0_points_undistorted,
            self.cam1_intrinsics,
            self.cam1_distortion_model,
            self.cam1_distortion_coeffs)
        cam1_points_copy = cam1_points.copy()

        # Track features using LK optical flow method.
        cam0_points = cam0_points.astype(np.float32)
        cam1_points = cam1_points.astype(np.float32)
        cam1_points, inlier_markers, _ = cv2.calcOpticalFlowPyrLK(
            self.curr_cam0_pyramid,
            self.curr_cam1_pyramid,
            cam0_points,
            cam1_points,
            **self.config.lk_params)

        cam0_points_, _, _ = cv2.calcOpticalFlowPyrLK(
            self.curr_cam1_pyramid,
            self.curr_cam0_pyramid,
            cam1_points,
            cam0_points.copy(),
            **self.config.lk_params)
        err = np.linalg.norm(cam0_points - cam0_points_, axis=1)

        disparity = np.abs(cam1_points_copy[:, 1] - cam1_points[:, 1])
        inlier_markers = np.logical_and.reduce([inlier_markers.reshape(-1), err < 3, disparity < 20])

        # Mark those tracked points out of the image region as untracked.
        img = self.cam1_curr_img_msg.image
        for i, point in enumerate(cam1_points):
            if not inlier_markers[i]:
                continue
            if (point[0] < 0 or point[0] > img.shape[1]-1 or 
                point[1] < 0 or point[1] > img.shape[0]-1):
                inlier_markers[i] = 0

        # Compute the relative rotation between the cam0 frame and cam1 frame.
        t_cam0_cam1 = self.R_cam1_imu.T @ (self.t_cam0_imu - self.t_cam1_imu)

        # Compute the essential matrix.
        E = skew(t_cam0_cam1) @ R_cam0_cam1

        # Further remove outliers based on the known essential matrix.
        cam0_points_undistorted = self.undistort_points(
            cam0_points,
            self.cam0_intrinsics,
            self.cam0_distortion_model,
            self.cam0_distortion_coeffs)

        cam1_points_undistorted = self.undistort_points(
            cam1_points,
            self.cam1_intrinsics,
            self.cam1_distortion_model,
            self.cam1_distortion_coeffs)

        norm_pixel_unit = 4.0 / (
            self.cam0_intrinsics[0] + self.cam0_intrinsics[1] +
            self.cam1_intrinsics[0] + self.cam1_intrinsics[1])

        for i in range(len(cam0_points_undistorted)):
            if not inlier_markers[i]:
                continue
            pt0 = np.array([*cam0_points_undistorted[i], 1.0])
            pt1 = np.array([*cam1_points_undistorted[i], 1.0])
            epipolar_line = E @ pt0
            error = np.abs((pt1 * epipolar_line)[0]) / np.linalg.norm(epipolar_line[:2])

            if error > self.config.stereo_threshold * norm_pixel_unit:
                inlier_markers[i] = 0

        return cam1_points, inlier_markers

    def undistort_points(self, pts_in, intrinsics, distortion_model, 
        distortion_coeffs, rectification_matrix=np.identity(3),
        new_intrinsics=np.array([1, 1, 0, 0])):
        """
        Arguments:
            pts_in: points to be undistorted.
            intrinsics: intrinsics of the camera.
            distortion_model: distortion model of the camera.
            distortion_coeffs: distortion coefficients.
            rectification_matrix:
            new_intrinsics:

        Returns:
            pts_out: undistorted points.
        """
        if len(pts_in) == 0:
            return []
        
        pts_in = np.reshape(pts_in, (-1, 1, 2))
        K = np.array([
            [intrinsics[0], 0.0, intrinsics[2]],
            [0.0, intrinsics[1], intrinsics[3]],
            [0.0, 0.0, 1.0]])
        K_new = np.array([
            [new_intrinsics[0], 0.0, new_intrinsics[2]],
            [0.0, new_intrinsics[1], new_intrinsics[3]],
            [0.0, 0.0, 1.0]])

        if distortion_model == 'equidistant':
            pts_out = cv2.fisheye.undistortPoints(pts_in, K, distortion_coeffs, rectification_matrix, K_new)
        else:
            # default: 'radtan'
            pts_out = cv2.undistortPoints(pts_in, K, distortion_coeffs, None, rectification_matrix, K_new)
        return pts_out.reshape((-1, 2))

    def distort_points(self, pts_in, intrinsics, distortion_model, 
            distortion_coeffs):
        """
        Arguments:
            pts_in: points to be distorted.
            intrinsics: intrinsics of the camera.
            distortion_model: distortion model of the camera.
            distortion_coeffs: distortion coefficients.

        Returns:
            pts_out: distorted points. (N, 2)
        """
        if len(pts_in) == 0:
            return []

        K = np.array([
            [intrinsics[0], 0.0, intrinsics[2]],
            [0.0, intrinsics[1], intrinsics[3]],
            [0.0, 0.0, 1.0]])

        if distortion_model == 'equidistant':
            pts_out = cv2.fisheye.distortPoints(pts_in, K, distortion_coeffs)
        else:   # default: 'radtan'
            homogenous_pts = cv2.convertPointsToHomogeneous(pts_in)
            pts_out, _ = cv2.projectPoints(homogenous_pts, 
                np.zeros(3), np.zeros(3), K, distortion_coeffs)
        return pts_out.reshape((-1, 2))

    def draw_features_stereo(self):
        img0 = self.cam0_curr_img_msg.image
        img1 = self.cam1_curr_img_msg.image

        kps0 = []
        kps1 = []
        matches = []
        for feature in chain.from_iterable(self.curr_features):
            matches.append(cv2.DMatch(len(kps0), len(kps0), 0))
            kps0.append(cv2.KeyPoint(*feature.cam0_point, 1))
            kps1.append(cv2.KeyPoint(*feature.cam1_point, 1))

        img = cv2.drawMatches(img0, kps0, img1, kps1, matches, None, flags=2)
        cv2.imshow('stereo features', img)
        cv2.waitKey(1)


def skew(vec):
    x, y, z = vec
    return np.array([
        [0, -z, y],
        [z, 0, -x],
        [-y, x, 0]])


def select(data, selectors):
    return [d for d, s in zip(data, selectors) if s]


if __name__ == '__main__':
    from queue import Queue
    from threading import Thread
    
    from config import ConfigEuRoC
    from dataset import EuRoCDataset, DataPublisher

    img_queue = Queue()
    imu_queue = Queue()

    config = ConfigEuRoC()
    image_processor = ImageProcessor(config)

    path = 'path/to/your/EuRoC_MAV_dataset/MH_01_easy'
    dataset = EuRoCDataset(path)
    dataset.set_starttime(offset=0.)
    
    duration = 3.
    ratio = 0.5
    imu_publisher = DataPublisher(
        dataset.imu, imu_queue, duration, ratio)
    img_publisher = DataPublisher(
        dataset.stereo, img_queue, duration, ratio)

    now = time.time()
    imu_publisher.start(now)
    img_publisher.start(now)

    def process_imu(in_queue):
        while True:
            msg = in_queue.get()
            if msg is None:
                return
            print(msg.timestamp, 'imu')
            image_processor.imu_callback(msg)
    t2 = Thread(target=process_imu, args=(imu_queue,))
    t2.start()

    while True:
        msg = img_queue.get()
        if msg is None:
            break
        print(msg.timestamp, 'image')
        # cv2.imshow('left', np.hstack([x.cam0_image, x.cam1_image]))
        # cv2.waitKey(1)
        # timestamps.append(x.timestamp)
        image_processor.stereo_callback(msg)

    imu_publisher.stop()
    img_publisher.stop()
    t2.join()
