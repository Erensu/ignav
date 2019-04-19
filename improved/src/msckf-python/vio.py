from queue import Queue
from threading import Thread
from multiprocessing import Lock
from config import ConfigEuRoC
from image import ImageProcessor
from msckf import MSCKF


class VIO(object):
    def __init__(self, config, img_queue, imu_queue, pos_queue, viewer=None):
        self.config = config
        self.viewer = viewer

        self.img_queue = img_queue
        self.imu_queue = imu_queue
        self.pos_queue = pos_queue
        self.feature_queue = Queue()

        self.mutex = Lock()

        self.image_processor = ImageProcessor(config)
        self.msckf = MSCKF(config)

        self.img_thread = Thread(target=self.process_img)
        self.imu_thread = Thread(target=self.process_imu)
        self.vio_thread = Thread(target=self.process_feature)
        self.pos_init_thread = Thread(target=self.position_init_callback)
        self.pos_upda_thread = Thread(target=self.process_position_update)
        self.img_thread.start()
        self.imu_thread.start()
        self.vio_thread.start()
        self.pos_init_thread.start()
        self.pos_upda_thread.start()

    def process_img(self):
        while True:
            if not self.msckf.isready():
                continue

            img_msg = self.img_queue.get()
            if img_msg is None:
                self.feature_queue.put(None)
                continue

            # viewer
            if self.viewer is not None:
                self.viewer.update_image(img_msg.cam0_image)

            feature_msg = self.image_processor.stereo_callback(img_msg)

            if feature_msg is not None:
                self.feature_queue.put(feature_msg)

    def process_imu(self):
        while True:
            if self.msckf.is_block:
                continue

            imu_msg = self.imu_queue.get()
            if imu_msg is None:
                continue

            self.image_processor.imu_callback(imu_msg)
            self.msckf.imu_callback(imu_msg)

    def process_feature(self):
        while True:
            if not self.msckf.isready():
                continue

            if self.msckf.is_block:
                continue

            feature_msg = self.feature_queue.get()
            if feature_msg is None:
                continue
            # print('feature_msg', feature_msg.timestamp)

            self.mutex.acquire()
            if self.config.monocam:
                result = self.msckf.feature_callback_mono(feature_msg)
            else:
                result = self.msckf.feature_callback(feature_msg)
            self.mutex.release()

            # viewer
            if result is not None and self.viewer is not None:
                self.viewer.update_pose(result.cam0_pose)

    def position_init_callback(self):
        while True:
            if self.msckf.is_block:
                continue

            pos_msg = self.pos_queue.get()
            if pos_msg is None:
                continue
            # print('pos_msg', pos_msg.timestamp)
            self.mutex.acquire()
            self.msckf.position_init_callback(pos_msg)
            self.mutex.release()

    def process_position_update(self):
        nhz = 100
        count = 0
        while True:
            if self.msckf.is_block:
                continue
            count += 1
            if count < nhz:
                continue
            count = 0
            # print('pos_msg', pos_msg.timestamp)
            self.mutex.acquire()
            self.msckf.position_update_callback()
            self.mutex.release()


if __name__ == '__main__':
    import time
    import argparse

    from dataset import EuRoCDataset, DataPublisher
    from viewer import Viewer

    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, default='path/to/your/EuRoC_MAV_dataset/MH_01_easy', help='Path of EuRoC MAV dataset.')
    parser.add_argument('--view', action='store_true', help='Show trajectory.')
    args = parser.parse_args()

    if args.view:
        viewer = Viewer()
    else:
        viewer = None

    img_queue = Queue()
    imu_queue = Queue()
    pos_queue = Queue()
    gt_queue = Queue()

    config = ConfigEuRoC()
    msckf_vio = VIO(config, img_queue, imu_queue, pos_queue, viewer=viewer)
    dataset = EuRoCDataset(args.path)

    if not config.init_pose:
        # start from static state
        dataset.set_starttime(offset=40.0)
    else:
        # start from first record
        dataset.set_starttime(offset=0.0)

    duration = float('inf')

    dataset.pos.readall(msckf_vio.msckf.all_position_data)

    # make it smaller if image processing and MSCKF computation is slow
    ratio = 0.4
    imu_publisher = DataPublisher(dataset.imu, imu_queue, duration, ratio)
    img_publisher = DataPublisher(dataset.stereo, img_queue, duration, ratio)
    pos_publisher = DataPublisher(dataset.pos, pos_queue, duration, ratio)

    now = time.time()
    imu_publisher.start(now)
    img_publisher.start(now)
    pos_publisher.start(now)
