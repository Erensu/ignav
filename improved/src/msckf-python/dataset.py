import numpy as np
import cv2
import os
import time

from collections import defaultdict, namedtuple
from threading import Thread


class GroundTruthReader(object):
    def __init__(self, path, scaler, starttime=-float('inf')):
        self.scaler = scaler   # convert timestamp from ns to second
        self.path = path
        self.starttime = starttime
        self.field = namedtuple('gt_msg', ['timestamp', 'p', 'q', 'v', 'bw', 'ba'])

    def parse(self, line):
        """
        line: (timestamp, p_RS_R_x [m], p_RS_R_y [m], p_RS_R_z [m], 
        q_RS_w [], q_RS_x [], q_RS_y [], q_RS_z [], 
        v_RS_R_x [m s^-1], v_RS_R_y [m s^-1], v_RS_R_z [m s^-1], 
        b_w_RS_S_x [rad s^-1], b_w_RS_S_y [rad s^-1], b_w_RS_S_z [rad s^-1], 
        b_a_RS_S_x [m s^-2], b_a_RS_S_y [m s^-2], b_a_RS_S_z [m s^-2])
        """
        line = [float(_) for _ in line.strip().split(',')]

        timestamp = line[0] * self.scaler
        p = np.array(line[1:4])
        q = np.array(line[4:8])
        v = np.array(line[8:11])
        bw = np.array(line[11:14])
        ba = np.array(line[14:17])
        return self.field(timestamp, p, q, v, bw, ba)

    def set_starttime(self, starttime):
        self.starttime = starttime

    def __iter__(self):
        with open(self.path, 'r') as f:
            next(f)
            for line in f:
                data = self.parse(line)
                if data.timestamp < self.starttime:
                    continue
                yield data


class IMUDataReader(object):
    def __init__(self, path, scaler, starttime=-float('inf')):
        self.scaler = scaler
        self.path = path
        self.starttime = starttime
        self.field = namedtuple('imu_msg', 
            ['timestamp', 'angular_velocity', 'linear_acceleration'])

    def parse(self, line):
        """
        line: (timestamp [ns],
        w_RS_S_x [rad s^-1], w_RS_S_y [rad s^-1], w_RS_S_z [rad s^-1],  
        a_RS_S_x [m s^-2], a_RS_S_y [m s^-2], a_RS_S_z [m s^-2])
        """
        line = [float(_) for _ in line.strip().split(',')]

        timestamp = line[0] * self.scaler
        wm = np.array(line[1:4])
        am = np.array(line[4:7])
        return self.field(timestamp, wm, am)

    def __iter__(self):
        with open(self.path, 'r') as f:
            next(f)
            for line in f:
                data = self.parse(line)
                if data.timestamp < self.starttime:
                    continue
                yield data

    def start_time(self):
        # return next(self).timestamp
        with open(self.path, 'r') as f:
            next(f)
            for line in f:
                return self.parse(line).timestamp

    def set_starttime(self, starttime):
        self.starttime = starttime


class PositionReader(object):
    def __init__(self, path, scaler, starttime=-float('inf')):
        self.scaler = scaler
        self.path = path
        self.starttime = starttime
        self.nhz = 100
        self.count = 0
        self.field = namedtuple('pos_msg', ['timestamp', 'position', 'variance', 'quaternion', 'velocity'])

    def parse(self, line):
        """
        line: (timestamp [ns],
        p_RS_R_x [m], p_RS_R_y [m], p_RS_R_z [m], q_RS_w [], q_RS_x [],
        q_RS_y [], q_RS_z [], v_RS_R_x [m s^-1], v_RS_R_y [m s^-1],
        v_RS_R_z [m s^-1], b_w_RS_S_x [rad s^-1], b_w_RS_S_y [rad s^-1],
        b_w_RS_S_z [rad s^-1], b_a_RS_S_x [m s^-2], b_a_RS_S_y [m s^-2],
        b_a_RS_S_z [m s^-2])
        """
        line = [float(_) for _ in line.strip().split(',')]

        timestamp = line[0] * self.scaler
        pos = np.array(line[1:4])

        white_noise = np.random.standard_normal(size=3)
        var = white_noise * 0.05

        quaternion = np.array(line[4:8])
        velocity = np.array(line[9:12])
        return self.field(timestamp, pos, var, quaternion, velocity)

    def __iter__(self):
        with open(self.path, 'r') as f:
            next(f)
            for line in f:
                data = self.parse(line)
                self.count += 1
                if data.timestamp < self.starttime or self.count < self.nhz:
                    continue
                self.count = 0
                yield data

    def start_time(self):
        # return next(self).timestamp
        with open(self.path, 'r') as f:
            next(f)
            for line in f:
                return self.parse(line).timestamp

    def set_starttime(self, starttime):
        self.starttime = starttime

    def readall(self, alldata):
        with open(self.path, 'r') as f:
            next(f)
            for line in f:
                alldata.append(self.parse(line))


class ImageReader(object):
    def __init__(self, ids, timestamps, starttime=-float('inf')):
        self.ids = ids
        self.timestamps = timestamps
        self.starttime = starttime
        self.cache = dict()
        self.idx = 0

        self.field = namedtuple('img_msg', ['timestamp', 'image'])

        self.ahead = 10   # 10 images ahead of current index
        self.wait = 1.5   # waiting time

        self.preload_thread = Thread(target=self.preload)
        self.thread_started = False

    def read(self, path):
        return cv2.imread(path, -1)
        
    def preload(self):
        idx = self.idx
        t = float('inf')
        while True:
            if time.time() - t > self.wait:
                return
            if self.idx == idx:
                time.sleep(1e-2)
                continue
            
            for i in range(self.idx, self.idx + self.ahead):
                if self.timestamps[i] < self.starttime:
                    continue
                if i not in self.cache and i < len(self.ids):
                    self.cache[i] = self.read(self.ids[i])
            if self.idx + self.ahead > len(self.ids):
                return
            idx = self.idx
            t = time.time()
    
    def __len__(self):
        return len(self.ids)

    def __getitem__(self, idx):
        self.idx = idx

        if idx in self.cache:
            img = self.cache[idx]
            del self.cache[idx]
        else:   
            img = self.read(self.ids[idx])
        return img

    def __iter__(self):
        for i, timestamp in enumerate(self.timestamps):
            if timestamp < self.starttime:
                continue
            yield self.field(timestamp, self[i])

    def start_time(self):
        return self.timestamps[0]

    def set_starttime(self, starttime):
        self.starttime = starttime


class Stereo(object):
    def __init__(self, cam0, cam1):
        assert len(cam0) == len(cam1)
        self.cam0 = cam0
        self.cam1 = cam1
        self.timestamps = cam0.timestamps

        self.field = namedtuple('stereo_msg', 
            ['timestamp', 'cam0_image', 'cam1_image', 'cam0_msg', 'cam1_msg'])

    def __iter__(self):
        for l, r in zip(self.cam0, self.cam1):
            assert abs(l.timestamp - r.timestamp) < 0.01, 'unsynced stereo pair'
            yield self.field(l.timestamp, l.image, r.image, l, r)

    def __len__(self):
        return len(self.cam0)

    def start_time(self):
        return self.cam0.starttime

    def set_starttime(self, starttime):
        self.starttime = starttime
        self.cam0.set_starttime(starttime)
        self.cam1.set_starttime(starttime)
        

class EuRoCDataset(object):   # Stereo + IMU
    '''
    path example: 'path/to/your/EuRoC Mav Dataset/MH_01_easy'
    '''
    def __init__(self, path):
        self.groundtruth = GroundTruthReader(os.path.join(path, 'mav0', 'state_groundtruth_estimate0', 'data.csv'), 1e-9)
        self.imu = IMUDataReader(os.path.join(path, 'mav0', 'imu0', 'data.csv'), 1e-9)
        self.cam0 = ImageReader(*self.list_imgs(os.path.join(path, 'mav0', 'cam0', 'data')))
        self.cam1 = ImageReader(*self.list_imgs(os.path.join(path, 'mav0', 'cam1', 'data')))
        self.pos = PositionReader(os.path.join(path, 'mav0', 'state_groundtruth_estimate0', 'data.csv'), 1e-9)

        self.stereo = Stereo(self.cam0, self.cam1)
        self.timestamps = self.cam0.timestamps

        self.starttime = max(self.imu.start_time(), self.stereo.start_time())
        self.starttime = max(self.pos.start_time(), self.starttime)
        self.set_starttime(0)

    def set_starttime(self, offset):
        self.groundtruth.set_starttime(self.starttime + offset)
        self.imu.set_starttime(self.starttime + offset)
        self.cam0.set_starttime(self.starttime + offset)
        self.cam1.set_starttime(self.starttime + offset)
        self.stereo.set_starttime(self.starttime + offset)
        self.pos.set_starttime(self.starttime + offset)

    def list_imgs(self, dir):
        xs = [_ for _ in os.listdir(dir) if _.endswith('.png')]
        xs = sorted(xs, key=lambda x:float(x[:-4]))
        timestamps = [float(_[:-4]) * 1e-9 for _ in xs]
        return [os.path.join(dir, _) for _ in xs], timestamps


# simulate the online environment
class DataPublisher(object):
    def __init__(self, dataset, out_queue, duration=float('inf'), ratio=1.0):
        self.dataset = dataset
        self.dataset_starttime = dataset.starttime
        self.out_queue = out_queue
        self.duration = duration
        self.ratio = ratio
        self.starttime = None
        self.started = False
        self.stopped = False

        self.publish_thread = Thread(target=self.publish)
        
    def start(self, starttime):
        self.started = True
        self.starttime = starttime
        self.publish_thread.start()

    def stop(self):
        self.stopped = True
        if self.started:
            self.publish_thread.join()
        self.out_queue.put(None)

    def publish(self):
        dataset = iter(self.dataset)
        while not self.stopped:
            try:
                data = next(dataset)
            except StopIteration:
                self.out_queue.put(None)
                return

            interval = data.timestamp - self.dataset_starttime
            if interval < 0:
                continue
            while (time.time() - self.starttime) * self.ratio < interval + 1e-3:
                time.sleep(1e-3)   # assumption: data frequency < 1000hz
                if self.stopped:
                    return

            if interval <= self.duration + 1e-3:
                self.out_queue.put(data)
            else:
                self.out_queue.put(None)
                return


if __name__ == '__main__':
    from queue import Queue

    path = '/home/sjl2019/Downloads/MH_01_easy'
    dataset = EuRoCDataset(path)
    dataset.set_starttime(offset=40)

    img_queue = Queue()
    imu_queue = Queue()
    gt_queue = Queue()
    pos_queue = Queue()

    duration = 1.0
    imu_publisher = DataPublisher(dataset.imu, imu_queue, duration)
    img_publisher = DataPublisher(dataset.stereo, img_queue, duration)
    gt_publisher = DataPublisher(dataset.groundtruth, gt_queue, duration)
    pos_publisher = DataPublisher(dataset.pos, pos_queue, duration)

    now = time.time()
    imu_publisher.start(now)
    img_publisher.start(now)
    gt_publisher.start(now)
    pos_publisher.start(now)

    def print_msg(in_queue, source):
        while True:
            x = in_queue.get()
            if x is None:
                return
            print(x.timestamp, source)

    t2 = Thread(target=print_msg, args=(imu_queue, 'imu'))
    t3 = Thread(target=print_msg, args=(gt_queue, 'groundtruth'))
    t4 = Thread(target=print_msg, args=(pos_queue, 'position'))
    t2.start()
    t3.start()
    t4.start()

    timestamps = []
    while True:
        x = img_queue.get()
        if x is None:
            break
        print(x.timestamp, 'image')
        cv2.imshow('left', np.hstack([x.cam0_image, x.cam1_image]))
        cv2.waitKey(1)
        timestamps.append(x.timestamp)

    imu_publisher.stop()
    img_publisher.stop()
    gt_publisher.stop()
    pos_publisher.stop()
    t2.join()
    t3.join()
    t4.join()

    print(f'\nelapsed time: {time.time() - now}s')
    print(f'dataset time interval: {timestamps[-1]} -> {timestamps[0]}'
          f'  ({timestamps[-1]-timestamps[0]}s)\n')
    print('Please check if IMU and image are synced')