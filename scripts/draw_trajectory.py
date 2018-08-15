#!/usr/bin/python
# encoding: utf-8

from PIL import Image
# from skimage import io,data,color
import os.path
import matplotlib.pyplot as plt
import matplotlib.image as img
import math as m
import numpy as np

class YogoMap:
    def __init__(self):
        self.log_occ = -2
        self.log_free = 2
        self.width = 500
        self.height = 500
        self.resolution = 0.04
        self.log_max = 100.
        self.log_min = 0.
        self.origin_x = 0.
        self.origin_y = 0.
        self.offset_x = int(self.width / 2)
        self.offset_y = int(self.height / 2)
        self.map = []
        for i in range(self.width * self.height):
            self.map.append(200)

    def convertWWorld2GridIndex(self, x, y):
        index = []
        index.append(int(m.ceil((x - self.origin_x) / self.resolution + self.offset_x) ) )
        index.append(int(m.ceil((y - self.origin_y) / self.resolution + self.offset_y) ) )

        return index

    def gridIndex2LinearIndex(self, grid_index):
        linear_index = grid_index[1] + grid_index[0] * self.width

        return int(linear_index)

    def isValidGridIndex(self, grid_index):
        if grid_index[0] >=0 and grid_index[0] < self.width and grid_index[1] >= 0 and grid_index[1] < self.height:
            return True
        else:
            return False

    def occupanyMapping(self, scans, angles, flags, robot_poses):
        print("Start mapping! This will take a while...")
        if len(scans) == len(robot_poses):
            print "length equal to each other. all is well."
        else:
            print "length not equal to each other! scans:{0}, poses:{1}".format(len(scans),len(robot_poses))
            return
        
        frames = 0
        for i in range(len(scans)):
            frames += 1
            robot_pose = robot_poses[i]
            scan = scans[i]
            flag = flags[i]
            print "process {0}-th frame...".format(i)
            # robot index in map
            robot_index = self.convertWWorld2GridIndex(robot_pose[0], robot_pose[1])

            for id in range(len(scan)):
                dist = scan[id]
                angle = angles[id]
                valid = flag[id]
                if not valid:
                    continue
                if dist < 0. or dist > 10.:
                    continue

                theta = robot_pose[2]
                laser_x = dist * m.cos(angle)
                laser_y = dist * m.sin(angle)
                world_x = m.cos(theta) * laser_x - m.sin(theta) * laser_y + robot_pose[0]
                world_y = m.sin(theta) * laser_x + m.cos(theta) * laser_y + robot_pose[1]

                # update map cells
                scan_index = self.convertWWorld2GridIndex(world_x, world_y)
                # update hits
                if self.isValidGridIndex(scan_index):                    
                    hit_index = self.gridIndex2LinearIndex(scan_index)
                    self.map[hit_index] += self.log_occ
                    # print "scan hit {0} and update to {1}".format(hit_index,self.map[hit_index])
                # update misses
                line_indexs = self.traceLine(robot_index[0], robot_index[1], scan_index[0], scan_index[1])
                # print line_indexs

                for line_index in line_indexs:
                    li = self.gridIndex2LinearIndex(line_index)
                    self.map[li] += self.log_free
                    # print "scan miss {0} and update to {1}".format(li, self.map[li])
            
            # if (frames >= 20):
            #     break

        print "End occupany mapping! You can check the map."

    def swapTwoValue(self, x, y):
        x_tmp = x
        x = y
        y = x_tmp


    def traceLine(self, x0, y0, x1, y1):
        gridIndexList = []
        
        range_x = abs(x1 - x0) + 1
        range_y = abs(y1 - y0) + 1
        delta_x = -1
        delta_y = -1
        if x1 - x0 > 0:
            delta_x = 1
        if y1 - y0 >= 0:
            delta_y = 1

        if x1 != x0:
            # y = a*x+b
            a = float((y1 - y0))/(x1 - x0)
            b = float(y0 - a * x0)
            theta = m.atan2(y1 - y0, x1 - x0)
            
            x = x0
            y = y0
            for i in range(range_x):
                for j in range(range_y):
                    x = x0 + i * delta_x
                    y = y0 + j * delta_y
                    d = abs(a * x + b - y) * m.cos(theta) 
                    if d < m.sqrt(2) * 0.5 and (x != x1 or y != y1):
                        gridIndexList.append([x, y])
        else: # x1 = x0
            y = y0
            for j in range(range_y):
                y = y0 + j * delta_y
                gridIndexList.append([x0, y])

        return gridIndexList


def read_trajectory(trajectory_file, poses):
    file = open(trajectory_file, "r")
    poses.append([0., 0., 0.]) # the first pose set to [0,0,0]
    try:
        while (True):
            line = file.readline()
            if not line:
                break
            pose = [ float(i) for i in line.split()]
            poses.append(pose)
    except:
        raise
    finally:
        file.close()

def read_g2o_trajectory(trajectory_file, poses, nlines):
    file = open(trajectory_file, "r")
    try:
        l = -1
        while (True):
            l += 1
            line = file.readline()
            if l >= nlines:
                break
            p = [ i for i in line.split()]
            pp = p[2:]
            pose = list(( float(j) for j in p[2:]))
            poses.append(pose)
    except:
        raise
    finally:
        file.close()

def read_ranges(range_file_path, ranges, flags):
    path_dir = os.listdir(range_file_path)
    file_num = len(path_dir)
    for i in range(1,file_num+1):
        pd = range_file_path + str(i) + ".txt"
        file = open(pd, "r")
        ranges_single_file = []
        flags_single_file = []
        try:            
            line = file.readline() # skpi the first line
            while (True):
                line = file.readline()
                if not line:
                    break
                line_data = [int(i) for i in line.split()]
                ranges_single_file.append(line_data[0] * 0.001)
                flags_single_file.append(line_data[1])
        except:
            raise
        finally:
            file.close()
            
        ranges.append(ranges_single_file)
        flags.append(flags_single_file)

    print("read {0} files.".format(len(ranges)))


def setAngles(angles, n):
    delta_theta = m.pi / 720
    for i in range(n):
        angles.append((i - (len(angles)-1)/2) * delta_theta)


# def draw_map():
    # draw map
    # map = YogoMap()
    # map.occupanyMapping(ranges, angles, flags, poses)
    # draw_map = [[0 for i in range(map.width)] for i in range(map.height)]
    # for i in range(map.width):
    #     for j in range(map.height):
    #         index = map.gridIndex2LinearIndex([i, j])
    #         scale = map.map[index]
    #         if scale > 255:
    #             scale = 255
    #         elif scale < 0:
    #             scale = 0
    #         draw_map[i][j] = scale
    
    # add trajectory to map
    # gp = [0., 0., 0.]
    # global_poses_x = []
    # global_poses_y = []
    # plt.figure()
    # for pose in poses:
        # gp[0] += pose[0] * m.cos(gp[2]) - pose[1] * m.sin(gp[2])
        # gp[1] += pose[0] * m.sin(gp[2]) + pose[1] * m.cos(gp[2])
        # gp[2] += pose[2]
        # global_poses_x.append(pose[0])
        # global_poses_y.append(pose[1])
        # plt.arrow(pose[0], pose[1], m.cos(pose[2]), m.sin(pose[2]), width=0.01)
        # plt.hold
        # x, y = map.convertWWorld2GridIndex(gp[0], gp[1])
        # draw_map[x][y] = 50
    
    # dm = np.matrix(draw_map, dtype='float')
    # image_map = Image.fromarray(dm.astype(np.uint8))
    # plt.figure("map")
    # plt.imshow(image_map, cmap='Greys_r')
    # plt.show()


if __name__ =='__main__':
    point_num = 937
    trajectory_file1 = "/home/vance/slam_ws/csm_yogo/bin/080920-5.txt"
    trajectory_file2 = "/home/vance/slam_ws/csm_yogo/bin/080920-5-kf.g2o"

    poses1 = []
    poses2 = []
    gp1_x = []
    gp1_y = []
    gp2_x = []
    gp2_y = []

    read_trajectory(trajectory_file1, poses1)
    read_g2o_trajectory(trajectory_file2, poses2, 818)
    for pose in poses1:
        gp1_x.append(pose[0])
        gp1_y.append(pose[1])
    for pose in poses2:
        gp2_x.append(pose[0])
        gp2_y.append(pose[1])

    plt.xlim(-40, 40)
    plt.ylim(-15, 15)
    # plt.scatter(gp1_x, gp1_y, c='r', s=1)
    # plt.scatter(gp2_x, gp2_y, c='g', s=1)
    #plt.plot(gp1_x, gp1_y, '.r' )
    plt.plot(gp2_x, gp2_y, '-g')
    plt.show()
    

