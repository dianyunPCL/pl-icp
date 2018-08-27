#!/usr/bin/python
# encoding: utf-8

from PIL import Image
# from skimage import io,data,color
import os.path
import matplotlib.pyplot as plt
import matplotlib.image as img
import math as m
import numpy as np


def read_trajectory(trajectory_file, poses):
    file = open(trajectory_file, "r")
    # poses.append([0., 0., 0.]) # the first pose set to [0,0,0]
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
    trajectory_file1 = "/home/vance/slam_ws/pl-icp/bin/080920-2-raw.txt"
    trajectory_file2 = "/home/vance/slam_ws/pl-icp/bin/080920-2-opt.txt"

    poses1 = []
    poses2 = []
    gp1_x = []
    gp1_y = []
    gp2_x = []
    gp2_y = []

    read_trajectory(trajectory_file1, poses1)
    read_trajectory(trajectory_file2, poses2)

    for pose in poses1:
        gp1_x.append(pose[0])
        gp1_y.append(pose[1])
    for pose in poses2:
        gp2_x.append(pose[0])
        gp2_y.append(pose[1])

    # plt.xlim(-6, 6)
    #plt.ylim(-2, 2)
    # plt.scatter(gp1_x, gp1_y, c='r', s=1)
    # plt.scatter(gp2_x, gp2_y, c='g', s=1)
    # plt.plot(gp1_x, gp1_y, '.r' )
    plt.plot(gp1_x, gp1_y, '.r', gp2_x, gp2_y, '.g')
    plt.show()

