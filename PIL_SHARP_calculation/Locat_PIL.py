# This is the main programme for calculating Polarity Inversion Line with built wheel.
# sunpy, sklearn are the two main packages.
"""
Brief introduction: This is the main programme for locating the area of polarity inversion line of a series of active regions with multiple cores.

Last modified: 2021.07.12 by RanHao

file structure:
             > PIL (main folder)
               Locate_PIL.py (This code)
               > data
               > res
                 > Br
                 > bitmap_positive
                 > bitmap_negative
                 > map_pil
                 > coord_pil

input: Br.fits files from JSOC.

output: .txt file containing coordinates of the polarity inversion line in 'coord_pil' folder
        Maps of Br, bitmap_positive, bitmap_negative, map_pil in the rest of the folders in 'res'

principle: 1. Read the Br.fits files and get the map.data, which is supposed to be a two-dimension list
           2. Set sigma = 100G (Hoekesema et al., 2014), we use 2*sigma as the threshold to get bitmap_positive and bitmap_negative.
              In bitmap_positive, any pixel with value mt 200G is set to be 1, the rest are 0.
              In bitmap_negative, any pixel with value lt -200G is set to be -1, the rest are 0.
           3. We use DBSCAN algorithm on the bitmaps, which is for removing small clusters that may influence the result.
           4. We set up a Gaussian kernel with width as 10 pixels(width is adjustable).
              Then convolve the two bitmaps with the Gaussian kernel.
           5. Multiply the two bitmaps that is Gaussian-kernel processed to get the needed PIL-Inversion_Line map.

==================================================================================================================
"""

import matplotlib.pyplot as plt
import numpy as np, sunpy, os, pandas as pd, sys, sunpy.map, sunpy.coordinates, concurrent.futures
from sklearn.cluster import DBSCAN
from scipy.ndimage import gaussian_filter
from astropy.coordinates import SkyCoord
from sunpy.io.special import srs


def get_coordinates(B_map_data, threshold, sign):
    """
    :param map: Directly get from the Br.fits.
    :param threshold: According to Hoeksema, this value is often chosen to be 200.
    :param sign: string, 'positive' or 'negative'
    :return: a list of coordinates, which show the position of pixels that meet the threshold.
    """
    coordinates = []
    if sign == 'positive':
        for i in range(B_map_data.shape[0]):
            for j in range(B_map_data.shape[1]):
                if B_map_data[i][j] > threshold:
                    coordinates.append([i,j])
    else:
        for i in range(B_map_data.shape[0]):
            for j in range(B_map_data.shape[1]):
                if B_map_data[i][j] < -threshold:
                    coordinates.append([i,j])
    return coordinates

def read_data(file_Br):
    """
    Read fits data and return any needed parameters.
    :param file_Br: Br.fits
    :return: br as a map
    """
    try:
        br = sunpy.map.Map(file_Br, autoalign = True)
    except:
        print("Fail to read the br.fits data")
        sys.exit(1)
    return br

def find_biggest_cluster(cluster,n):
    """
    :param cluster: cluster is directly get from DBSCAN algorithm
    :param n: the number of clusters one wants to select, we give out the top n biggest clusters
    :return: a list of n two-element list in which the two elements are the cluster index and the number of points in the cluster.
    """
    cluster_index = cluster.labels_.tolist() # Turn np.array to list for convenience
    count_n = []
    count_n_tem = []
    for i in range(max(cluster_index)):
        count_n.append(cluster_index.count(i))
        count_n_tem .append([i, cluster_index.count(i)])
    count_n.sort()
    count_n.reverse()
    count_n = count_n[0:n]

    return [count_n_tem[i] for i in range(len(count_n_tem)) if count_n_tem[i][1] in count_n]

def coordinates_of_clusters(cluster_dbscan,cluster_fbc,n_cluster):
    """
    This function will return the coordinates of the clusters one selected.
    :param cluster_dbscan: cluster variables directly from DBSCAN
    :param cluster_fbc: cluster variables from find_biggest_cluster(fbc)
    :return: a n-element list, in which the elements are coordinates of the chosen clusters.
    """
    coordinates_nc = []
    for times in range(n_cluster):
        coordinates_nc.append([])
    labels = [x for x in cluster_dbscan.labels_ if x != -1]
    for i in range(len(labels)):
        for j in range(n_cluster):
            if labels[i] == cluster_fbc[j][0]:
                coordinates_nc[j].append(cluster_dbscan.components_[i].tolist())

    return coordinates_nc

def coordinates2map(map_size, coordinates, sign):
    """
    :param map_size: Size of the map data, list
    :param coordinates: The coordinates that will be set non-zero.
    :param sign: 'positive' or 'negative', wihch decides the sign of the bitmap.
    :return: a bitmap data in which the chosen coordinates are non-zero which others are zero.
    """
    if sign == 'positive':
        a = 1
    else:
        a = -1

    map = np.zeros(map_size).tolist()
    for i in range(len(coordinates)):
        for j in range(len(coordinates[i])):
            map[int(coordinates[i][j][0])][int(coordinates[i][j][1])] = a

    return map

def map_PIL(bitmapp, bitmapn):
    """
    using "dilate" algorithm to expand the two kinds of bitmaps to get the pil area.
    :param bitmapp: the positive bitmap
    :param bitmapn: the negative bitmap
    :return: a list containing the map data of the PIL and the coordinates of the PIL.
    """
    size = [len(bitmapp), len(bitmapp[0])]
    map = np.zeros(size).tolist()
    coordinates = []
    for i in range(size[0]):
        for j in range(size[1]):
            if bitmapp[i][j] != 0 and bitmapn[i][j] != 0:
                map[i][j] = bitmapp[i][j]
                coordinates.append([i, j])
    print(len(coordinates))
    return map, coordinates

def change_mapdata(br, new_mapdata):
    # Turn br.data into new_mapdata, for that we can't update the data directly.
    shape = br.data.shape
    for i in range(shape[0]):
        for j in range(shape[1]):
            br.data[i][j] = new_mapdata[i][j]

    return br


def main_run(input):
    file_br = input[0]
    index = input[1]
    br = read_data(file_br)

    # Get bitmap coordinates step:
    br_mapdata = br.data
    coordinates_posi = get_coordinates(br_mapdata, 200, 'positive')
    coordinates_nega = get_coordinates(br_mapdata, 200, 'negative')

    # DBSCAN to get clusters step:
    cluster_positive = DBSCAN(eps=1., min_samples=2).fit(coordinates_posi)
    cluster_negative = DBSCAN(eps=1., min_samples=2).fit(coordinates_nega)

    # Note: DBSCAN returns a .label_, which tells which cluster every single point belongs to.
    n_cluster = 5 # n_cluster is the number of clusters one selects
    clusterp = find_biggest_cluster(cluster_positive, n_cluster) # Return top 3 biggest clusters of positive clusters and negative clusters.
    clustern = find_biggest_cluster(cluster_negative, n_cluster)

    # The following codes tells the coordinates of the top n biggest clusters
    # For positive and negative, it both returns a 3-element list, in which every element is the list of the coordinates of the cluster points.
    coordinates_ncp = coordinates_of_clusters(cluster_positive, clusterp, n_cluster)
    coordinates_ncn = coordinates_of_clusters(cluster_negative, clustern, n_cluster)
    # "coordinates_ncp" means coordinates of the n positive clusters. Construction please refer to the function 'coordinates_of_clusters'
    # Following is using Gaussian kernel to process the two bitmaps which contains the top n biggest clusters.
    # 1. Turn the known coordinates to maps.
    size_map = [len(br_mapdata), len(br_mapdata[0])] # list object has no attribute 'shape'
    bitmap_positive = coordinates2map(size_map, coordinates_ncp, 'positive')
    bitmap_negative = coordinates2map(size_map, coordinates_ncn, 'negative')

    # Using Gaussian kernel to convolve the two bitmaps.
    bitmapp_gs = gaussian_filter(bitmap_positive, sigma = 2)
    bitmapn_gs = gaussian_filter(bitmap_negative, sigma = 2)  #option1: using gaussian kernel, not very good result

    # Using dilate on bitmaps directly to find the overlap area.
    map_pil, coordinates_pil = map_PIL(bitmapp_gs, bitmapn_gs)

    # Draw the maps in the following stop:
    br.plot()  # Original Br pic
    plt.savefig(r'res/Br/' + str('%03d' % index) + '.jpg')

    br_temp = br
    change_mapdata(br_temp, bitmap_positive)
    br_temp.plot()
    plt.savefig(r'res/bitmap_positive/' + str('%03d' % index) + '.jpg')

    change_mapdata(br_temp, bitmap_negative)
    br_temp.plot()
    plt.savefig(r'res/bitmap_negative/' + str('%03d' % index) + '.jpg')

    change_mapdata(br_temp, map_pil)
    br_temp.plot()
    plt.savefig(r'res/map_pil/' + str('%03d' % index) + '.jpg')

    coor_path = open(r'res/coord_pil/' + str('%03d' % index) + '.txt', 'w')
    for i in range(len(coordinates_pil)):
        coor_path.write(str(coordinates_pil[i][0]))
        coor_path.write(' ')
        coor_path.write(str(coordinates_pil[i][1]))
        coor_path.write('\n')
    coor_path.close()

#    print('break here')

    return 0

def main():
    path_br = r'data/JSOC'
    list_br = os.listdir(path_br)
    list_br = [x for x in list_br if 'Br' in x and 'Br_err' not in x]
    list_br.sort()
    list_br = list_br[0:780]  #Rest of the data is broken
    input = []
    for index in range(len(list_br)):
        file_br = r'data/JSOC/' + list_br[index]
        temp = []
        temp.append(file_br)
        temp.append(index)
        input.append(temp)
    with concurrent.futures.ProcessPoolExecutor(max_workers = 20) as pool:
        pool.map(main_run, input)

    return 0

if __name__ == "__main__":
    main()

# Happy Coding!
# __author__ = RanHao



















