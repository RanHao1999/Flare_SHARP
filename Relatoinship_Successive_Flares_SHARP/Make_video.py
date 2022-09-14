"""
Turn pictures into video
"""

import cv2, os

list_path = ['map_pil', 'bitmap_positive', 'bitmap_negative', 'Br']
for j in range(len(list_path)):
    path = r'res/'+list_path[j]
    files = os.listdir(path)
    files_use = [x for x in files if '.jpg' in x]
    count = len(files_use)

    fps = 30
    fourcc = cv2.VideoWriter_fourcc('M', 'J', 'P', 'G')
    videoWriter = cv2.VideoWriter(r'res/'+list_path[j] +'/'+list_path[j]+'.avi', fourcc, fps, (640,480))
    for i in range(1,count):
        img12 = cv2.imread(r'res/'+list_path[j]+'/'+str('%03d' % i)+'.jpg')
        videoWriter.write(img12)
    videoWriter.release()
