"""
This code is for inframe the pictures with python.
"""
import os
from PIL import Image

width_i = 320
height_i = 240

row_max = 4
line_max = 4

pic_max = row_max*line_max

list_cat = ['original', 'pilmasked', 'original_calculated']

for k in range(len(list_cat)):
    dir_pics = 'res/SHARP_Goes_pic/' + list_cat[k]
    dir_names = os.listdir(dir_pics)

    coordinates_leftup_corner = []

    # Derive all the names of the wanted pictures
    pic_names = [x for x in dir_names if '.jpg' in x]
    pic_names.sort()

    # total image to save the output picture
    totImage = Image.new('RGBA', (width_i*line_max, height_i*row_max))
    num = 0

    for i in range(row_max):
        for j in range(line_max):
            # 1. Get the path of each picture.
            pic = Image.open(dir_pics + '/' + pic_names[num])
            width, height = pic.size
            resized_pic = pic.resize((width_i, height_i))

            # Calculate the coordinates of the left up corners of each picture
            loc = (int(j % row_max * width_i), int(i % line_max * height_i))
            print('The left uo corner of the no.'+ str(num) + ' is ' + str(loc))
            totImage.paste(resized_pic, loc)
            num += 1

            if num > len(pic_names):
                print('break here')
                break
        if num > pic_max:
            break
    totImage.save('res/SHARP_Goes_pic/'+ list_cat[k] +'/SHARP_goes.png')
