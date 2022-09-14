import os

pathes = [
    'data', 'data/JSOC',
    'res', 'res/bitmap_negative', 'res/bitmap_positive',
    'res/Br', 'res/coord_pil', 'res/map_pil',
    'res/SHARP_csv', 'res/SHARP_csv/csv_all', 'res/SHARP_csv/csv_correlation', 'res/SHARP_csv/csv_Goes', 'res/SHARP_CSV/csv_SF',
    'res/SHARP_Goes_pic', 'res/SHARP_Goes_pic/original', 'res/SHARP_Goes_pic/original_calculated', 'res/SHARP_Goes_pic/pilmasked',
    'res/SHARP_pic'
]
for path in pathes:
    if os.path.exists(path) == False:
        os.mkdir(path)