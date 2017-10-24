import cv2
import os
import csv
import ast
import numpy as np
import sys

csv.field_size_limit(sys.maxsize)  # let's rock


class BoundingBox:
    '''
    We Use this bounding box extraction regularly so lets extract and abstract
        bool_array is 2d array where region to be bound
        (edges determined) is true and everything else is false
    '''
    @staticmethod
    def get_edges(bool_twodarray):
        xmin = 0
        ymin = 0
        xmax = bool_twodarray.shape[1] - 1
        ymax = bool_twodarray.shape[0] - 1
        # find xmin
        looking = True
        x = xmin
        while looking:
            col = bool_twodarray[:, x]
            if col.any():
                looking = False
                xmin = x
            x += 1
            if x >= xmax:
                xmin = 0
                looking = False
        # find xmax
        looking = True
        x = xmax
        while looking:
            col = bool_twodarray[:, x]
            if col.any():
                looking = False
                xmax = x
            x -= 1
            if x <= 0:
                xmax = bool_twodarray.shape[1] - 1
                looking = False
        # find ymin
        looking = True
        y = ymin
        while looking:
            row = bool_twodarray[y, :]
            if row.any():
                looking = False
                ymin = y
            y += 1
            if y >= ymax:
                ymin = 0
                looking = False
        # find xmax
        looking = True
        y = ymax
        while looking:
            row = bool_twodarray[y, :]
            if row.any():
                looking = False
                ymax = y
            y -= 1
            if y <= 0:
                ymax = bool_twodarray.shape[0] - 1
                looking = False
        return (xmin, xmax, ymin, ymax)


# atlas - rgb atlas file
# t - ndarray with dtype = bool
#     true for any non-white pixel element in rgb atlas file
# returns 4 element tuple (
# xmin - leftmost non-white atlas pixel (region)
# xmax - rightmost non-white atlas pixel (region)
# ymin - topmost non-while atlas pixel (region)
# ymax - bottommost non-white atlas pixel (region)
def get_edges(rgb_code):
    tb = (rgb_code[:, :, 0] != 255)
    tg = (rgb_code[:, :, 1] != 255)
    tr = (rgb_code[:, :, 2] != 255)
    t = tb | tg | tr
    return BoundingBox.get_edges(t)


#  return sorted list of frozensets corresponding to communities in consensus
#   only of images with threshold level
#   note this method assumes cmt structure is made up of grid cells
def cons_cmt_str(cons_cmt_csv_path, lvl):
    assert os.path.isfile(cons_cmt_csv_path), "No csv {}".format(
        cons_cmt_csv_path)
    cons_cmt_str = None

    with open(cons_cmt_csv_path) as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:
            if 'consensus com str' in row[0]:
                cons_cmt_str_lst_lst = ast.literal_eval(row[1])
                cons_cmt_str = []
                for lst in cons_cmt_str_lst_lst:
                    lvl_lst = []
                    for cell in lst:
                        if cell.startswith('({}:'.format(lvl)):
                            lvl_lst.append(cell)
                    cons_cmt_str.append(frozenset(lvl_lst))
    return cons_cmt_str


# return community index of key defined by lvel, hemisphere, col and row
#  or return None if no community found
def cmt_idx(cons_cmt_str, lvl, hemi, col, row):
    key = '({}:{}:{}:{})'.format(lvl, hemi, col, row)
    for idx, cons_cmt in enumerate(cons_cmt_str):
        if key in cons_cmt:
            return idx
    return None


# assumes and returns grayscale
def thresh_tif(thresh_tif_path):
    assert os.path.isfile(thresh_tif_path), \
        "No tif {}".format(thresh_tif_path)
    return cv2.imread(thresh_tif_path, cv2.IMREAD_GRAYSCALE)


# does not change color channels
def atlas_tif(atlas_tif_path):
    assert os.path.isfile(atlas_tif_path), \
        "No tif {}".format(atlas_tif_path)
    return cv2.imread(atlas_tif_path, cv2.IMREAD_UNCHANGED)


# get ending cell boundaries
def stops(grid_thresh_img, y, x, gcs, hemi, edges_tup):
    (xmin, xmax, ymin, ymax) = edges_tup
    midx = grid_thresh_img.shape[1] / 2
    y_stop = min(y + gcs, ymax)
    if hemi == 'l':
        x_stop = min(x + gcs, midx)
    if hemi == 'r':
        offset_x = (midx / gcs + 1) * gcs + 1
        if x < offset_x and midx % gcs > 0:
            x_stop = min(x + gcs, offset_x)
        else:
            x_stop = min(x + gcs, xmax)
    return (x_stop, y_stop)


# return img corresponding to gcs size squared grid cell at row, col
#  edges_tup - (xmin, xmax, ymin, ymax)
#  row - row in pixels (not grid cells)
#  col - col in pixels (not grid cells)
#  gcs - grid cell size
def cell_img(grid_thresh_img, y, x, gcs, hemi, edges_tup):
    (x_stop, y_stop) = stops(grid_thresh_img=grid_thresh_img,
                             y=y, x=x, gcs=gcs, hemi=hemi,
                             edges_tup=edges_tup)
    cell_img = grid_thresh_img[y:y_stop, x:x_stop]
    return cell_img


# sort of complicated to determining number of channels
def num_channels(img):
    if len(img.shape) <= 2:
        return 1
    else:
        return img.shape[2]


# similar to cell_img(), but paste cell_img at location
#  idx_lst is the index into the color channel... hmmm, could be confusing
#  does not return anything
def paste_cell_img(cell_img, y, x, gcs, hemi, edges_tup, grid_thresh_img):
    assert num_channels(grid_thresh_img) == num_channels(cell_img), \
        "grid_thresh_img {} channels, cell_img {}".format(
            num_channels(grid_thresh_img), num_channels(cell_img))
    (x_stop, y_stop) = stops(grid_thresh_img=grid_thresh_img,
                             y=y, x=x, gcs=gcs, hemi=hemi,
                             edges_tup=edges_tup)
    grid_thresh_img[y:y_stop, x:x_stop] = cell_img


def has_thresh(cell_img):
    # do bitwise_not since are actually testing for black since thresh
    #  represented as black... bitwise_not converts black to white, then
    #  we test for white
    # any() checks if anything evaluates to True i.e. not zero
    return (cv2.bitwise_not(cell_img).any())


# does not change color channels
def gray2bgra_tif(tif_path):
    assert os.path.isfile(tif_path), \
        "No tif {}".format(tif_path)
    gray_img = cv2.imread(tif_path, cv2.IMREAD_GRAYSCALE)
    return cv2.cvtColor(gray_img, cv2.COLOR_GRAY2BGRA)


def visual_img(grid_ref_tif_path, output_img_path, to_erode_compose_img):
    # read and convert ref img as BGR
    grid_ref_img_bgra = cv2.imread(grid_ref_tif_path, cv2.IMREAD_UNCHANGED)
    grid_ref_img = grid_ref_img_bgra[:, :, :3]  # just use the bgr part (if a)

    if "degenerate" in output_img_path:
        eroded_img = erode(img=to_erode_compose_img)
        visual_img = compose(thresh_img=eroded_img,
                             ref_img=grid_ref_img)
    else:
        visual_img = compose(thresh_img=to_erode_compose_img,
                             ref_img=grid_ref_img)
    return visual_img


# expects BGRA thresh_img (or will convert) and BGR ref_img (will not convert)
def compose(thresh_img, ref_img):
    assert num_channels(ref_img) == 3, \
        "expected 3 channels in ref_img found {}".format(num_channels(ref_img))

    # convert to BGRA thresh img if needed
    if num_channels(thresh_img) < 3:
        bgra_thresh_img = cv2.cvtColor(thresh_img, cv2.COLOR_GRAY2BGRA)
    elif num_channels(thresh_img) < 4:
        bgra_thresh_img = cv2.cvtColor(thresh_img, cv2.COLOR_BGR2BGRA)
    else:
        bgra_thresh_img = thresh_img

    # convert all white of threshold image to transparent
    bgra_thresh_img[np.where(
        (bgra_thresh_img == [255, 255, 255, 255]).all(axis=2))] = \
        [255, 255, 255, 0]

    # now create masks and return blend
    overlay_img = bgra_thresh_img[:, :, :3]
    overlay_msk = bgra_thresh_img[:, :, 3:]
    background_msk = 255 - overlay_msk  # everything background is 255

    # convert masks to 3 channel
    overlay_msk = cv2.cvtColor(overlay_msk, cv2.COLOR_GRAY2BGR)
    background_msk = cv2.cvtColor(background_msk, cv2.COLOR_GRAY2BGR)

    # now the magical part, convert to 0 - 1 and multiply with mask
    msked_ref_img = \
        (ref_img * (1 / 255.0)) * (background_msk * (1 / 255.0))
    msked_thresh_img = \
        (overlay_img * (1 / 255.0)) * (overlay_msk * (1 / 255.0))

    return np.uint8(cv2.addWeighted(msked_ref_img, 255.0,
                                    msked_thresh_img, 255.0,
                                    0.0))


# if cell_img is one channel then convert
#  returns BGRA image
def clr_thresh(cell_img, clr_idx):
    if(num_channels(cell_img) != 4):
        # got this technique from
        #  https://stackoverflow.com/questions/14786179/
        #  how-to-convert-a-1-channel-image-into-a-3-channel-with-opencv2
        new_img = cv2.cvtColor(cell_img, cv2.COLOR_GRAY2BGRA)

    else:
        new_img = cell_img

    assert num_channels(new_img) == 4, "cell image only has {} channels".\
        format(num_channels(new_img))

    clr_arr = [[0,   0,   255, 255],
               [0,   255,   0, 255],
               [255,    0,  0, 255]]
    assert clr_idx < len(clr_arr), "Only {} colors supported".format(
        len(clr_arr))

    new_img[np.where((new_img == [0, 0, 0, 255]).all(axis=2))] = \
        clr_arr[clr_idx]
    return new_img


def erode(img):
    kernel = np.ones((50, 50), np.uint8)
    return cv2.erode(img, kernel, iterations=1)
