from skimage.io import imread
import png
import os
import sys
import uuid
from queue import Queue
from PIL import Image, ImageDraw
import copy

seg_q = Queue()

def gain_image_pixel_simple_matrix(image_path, num_bins=10):
    pixel_matrix = imread(image_path)
    pixel_max = 0
    pixel_min = 255

    # gain the floor and cieling of pixel values ( in theory you can assume a floor of 0 and a cieling of 255 and skip this step but i didnt want empty bins)
    for row in pixel_matrix:
        for col in row:
            try:
                rgb_ave = (float(int(col[0]) + int(col[1]) + int(col[2])) / 3)
                if rgb_ave >= pixel_max:
                    pixel_max = rgb_ave
                if rgb_ave <= pixel_min:
                    pixel_min = rgb_ave
            except Exception as ex:
                rgb_ave = col[0]

    print(f"DEBUG - pixel min : {pixel_min} - {pixel_max}")

    # set our bin size based on the pixel floor and ceiling and the number of bins we want and generate bins
    bin_size = (pixel_max-pixel_min)/num_bins
    bins = []
    current_bin_start = pixel_min
    for bin_ix in range(1, num_bins):
        bins.append({"bin_min":current_bin_start, "bin_max":(current_bin_start+bin_size), "pixels":[], "ave":0})
        current_bin_start = (current_bin_start+bin_size)
    bins.append({"bin_min":current_bin_start, "bin_max":pixel_max, "pixels":[], "ave":0})    
    print(f"DEBUG - bins : {bins}")


    # generate a matrix of values based on the index of the bin that the pixel average falls in
    # this is needed so that we can then assign the average value of the bin to each pixel
    indexed_pixel_matrix = []
    for row in pixel_matrix:
        indexed_row = []
        for col in row:
            try:
                rgb_ave = (float(int(col[0]) + int(col[1]) + int(col[2])) / 3)
                [pbin for pbin in bins if rgb_ave>=pbin['bin_min'] and rgb_ave<=pbin['bin_max']][0]['pixels'].append(rgb_ave)
                indexed_row.append([pbix for pbix, pbin in enumerate(bins) if rgb_ave>=pbin['bin_min'] and rgb_ave<=pbin['bin_max']][0])
            except Exception as ex:
                # print(f"ERROR - adding pixel : {ex}")
                [pbin for pbin in bins if rgb_ave>=pb['bin_min'] and rgb_ave<=pb['bin_max']][0]['pixels'].append(0)
                indexed_row.append(0)
        indexed_pixel_matrix.append(indexed_row)
    
    # average the bins
    for pbin in bins:
        try:
            pbin['ave'] = int(round(sum(pbin['pixels'])/len(pbin['pixels']),0))
        except:
            pbin['ave'] = 0
    print(f"DEBUG - bins averaged...")


    # tranform the orginal pixel matrix into a simple matrix of 1 and 0 values
    # a pixel is assigned a 1 if the pixel next to it to the right is in a different bin
    # this is technically a very naive MSE but rather than look at more surrounding pixels we are looking at just the next.
    # we can afford to do this on documents because typically we are operating on monochrome OR a small range of colours/likely bins.
    # a more traditional Sobel xformation in a domain where there are a majority white pixels can mean that we loose singular values like "1" or "I" etc.
    # this approach accounts for that by not performing a 'mean' of local pixels.    
    image_pixel_png_matrix = [] # only used for testing output 
    image_pixel_edge_matrix = [] # only used for testing output 
    image_pixel_simple_matrix = []
    for row in indexed_pixel_matrix:
        png_row = []
        png_row_edge = []
        simple_row = []
        for col_ix,col in enumerate(row):
            png_row.append([bins[col]['ave'],bins[col]['ave'],bins[col]['ave']])

            if col_ix <= len(row)-2 and row[col_ix+1] != col:
                png_row_edge.append([0, 0, 0])
                simple_row.append(1)
            else:
                png_row_edge.append([255, 255, 255])
                simple_row.append(0)

        image_pixel_png_matrix.append(png_row)
        image_pixel_edge_matrix.append(png_row_edge)
        image_pixel_simple_matrix.append(simple_row)

    # output for testing not required
    # costly file system operation
    name = os.path.join("tmp", "local_binned_image_"+str(uuid.uuid4())+".png")
    png.from_array(image_png_matrix, 'RGB').save(name)

    name = os.path.join("tmp", "local_binned_image_edge_"+str(uuid.uuid4())+".png")
    png.from_array(image_png_array_edge, 'RGB').save(name)

    # the next step, and this is because i hate leaving data on the table, is to gather as much statistically important data as possible.
    # for example:
    # - the average (this may be better as a modal average, but i hate doing modals) horizotal distance between two black pixels
    # can be used to indicate how far you should look to join/split bounding boxes. this is the horizontal_tolerance of the data on the page
    # - the average length of contiguous black pixels can be used to determine both the font size and also any harsh vertical lines that you want to remove 
    # (vertical lines are an issue because they will force a bounding box to breach multiple lines when you perform a line segment trace).

    horizontal_space_between_black_pixels = []
    for row in image_pixel_simple_matrix:
        if row.count(1) >= len(row)*0.05:
            for col_ix, col in enumerate(row):
                if col == 1 and row[col_ix+1] == 0:
                    for i in range(col_ix+1, len(row)-1):
                        if row[i] == 1:
                            horizontal_space_between_black_pixels.append((i-col_ix))
                            break

    horizontal_space_between_black_pixels = list(set(horizontal_space_between_black_pixels))
    horizontal_space_between_black_pixels.sort()
    print(f"DEBUG - white spaces : {horizontal_space_between_black_pixels}")
    horizontal_tolerance = int(round(sum(horizontal_space_between_black_pixels)/len(horizontal_space_between_black_pixels)*0.035,0))
    print(f"DEBUG - horizontal tolerance : {horizontal_tolerance}")

    # generate the vertical tolerance
    # this may be more efficient in numpy as this is just an array transformation / rotation and i know you can do a .rot()
    # but i am still stuck in the late 90s and this dosent feel to costly... 
    vertical_segments = []
    for col_ix in range(0, len(image_pixel_simple_matrix[0])):
        # check the column has at least two black pixels, i.e. there is statistical line value in the column
        if sum([col[col_ix] for col in [row for row in image_pixel_simple_matrix]]) >= 2:            
            # gain all pixels in that column and walk them from each black pixel (1) until you hit a 0
            pxl_pos = -1            
            xformed_pxls = [col[col_ix] for col in [row for row in image_pixel_simple_matrix]]
            for pxl_ix, pxl in enumerate(xformed_pxls):
                # basically ignore any pixel you already included as part of a previous line segment.
                if pxl_ix >= pxl_pos and (pxl_ix+1) <= (len(xformed_pxls)-1) and pxl == 1 and xformed_pxls[pxl_ix+1] == 1:
                    pxl_pos = pxl_ix+1                    
                    while xformed_pxls[pxl_pos] == 1:
                        pxl_pos += 1
                        if pxl_pos >= len(xformed_pxls)-1:
                            pxl_pos = len(xformed_pxls)-1
                            break

                    # 10 is entirely arbitary here, this number should be a function of the height of the page, but 10 feels like a good floor.
                    if pxl_pos - pxl_ix >= 10:
                        # we are storing the start row (pix_ix - pixel index) the start column (col_ix) and the line segment length               
                        vertical_segments.append([pxl_ix, col_ix, pxl_pos-pxl_ix])
                        
    print(f"DEBUG - vertical segments : {len(vertical_segments)}")    

    # this final step uses the derived vertical tolerance to remove harsh lines from the simple pixel matrix.
    # in essence any vertical segment that looks like it breaches the average by x amount is overwritten with 0 values.
    vertical_tolerance = sum([v[2] for v in vertical_segments])/len(vertical_segments)
    for vs in vertical_segments:
        if vs[2] >= vertical_tolerance*1.5:
            print(f"DEBUG - likely hard line to remove at : {vs} : {vertical_tolerance*1.5}")
            for i in range(vs[0], vs[0]+vs[2]):
                image_pixel_simple_matrix[i][vs[1]] = 0

    return image_pixel_simple_matrix, horizontal_tolerance, vertical_tolerance

def gain_bounding_boxes(image_pixel_simple_matrix, horizontal_tolerance, vertical_tolerance, type=None):
    # the purpose of this function is to walk/trace line segments (pixels with a value of 1) until they are fully enumerated and then 
    # hallucinate a box around the last pixel and IF a black pixel exists in that projection hop to that pixel and continue.
    # we use the horizontal and vertical tolerances to project that box.

    bounding_boxes = []
    for row_ix, row in enumerate(image_pixel_simple_matrix):
        for col_ix, col in enumerate(row):
            if col == 1:
                # pixel value is 1, start enumerating.                

                # set some min/max values for the current bounding box
                sx = col_ix
                sy = row_ix
                ex = col_ix+1
                ey = row_ix
            

                # generate a projected box around the current pixel
                pixel_x_min = int(round(col_ix-(horizontal_tolerance*0.5),0))
                pixel_x_max = col_ix+horizontal_tolerance
                pixel_y_min = int(round(row_ix-(horizontal_tolerance*0.25),0))
                pixel_y_max = int(round(row_ix+(horizontal_tolerance*0.25),0))       

                if type == "text":
                    # constrains the vertical window
                    # the purpose of this is to essentially limit the growth of the projected box in the y axis because we have an understanding of expected line height.
                    abs_min_y = int(round(row_ix-vertical_tolerance*3,0))
                    abs_max_y = int(round(row_ix+vertical_tolerance*3,0))
                    if pixel_y_min <= abs_min_y:
                        pixel_y_min = abs_min_y
                    if pixel_y_max >= abs_max_y:
                        pixel_y_max = abs_max_y

                # this may be able to be covered within list comprehension but for the sake of explanation i have kept it long form
                # walk the collected extants of the box and if the pixel is 1 add it to the queue to be processed.
                for yy in range(pixel_y_min, pixel_y_max):
                    for xx in range(pixel_x_min, pixel_x_max):
                        try:
                            if image_pixel_simple_matrix[yy][xx] == 1:
                                # seg q is a queue that is defined to house/contain all the pixels in the current window that need to be investigated.
                                # the reason we do this rather than hop from one pixel to the next linearly is because that means we only trace a singular path.
                                # non linear exploration like this means that we investigate each pixel within scope.
                                seg_q.put([xx, yy])
                        except:
                            pass
                
                # until the queue is empty keep enumerating and expanding the box
                while seg_q.qsize() != 0:
                    # get next pixel
                    pixel_x_y = seg_q.get()                    
                    if image_pixel_simple_matrix[pixel_x_y[1]][pixel_x_y[0]] == 1:
                        image_pixel_simple_matrix[pixel_x_y[1]][pixel_x_y[0]] = 2 # means that each pixel can only be processed once.

                        if pixel_x_y[0] <= sx:
                            sx = pixel_x_y[0]
                        if pixel_x_y[0] >= ex:
                            ex = pixel_x_y[0]
                        if pixel_x_y[1] <= sy:
                            sy = pixel_x_y[1]
                        if pixel_x_y[1] >= ey:
                            ey = pixel_x_y[1]

                        pixel_x_min = int(round(pixel_x_y[0]-(horizontal_tolerance*0.5),0))
                        pixel_x_max = pixel_x_y[0]+horizontal_tolerance
                        pixel_y_min = int(round(pixel_x_y[1]-(horizontal_tolerance*0.3),0))
                        pixel_y_max = int(round(pixel_x_y[1]+(horizontal_tolerance*0.3),0))

                        if type == "text":
                            if pixel_y_min <= abs_min_y:
                                pixel_y_min = abs_min_y
                            if pixel_y_max >= abs_max_y:
                                pixel_y_max = abs_max_y

                        for yy in range(pixel_y_min, pixel_y_max):
                            for xx in range(pixel_x_min, pixel_x_max):
                                try:
                                    if image_pixel_simple_matrix[yy][xx] == 1:
                                        seg_q.put([xx, yy])                
                                except:
                                    pass

                # once complete generate your bounding box.
                bounding_boxes.append([sx,sy,ex,ey,ex-sx,ey-sy,(ex-sx)*(ey-sy), str(uuid.uuid4())])
    return bounding_boxes

def order_boxes_on_key(boxes, key, reverse):
    ordered_boxes = []
    key_vals = list(set([b[key] for b in boxes]))
    key_vals.sort()
    if reverse:
        key_vals.reverse()
    for k in key_vals:
        for sb in [b for b in boxes if b[key] == k]:
            ordered_boxes.append(sb)
    return ordered_boxes

# overlap of x, then over lap of y

def bbox_overlap_percent(box1, box2):
    # Each box is (sx, sy, ex, ey)
    x1_min, y1_min, x1_max, y1_max = box1
    x2_min, y2_min, x2_max, y2_max = box2
    
    # Calculate overlap rectangle
    x_overlap = max(0, min(x1_max, x2_max) - max(x1_min, x2_min))
    y_overlap = max(0, min(y1_max, y2_max) - max(y1_min, y2_min))
    intersection = x_overlap * y_overlap
    
    # Areas of each box
    area1 = (x1_max - x1_min) * (y1_max - y1_min)
    area2 = (x2_max - x2_min) * (y2_max - y2_min)
    
    # Union area
    union = area1 + area2 - intersection
    
    # IoU (intersection over union) as percent
    if union == 0:
        return 0.0  # avoid division by zero
    
    percent_overlap = (intersection / union) * 100
    return percent_overlap

def collide_boxes(bounding_boxes):
    print(f"DEBUG - starting collision , {len(bounding_boxes)} boxes")
    # order the boxes on their area, this basically means that we process the larger boxes first
    # the theory here is that the larger boxes are more likely to have more interactions with other boxes.

    ordered_bounding_boxes = order_boxes_on_key(bounding_boxes, 6, True)
    # maintain a list of boxes that have already been included within at least one collision
    # maintain a list of boxes that we wish to remove.
    used_boxes = []
    boxes_to_remove = []
    for ob in ordered_bounding_boxes:
        if ob[7] not in used_boxes:
            used_boxes.append(ob[7])
            # check every other box that is not this box (can be further restricted to 'local' boxes)
            for obb in [bb for bb in ordered_bounding_boxes if bb[7] != ob[7]]:
                if bbox_overlap_percent((ob[0], ob[1], ob[2], ob[3]), (obb[0], obb[1], obb[2], obb[3])) >= 5:
                    used_boxes.append(obb[7])
                    boxes_to_remove.append(obb[7])
                    # update box sx, sy, ex, ey based on the extents of the overlapped box e.g grow box
                    if obb[0] <= ob[0]:
                        ob[0] = obb[0]
                    if obb[1] <= ob[1]:
                        ob[1] = obb[1]
                    if obb[2] >= ob[2]:
                        ob[2] = obb[2]
                    if obb[3] >= ob[3]:
                        ob[3] = obb[3]
                    print(f"DEBUG - box {obb[0:5]} collides with : {ob[0:5]} : {bbox_overlap_percent((ob[0], ob[1], ob[2], ob[3]), (obb[0], obb[1], obb[2], obb[3]))}")

    bounding_boxes = [ob for ob in ordered_bounding_boxes if ob[7] not in boxes_to_remove]
    print(f"DEBUG - ending collision , {len(bounding_boxes)} boxes, {len(boxes_to_remove)}")
    return bounding_boxes, len(boxes_to_remove)

def draw_entities_on_image(image_path, entities, image_name): 
    im_res_iv = Image.open(image_path)
    draw_iv = ImageDraw.Draw(im_res_iv)
    for et in entities:        
        try:
            draw_iv.rectangle([(et[0]), (et[1]), (et[2]), (et[3])], outline="blue")
        except:
            pass

    new_image_path_iv = os.path.join("tmp", image_name)
    del draw_iv
    im_res_iv.save(new_image_path_iv, "PNG")

    return new_image_path_iv
    
image_path = sys.argv[1]
print(f"DEBUG - binning image : {image_path}")
# gain the simplified pixel array for the image and horizontal and vertical tolerances.
image_pixel_simple_matrix, horizontal_tolerance, vertical_tolerance = gain_image_pixel_simple_matrix(image_path, 5)
# im deepcopying here just because i use the unmolested simple pixel array to perform modal analysis.
base_doc_matrix = copy.deepcopy(image_pixel_simple_matrix)

# get the bounding boxes most likely to be 'text'
bounding_boxes_text = gain_bounding_boxes(image_pixel_simple_matrix, horizontal_tolerance, vertical_tolerance, "text")
# get all the bounding boxes without restriction
bounding_boxes = gain_line_segments(base_doc_matrix, horizontal_tolerance, vertical_tolerance)

print(f"DEBUG - gained : {len(bounding_boxes_text)} text bounding boxes")
print(f"DEBUG - gained : {len(bounding_boxes)} text bounding boxes")

# this is purely for debug but this lets you plot the pixels back out over the orignal image.
draw_entities_on_image(image_path, [b for b in bounding_boxes_text if (b[2]-b[0] >= 10) and (b[3]-b[1] >= 10)], f"initial_blind_boxes_{str(uuid.uuid4())}.png")
draw_entities_on_image(image_path, bounding_boxes, f"initial_blind_boxes_raw_{str(uuid.uuid4())}.png")

# the final step (AND THIS NEEDS work) is to regressively collide these boxes.
# technically this could probably be a use for centroid based Lloyd (K-Means) but for our sake something a little more naive will probably work and its cheaper.
# the goal of this is to join boxes based on their interference / overlap, simple regression down to 0 for the number of collisions detected.
bounding_boxes_text, collisions = collide_boxes([b for b in bounding_boxes_text if b[6] >= 1], horizontal_tolerance, vertical_tolerance)
while collisions != 0:
    bounding_boxes_text, collisions = collide_boxes([b for b in bounding_boxes_text if b[6] >= 1], horizontal_tolerance, vertical_tolerance)

draw_entities_on_image(image_path, [b for b in bounding_boxes_text if (b[2]-b[0] >= 10) and (b[3]-b[1] >= 10)], f"collided_blind_boxes_{str(uuid.uuid4())}.png")
