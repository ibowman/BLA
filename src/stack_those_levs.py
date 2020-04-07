#!/usr/bin/env python
from __future__ import print_function
import argparse
from cic_dis import cic_utils
from cic_dis import cic_plot
from cic_dis import cic_overlap
from cic_dis import cic_outspector
import time
import chart_studio.plotly as py
import plotly.graph_objs as go
from collections import defaultdict
import os
import cv2
import cPickle as pickle


def main():
    parser = argparse.ArgumentParser(
        description="Render stacked community connectivity matrix, "
        "same as stack_cmt_mat.py except stack ARA lvls and group regions"
    )
    parser.add_argument('-cmc', '--ctx_mat_csv',
                        help='Input, labeled row source, column destination '
                        'annotation CSV',
                        required=True)
    parser.add_argument('-ccsc', '--char_com_str_csv',
                        help='Path to CSV containing consensus community '
                        'structure ',
                        required=True)
    parser.add_argument('-aoc', '--agg_overlap_csv',
                        help='Input aggregated overlap csv',
                        required=True)
    parser.add_argument('-lvls', '--levels',
                        help='List of levels to visualize',
                        required=True,
                        nargs='+')
    parser.add_argument('-vimgs', '--visual_images',
                        help='List of images for visualizing split ROIs',
                        nargs='+')
    parser.add_argument('-isc', '--injection_site_colors',
                        help='List of injection sites color tuples e.g. {}'
                        .format('-iso BLA_am::228:26:28 BLA_al::255:127:0'),
                        required=True,
                        nargs='+')
    parser.add_argument('-srs', '--split_rois',
                        help="List of ROIs to split into quadrants, can specify mode with e.g. 'ACB:mdlvl', if no mode specified then 'quad' is default",  # noqa
                        nargs='+',
                        default=[])
    parser.add_argument('-fcs', '--focus_roi_lst',
                        help='ROIs to plot exclusively',
                        nargs='+')
    parser.add_argument('-gcs', '--grid_cell_size',
                        help='Grid cell size to use for visual image',
                        default=175,
                        required=False)
    parser.add_argument('-norm', '--norm_mode',
                        help='Indicates that matrix values are normalized',
                        action='store_true')
    parser.add_argument('-v', '--verbose',
                        help='Print relevant but optional output',
                        action='store_true')

    args = parser.parse_args()
    ctx_mat_csv = args.ctx_mat_csv
    char_com_str_csv = args.char_com_str_csv
    agg_overlap_csv = args.agg_overlap_csv
    levels = args.levels
    verbose = args.verbose

    # get, check and parse injecction site color list
    inj_site_clr_lst = args.injection_site_colors
    assert len(inj_site_clr_lst) > 0, "Invalid format for injection_site_colors {}".format(args['injection_site_colors'])  # NOQA
    for inj_site_tup in inj_site_clr_lst:
        assert '::' in inj_site_tup, \
            "Invalid injection site color tup {}".format(inj_site_tup)
        assert ':' in inj_site_tup, \
            "Invalid injection site color tup {}".format(inj_site_tup)
    inj_site_clr_map, inj_site_lst = cic_plot.parse_inj_site_clr_lst(
        inj_site_clr_lst=inj_site_clr_lst)

    split_rois = args.split_rois
    # create dct for roi mode lookup later
    roi_mode_dct = {}
    for roi_mode in split_rois:
        roi_mode_lst = roi_mode.split(':')
        if len(roi_mode_lst) > 1:
            roi_mode_dct[roi_mode_lst[0]] = roi_mode_lst[1]
        else:
            roi_mode_dct[roi_mode_lst[0]] = 'quad'  # quad is default
    focus_roi_lst = args.focus_roi_lst
    visual_images = args.visual_images
    grid_cell_size = int(args.grid_cell_size)
    norm_mode = args.norm_mode

    # what we want to do here:
    #  1) Open all the freakin CSVs
    #  2) Create a DICT from community to roi list of rois with overlap amount
    #     { lvl: { cmt_inj_site : { roi : ovlp } } }
    #  3) Plot DICT as a bar chart or matrix or something

    # 1) opening all the freakin CSVs
    # OPEN CONNECTIVITY MATRIX
    if verbose:
        print("Opening ctx matrix from {}...".format(ctx_mat_csv))
        start = time.time()
    (row_roi_name_npa, col_roi_name_npa, ctx_mat_npa) = \
        cic_utils.read_ctx_mat(ctx_mat_csv)

    # check for retrograde which has many rows going to few cols
    if 'ret' in ctx_mat_csv:
        assert len(row_roi_name_npa) > len(col_roi_name_npa)
        col_roi_name_npa_tmp = col_roi_name_npa
        col_roi_name_npa = row_roi_name_npa
        row_roi_name_npa = col_roi_name_npa_tmp
        ctx_mat_npa = ctx_mat_npa.transpose()
        if verbose:
            print('Detected retrograde in csv name, swapped row and cols')
            print('and transposed matrix')

    # check for capitlization/duplicates
    cic_utils.dup_check_container(dup_check_roi_container=row_roi_name_npa,
                                  input_csv_path=ctx_mat_csv)

    cic_utils.dup_check_container(dup_check_roi_container=col_roi_name_npa,
                                  input_csv_path=ctx_mat_csv)
    if verbose:
        print("opened ctx matrix in {:.04}s".format(time.time() - start))
        print("(num rows, num cols) {}".format(ctx_mat_npa.shape))
        print("top left val {}".format(ctx_mat_npa[0][0]))

    # OPEN CHAR COMMUNITY STRUCTURE
    if verbose:
        print("Opening char com str csv {}...".format(char_com_str_csv))
        start = time.time()
    cons_cmt_str = cic_plot.cons_cmt_str(char_cmt_csv_path=char_com_str_csv)

    if verbose:
        print("opened char com str csv in {:.04}s".format(time.time() - start))
        print("num communities {}".format(len(cons_cmt_str)))

    # OPEN AGG OVERLAP CSV
    if verbose:
        print("Opening aggregated overlap csv {}...".format(agg_overlap_csv))
        start = time.time()
    (agg_overlap_csv_header, agg_overlap_rows) = \
        cic_overlap.read_agg_overlap_csv(input_csv_path=agg_overlap_csv)
    if verbose:
        print("opened agg overlap csv in {:.04}s".format(time.time() - start))
        print("num rows {}".format(len(agg_overlap_rows)))

    # 2) Creating a DICT {cmt : frozenset([roi])}
    all_cmt_inj_sites_set = frozenset()
    if verbose:
        print("Creating community -> roi dictionary for levels {}...".
              format(levels))
        start = time.time()

    # first populate split roi dct
    lvl_split_roi_dct = {}  # { lvl: { roi: { min_col: ..., max_col,... } } }
    for idx, grid_tup_str in enumerate(col_roi_name_npa):
        if cic_plot.in_lvl_lst(grid_tup_str=grid_tup_str,
                               lvl_lst=levels):
            (lvl, hemi, col_str, row_str) = \
                cic_plot.grid_tup_str_to_lvl_hemi_col_row(grid_tup_str)
            col = int(col_str)
            row = int(row_str)
            # populate split_roi_dct either default, of from lvl
            # get roi and check if it is in the split_rois list
            unsplit_roi_str = cic_plot.col_val(
                col_hdr_str='REGION(S)',
                agg_overlap_csv_header=agg_overlap_csv_header,
                agg_overlap_csv_rows=agg_overlap_rows,
                grid_tup_str=grid_tup_str)

            # now we need to populate split dict

            if lvl not in lvl_split_roi_dct:
                # { roi: { min_col: ..., max_col,... } }
                split_roi_dct = {}  # create empty dct on every lvl by default
                if verbose:
                    print("creating split roi dct vals, level {}...".format(
                        lvl))
            else:
                split_roi_dct = lvl_split_roi_dct[lvl]

            roi_str = cic_plot.in_roi_lst(
                unsplit_roi_str=unsplit_roi_str,
                roi_str_lst=[x.split(':')[0] for x in split_rois])

            if (roi_str is not None and
                    (focus_roi_lst is None or cic_plot.in_focus_roi_lst(
                        roi_str=roi_str,
                        focus_roi_lst=focus_roi_lst))):
                if roi_str not in split_roi_dct:
                    roi_dct = {'min_col': col,
                               'max_col': col,
                               'min_row': row,
                               'max_row': row,
                               'mode': roi_mode_dct[roi_str]}
                    if verbose:
                        print("populating {} split {} vals".format(
                            roi_str, roi_mode_dct[roi_str]))
                else:
                    roi_dct = split_roi_dct[roi_str]

                # update split_roi_dct
                roi_dct['min_col'] = min(roi_dct['min_col'], col)
                roi_dct['max_col'] = max(roi_dct['max_col'], col)
                roi_dct['min_row'] = min(roi_dct['min_row'], row)
                roi_dct['max_row'] = max(roi_dct['max_row'], row)
                split_roi_dct[roi_str] = roi_dct
            lvl_split_roi_dct[lvl] = split_roi_dct

    # now check if visual images list specified, and if so then write out imgs
    if visual_images is not None:
        assert len(visual_images) == len(levels)  # must be one image per lvl
        for lvl_idx, lvl in enumerate(levels):
            img_path = visual_images[lvl_idx]
            assert lvl in img_path, "level {} not in {}... \nvisual image order {} does not match order level order {}".format(lvl, img_path, visual_images, levels)  # noqa
            assert os.path.isfile(img_path)
            if verbose:
                start = time.time()
                print("Making split roi image with {} and gcs {}...".format(
                    img_path, grid_cell_size))
            split_roi_atlas_img = cic_outspector.make_split_roi_atlas_img(
                ann_atlas_tif_path=img_path,
                gcs=grid_cell_size,
                lvl=lvl,
                split_roi_dct=lvl_split_roi_dct[lvl],
                hemi='r',
                agg_overlap_csv_header=agg_overlap_csv_header,
                agg_overlap_csv_rows=agg_overlap_rows,
                verbose=verbose)

            split_roi_atlas_img_path = \
                cic_outspector.split_roi_atlas_img_path(
                    out_dir_path=os.path.dirname(img_path),
                    atlas_tif_path=img_path,
                    split_rois=split_rois,
                    focus_roi_lst=focus_roi_lst)

            cv2.imwrite(split_roi_atlas_img_path, split_roi_atlas_img)
            if verbose:
                print("done making split roi image in {:.04}s".format(
                    time.time() - start))

    # walk through each non-zero entry of ctx mat csv, get roi and bin by cmt
    #  cmt_inj_site_lvl_dct  { cmt_inj_site { lvl { roi : ovlp } } }
    cmt_inj_site_lvl_dct = defaultdict(dict)
    for idx, grid_tup_str in enumerate(col_roi_name_npa):
        for row_inj_site in row_roi_name_npa:
            # lvl_roi_ovlp_dct  { lvl { roi : ovlp } }
            lvl_roi_ovlp_dct = cmt_inj_site_lvl_dct[row_inj_site]
            if cic_plot.in_lvl_lst(grid_tup_str=grid_tup_str,
                                   lvl_lst=levels):
                lvl = cic_plot.grid_tup_str_to_lvl_hemi_col_row(
                    grid_tup_str)[0]
                # roi_ovlp_dct { roi : ovlp }
                roi_ovlp_dct = lvl_roi_ovlp_dct.get(lvl, {})
                cmt_inj_site_lst = cic_plot.cmt_inj_site_lst_from_grid_tup_str(
                    inj_site_lst=inj_site_lst,
                    cons_cmt_str=cons_cmt_str,
                    grid_tup_str=grid_tup_str)
                assert len(cmt_inj_site_lst) == 1, print(cmt_inj_site_lst)
                all_cmt_inj_sites_set = all_cmt_inj_sites_set.union(
                    cmt_inj_site_lst)
                cmt_inj_site = cmt_inj_site_lst[0]
                if cmt_inj_site == row_inj_site:
                    unsplit_roi_str = cic_plot.col_val(
                        col_hdr_str='REGION(S)',
                        agg_overlap_csv_header=agg_overlap_csv_header,
                        agg_overlap_csv_rows=agg_overlap_rows,
                        grid_tup_str=grid_tup_str)
                    assert unsplit_roi_str is not None
                    if focus_roi_lst is None or cic_plot.in_roi_lst(
                            unsplit_roi_str,
                            roi_str_lst=focus_roi_lst):
                        if verbose:
                            if lvl not in lvl_roi_ovlp_dct:
                                print(
                                    "populating {} lvl {} roi overlap values".
                                    format(row_inj_site, lvl))

                        # split rois if they are in lvl_split_roi_dct
                        roi_str = cic_plot.split_roi(
                            unsplit_roi_str=unsplit_roi_str,
                            grid_tup_str=grid_tup_str,
                            split_roi_dct=lvl_split_roi_dct[lvl])
                        # update overlap val in roi overlap dict for this level
                        row_idx = row_roi_name_npa.tolist().index(row_inj_site)

                        # do sanity checks before setting roi_ovlp_dct val
                        grid_only_str = cic_plot.col_val(
                            col_hdr_str='GRID ONLY',
                            agg_overlap_csv_header=agg_overlap_csv_header,
                            agg_overlap_csv_rows=agg_overlap_rows,
                            grid_tup_str=grid_tup_str)
                        overlap_str = cic_plot.col_val(
                            col_hdr_str='OVERLAP',
                            agg_overlap_csv_header=agg_overlap_csv_header,
                            agg_overlap_csv_rows=agg_overlap_rows,
                            grid_tup_str=grid_tup_str)
                        grid_only = float(grid_only_str)
                        overlap = int(overlap_str)
                        fraction_overlap = overlap / (grid_only + overlap)
                        assert fraction_overlap <= 1.0
                        assert ctx_mat_npa[row_idx][idx] <= 1.0
                        # use regular overlap pixel count since more intuitive
                        dct_overlap = overlap

                        roi_ovlp_dct[roi_str] = \
                            roi_ovlp_dct.get(roi_str, 0) + dct_overlap

                        # point to all the potentially new dict vals
                        lvl_roi_ovlp_dct[lvl] = roi_ovlp_dct

    # VISUALIZATION CODE BEGIN
    if verbose:
        print("Finished creating dict in in {:.04}s".format(
            time.time() - start))

    # make new all_cmt_inj_sites that matches order of inj site color list
    assert len(all_cmt_inj_sites_set) == len(inj_site_lst), "{} != {}".format(
        sorted(all_cmt_inj_sites_set), sorted(inj_site_lst))
    all_cmt_inj_sites = inj_site_lst

    # 3) Plot DICT as a bar chart or matrix or something
    if verbose:
        print("Plotting a bar chart or matrix or something for {} inj_sites".
              format(len(cmt_inj_site_lvl_dct)))
        start = time.time()
        print("injection site order is {}".format(all_cmt_inj_sites))

    # create top roi lst
    num_top = 2
    if verbose:
        print("creating roi list of top {} rois ...".format(num_top))
    # create { lvl : [ "cmt_a top roi str", ... ] } } list of most proj. roi
    lvl_top_roi_lst_dct = defaultdict(list)
    for inj_site in all_cmt_inj_sites:
        lvl_roi_ovlp_dct = cmt_inj_site_lvl_dct[inj_site]
        for lvl in sorted(args.levels):
            if lvl in lvl_roi_ovlp_dct:  # level not in every community
                roi_ovlp_dct = lvl_roi_ovlp_dct[lvl]

                # lst of tups will be sorted by overlap value
                lst = sorted(roi_ovlp_dct.iteritems(),
                             key=lambda (k, v): (v, k), reverse=True)
                roi_str = "len lst {} ".format(len(lst))
                roi_str = "({}) ".format(lst[0][0]).replace("|", " ")
                for tup in lst[1:num_top]:
                    roi_str += " ({})".format(tup[0]).replace("|", " ")
                    # create list of strings
                    lvl_top_roi_lst_dct[lvl].append(roi_str)
                if verbose:
                    print("{} {} {}".format(inj_site, lvl,
                                            lvl_top_roi_lst_dct[lvl]))
            else:
                lvl_top_roi_lst_dct[lvl].append('')

    if verbose:
        for inj_site in all_cmt_inj_sites:
            lvl_roi_ovlp_dct = cmt_inj_site_lvl_dct[inj_site]
            print("{}".format(inj_site))
            for lvl in sorted(args.levels):
                roi_ovlp_dct = lvl_roi_ovlp_dct[lvl]
                print("-- lvl {} roi count {} total ovlp {} {}".
                      format(lvl, len(roi_ovlp_dct),
                             sum([roi_ovlp_dct[roi] for roi in roi_ovlp_dct]),
                             lvl_top_roi_lst_dct[lvl]))

    if verbose:
        print("making traces for each level...")
    # create traces, or levels that will be plotted as bars
    # first create the values
    # { lvl : [ cmt trace val, ...] }
    lvl_trace_y_dct = defaultdict(list)
    # populate total overlap per cmt for normalization
    cmt_trace_total_dct = defaultdict(float)
    for cmt_inj_site in all_cmt_inj_sites:
        lvl_roi_ovlp_dct = cmt_inj_site_lvl_dct[cmt_inj_site]

        # sort trace values to match across bars
        for lvl in sorted(args.levels):
            if lvl in lvl_roi_ovlp_dct:  # level not in every community
                roi_ovlp_dct = lvl_roi_ovlp_dct[lvl]
                # append num rois to lvl_trace_y_dct
                lvl_trace_y_dct[lvl].append(
                    sum([float(roi_ovlp_dct[roi]) for roi in roi_ovlp_dct]))
                cmt_trace_total_dct[cmt_inj_site] = (
                    cmt_trace_total_dct[cmt_inj_site] + 
                    sum([float(roi_ovlp_dct[roi]) for roi in roi_ovlp_dct]))
            else:
                # append 0 amount of trace
                lvl_trace_y_dct[lvl].append(0)

    if verbose:
        print("normalizing traces to 100% ...")
    # normalize traces to 100 %
    # get community from the first level
    # need to sort this too since cmt_trace_y_dct points to a list
    for cmt_idx, cmt_inj_site in enumerate(all_cmt_inj_sites):
        lvl_roi_ovlp_dct = cmt_inj_site_lvl_dct[cmt_inj_site]
        for lvl in sorted(args.levels):
            if lvl in lvl_roi_ovlp_dct:  # level not in every community
                # normalize to 100% using level values
                lvl_trace_y_dct[lvl][cmt_idx] = \
                    (float(lvl_trace_y_dct[lvl][cmt_idx]) /
                     cmt_trace_total_dct[cmt_inj_site]) * 100.0

    if verbose:
        print("populating graph object and writing file...")
    # now populate graph object
    # not all lvels have all communities, so use global max to cover all
    traces = []

    for lvl in sorted(args.levels):
        traces.append(
            go.Bar(
                # don't forget to sort here bro
                x=["{}".format(x) for x in all_cmt_inj_sites],
                y=lvl_trace_y_dct[lvl],  # contains [ cmt trace val, ... ]
                text=lvl_top_roi_lst_dct[lvl],
                textposition='auto',
                textfont=dict(
                    color='rgb(255, 255, 255)'
                ),
                #                marker=dict(
                #                    color='rgb{}'.format(inj_site_clr_map[lvl])
                #                ),
                name="lvl {}".format(lvl))
        )

    data = traces
    layout = go.Layout(
        barmode='stack',
        width=1500,
        height=1500
    )

    fig = go.Figure(data=data, layout=layout)
    if norm_mode:
        base_agg_overlap_csv = 'norm_stacked-' + os.path.basename(
            agg_overlap_csv)
    else:
        base_agg_overlap_csv = 'stcklvl-' + os.path.basename(
            agg_overlap_csv)
    out_img_path = base_agg_overlap_csv.replace(
        '.csv', '_lvls-{}.png'.format(levels))

    if focus_roi_lst is not None:
        out_img_path = out_img_path.replace(
            '.png',
            '_fcs-{}.png'.format(
                [roi_name[0] for roi_name in focus_roi_lst]))
    out_img_path = out_img_path.replace(",", "")
    out_img_path = out_img_path.replace("'", "")
    out_img_path = out_img_path.replace("[", "")
    out_img_path = out_img_path.replace("]", "")
    out_img_path = out_img_path.replace(" ", "")
    py.image.save_as(fig, out_img_path)
    if verbose:
        print("Finished plotting {} in {:.04}s".format(
            out_img_path,
            time.time() - start))

    output_pickle_path = cic_utils.pickle_path(out_img_path)
    pickle.dump(args, open(output_pickle_path, "wb"))
    if verbose:
        print("Wrote pickle args to {}".format(output_pickle_path))


if __name__ == '__main__':
    main()
