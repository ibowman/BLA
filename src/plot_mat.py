#!/usr/bin/env python
from __future__ import print_function
import matplotlib
matplotlib.use('Agg')  # noqa
import matplotlib.pyplot as plt
import argparse
import os
import numpy as np
from cic_dis import cic_utils
from cic_dis import cic_plot
import ast
import sys
import bct
import cPickle as pickle
import pydot
import math


def main():
    parser = argparse.ArgumentParser(
        description="Order nodes of matrix by modules defined in "
        "characteristic community structure csv, and plot.")
    parser.add_argument('-c', '--char_cmt_str_csv',
                        help='Characteristic community structure CSV')

    parser.add_argument('-m', '--matrix_csv',
                        help='Row/column based connectivity matrix',
                        required=True)

    parser.add_argument('-v', '--verbose',
                        help='Print matrix, cluster related output',
                        action='store_true')

    parser.add_argument('-p', '--partition',
                        help="Double quote delimited partition e.g. "
                        "\"[['A', 'C'], ['B']\"")

    parser.add_argument('-dl', '--dictionary_list',
                        help="Double quoted nested list of modules, with "
                        "subnetworks defined at the top level e.g. "
                        "\"[{'dorsal' : ['dorsal_DG', 'dorsal_CA3', "
                        "'dorsal_CA1'], 'ventral' : ['ventral_DG', "
                        "'ventral_CA3', 'ventral_CA1']}, {'dorsal_DG' : "
                        "[['DGd', 'DGi', 'DGpod']], 'dorsal_CA3' : "
                        "[['CA3d', 'CA3dd'], ['CA3ic', 'CA3id']], 'dorsal_CA1'"
                        ": [['CA1dc', 'CA2', 'ProSUB'], ['CA1dr', 'SUBdd']], "
                        "'ventral_DG': [['DGpov', 'DGv']], 'ventral_CA1': "
                        "[['CA1i', 'CA1v', 'CA3v', 'SUBdv', 'SUBv']], "
                        "'ventral_CA3': [['CA1vv', 'CA3vv', 'SUBvv']]},]\"")

    parser.add_argument('-pt', '--plot_type',
                        help='Type of plot to generate: modularized matrix, '
                        'sparse matrix, schematic graph, or all',
                        choices=['mod_mat', 'schematic', 'exp_prof',
                                 'all'],
                        default='all')

    parser.add_argument('-pts', '--pad_to_square',
                        help='Pad to square rather than convert',
                        action='store_true')

    parser.add_argument('-ps', '--pretend_square',
                        help='Pretend square',
                        action='store_true')

    parser.add_argument('-fs', '--fontsize',
                        help='Font size of tick labels in x and y',
                        type=int,
                        default=2)

    parser.add_argument('-fmt', '--format',
                        help='Format of plot (e.g. svg, png,...)',
                        default='svg')

    parser.add_argument('-cmt_ln_wt', '--community_line_weight',
                        help='Weight of edges between intra community nodes',
                        type=float,
                        default=2.0)

    parser.add_argument('-mod_ln_wt', '--module_line_weight',
                        help='Weight of edges between intra module nodes',
                        type=float,
                        default=0.5)

    parser.add_argument('-inter_ln_wt', '--inter_line_weight',
                        help='Weight of edges between non intra nodes',
                        type=float,
                        default=0.1)

    parser.add_argument('-rd', '--rankdir',
                        help='Rank direction of graph e.g. TB, RL',
                        default='TB')
    parser.add_argument('-g', '--gamma', help='Gamma value use in filename',
                        type=float, default=1.0)
    parser.add_argument('-nr', '--num_runs',
                        help='Num runs to use in filename',
                        type=int, default=100)
    parser.add_argument('-mlw', '--mat_line_width',
                        type=float, default=1.5)

    parser.add_argument(
        '-rvcs', '--row_value_community_sort',
        help='Sort communities by average value of member rows',
        action='store_true')

    parser.add_argument(
        '-abc', '--alpha_community_rows',
        help='Sort community rows rois alphabetically',
        action='store_true')

    parser.add_argument('-dsn', '--draw_subnetwork',
                        help='Draw only the subnetwork specified')

    args = vars(parser.parse_args())

    if args['char_cmt_str_csv'] is None and \
       args['partition'] is None:
        parser.error("You need to specify either a --char_cmt_str_csv or a --partition")  # noqa

    matrix_csv_path = args['matrix_csv']
    assert os.path.isfile(matrix_csv_path),\
        "can't find row column matrix csv file {}".format(matrix_csv_path)

    # PARSE ARGMENTS
    verbose = args['verbose']

    gamma = args['gamma']
    num_runs = args['num_runs']

    mlw = args['mat_line_width']

    pad_to_sq = args['pad_to_square']
    pretend_sq = args['pretend_square']

    # determine plot type
    plot_type = args['plot_type']

    fontsize = args['fontsize']

    fmt = "{}".format(args['format'])

    cmt_ln_wt = args['community_line_weight']
    mod_ln_wt = args['module_line_weight']
    inter_ln_wt = args['inter_line_weight']

    rankdir = args['rankdir']

    row_value_community_sort = args['row_value_community_sort']

    alphabetize_rows = args['alpha_community_rows']

    draw_subnetwork = args['draw_subnetwork']

    if plot_type == 'all':
        mod_mat = schematic = True

    else:
        mod_mat = plot_type == 'mod_mat'
        schematic = plot_type == 'schematic'

    partition = args['partition']

    dictionary_list = args['dictionary_list']
    submod_roi_name_lst_dct = []
    if dictionary_list:
        submod_roi_name_lst_dct = ast.literal_eval(dictionary_list)

    char_cmt_str_csv_path = args['char_cmt_str_csv']
    if partition is not None:
        # convert get rid of sets and frozensets and just read in list
        partition = partition.replace('frozenset([', '[').replace('])', ']')
        cmt_str_lst = ast.literal_eval(partition)

    if char_cmt_str_csv_path:
        assert os.path.isfile(char_cmt_str_csv_path),\
            "can't find char cmt str csv file {}".format(char_cmt_str_csv_path)

        cmt_str_lst = cic_plot.cons_cmt_str(char_cmt_str_csv_path)

        assert len(cmt_str_lst) > 0
        if verbose:
            print("read {}\nwith value\n{}".
                  format(type(cmt_str_lst),
                         cic_utils.rm_startstr_sublist('NONE', cmt_str_lst)))

    # OPEN, READ INPUT CSV, PAD TO SQUARE IF NECESSARY
    (row_roi_name_npa, col_roi_name_npa, ctx_mat_npa) = \
        cic_utils.read_ctx_mat(matrix_csv_path)

    if verbose:
        print("Read row ROIs: \n{}\n".format(row_roi_name_npa))
        print("Read col ROIs: \n{}\n".format(col_roi_name_npa))

    # double check formatting, shape of array
    if not cic_utils.is_sq(row_roi_name_npa=row_roi_name_npa,
                           col_roi_name_npa=col_roi_name_npa,
                           ctx_mat_npa=ctx_mat_npa):
        if verbose:
            shape = ctx_mat_npa.shape
            print("Matrix in {} is not square\nrows {}, columns {}".
                  format(matrix_csv_path, shape[0], shape[1]))
            print("and/or row ROIs \n{}\n don't match col ROIs\n {}".format(
                row_roi_name_npa, col_roi_name_npa))

        if pad_to_sq:
            if verbose:
                print("padding to square...")
            (pad_row_roi_name_npa, pad_col_roi_name_npa, sq_ctx_mat_npa) = \
                cic_utils.pad_rect_ctx_mat_to_sq(row_roi_name_npa,
                                                 col_roi_name_npa,
                                                 ctx_mat_npa)
            row_roi_name_npa = pad_row_roi_name_npa
            col_roi_name_npa = pad_col_roi_name_npa

        elif pretend_sq:  # it turns out being a square matrix not necessary
            if verbose:
                print("pretending square...")
            sq_ctx_mat_npa = ctx_mat_npa

        # else conv to square
        else:
            if verbose:
                print("converting to square...")
            (sq_roi_name_npa, sq_ctx_mat_npa) = \
                cic_utils.conv_rect_ctx_mat_to_sq(row_roi_name_npa,
                                                  col_roi_name_npa,
                                                  ctx_mat_npa)
            row_roi_name_npa = sq_roi_name_npa
            col_roi_name_npa = sq_roi_name_npa

        if verbose:
            shape = sq_ctx_mat_npa.shape
            print("with rows/cols {}/{} done".format(shape[0], shape[1]))

        ctx_mat_npa = sq_ctx_mat_npa

    # row_roi_name_npa and col_roi_name_npa as cmt_str_lst should be defined
    # build initial row_ci and col_ci
    flat_cmt_str_lst = cic_utils.flatten(cmt_str_lst)
    row_ci = np.array(cic_plot.build_ci_lst(roi_name_npa=row_roi_name_npa,
                                            cmt_str_lst=cmt_str_lst))
    col_ci = np.array(cic_plot.build_ci_lst(roi_name_npa=col_roi_name_npa,
                                            cmt_str_lst=cmt_str_lst))

    if mod_mat:
        # NOW, USE COL ROIS FOR COMMUNITY REORDER
        # don't check ROIs if pad to sq
        if not (pad_to_sq or pretend_sq):
            if verbose:
                print("Checking ROIs...")
            # check rois
            roi_set1 = frozenset(flat_cmt_str_lst)
            roi_set2 = frozenset(col_roi_name_npa)

            not_in2 = roi_set1 - roi_set2
            not_in1 = roi_set2 - roi_set1

            if(len(not_in1) > 0 or len(not_in2) > 0):
                if args['char_cmt_str_csv'] is None:
                    char_cmt_str_csv_path = "user provided community partition"
                    print(
                        "ERROR: following ROIs are in {} but not in {}: {}".
                        format(
                            char_cmt_str_csv_path, matrix_csv_path, not_in2))

                    print(
                        "ERROR: Following ROIs are in {} but not in {}: {}".
                        format(
                          matrix_csv_path, char_cmt_str_csv_path, not_in1))

                    assert 0, \
                        "Can not reorder by module. Fix ERRORs above, re-run."

        # PLOT REORDERED MATRIX
        #  reorder matrix using cmt str rois and then plot that
        #  make community indices

        # set row cmt str list differently if row value sort option selected
        row_cmt_str_lst = cmt_str_lst
        col_cmt_str_lst = cmt_str_lst

        # if option specificied, sort both rows and cols... starting with cols
        if row_value_community_sort:
            # first sort cols and create temp sorted matrix for row sort
            col_sorted_cmt_str_lst = cic_plot.sort_cmt_str_lst_by_mat_col_vals(
                roi_name_npa=col_roi_name_npa,
                ctx_mat_npa=ctx_mat_npa,
                cmt_str_lst=cmt_str_lst)
            col_cmt_str_lst = col_sorted_cmt_str_lst

            # all vars are temporary for row sort
            temp_new_col_roi_indices = cic_plot.new_indices(
                roi_name_npa=col_roi_name_npa,
                sorted_cmt_str_lst=col_cmt_str_lst)
            temp_reorder_ctx_mat_npa_cols = \
                ctx_mat_npa[:, temp_new_col_roi_indices]

            row_sorted_cmt_str_lst = cic_plot.sort_cmt_str_lst_by_mat_row_vals(
                roi_name_npa=row_roi_name_npa,
                ctx_mat_npa=temp_reorder_ctx_mat_npa_cols,
                cmt_str_lst=cmt_str_lst)
            row_cmt_str_lst = row_sorted_cmt_str_lst

        row_ci = np.array(cic_plot.build_ci_lst(roi_name_npa=row_roi_name_npa,
                                                cmt_str_lst=row_cmt_str_lst))
        col_ci = np.array(cic_plot.build_ci_lst(roi_name_npa=col_roi_name_npa,
                                                cmt_str_lst=col_cmt_str_lst))

        # calculated indices differently based on row sorted cmt option val
        if row_value_community_sort:
            y_bounds = bct.grid_communities(row_ci)[0]
            new_row_roi_indices = cic_plot.new_indices(
                roi_name_npa=row_roi_name_npa,
                sorted_cmt_str_lst=row_cmt_str_lst)
            x_bounds = bct.grid_communities(col_ci)[0]
            new_col_roi_indices = cic_plot.new_indices(
                roi_name_npa=col_roi_name_npa,
                sorted_cmt_str_lst=col_cmt_str_lst)

        else:
            (y_bounds, new_row_roi_indices) = bct.grid_communities(row_ci)
            (x_bounds, new_col_roi_indices) = bct.grid_communities(col_ci)

        if alphabetize_rows:
            new_row_roi_indices = \
                cic_plot.alpha_row_roi_indices(
                    row_roi_name_npa, new_row_roi_indices, row_ci)

        new_col_roi_name_npa = np.array(col_roi_name_npa[new_col_roi_indices])
        new_row_roi_name_npa = np.array(row_roi_name_npa[new_row_roi_indices])
        temp_rs = ctx_mat_npa[new_row_roi_indices, :]
        reorder_ctx_mat_npa = temp_rs[:, new_col_roi_indices]

    # if plot expression profile mode, then use reordered ctx mat, as well as
    # y bounds to calculate expression profiles a) total b) per cluster, then
    # create expression profile taking ration of expression per gene vs total
    if plot_type == 'exp_prof':
        if verbose:
            print('y bounds {} rows {}'.format(
                y_bounds, reorder_ctx_mat_npa.shape[0]))

        # get col totals
        mat_col_total = cic_plot.col_totals(reorder_ctx_mat_npa)
        if verbose:
            print('mat col total {}'.format(mat_col_total))

        lower_bound = -1
        submat_col_totals = []
        row_roi_name_col_totals = []
        for idx, y_bound in enumerate(y_bounds):
            upper_bound = int(math.ceil(y_bound))
            if upper_bound > 0:
                assert(lower_bound >= 0)
                submat_col_totals.append(
                    cic_plot.col_totals(
                        reorder_ctx_mat_npa[lower_bound:upper_bound, :]))
                row_roi_name_col_totals.append(
                    new_row_roi_name_npa[lower_bound:upper_bound])
            lower_bound = upper_bound
        submat_col_totals_ratio = []

        for idx, submat_col_total in enumerate(submat_col_totals):
            submat_col_totals_ratio.append(
                np.array(submat_col_totals[idx])/np.array(mat_col_total))

        for idx, genes in enumerate(row_roi_name_col_totals):
            if verbose:
                print('genes in cluster {}: {}'.format(idx+1, genes))
            plot_exp_prof(col_totals=submat_col_totals_ratio[idx],
                          legend='Expression Profile Cluster {}'.format(idx+1),
                          legend_loc='upper left',
                          color='blue',
                          out_path=matrix_csv_path.replace(
                              '.csv',
                              '_{}-{}_runs-{}_clust-{}_exp_prof.{}'.
                              format(gamma,
                                     len(np.unique(row_ci)),
                                     num_runs, idx+1, fmt)),
                          x_roi_name_npa=new_col_roi_name_npa,
                          fontsize=fontsize,
                          args=args)

    if mod_mat:
        if verbose:
            print("plotting reordered matrix... ", end='')
            sys.stdout.flush()

        if verbose:
            shape = reorder_ctx_mat_npa.shape
            print("reorder pre rm none ctx rows/cols {}/{} done".format(
                shape[0], shape[1]))

        (plot_new_col_roi_name_npa, plot_ctx_mat_npa) =\
            rm_NONE_cols(col_roi_name_npa=new_col_roi_name_npa,
                         ctx_mat_npa=reorder_ctx_mat_npa)

        if verbose:
            shape = plot_ctx_mat_npa.shape
            print("reorder post rm none ctx rows/cols {}/{} done".format(
                shape[0], shape[1]))

        plot_matrix(connectivity_matrix=plot_ctx_mat_npa,
                    cmap=plt.cm.gray_r,
                    out_path=matrix_csv_path.replace(
                        '.csv',
                        '_gamma-{}_{}-runs_mod_reorder.{}'.
                        format(gamma, num_runs, fmt)),
                    x_roi_name_npa=plot_new_col_roi_name_npa,
                    y_roi_name_npa=new_row_roi_name_npa,
                    fontsize=fontsize,
                    args=args,
                    x_bounds=x_bounds,
                    y_bounds=y_bounds,
                    pad_to_sq=pad_to_sq,
                    linewidth=mlw)

        if verbose:
            print("done")

    if schematic:
        if verbose:
            print("Creating directed graph from row col matrix...")
        dg = cic_plot.dir_graph_from_ctx_mat(
            row_roi_name_npa=row_roi_name_npa,
            col_roi_name_npa=col_roi_name_npa,
            cmt_str_lst=cmt_str_lst,
            row_ci=row_ci,
            col_ci=col_ci,
            ctx_mat_npa=ctx_mat_npa,
            submod_roi_name_lst_dct=submod_roi_name_lst_dct,
            cmt_ln_wt=cmt_ln_wt,
            mod_ln_wt=mod_ln_wt,
            inter_ln_wt=inter_ln_wt,
            rankdir=rankdir,
            verbose=verbose)

        if verbose:
            print("done")
        if verbose:
            print("Drawing...")
        sch_path = matrix_csv_path. \
            replace(
                '.csv', '_{}_runs-{}_schematic.{}'.
                format(gamma, num_runs, fmt))

        draw_graph = None
        if draw_subnetwork:
            sch_path = matrix_csv_path.\
                replace('.csv', '_sn_{}-{}_runs-{}_schematic.{}'.
                        format(
                            draw_subnetwork, gamma, num_runs,
                            fmt))
            if verbose:
                print("Will draw {} subnetwork exclusively".
                      format(draw_subnetwork))
            cluster_name = 'cluster_{}'.format(draw_subnetwork)
            for subgraph in dg.get_subgraphs():
                if subgraph.get_name() == cluster_name:
                    draw_graph = pydot.Dot(graph_type='digraph',
                                           rankdir=rankdir)
                    draw_graph.add_subgraph(subgraph)

            if not draw_graph:
                print("Uh-oh, subnetwork {} not found in {}".
                      format(
                          draw_subnetwork,
                          [sg.get_name().replace('cluster_', '')
                              for sg in dg.get_subgraphs()]))
                exit(1)
        else:
            draw_graph = dg

        draw_graph.write(sch_path, prog='dot', format=fmt)
        pickle_path = cic_utils.pickle_path(sch_path)
        pickle.dump(args, open(pickle_path, "wb"))

        if verbose:
            print("done")


# plot expression profile
def plot_exp_prof(col_totals, legend, legend_loc, color, out_path,
                  x_roi_name_npa, fontsize, args):

    p, = plt.plot(col_totals, color=color)
    plt.legend([p], [legend], legend_loc)
    plt.xticks(np.arange(0, len(x_roi_name_npa)), x_roi_name_npa,
               rotation='vertical')
    yticks = np.arange(0, 1.1, 0.1)
    plt.yticks(yticks)

    plt.savefig(out_path, dpi=600, bbox_inches='tight')
    plt.close()
    pickle_path = cic_utils.pickle_path(out_path)
    pickle.dump(args, open(pickle_path, "wb"))


def plot_matrix(connectivity_matrix, cmap, out_path, x_roi_name_npa,
                y_roi_name_npa, fontsize, args, x_bounds=[], y_bounds=[],
                pad_to_sq=False, linewidth=1.5):
    plt.matshow(connectivity_matrix, cmap=cmap)
    plt.xticks(range(len(x_roi_name_npa)), x_roi_name_npa,
               rotation='vertical',
               fontsize=fontsize)
    plt.yticks(range(len(y_roi_name_npa)), y_roi_name_npa,
               fontsize=fontsize)
    for x_b in x_bounds:
        if not pad_to_sq:
            plt.axvline(x=x_b, color='blue', linewidth=linewidth)
    for y_b in y_bounds:
        plt.axhline(y=y_b, color='blue', linewidth=linewidth)

    plt.savefig(out_path, dpi=600)
    pickle_path = cic_utils.pickle_path(out_path)
    pickle.dump(args, open(pickle_path, "wb"))


# returns tuple:
#  (plot_col_roi_name_npa, plot_ctx_mat_npa)
def rm_NONE_cols(col_roi_name_npa, ctx_mat_npa):
    # remove all columns with NONE in ROI name
    rm_index_arr = []
    for col_index, col_roi in enumerate(col_roi_name_npa):
        if col_roi.startswith('NONE'):
            rm_index_arr.append(col_index)

    plot_col_roi_name_npa = np.delete(col_roi_name_npa, rm_index_arr)
    plot_ctx_mat_npa = np.delete(ctx_mat_npa, rm_index_arr, 1)  # 1 -> y axis
    return (plot_col_roi_name_npa, plot_ctx_mat_npa)


if __name__ == '__main__':
    main()
