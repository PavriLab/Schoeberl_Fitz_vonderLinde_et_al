# setting number threads environment variable
import os

from numpy.lib.function_base import flip
os.environ['NUMEXPR_MAX_THREADS'] = '8'

import numpy as np
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import matplotlib.gridspec as gs
import pandas as pd
import pyBigWig as pbw
import logging
import argparse as ap
import regex as re

# setting rcParams to enable editable text in Adobe Illustrator
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.size'] = 20
plt.rcParams['lines.linewidth'] = 1.25
def latex(bool):
    if bool:
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = 'Helvetica'
        plt.rcParams['text.usetex'] = True
        plt.rcParams['text.latex.preamble'] = "\n".join([r'\usepackage[Symbol]{upgreek}', r'\usepackage{helvet}', r'\renewcommand{\familydefault}{\sfdefault}'])
    else: 
        plt.rcParams['text.usetex'] = False
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = 'Arial'

def get_regional_matrix(contactMatrix, intervalstarts, start, end, csr=False):
    indices, include = [], []
    for i, loc in enumerate(intervalstarts):
        if loc >= start and loc < end:
            indices.append(i)
            include.append(loc)

    contactMatrix = contactMatrix[indices, :][:, indices]

    return contactMatrix if not csr else contactMatrix.toarray()


def load_profile_table(filename):
    tab = pd.read_csv(filename, sep='\t')
    header = [s.strip("'") for s in open(filename, 'r').readline()[1:].rstrip().split('\t')]
    tab.columns = header
    tab.sort_values(by=['chr', 'start'], inplace=True)

    return tab


def add_annotation_marker(ax, annotation, increment, xmin, xmax, flipped):
    x = []
    tab = pd.read_csv(annotation, sep='\t')
    subset = tab.loc[(tab.start > xmin) & (tab.end < xmax), :]

    for i, locus in subset.iterrows():
        if flipped:
            start, end = locus['start'], locus['end']
            x1 = xmin if end < xmin else (xmax - end) * increment
            x2 = xmax if end > xmax else (xmax - start) * increment
        else:
            start, end = locus['start'], locus['end']
            x1 = xmin if start < xmin else (start - xmin) * increment
            x2 = xmax if end > xmax else (end - xmin) * increment
        x.append((x1 + x2) / 2)

    ax.plot(x, [0.6] * len(x), '|', color='black', markeredgewidth=0.1)


def add_annotation_line2D(ax, annotation, increment, xmin, xmax, flipped, alternating=False, mirror_horizontal=False):
    tab = pd.read_csv(annotation, sep='\t')
    subset = tab.loc[(tab.start > xmin) & (tab.end < xmax), :]
    
    latex(True)

    if flipped:
        for i, locus in subset.iterrows():
            start, end = locus['start'], locus['end']
            x1 = xmin if end < xmin else (xmax - end) * increment
            x2 = xmax if end > xmax else (xmax - start) * increment

            color = 'black' if pd.isna(locus['color']) else locus['color']
            alpha = 1 if pd.isna(locus['alpha']) else float(locus['alpha'])

            if alternating:
                if i % 2 == 0:
                    ax.add_line(Line2D([x1, x2], [0.425, 0.425], lw=5, solid_capstyle='butt', color=color, alpha=alpha))

                else:
                    ax.add_line(Line2D([x1, x2], [0.575, 0.575], lw=5, solid_capstyle='butt', color=color, alpha=alpha))

            else:
                if not mirror_horizontal:
                    ax.add_line(Line2D([x1, x2], [0.75, 0.75], lw=5, solid_capstyle='butt', color=color, alpha=alpha))
                else:
                    ax.add_line(Line2D([x1, x2], [0.25, 0.25], lw=5, solid_capstyle='butt', color=color, alpha=alpha))


            if alternating:
                if i % 2 == 0:
                    ax.text((x2 - x1) / 2 + x1, 0.325, 
                    '' if pd.isna(locus['display_name']) else locus['display_name'],
                    ha='center' if pd.isna(locus['pos']) else locus['pos'],
                    va='top')

                else:
                    ax.text((x2 - x1) / 2 + x1, 0.675, 
                    '' if pd.isna(locus['display_name']) else locus['display_name'],
                    ha='center' if pd.isna(locus['pos']) else locus['pos'],
                    va='bottom')

            else:
                if not mirror_horizontal:
                    ax.text((x2 - x1) / 2 + x1, 0.5, 
                    '' if pd.isna(locus['display_name']) else locus['display_name'],
                    ha='center' if pd.isna(locus['pos']) else locus['pos'],
                    va='top')
                else:
                    ax.text((x2 - x1) / 2 + x1, 0.5, 
                    '' if pd.isna(locus['display_name']) else locus['display_name'],
                    ha='center' if pd.isna(locus['pos']) else locus['pos'],
                    va='bottom')

    else:
        for i, locus in subset.iterrows():
            start, end = locus['start'], locus['end']
            x1 = xmin if start < xmin else (start - xmin) * increment
            x2 = xmax if end > xmax else (end - xmin) * increment

            color = 'black' if pd.isna(locus['color']) else locus['color']
            alpha = 1 if pd.isna(locus['alpha']) else float(locus['alpha'])

            if alternating:
                if i % 2 == 0:
                    ax.add_line(Line2D([x1, x2], [0.425, 0.425], lw=5, solid_capstyle='butt', color=color, alpha=alpha))

                else:
                    ax.add_line(Line2D([x1, x2], [0.575, 0.575], lw=5, solid_capstyle='butt', color=color, alpha=alpha))

            else:
                if not mirror_horizontal:
                    ax.add_line(Line2D([x1, x2], [0.75, 0.75], lw=5, solid_capstyle='butt', color=color, alpha=alpha))
                else:
                    ax.add_line(Line2D([x1, x2], [0.25, 0.25], lw=5, solid_capstyle='butt', color=color, alpha=alpha))


            if alternating:
                if i % 2 == 0:
                    ax.text((x2 - x1) / 2 + x1, 0.325, 
                    '' if pd.isna(locus['display_name']) else locus['display_name'],
                    ha='center' if pd.isna(locus['pos']) else locus['pos'],
                    va='top')

                else:
                    ax.text((x2 - x1) / 2 + x1, 0.675, 
                    '' if pd.isna(locus['display_name']) else locus['display_name'],
                    ha='center' if pd.isna(locus['pos']) else locus['pos'],
                    va='bottom')

            else:
                if not mirror_horizontal:
                    ax.text((x2 - x1) / 2 + x1, 0.5, 
                    '' if pd.isna(locus['display_name']) else locus['display_name'],
                    ha='center' if pd.isna(locus['pos']) else locus['pos'],
                    va='top')
                else:
                    ax.text((x2 - x1) / 2 + x1, 0.5, 
                    '' if pd.isna(locus['display_name']) else locus['display_name'],
                    ha='center' if pd.isna(locus['pos']) else locus['pos'],
                    va='bottom')
    latex(False)

def smooth(values, smoothwindow):
    if len(values) % smoothwindow:
        rows = len(values) // smoothwindow + 1

    else:
        rows = len(values) // smoothwindow

    return np.nanmean(values.reshape(rows, smoothwindow), axis=1)


def add_bigwig_track(ax, bigwig, chrom, start, end, xmin, xmax, flipped, smooth_track=True, smoothwindow=100):
    bw = pbw.open(bigwig)
    values = np.array(bw.values(chrom, int(start), int(end)))
    if smooth_track:
        values = smooth(values, smoothwindow)

    if flipped:
        values = np.flip(values)

    ax.fill_between(np.linspace(xmin, xmax, len(values)), values, color='grey', linewidth=0.05)


def plot_annotation(ax,
                    track,
                    ylabel,
                    alternating,
                    plottype,
                    ylim,
                    number_of_bins,
                    start,
                    end,
                    chrom,
                    scaling=1000,
                    xticknum=None,
                    mirror_horizontal=False,
                    flipped=True):
    increment = number_of_bins / (end - start)
    ax.set_ylim(bottom=ylim[0], top=ylim[1])
    ax.set_yticks([])
    ax.set_ylabel(ylabel, fontsize=15)

    if plottype == 'Line2D':
        add_annotation_line2D(ax, track, increment, start, end, flipped, alternating, mirror_horizontal)

    elif plottype == 'Marker':
        add_annotation_marker(ax, track, increment, start, end, flipped)
        ax.set_ylabel(ylabel, ha='right', rotation='horizontal', va='center')


    elif plottype == 'bigwig':
        add_bigwig_track(ax,
                         track,
                         chrom,
                         start,
                         end,
                         0,
                         number_of_bins,
                         flipped)
        # ax.set_yticks(ylim)
        ax.set_ylabel(ylabel, ha='right', rotation='horizontal', va='center')

    else:
        raise Exception("plottype not supported")

    if not xticknum:
        for loc in ['left', 'top', 'right', 'bottom']:
            ax.spines[loc].set_visible(False)
        ax.set_xticks([])

    else:
        for loc in ['left', 'top', 'right']:
            ax.spines[loc].set_visible(False)

        ax.set_xlabel(chrom)
        ax.set_xticks(np.linspace(0, number_of_bins, xticknum))
        ax.set_xticklabels(['{val:,}'.format(val=i) for i in np.linspace(start/scaling, end/scaling, xticknum, dtype=int)])

    ax.set_xlim((0, number_of_bins))
    return ax


def plot_matrix(ax,
                mat,
                cmap,
                xrange,
                chrom,
                flipped,
                scaling=1000,
                capturebins=None,
                highlightbins=None,
                xlabel=None,
                xticknum=0,
                cbarwidth=0.025,
                vmin=0,
                vmax=50,
                mirror_horizontal=False,
                subplot_label=None,
                colorbar_ticks = True,
                colorbar_label='Norm. Interactions',
                colorbar_range=None):
    '''
    plotting function for triC results

    :param ax:                  plt.Axes object to generate the plot in
    :param mat:                 TriC matrix as generated by the CCseq pipeline
    :param cmap:                colormap to use for plotting
    :param xrange:              tuple of integer coordinates in bp of the genomic region plotted in the matrix
    :param chrom:               string naming the chromosome of the locus
    :param chrom:               boolean indicating whether the matrix should be displayed in 5 prime 3 prime 
                                direction (flipped=True)
    :param scaling:             divisor for scaling the xrange
    :param capturebins:         bins containing the capture probes (optional)
    :param highlightbins:       list of tuples of startbin, endbin and highlightcolor for bins to highlight
                                if endbin == startbin only the startbin will be highlighted
    :param xlabel:              xaxis label
                                the annotation given at index -1 will be the bottom most annotation
    :param xticknum:            number of xticks to plot. xticklabels will be interpolated with np.linspace
    :param cbarwidth:           width of the colorbar as fraction of x-axis length
    :param vmin:                minimum value of the colormap
    :param vmax:                maximum value of the colormap
    :param mirror_horizontal:   indicates if generated matrix plot should be mirrored at a horizontal line
    :param subplot_label:       label of the matrix, None = no label
    :param colorbar_label:      label of the colorbar, Default = 'Normalized Interaction counts'
    :param colorbar_range:      colorbar labels at the top and bottom, describing the range, Default: None = no labels

    :return:                    plt.Axes
    '''


    xrange = (xrange[0] / scaling, xrange[1] / scaling)

    N = mat.shape[0] # spacing made with the help of N always has to be done relative, as N can change drastically.
    # Get the lower triangle of the matrix.
    C = np.triu(mat)
    # Mask the upper triangle
    C = np.ma.masked_array(C, C == 0)

    # Check for flipped locus:
    if flipped:
        C = np.rot90(np.fliplr(C), 3)

    # Transformation matrix for rotating the heatmap.
    A = np.array([(y, x) for x in range(N, -1, -1) for y in range(N + 1)])
    t = np.array([[1, 0.5], [-1, 0.5]]) if not mirror_horizontal else np.array([[-1, 0.5], [1, 0.5]])
    A = np.dot(A, t)

    # Plot the heatmap triangle.
    X = A[:, 1].reshape(N + 1, N + 1)
    Y = A[:, 0].reshape(N + 1, N + 1)
    mesh = ax.pcolormesh(X, Y, np.flipud(C), cmap=cmap, vmin=vmin, vmax=vmax, zorder=2)
    # mesh = ax.pcolormesh(X, Y, C, cmap=cmap, vmin=vmin, vmax=vmax, zorder=2)

    # mark bin containing capturesite in grey
    if capturebins:
        greymap = clr.LinearSegmentedColormap.from_list('greymap', ['Grey', 'Grey'], N=256)
        capturemat = np.zeros(shape=mat.shape)
        
        # flip if necessary
        if flipped:
            for capturebin in capturebins:
                if capturebin is not None:
                    capturemat[N-capturebin, :] = 1
                    capturemat[:, N-capturebin] = 1
        else:
            for capturebin in capturebins:
                if capturebin is not None:
                    capturemat[capturebin, :] = 1
                    capturemat[:, capturebin] = 1

        capturemat = np.triu(capturemat)
        capturemat = np.ma.masked_array(capturemat, capturemat == 0)
        capturemesh = ax.pcolormesh(X, Y, np.flipud(capturemat), cmap=greymap, vmin=0, vmax=1, zorder=3)

    if highlightbins:
        for hlstartbin, hlendbin, hlclr in highlightbins:
            hlmap = clr.LinearSegmentedColormap.from_list(hlclr, [hlclr, hlclr], N=256)
            hlmat = np.zeros(shape=mat.shape)

            if flipped:
                if not hlendbin == None:
                    hlidx = np.arange(N-hlendbin, N-hlstartbin + 1)
                else:
                    hlidx = N-hlstartbin
            else:
                if not hlendbin == None:
                    hlidx = np.arange(hlstartbin, hlendbin + 1)
                else:
                    hlidx = hlstartbin

            hlmat[hlidx, :] = 1
            hlmat[:, hlidx] = 1
            hlmat = np.triu(hlmat)
            hlmat = np.ma.masked_array(hlmat, hlmat == 0)
            hlmesh = ax.pcolormesh(X, Y, np.flipud(hlmat), cmap=hlmap, vmin=0, vmax=1, alpha=0.2, zorder=1)

    # draw outlines of triangle plot
    vertices = np.array([[0, 0], [N / 2, N], [N, 0]]) if not mirror_horizontal else np.array(
        [[0, 0], [N / 2, -N], [N, 0]])
    triangle = patches.Polygon(vertices, fill=False, edgecolor='black')
    ax.add_patch(triangle)

    ax.set_xlim(left=0, right=N)

    if not mirror_horizontal:
        ax.set_ylim(bottom=0, top=N)

    else:
        ax.set_ylim(bottom=-N, top=0)

    if not mirror_horizontal:
        for loc in ['left', 'right', 'top']:
            ax.spines[loc].set_visible(False)

    else:
        for loc in ['left', 'right', 'bottom']:
            ax.spines[loc].set_visible(False)

        ax.xaxis.tick_top()

    ax.set_yticks([])
    ax.set_xticks([])

    if xticknum:
        ax.set_xticks(np.linspace(0, N, xticknum))
        ax.tick_params(labelbottom=False, labeltop=False)    

    if xlabel:
        ax.set_xlabel(xlabel)

    # plot colorbar
    if not mirror_horizontal:

        if colorbar_ticks:

            rect = patches.Rectangle((N * cbarwidth, N * 3/8), N * cbarwidth, N / 2, fill=False, edgecolor='white') 

            ax.add_patch(rect)

            cbarY = np.tile(np.linspace(N * 3/8, N * 7/8, cmap.N).reshape(-1, 1), 2)
            cbarX = np.tile(np.array([N * cbarwidth, N * cbarwidth * 2]), (cbarY.shape[0], 1))
            cbarmesh = ax.pcolormesh(cbarX, cbarY, np.linspace(0, 1, cmap.N - 1).reshape(-1, 1), cmap=cmap, vmin=0, vmax=1)

            ys = np.linspace(N * 3/8, N * 7/8, 5)
            for y, cmapval in zip(ys, np.linspace(vmin, vmax, 5)):
                ax.add_line(
                    Line2D([N * cbarwidth * 2, N * cbarwidth * 2 + N * 0.005], [y, y], color='black', lw=plt.rcParams['patch.linewidth']))
                ax.text(N * cbarwidth * 2 + N * 0.0075, y, '{:.01f}'.format(cmapval), ha='left', va='center', fontsize=15)

            ax.text(N * cbarwidth * 0.9, N * 5/8 , colorbar_label, ha='right', va='center', rotation=90, fontsize=15)

        elif colorbar_range:
            rect = patches.Rectangle((N - N * cbarwidth, N / 2), N * cbarwidth, N / 2, fill=False, edgecolor='white') 

            ax.add_patch(rect)

            cbarY = np.tile(np.linspace(N / 2, N, cmap.N).reshape(-1, 1), 2)
            cbarX = np.tile(np.array([N - N * cbarwidth, N]), (cbarY.shape[0], 1))
            cbarmesh = ax.pcolormesh(cbarX, cbarY, np.linspace(0, 1, cmap.N - 1).reshape(-1, 1), cmap=cmap, vmin=0, vmax=1)

            ys = np.linspace(N / 2, N, 5)

            ax.text(N - N * cbarwidth - N * 0.0075, ys[0], colorbar_range[0], ha='right', fontsize=15)
            ax.text(N - N * cbarwidth - N * 0.0075, ys[-1], colorbar_range[1], ha='right', fontsize=15)

            ax.text(N + 1, 3 * N / 4 , colorbar_label, ha='left', va='center', rotation=90, fontsize=15)

    if subplot_label:
        ax.text(0, N + N * 0.02 if not mirror_horizontal else -N - N * 0.06, subplot_label,
                ha='left')

    # add chromosome location
    if mirror_horizontal: 
        ax.text(N, -N - N * 0.07, f"{chrom}: {'{val:,}'.format(val=int(xrange[0]))},000-{'{val:,}'.format(val=int(xrange[1]))},000", ha='right')

    return ax


def plot_profile_overlay(ax,
                         profiledict,
                         number_of_bins,
                         xrange,
                         colors,
                         flipped,
                         yrange=None,
                         scaling=1000,
                         capturebins=None,
                         ylabel=None,
                         xlabel=None,
                         xticknum=0):
    xrange = (xrange[0] / scaling, xrange[1] / scaling)

    if ylabel:
        ax.set_ylabel(ylabel)

    ax.set_xticks([])
    ax.set_xlim((0, number_of_bins))

    if xticknum:
        ax.set_xticks(np.linspace(0, number_of_bins, xticknum))
        ax.set_xticklabels(['{val:,}'.format(val=i) for i in np.linspace(xrange[0], xrange[1], xticknum, dtype=int)])

    if yrange:
        ax.set_ylim(yrange)
        ax.set_yticks(yrange)
        ax.set_yticklabels(['', str(yrange[1])], fontsize = 15)

    if xlabel:
        ax.set_xlabel(xlabel)

    for loc in ['left', 'top', 'right']:
        ax.spines[loc].set_visible(False)

    # iterate through profiles for plotting, flip the dict if necessary
    for color, (profilename, profile) in zip(colors, profiledict.items()):
        
        if flipped:
            profile = list(reversed(profile))

        ax.plot(np.arange(number_of_bins) + 0.5, profile, color=color, label=profilename)

        if capturebins:
            for capturebin in capturebins:
                if capturebin is not None:
                    if flipped:
                        capturebin = len(profile) - capturebin

                    ax.bar(capturebin + 0.5, ax.get_ylim()[1], align='center', width=0.75, color='black')

    if flipped:
        # ax.legend(loc="upper right", frameon=False, handlelength=1, ncol = 2)
        ax.legend(loc=(0.9, 0.2), frameon=False, handlelength=1, fontsize=15)
    else:
        ax.legend(loc=(0, 0.5), frameon=False, handlelength=1, fontsize=15)
    return ax


def compute_average_matrix(matrices):
    summed = np.zeros(shape=matrices[0].shape)
    for mat in matrices:
        summed += mat

    return summed / len(matrices)


def make_difference_matrix(mat1, mat2):
    # calculate the % of counts from mat1 if sum(mat1, mat2), throws unnecessary warnings if both bins are 0
    with np.errstate(divide='ignore',invalid='ignore'): percentual_diff = (mat1 / (mat1 + mat2)) * 100
    percentual_diff = np.nan_to_num(percentual_diff, nan=50) - 50

    return percentual_diff


def get_bin_index(captureSiteStart, leftBound, rightBound, binsize):
    binbounds = np.arange(leftBound, rightBound, binsize)
    return len(np.where(binbounds < captureSiteStart)[0]) \
           if not (captureSiteStart < binbounds[0] or captureSiteStart > binbounds[-1]) \
           else None


def get_highlight_bin_argument_from_annotation(annotation_name, features, leftBound, rightBound, binsize, hlcolor='cyan'):
    '''
    takes an annotation table, a list of names of features in the table, the boundaries and binsize of the plot region
    and returns the argumentlist that can be passed to highlightbins in plotMatrix

    :param annotation_name:     annotation file name
    :param features:            names of the features in annotation that should be highlighted
    :param leftBound:           left boundary position of the plotting range
    :param rightBound:          right boundary position of the plotting range
    :param binsize:             size of bins plotted
    :param hlcolor:             color to use to highlight or list of colors if list has to be the same length as features

    :return:            highlightbins argument list
    '''

    annotation_df = pd.read_csv(annotation_name, sep = '\t')
    hlargument = []
    features_df = annotation_df.loc[annotation_df.name.isin(features), :] \
        .reset_index()

    if isinstance(hlcolor, list):
        if len(hlcolor) == len(features):
            for i, feature in features_df.iterrows():
                hlstartbin = get_bin_index(feature['start'], leftBound, rightBound, binsize)
                hlendbin = get_bin_index(feature['end'], leftBound, rightBound, binsize)
                if not (hlstartbin is None or hlendbin is None):
                    hlargument.append((hlstartbin, hlendbin, hlcolor[i]))
        else:
            raise Exception(
                'numbers of color has to match the number of features to highlight if multiple colors are passed')

    else:
        for i, feature in features_df.iterrows():
            hlstartbin = get_bin_index(feature['start'], leftBound, rightBound, binsize)
            hlendbin = get_bin_index(feature['end'], leftBound, rightBound, binsize)
            if not (hlstartbin is None or hlendbin is None):
                hlargument.append((hlstartbin, hlendbin, hlcolor))

    return hlargument


def same_length(*args):
    l = len(args[0])
    for arg in args[1:]:
        if len(arg) != l:
            return False

    return True


def get_region(region):
    chrom, bounds = region.split(':')
    start, end = [int(i) for i in bounds.split('-')]
    return chrom, start, end


def get_zoom_matrix(mat, region, subregion, binsize):
    r_chrom, r_start, r_end = get_region(region)
    s_chrom, s_start, s_end = get_region(subregion)

    r_bins = (r_end - r_start) // binsize
    r_bin_bounds = np.linspace(r_start, r_end, r_bins + 1)

    s_start_bin_in_r = np.where(r_bin_bounds < s_start)[0][-1] + 1
    s_end_bin_in_r = np.where(r_bin_bounds < s_end)[0][-1] + 1

    return mat[s_start_bin_in_r: s_end_bin_in_r, s_start_bin_in_r: s_end_bin_in_r]


def load_profiles(treatment_profile, control_profile, treatment_label, control_label, leftBound, rightBound,
                  capturebins):
    profiles = {}
    for k, file in zip([treatment_label, control_label], [treatment_profile, control_profile]):
        profiletab = load_profile_table(file)
        
        meanprofile = profiletab.loc[:, ~profiletab.columns.isin(['chr', 'start', 'end'])].mean(axis=1)
        if meanprofile.sum():
            totalnorm = 100000 / meanprofile.sum()

        else:
            totalnorm = 1

        # binnorm = 1000 / (meanprofile * totalnorm).max()
        for capturebin in capturebins:
            if capturebin is not None:
                meanprofile.loc[capturebin + 1] = 0  # setting capture site counts to 0, +1 because of count start at 0

        profiletab['meanprofile'] = meanprofile
        profiletab = profiletab \
                          .loc[(profiletab['start'] >= leftBound) & (profiletab['start'] < rightBound), :] \
                          .reset_index(drop=True)
        profiles[k] = profiletab['meanprofile'] * totalnorm  # * binnorm

    return profiles


def get_colormap(colors, N = 256):
    return clr.LinearSegmentedColormap.from_list('custom', colors, N=N) if len(colors) > 1 else plt.get_cmap(*colors)

#wyorb = clr.LinearSegmentedColormap.from_list('wyorb', ['White', 'Yellow', 'Orange', 'Red', 'Black'], N=256)
#gyorb = clr.LinearSegmentedColormap.from_list('gorb', ['whitesmoke', 'Yellow', 'Orange', 'Red', 'Black'], N=256)
#gorb = clr.LinearSegmentedColormap.from_list('gorb', ['whitesmoke', 'Orange', 'Red', 'Black'], N=256)
#bwr = plt.get_cmap('bwr')
# chrom, leftBound, rightBound, binsize = 'chr12', 114435000, 114669000, 1000

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('--treatment', '-t', nargs='+',
                    help='one or more TriC interaction files for treatment condition'
                         'if more than one are given a mean is computed over all samples')
parser.add_argument('--treatment_label', '-tl',
                    help='name of the treatment condition to use for plot naming etc.'
                         'if not given the script tries to infer a name from the filename')
parser.add_argument('--control', '-c', nargs='+',
                    help='one or more TriC interaction files for control condition'
                         'if more than one are given a mean is computed over all samples')
parser.add_argument('--control_label', '-cl',
                    help='name of the control condition to use for plot naming etc.'
                         'if not given the script tries to infer a name from the filename')
parser.add_argument('--region', '-r', required=True,
                    help='region the input matrices span and to consider for plotting in the format of chr:start-stop')
parser.add_argument('--capture_bins', required=True,
                    help='oligo file used for processing with CCseqBasic denoting the bins in the matrix and profile to mask')
parser.add_argument('--binsize', '-bs', default=1000, type=int,
                    help='size of the bins the region is subdivided into')
parser.add_argument('--annotation', '-a', nargs='*',
                    help='space-separated list of filenames that should be added to the figure'
                         'currently only bed and bigwig files are supported. Bed files can be in BED3 or BED4 format'
                         'which means chr, start, end or chr, start, end, name columns present with header.'
                         'if name is present the plotted annotations will be named')
parser.add_argument('--annotation_drawstyle', '-ad', nargs='*',
                    help='space-separated list of drawstyles for each annotation file'
                         'in the same order as the files (has to be one of Line2D, Marker or bigwig'
                         'Line2D is usually use for bedfiles with small numbers of regions to display'
                         'e.g. genes in the locus, Marker is used for bedfiles with a large number of regions'
                         'e.g. NlaIII restriction sites and bigwig is reserved for bigwig files)')
parser.add_argument('--annotation_labels', '-al', nargs='*',
                    help='names of the annotations used to name the generated plots'
                         'in the same order as the files given')
parser.add_argument('--annotation_yMin', nargs='*',
                    help='minimum value of y axis of the annotation plots')
parser.add_argument('--annotation_yMax', nargs='*',
                    help='maximum value of y axis of the annotation plots')
parser.add_argument('--alternating', nargs='*',
                    help='space-separated list of 0 and 1 indicating if the the given annotation'
                         'should be drawn in an alternating fashing. Most suitable for dense annotations'
                         'to avoid label overlap')
parser.add_argument('--highlight_annotation', type=int,
                    help='the 1-based --annotation list index of the annotation to use for highlighting')
parser.add_argument('--highlight_features', nargs='*',
                    help='space-separated list of features in the --highlight_annotation to use for highlighting')
parser.add_argument('--treatment_3plus',
                    help='tab-separated table containing counts of reads with 3 or more valid restriction fragments in the treatment condition')
parser.add_argument('--control_3plus',
                    help='tab-separated table containing counts of reads with 3 or more valid restriction fragments in the control condition')
parser.add_argument('--profile_yMax', type = float, default = 200,
                    help='maximum value of the y axis of the profile plot')
parser.add_argument('--profile_labels', nargs='*',
                    help='space-separated list of the labels of the profile plot')
parser.add_argument('--compare_vMin', default = 0, type = float,
                    help = 'minimum value of colorbars in the compare matrix plots')
parser.add_argument('--compare_vMax', default = 50, type = float,
                    help = 'maximum value of colorbars in the compare matrix plots')
parser.add_argument('--diff_vMin', default = -50, type = float,
                    help = 'minimum value of colorbars in the difference matrix plots')
parser.add_argument('--diff_vMax', default = 50, type = float,
                    help = 'maximum value of colorbars in the difference matrix plots')
parser.add_argument('--figwidth', default=10, type=float,
                    help='width of the generated figure. Height is computed accordingly.')
parser.add_argument('--compare_colormap', default = 'whitesmoke,orange,red,black',
                    help = 'either a name of colormap predefined in matplotlib or a comma-separated list of colors'
                           'where position of the color in the list corresponds to the value it represents'
                           'with first = smallest, last = highest')
parser.add_argument('--diff_colormap', default = 'bwr',
                    help = 'either a name of colormap predefined in matplotlib or a comma-separated list of colors'
                           'where position of the color in the list corresponds to the value it represents'
                           'with first = smallest, last = highest')
parser.add_argument('--flipped', action='store_true',
                    help="flipping of the matirces and annotation so that the locus is displayed in 5' to 3' direction. Default is False.")
parser.add_argument('--outputFilePrefix', '-o', required=True,
                    help='prefix to use for output files')
args = parser.parse_args()

chrom, leftBound, rightBound = get_region(args.region)
n_bins = (rightBound - leftBound) // args.binsize

compare_cmap = get_colormap(args.compare_colormap.split(','))
diff_cmap = get_colormap(args.diff_colormap.split(','))

profile_treatment_label = args.treatment_label if not args.profile_labels else args.profile_labels[0]
profile_control_label = args.control_label if not args.profile_labels else args.profile_labels[1]

annotations = []
number_of_annotation_axes = 0
if args.annotation:
    assert same_length(args.annotation,
                       args.annotation_drawstyle,
                       args.annotation_labels,
                       args.annotation_yMin,
                       args.annotation_yMax,
                       args.alternating), \
        'If annotations are given, all arguments must have the same number of parameters'

    annotations = list(zip(args.annotation,
                           args.annotation_labels,
                           [int(i) for i in args.alternating],
                           args.annotation_drawstyle,
                           zip([float(i) for i in args.annotation_yMin],
                               [float(i) for i in args.annotation_yMax]),
                           [0] * (len(args.annotation) - 1) + [10]))  # xtick visibility

    number_of_annotation_axes = len(annotations)

treatments = [np.loadtxt(file) for file in args.treatment]
treatment_avg = compute_average_matrix(treatments)

controls = [np.loadtxt(file) for file in args.control]
control_avg = compute_average_matrix(controls)

# setting up figure layout
#plt.rcParams['figure.constrained_layout.use'] = True
fig1 = plt.figure(dpi=300)
matrix_subplot_height = args.figwidth / 2
annotation_height = 0.3
profile_height = 1
hspace = 0.2
hspaceFig2 = hspace + 0.05

profile_args = [args.treatment_3plus, args.control_3plus, args.profile_yMax]
if any(profile_args):
    assert all(profile_args), 'all profile arguments have to be set if one is used'
    number_of_annotation_axes += 1
    fig1_height = 2 * matrix_subplot_height + \
                  annotation_height * (number_of_annotation_axes - 1) + \
                  (annotation_height if annotations else 0) + \
                  profile_height + \
                  hspace * (number_of_annotation_axes + 2)
    number_of_axes_fig1 = number_of_annotation_axes + 3
    height_ratios_fig1 = [matrix_subplot_height / fig1_height] + \
                         [annotation_height / fig1_height] * (number_of_annotation_axes - 1) + \
                         [profile_height / fig1_height] + \
                         ([annotation_height / fig1_height] if annotations else []) + \
                         [matrix_subplot_height / fig1_height]

    number_of_axes_fig2 = number_of_annotation_axes + 1
    fig2_height = matrix_subplot_height + \
                  annotation_height * (number_of_annotation_axes) + \
                  profile_height + \
                  hspaceFig2 * number_of_annotation_axes
    height_ratios_fig2 = [matrix_subplot_height / fig2_height] + \
                         [annotation_height / fig2_height] * (number_of_annotation_axes - 1) + \
                         [profile_height / fig1_height]

else:
    fig1_height = 2 * matrix_subplot_height + \
                  annotation_height * (number_of_annotation_axes) + \
                  (annotation_height if annotations else 0) + \
                  hspace * (number_of_annotation_axes + 2)
    number_of_axes_fig1 = number_of_annotation_axes + 3
    height_ratios_fig1 = [matrix_subplot_height / fig1_height] + \
                         [annotation_height / fig1_height] * (number_of_annotation_axes) + \
                         ([annotation_height / fig1_height] if annotations else []) + \
                         [matrix_subplot_height / fig1_height]

    number_of_axes_fig2 = number_of_annotation_axes + 1
    fig2_height = matrix_subplot_height + \
                  annotation_height * (number_of_annotation_axes) + \
                  hspaceFig2 * number_of_annotation_axes
    height_ratios_fig2 = [matrix_subplot_height / fig2_height] + \
                         [annotation_height / fig2_height] * (number_of_annotation_axes)


# compare figure
fig1.set_figwidth(args.figwidth)
fig1.set_figheight(fig1_height)
gridspec1 = gs.GridSpec(number_of_axes_fig1,
                        1,
                        height_ratios=height_ratios_fig1,
                        figure=fig1,
                        hspace=hspace)
treatment_ax = fig1.add_subplot(gridspec1[0])
control_ax = fig1.add_subplot(gridspec1[-1])
profile_ax1 = fig1.add_subplot(gridspec1[number_of_annotation_axes]) if any(profile_args) else None
annotation_axs1 = [fig1.add_subplot(gridspec1[i + 1]) for i in range(len(annotations))] if annotations else []
secondLine2d_ax1 = fig1.add_subplot(gridspec1[-2]) if annotations else None

# diff figure
fig2 = plt.figure(dpi=300)
fig2.set_figwidth(args.figwidth)
fig2.set_figheight(fig2_height)
gridspec2 = gs.GridSpec(number_of_axes_fig2,
                        1,
                        height_ratios=height_ratios_fig2,
                        figure=fig2,
                        hspace=hspaceFig2)
diff_ax = fig2.add_subplot(gridspec2[0])
profile_ax2 = fig2.add_subplot(gridspec2[-1]) if any(profile_args) else None
annotation_axs2 = [fig2.add_subplot(gridspec2[i + 1]) for i in range(len(annotations))] if annotations else []

oligotab = pd.read_csv(args.capture_bins,
                       sep='\t',
                       header=None,
                       names=['name', 'chrom', 'start', 'end'],
                       usecols=[0, 1, 2, 3])

capturebins = [get_bin_index(r['start'], leftBound, rightBound, args.binsize) for i, r in oligotab.iterrows()]
highlightbins = []
if args.highlight_annotation:
    highlightbins = get_highlight_bin_argument_from_annotation(annotations[args.highlight_annotation - 1][0],
                                                               args.highlight_features,
                                                               leftBound,
                                                               rightBound,
                                                               args.binsize,
                                                               'cyan')

treatment_ax = plot_matrix(treatment_ax,
                           treatment_avg,
                           compare_cmap,
                           (leftBound, rightBound),
                           chrom,
                           args.flipped,
                           capturebins=capturebins,
                           highlightbins=highlightbins,
                           vmin=args.compare_vMin,
                           vmax=args.compare_vMax,
                           xticknum=11,
                           subplot_label=args.treatment_label)

control_ax = plot_matrix(control_ax,
                         control_avg,
                         compare_cmap,
                         (leftBound, rightBound),
                         chrom,
                         args.flipped,
                         capturebins=capturebins,
                         highlightbins=highlightbins,
                         vmin=args.compare_vMin,
                         vmax=args.compare_vMax,
                         xticknum=11,
                         mirror_horizontal=True,
                         subplot_label=args.control_label)

diff_ax = plot_matrix(diff_ax,
                      make_difference_matrix(treatment_avg, control_avg),
                      diff_cmap,
                      (leftBound, rightBound),
                      chrom,
                      args.flipped,
                      capturebins=capturebins,
                      highlightbins=highlightbins,
                      vmin=args.diff_vMin,
                      vmax=args.diff_vMax,
                      xticknum=11,
                      subplot_label=' - '.join((args.treatment_label, args.control_label)),
                      colorbar_ticks=False,
                      colorbar_label='Enrichment',
                      colorbar_range=[profile_control_label, profile_treatment_label])

if any(profile_args):
    profiles = load_profiles(args.treatment_3plus, args.control_3plus,
                             profile_treatment_label, 
                             profile_control_label,
                             leftBound, rightBound,
                             capturebins)

    for profile_ax in [profile_ax1, profile_ax2]:
        ax = plot_profile_overlay(profile_ax,
                                  profiles,
                                  n_bins,
                                  (leftBound, rightBound),
                                  yrange=(0, args.profile_yMax),
                                  capturebins=capturebins,
                                  colors=('black', 'red'),
                                  flipped=args.flipped)

if annotations:
    for annotation_axs in [annotation_axs1, annotation_axs2]:
        for (annoax, (track, label, alternating, plottype, ylim, number_of_xticks)) in zip(annotation_axs, annotations):
            ax = plot_annotation(annoax,
                                 track,
                                 label,
                                 alternating,
                                 plottype,
                                 ylim,
                                 n_bins,
                                 leftBound,
                                 rightBound,
                                 chrom,
                                 flipped=args.flipped)

    # plot only the Line2d a second time for bottom matrix
    ax = plot_annotation(secondLine2d_ax1, 
                         annotations[0][0], 
                         annotations[0][1], 
                         annotations[0][2], 
                         annotations[0][3], 
                         annotations[0][4], 
                         n_bins, 
                         leftBound, 
                         rightBound, 
                         chrom,
                         mirror_horizontal=True,
                         flipped=args.flipped)

# fig1.tight_layout(pad = 3, h_pad = hspace)
# fig2.tight_layout(pad = 3, h_pad = hspace)
fig1.savefig(args.outputFilePrefix + '_compare.pdf')
fig2.savefig(args.outputFilePrefix + '_difference.pdf')
