import numpy as np
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import matplotlib.gridspec as gs
import matplotlib as mpl
import pandas as pd
import pyBigWig as pbw
import itertools as it
import tables
from scipy.sparse.csr import csr_matrix
import sys
import logging
mpl.rcParams['pdf.fonttype'] = 42

def toString(s):
    """
    This takes care of python2/3 differences
    """
    if isinstance(s, str):
        return s

    if isinstance(s, bytes):  # or isinstance(s, np.bytes_):
        if sys.version_info[0] == 2:
            return str(s)
        return s.decode('ascii')

    if isinstance(s, list):
        return [toString(x) for x in s]

    if isinstance(s, np.ndarray):
        return s.astype(str)

    return s


def loadH5subset(filename, includechroms=None, csr = True):
    '''
    loadH5(filename, includechroms=None, csr=True)

    loads an *.h5 hic matrix as created by hicexplorer

    :param filename:        name of the *.h5 file containing the matrix
    :param includechroms:   list of chromosomes to include in the returned objects
                            if not given all chromosomes in the *.h5 file are included
    :param csr:             if True returns a csr_matrix object else a full numpy.array

    :return:                csr_matrix containing the data in the matrix
    '''
    with tables.open_file(filename) as f:
        parts = {}
        try:
            for matrix_part in ('data', 'indices', 'indptr', 'shape'):
                parts[matrix_part] = getattr(f.root.matrix, matrix_part).read()
        except Exception:
            logging.info('No h5 file. Please check parameters concerning the file type!')
            exit(1)

        matrix = csr_matrix(tuple([parts['data'], parts['indices'], parts['indptr']]),
                            shape=parts['shape'], dtype = int)

        intervals = {}
        for interval_part in ('chr_list', 'start_list', 'end_list', 'extra_list'):
            if toString(interval_part) == toString('chr_list'):
                chrom_list = getattr(f.root.intervals, interval_part).read()
                intervals[interval_part] = toString(chrom_list)
            else:
                intervals[interval_part] = getattr(f.root.intervals, interval_part).read()

        cut_intervals = list(
            zip(intervals['chr_list'], intervals['start_list'], intervals['end_list'], intervals['extra_list']))

        assert len(cut_intervals) == matrix.shape[0], \
            "Error loading matrix. Length of bin intervals ({}) is different than the " \
            "size of the matrix ({})".format(len(cut_intervals), matrix.shape[0])

        # compute index array and chromosome list
        inds, chr_list, chroms = [], [], set()
        for i, (chr, start, end, extra) in enumerate(cut_intervals):
            if chr not in chroms:
                chroms.add(chr)
                inds.append(i)
                chr_list.append(chr)

        # if includechroms is given we filter the output for the chromosomes listed
        # and recompute indices of chromosome boundaries in the resulting matrix
        if includechroms:
            includechroms = set(includechroms)

            intervals = {k: [] for k in ['chr_list', 'start_list', 'end_list', 'extra_list']}
            for i, vals in enumerate(cut_intervals):
                if vals[0] in includechroms:
                    for k, v in zip(['chr_list', 'start_list', 'end_list', 'extra_list'], vals):
                        intervals[k].append(v)

            filterinds, filterchrs = [], []
            for i, chr in zip(range(len(inds)), chr_list):
                if chr in includechroms:
                    filterinds.append([inds[i], inds[i + 1] if i + 1 != len(inds) else matrix.shape[0]])
                    filterchrs.append(chr)

            matrixinds = np.zeros(shape=matrix.shape[0], dtype=bool)
            ncuts, tmpe = [], 0
            for s, e in filterinds:
                matrixinds[s: e] = True

                if s == tmpe:
                    ncuts.append(s)
                    tmpe = e

                else:
                    ncuts.append(tmpe)
                    tmpe = e - s + tmpe

            matrix = matrix[matrixinds, :][:, matrixinds]

    if not csr:
        x = matrix.toarray()
        xi, yi = np.triu_indices(x.shape[0], k=1)
        x[yi, xi] = x[xi, yi]
        matrix = x

    return matrix, intervals


def getRegionalMatrix(contactMatrix, intervalstarts, start, end, csr = False):
    indices, include = [], []
    for i, loc in enumerate(intervalstarts):
        if loc >= start and loc < end:
            indices.append(i)
            include.append(loc)

    contactMatrix = contactMatrix[indices,:][:, indices]

    return contactMatrix if not csr else contactMatrix.toarray()


def loadProfileTab(filename, interval = None):
    tab = pd.read_csv(filename, sep = '\t')
    header = [s.strip("'") for s in open(filename, 'r').readline()[1:].rstrip().split('\t')]
    tab.columns = header
    tab.sort_values(by = ['chr', 'start'], inplace = True)

    if interval:
        tab = tab.loc[(tab['start'] >= interval[0]) & (tab['start'] < interval[1]), :].reset_index(drop = True)

    return tab


def addAnnotationMarker(ax, annotation, increment, scale, xmin, xmax):
    x = []
    for i, loci in annotation.iterrows():
        start, end = loci['start'] / scale, loci['end'] / scale
        x1 = xmin if start < xmin else (start - xmin) * increment
        x2 = xmax if end > xmax else (end - xmin) * increment
        x.append((x1 + x2)/2)

    ax.plot(x, [0.5] * len(x), '|')


def addAnnotationLine2D(ax, annotation, increment, scale, xmin, xmax, alternating = False):
    for i, loci in annotation.iterrows():
        start, end = loci['start'] / scale, loci['end'] / scale
        x1 = xmin if start < xmin else (start - xmin) * increment
        x2 = xmax if end > xmax else (end - xmin) * increment

        if alternating:
            if i % 2 == 0:
                ax.add_line(Line2D([x1, x2], [0.425, 0.425], lw=5, solid_capstyle = 'butt'))

            else:
                ax.add_line(Line2D([x1, x2], [0.575, 0.575], lw=5, solid_capstyle = 'butt'))

        else:
            ax.add_line(Line2D([x1, x2], [0.5, 0.5], lw=5, solid_capstyle = 'butt'))

        if loci['name']:
            if i % 2 == 0:
                ax.text((x2 - x1) / 2 + x1, 0.25 if alternating else 0.325, loci['name'], ha='center', va='top', fontsize = 10)

            else:
                ax.text((x2 - x1) / 2 + x1, 0.75 if alternating else 0.675, loci['name'], ha='center', va='bottom', fontsize = 10)


def smooth(values, smoothwindow):
    if len(values)%smoothwindow:
        rows = len(values) // smoothwindow + 1

    else:
        rows = len(values) // smoothwindow

    return np.nanmean(values.reshape(rows, smoothwindow), axis = 1)


def addBigWigTrack(ax, bigwig, chrom, start, end, xmin, xmax, smooth_track = True, smoothwindow = 100):
    bw = pbw.open(bigwig)
    values = np.array(bw.values(chrom, int(start), int(end)))
    if smooth_track:
        values = smooth(values, smoothwindow)

    ax.fill_between(np.linspace(xmin, xmax, len(values)), values)


def plotMatrix(mat,
               cmap,
               chrom,
               xrange,
               scale = 'kbp',
               capturebin = None,
               highlightbins = None,
               xlabel = None,
               annotations = None,
               xticknum = 0,
               cbarwidth = 5,
               vmin = 0,
               vmax = 50,
               unit = 0.05):
    '''
    plotting function for triC results

    :param mat:             TriC matrix as generated by the CCseq pipeline
    :param cmap:            colormap to use for plotting
    :param xrange:          tuple of integer coordinates in bp of the genomic region plotted in the matrix
    :param scale:           determines the scale used (either 'none', 'kbp' or 'Mbp')
    :param capturebin:      bin containing the capture site (optional)
    :param highlightbins:   list of tuples of startbin, endbin and highlightcolor for bins that should be highlighted
                            if endbin == startbin only the startbin will be highlighted
    :param xlabel:          xaxis label
    :param annotations:     list of tuples containing five items:
                            (0. pandas.Dataframe with coordinates of annotations to add to the plot (has to have 4 columns chr, start, end, name),
                             1. the label of the annotation axis,
                             2. bool indicating annotations should be plotted alternating,
                             3. string indicating the type of plotting function to use, currently 'Line2D', 'Marker' and 'bigwig' are supported
                                in case of 'Marker' 2. has no effect
                             4. tuple of y-axis limits (ymin, ymax))
                            the annotation given at index -1 will be the bottom most annotation
    :param xticknum:        number of xticks to plot. xticklabels will be interpolated with np.linspace
    :param cbarwidth:       width of the colorbar
    :param vmin:            minimum value of the colormap
    :param vmax:            maximum value of the colormap
    :param unit:            plotting unit used to determine the height and width of the resulting plot

    :return:                plt.figure and list(plt.axes)
    '''

    scales = {'kbp': 1000, 'Mbp': 1000000}
    xrange = (xrange[0]/scales[scale], xrange[1]/scales[scale])

    N = mat.shape[0]
    # Get the lower triangle of the matrix.
    C = np.triu(mat)
    # Mask the upper triangle
    C = np.ma.masked_array(C, C == 0)

    # Transformation matrix for rotating the heatmap.
    A = np.array([(y, x) for x in range(N, -1, -1) for y in range(N + 1)])
    t = np.array([[1, 0.5], [-1, 0.5]])
    A = np.dot(A, t)

    # Plot the correlation heatmap triangle.
    X = A[:, 1].reshape(N + 1, N + 1)
    Y = A[:, 0].reshape(N + 1, N + 1)

    if annotations:
        matrixheight = unit * N / 2
        annotationheight = matrixheight/(4.25 + len(annotations)) * len(annotations)
        fig = plt.figure(figsize = (unit * N, matrixheight + annotationheight))
        grid = gs.GridSpec(len(annotations) + 1, 1, height_ratios = [10] + [1] * len(annotations))
        axs = [fig.add_subplot(grid[i]) for i in range(len(annotations) + 1)]
        xax = axs[-1]

        increment = N / (xrange[1] - xrange[0])
        annoaxnum = len(annotations) - 1
        for annoax, (track, ylabel, alternating, plottype, ylim) in zip(axs[-len(annotations):], annotations):
            annoax.set_ylim(bottom=ylim[0], top=ylim[1])
            annoax.set_yticks([])
            annoax.set_ylabel(ylabel)

            if plottype == 'Line2D':
                addAnnotationLine2D(annoax, track, increment, scales[scale], xrange[0], xrange[1], alternating)

            elif plottype == 'Marker':
                addAnnotationMarker(annoax, track, increment, scales[scale], xrange[0], xrange[1])

            elif plottype == 'bigwig':
                addBigWigTrack(annoax,
                               track,
                               chrom,
                               xrange[0] * scales[scale],
                               xrange[1] * scales[scale],
                               0, N)

            else:
                raise Exception("plottype not supported")

            if annoaxnum > 0:
                for loc in ['left', 'top', 'right', 'bottom']:
                    annoax.spines[loc].set_visible(False)

                annoax.set_xticks([])
                annoax.set_xlim((0, N))

            else:
                for loc in ['left', 'top', 'right']:
                    annoax.spines[loc].set_visible(False)

            annoaxnum -= 1

    else:
        fig, axs = plt.subplots()
        axs = [axs]
        xax = axs[0]

    mesh = axs[0].pcolormesh(X, Y, np.flipud(C), cmap=cmap, vmin = vmin, vmax = vmax, zorder = 2)

    # mark bin containing capturesite in grey
    if capturebin:
        greymap = clr.LinearSegmentedColormap.from_list('greymap', ['Grey', 'Grey'], N = 256)
        capturemat = np.zeros(shape = mat.shape)
        capturemat[capturebin, :] = 1
        capturemat[:, capturebin] = 1
        capturemat = np.triu(capturemat)
        capturemat = np.ma.masked_array(capturemat, capturemat == 0)
        capturemesh = axs[0].pcolormesh(X, Y, np.flipud(capturemat), cmap = greymap, vmin=0, vmax=1, zorder = 3)

    if highlightbins:
        for hlstartbin, hlendbin, hlclr in highlightbins:
            hlmap = clr.LinearSegmentedColormap.from_list(hlclr, [hlclr, hlclr], N=256)
            hlmat = np.zeros(shape=mat.shape)
            if not hlendbin == None:
                hlidx = np.arange(hlstartbin, hlendbin + 1)
            else:
                hlidx = hlstartbin

            hlmat[hlidx, :] = 1
            hlmat[:, hlidx] = 1
            hlmat = np.triu(hlmat)
            hlmat = np.ma.masked_array(hlmat, hlmat == 0)
            hlmesh = axs[0].pcolormesh(X, Y, np.flipud(hlmat), cmap=hlmap, vmin=0, vmax=1, alpha = 0.2, zorder = 1)

    # draw outlines of triangle plot
    vertices = np.array([[0, 0], [N/2, N], [N, 0]])
    triangle = patches.Polygon(vertices, fill = False, edgecolor = 'black')
    axs[0].add_patch(triangle)

    xax.set_xlim(left = 0, right = N)
    axs[0].set_ylim(bottom = 0, top = N)

    for loc in ['left', 'right', 'top']:
        axs[0].spines[loc].set_visible(False)

    axs[0].set_yticks([])
    axs[0].set_xticks([])

    if not xticknum:
        xax.set_xticks([])

    else:
        xax.set_xticks(np.linspace(0, N, xticknum))
        xax.set_xticklabels(['{val:,}'.format(val = i) for i in np.linspace(xrange[0], xrange[1], xticknum, dtype = int)])

    if xlabel:
        xax.set_xlabel(xlabel)

    # plot colorbar
    rect = patches.Rectangle((N - cbarwidth, N/2), cbarwidth, N/2, fill = False, edgecolor = 'black')
    axs[0].add_patch(rect)

    cbarY = np.tile(np.linspace(N/2, N, cmap.N).reshape(-1, 1), 2)
    cbarX = np.tile(np.array([N - cbarwidth, N]), (cbarY.shape[0], 1))
    cbarmesh = axs[0].pcolormesh(cbarX, cbarY, np.linspace(0, 1, cmap.N).reshape(-1, 1), cmap = cmap, vmin = 0, vmax = 1)

    for y, cmapval in zip(np.linspace(N/2, N, 5), np.linspace(vmin, vmax, 5)):
        axs[0].add_line(Line2D([N - cbarwidth - 1, N - cbarwidth], [y, y], color = 'black', lw = mpl.rcParams['patch.linewidth']))
        axs[0].text(N - cbarwidth - 1.5, y, '{:.01f}'.format(cmapval), ha = 'right', va = 'center')

    axs[0].text(N + 1, 3*N/4, 'RPM', ha = 'left', va = 'center', rotation = 90)

    if not annotations:
        fig.set_figheight(unit * N / 2)
        fig.set_figwidth(unit * N)

    #axs[0].set_aspect(0.5)
    fig.tight_layout()

    return fig, axs


def plotProfile(profiledict,
                order,
                N,
                chrom,
                xrange,
                color,
                yranges = None,
                scale = 'kbp',
                capturebin = None,
                xlabel = None,
                annotations = None,
                xticknum = 0,
                profileheight = 2,
                figurewidth = 15):
    '''
    plotting function for profiles as bargraphs

    :param profiledict:     dictonary containing the profiles as iterable of values to plot {profname1: profiledata, profname2: profiledata, etc.}
    :param order:           list of ordered keys of the profiledict specifying the plotting order from top (listindex 0) to bottom
    :param N:               number of bins of the profile (i.e. how many entries each profiledata has)
    :param xrange:          tuple (start, end) of the region plotted
    :param color:           color of the bars
    :param yranges:         dictionary of tuples containing y-axis limits for the profiles in the form of {i: (low, high)}
                            where i corresponds to the listindex given in the order argument
    :param scale:           specifies the scale used for the x-axis annotation
    :param capturebin:      bin containing the capture site (if given data plotted there is masked and a black line is plotted for indication)
                            bin number is determined by the extent of the plotted loci divided by the bin size
    :param xlabel:          label of the x-axis
    :param annotations:     list of tuples containing five items:
                            (0. pandas.Dataframe with coordinates of annotations to add to the plot (has to have 4 columns chr, start, end, name),
                             1. the label of the annotation axis,
                             2. bool indicating annotations should be plotted alternating,
                             3. string indicating the type of plotting function to use, currently 'Line2D', 'Marker' and 'bigwig' are supported
                                in case of 'Marker' 2. has no effect
                             4. tuple of y-axis limits (ymin, ymax))
                            the annotation given at index -1 will be the bottom most annotation
    :param xticknum:        number of ticks to plot on the x-axis
    :param profileheigh:    height of a single profile in inch (determines the final figure height; if annotations are 1 inch is added per annotation axis)
    :param figurewidth:     width of the figure in inch

    :return:                figure, list of axes
    '''

    scales = {'kbp': 1000, 'Mbp': 1000000}
    xrange = (xrange[0]/scales[scale], xrange[1]/scales[scale])

    figureheight = profileheight * len(profiledict) + 1 * len(annotations) if annotations else profileheight * len(profiledict)
    fig = plt.figure(figsize = (figurewidth, figureheight))

    if annotations:
        grid = gs.GridSpec(len(profiledict) + len(annotations), 1, height_ratios=[3] * len(profiledict) + [1] * len(annotations))
        axs = [fig.add_subplot(grid[i]) for i in range(len(profiledict) + len(annotations))]
        increment = N / (xrange[1] - xrange[0])

        annoaxnum = len(annotations) - 1
        for annoax, (track, ylabel, alternating, plottype, ylim) in zip(axs[-len(annotations):], annotations):
            annoax.set_ylim(bottom=0, top=1)
            annoax.set_yticks([])
            annoax.set_ylabel(ylabel)

            if plottype == 'Line2D':
                addAnnotationLine2D(annoax, track, increment, scales[scale], xrange[0], xrange[1], alternating)

            elif plottype == 'Marker':
                addAnnotationMarker(annoax, track, increment, scales[scale], xrange[0], xrange[1])

            elif plottype == 'bigwig':
                addBigWigTrack(annoax,
                               track,
                               chrom,
                               xrange[0] * scales[scale],
                               xrange[1] * scales[scale],
                               0, N)

            else:
                raise Exception("plottype not supported")

            if annoaxnum > 0:
                annoax.spines['bottom'].set_visible(False)

            annoaxnum -= 1

    else:
        grid = gs.GridSpec(len(profiledict), 1, height_ratios= [3] * len(profiledict))
        axs = [fig.add_subplot(grid[i]) for i in range(len(profiledict))]

    for ax in axs:
        ax.set_xticks([])
        ax.set_xlim((0, N))

    if not xticknum:
        axs[-1].set_xticks([])

    else:
        axs[-1].set_xticks(np.linspace(0, N, xticknum))
        axs[-1].set_xticklabels(['{val:,}'.format(val = i) for i in np.linspace(xrange[0], xrange[1], xticknum, dtype = int)])

    if yranges:
        for i, yrange in yranges.items():
            axs[i].set_ylim(yrange)
            axs[i].set_yticks(yrange)

    else:
        for ax in axs[:-1]:
            ax.set_yticks(ax.get_ylim())

    if xlabel:
        axs[-1].set_xlabel(xlabel)

    for ax in axs:
        for loc in ['left', 'top', 'right']:
            ax.spines[loc].set_visible(False)

    for i, profilename in enumerate(order):
        axs[i].bar(range(N), profiledict[profilename], align = 'edge', color = color, width = 1)
        axs[i].set_ylabel(profilename)

        if capturebin:
            axs[i].bar(capturebin + 0.5, axs[i].get_ylim()[1], align = 'center', width = 0.75, color = 'black')

    fig.tight_layout()
    return fig, axs


def plotProfileOverlay(profiledict,
                       N,
                       chrom,
                       xrange,
                       colors,
                       yrange = None,
                       scale = 'kbp',
                       capturebin = None,
                       ylabel = None,
                       xlabel = None,
                       annotations = None,
                       xticknum = 0,
                       profileheight = 5,
                       figurewidth = 15):
    scales = {'kbp': 1000, 'Mbp': 1000000}
    xrange = (xrange[0]/scales[scale], xrange[1]/scales[scale])

    if annotations:
        figureheight = profileheight + 1 * len(annotations) if annotations else profileheight
        fig = plt.figure(figsize=(figurewidth, figureheight))
        grid = gs.GridSpec(1 + len(annotations), 1, height_ratios=[5] + [1] * len(annotations))
        axs = [fig.add_subplot(grid[i]) for i in range(len(annotations) + 1)]
        increment = N / (xrange[1] - xrange[0])

        annoaxnum = len(annotations) - 1
        for annoax, (track, ylabel, alternating, plottype, ylim) in zip(axs[-len(annotations):], annotations):
            annoax.set_ylim(bottom=0, top=1)
            annoax.set_yticks([])
            annoax.set_ylabel(ylabel)

            if plottype == 'Line2D':
                addAnnotationLine2D(annoax, track, increment, scales[scale], xrange[0], xrange[1], alternating)

            elif plottype == 'Marker':
                addAnnotationMarker(annoax, track, increment, scales[scale], xrange[0], xrange[1])

            elif plottype == 'bigwig':
                addBigWigTrack(annoax,
                               track,
                               chrom,
                               xrange[0] * scales[scale],
                               xrange[1] * scales[scale],
                               0, N)

            else:
                raise Exception("plottype not supported")

            if annoaxnum > 0:
                annoax.spines['bottom'].set_visible(False)

            annoaxnum -= 1

    else:
        fig, ax = plt.subplots()
        axs = [ax]

    if ylabel:
        axs[0].set_ylabel(ylabel)

    for ax in axs:
        ax.set_xticks([])
        ax.set_xlim((0, N))

    if not xticknum:
        axs[-1].set_xticks([])

    else:
        axs[-1].set_xticks(np.linspace(0, N, xticknum))
        axs[-1].set_xticklabels(['{val:,}'.format(val = i) for i in np.linspace(xrange[0], xrange[1], xticknum, dtype = int)])

    if yrange:
        axs[0].set_ylim(yrange)
        axs[0].set_yticks(yrange)

    else:
        for ax in axs[:-1]:
            ax.set_yticks(ax.get_ylim())

    if xlabel:
        axs[-1].set_xlabel(xlabel)

    for ax in axs:
        for loc in ['left', 'top', 'right']:
            ax.spines[loc].set_visible(False)

    for color, (profilename, profile) in zip(colors, profiledict.items()):
        axs[0].plot(np.arange(N) + 0.5, profile, color = color, label = profilename)

        if capturebin:
            axs[0].bar(capturebin + 0.5, axs[0].get_ylim()[1], align = 'center', width = 0.75, color = 'black')

    axs[0].legend(loc = 'upper right')
    fig.tight_layout()
    return fig, axs


def computeAverageMatrix(matrices):
    summed = np.zeros(shape = matrices[0].shape)
    for mat in matrices:
        summed += mat

    return summed/len(matrices)


def makeDiffMatrix(mat1, mat2):
    sum1 = mat1.sum()
    sum2 = mat2.sum()

    if sum1 > sum2:
        mat1 = mat1 * sum2/sum1

    else:
        mat2 = mat2 * sum1/sum2

    return mat1 - mat2


def getBinIndex(captureSiteStart, leftBound, rightBound, binsize):
    binbounds = np.arange(leftBound, rightBound, binsize)
    # -1 because 0-based indices
    return len(np.where(binbounds < captureSiteStart)[0]) - 1


def getHighlightBinArgumentFromAnnotation(annotation, features, leftBound, rightBound, binsize, hlcolor = 'cyan'):
    '''
    takes an annotation table, a list of names of features in the table, the boundaries and binsize of the plot region
    and returns the argumentlist that can be passed to highlightbins in plotMatrix

    :param annotation:  pd.DataFrame containing the annotations
    :param features:    names of the features in annotation that should be highlighted
    :param leftBound:   left boundary position of the plotting range
    :param rightBound:  right boundary position of the plotting range
    :param binsize:     size of bins plotted
    :param hlcolor:     color to use to highlight or list of colors if list has to be the same length as features

    :return:            highlightbins argument list
    '''

    hlargument = []
    features_df = annotation.loc[annotation.name.isin(features), :] \
                            .reset_index()

    if isinstance(hlcolor, list):
        if len(hlcolor) == len(features):
            for i, feature in features_df.iterrows():
                hlstartbin = getBinIndex(feature['start'], leftBound, rightBound, binsize)
                hlendbin = getBinIndex(feature['end'], leftBound, rightBound, binsize)
                hlargument.append((hlstartbin, hlendbin, hlcolor[i]))
        else:
            raise Exception('numbers of color has to match the number of features to highlight if multiple colors are passed')

    else:
        for i, feature in features_df.iterrows():
            hlstartbin = getBinIndex(feature['start'], leftBound, rightBound, binsize)
            hlendbin = getBinIndex(feature['end'], leftBound, rightBound, binsize)
            hlargument.append((hlstartbin, hlendbin, hlcolor))

    return hlargument


if __name__ == '__main__':

    wyorb = clr.LinearSegmentedColormap.from_list('wyorb', ['White', 'Yellow', 'Orange', 'Red', 'Black'], N=256)
    gyorb = clr.LinearSegmentedColormap.from_list('gyorb', ['lightgrey', 'Yellow', 'Orange', 'Red', 'Black'], N=256)
    bwr = plt.get_cmap('bwr')

    for prefix, capture, capturebin, capturebin3rr, features2highlight, profileymax in zip(['TriC_7', 'TriC_5', 'TriC_6', 'TriC_3', 'TriC_4'],
                                                                                           ['EMu_2', 'IgG1_2', 'HS4_2', 'EMu', 'IgG1'],
                                                                                           [229, 142, 31, 229, 142],
                                                                                           [None, None, 31, None, None],
                                                                                           [['3a', '3b', '1,2', 'Ig1', 'Em'],
                                                                                            ['4', '3a', '3b', '1,2', 'Ig1'],
                                                                                            ['4', '3a', '3b', '1,2', 'Em']],
                                                                                           [200, 200, 200, 200, 200]):

        cells = ['mESC', 'priB_d0', 'priB_d2']
        samplenames = ['_'.join([prefix, suffix]) for suffix in ['_'.join(t) for t in it.product(cells, ['1', '2', '3'])]]

        chrom, leftBound, rightBound, binsize = 'chr12', 114435000, 114669000, 1000
        n_bins = (rightBound - leftBound) // binsize
        geneannotation = pd.read_csv('../TriC/vdjanno_new.sort.bed', sep = '\t', header = None, names = ['chr', 'start', 'end', 'name'])
        REsiteannotation = pd.read_csv('../TriC/vdjREs.bed', sep = '\t', header = None, usecols = [0, 1, 2], names = ['chr', 'start', 'end'])
        REsiteannotation['name'] = ''
        annotations = [(geneannotation, 'Genes', False, 'Line2D', (0, 1)),
                       (REsiteannotation, 'NlaIII', False, 'Marker', (0, 1)),
                       ('../TriC/mm9_mappability.bw', 'mappability', None, 'bigwig', (0, 1))]
        highlightbins = getHighlightBinArgumentFromAnnotation(geneannotation, features2highlight, leftBound, rightBound, binsize, 'cyan')

        replicates = {name: np.loadtxt('../TriC/matrices/{0}_TriC_interactions_1000_RAW.tab'.format(name)) for name in samplenames}

        averageMats = {}
        for cell in cells:
            averageMats[cell] = computeAverageMatrix([replicates[f'{prefix}_{cell}_{i}'] for i in range(1, 4)])
            fig, axs = plotMatrix(averageMats[cell],
                                  gyorb,
                                  chrom,
                                  (leftBound, rightBound),
                                  capturebin = capturebin,
                                  highlightbins = highlightbins,
                                  xlabel = f'{chrom} [kbp]',
                                  annotations = annotations,
                                  xticknum = 10,
                                  vmax=50)
            fig.savefig(f'../TriC/TriC_{capture}_{cell}_contacts.pdf')
            plt.close(fig)

        differenceMats = {}
        for cell in cells:
            if cell == 'priB_d0':
                differenceMats[cell + 'mESC'] = makeDiffMatrix(averageMats[cell], averageMats['mESC'])

            elif cell == 'priB_d2':
                differenceMats[cell + 'mESC'] = makeDiffMatrix(averageMats[cell], averageMats['mESC'])
                differenceMats[cell + 'priB_d0'] = makeDiffMatrix(averageMats[cell], averageMats['priB_d0'])


        for diffcells in [('priB_d0', 'mESC'), ('priB_d2', 'mESC'), ('priB_d2', 'priB_d0')]:
            fig, axs = plotMatrix(differenceMats[''.join(diffcells)],
                                  bwr,
                                  chrom,
                                  (leftBound, rightBound),
                                  capturebin = capturebin,
                                  highlightbins = highlightbins,
                                  xlabel = f'{chrom} [kbp]',
                                  annotations = annotations,
                                  xticknum = 10,
                                  vmin = -50,
                                  vmax=50)
            fig.savefig(f'../TriC/TriC_{capture}_{diffcells[0]}_{diffcells[1]}_diff.pdf')
            plt.close(fig)

        # plot 3'RR zoom
        RRrightBound = leftBound + 61 * binsize

        for cell in cells:
            zoom3RR = averageMats[cell][:61, :][:, :61]
            fig, axs = plotMatrix(zoom3RR,
                                  gyorb,
                                  chrom,
                                  (leftBound, RRrightBound),
                                  capturebin = capturebin3rr,
                                  xlabel=f'{chrom} [kbp]',
                                  annotations = annotations,
                                  xticknum = 4,
                                  vmax=50,
                                  unit = 0.2,
                                  cbarwidth = 3)
            fig.savefig(f'../TriC/TriC_{capture}_{cell}_contacts_3RR.pdf')
            plt.close(fig)

        for diffcells in [('priB_d0', 'mESC'), ('priB_d2', 'mESC'), ('priB_d2', 'priB_d0')]:
            zoom3RR = differenceMats[''.join(diffcells)][:61, :][:, :61]
            fig, axs = plotMatrix(zoom3RR,
                                  bwr,
                                  chrom,
                                  (leftBound, RRrightBound),
                                  capturebin = capturebin3rr,
                                  xlabel = f'{chrom} [kbp]',
                                  annotations = annotations,
                                  xticknum = 4,
                                  vmin = -50,
                                  vmax=50,
                                  unit = 0.2,
                                  cbarwidth = 3)
            fig.savefig(f'../TriC/TriC_{capture}_{diffcells[0]}_{diffcells[1]}_diff_3RR.pdf')
            plt.close(fig)

        # 2way interaction matrices
    #    twoWayMatrices = {}
    #    for sample in samplenames:
    #        mat, intervals = loadH5subset('../TriC/matrices/{0}_2way.h5'.format(sample), includechroms = ['chr12'])
    #        twoWayMatrices[sample] = getRegionalMatrix(mat, intervals['start_list'], 114435000, 114669000, csr = True)
    #
    #    average2wayMats = {}
    #    for cell in cells:
    #        tmpsum = np.zeros(shape = ((114669000 - 114435000)//1000,) * 2)
    #        for i in range(1, 4):
    #            tmpsum += twoWayMatrices['TriC_{0}_{1}'.format(cell, i)]
    #
    #        average2wayMats[cell] = tmpsum * 1000000 / (3 * tmpsum.sum())
    #        fig, axs = plotMatrix(average2wayMats[cell], gyorb, (114435000, 114669000), xlabel='chr12 [kbp]',
    #                              annotations = annotations, xticknum=10, vmax=100)
    #        fig.savefig('../TriC/{0}_2waycontacts.pdf'.format(cell))
    #        plt.close(fig)

        # plotting profiles
        profiles = {cell: {} for cell in cells}
        allprofile = loadProfileTab(f'../TriC/profiles/{prefix}_all.tsv', interval = (leftBound, rightBound))
        #allprofile.loc[capturebin, ~allprofile.columns.isin(['chr', 'start', 'end'])] = 0 #setting capture site counts to 0

        for origin in ['2way', '3plus', 'all']:
            profiletab = loadProfileTab(f'../TriC/profiles/{prefix}_{origin}.tsv', interval = (leftBound, rightBound))
            for cell in cells:
                norm = allprofile.filter(regex = f'{prefix}_{cell}_*').mean(axis = 1)
                totalnorm = 100000 / norm.sum()
                binnorm = 1000 / (norm * totalnorm).max()

                tmp = profiletab.filter(regex = f'{prefix}_{cell}_*').mean(axis = 1)
                tmp.loc[capturebin] = 0 #setting capture site counts to 0

                profiles[cell][origin] = tmp * totalnorm #* binnorm


        for cell in cells:
            print(cell, profiles[cell]['3plus'].max())
            fig, axs = plotProfile(profiles[cell],
                                   ['2way', '3plus', 'all'],
                                   n_bins,
                                   chrom,
                                   (leftBound, rightBound),
                                   'grey',
                                   yranges = {0: (0, profileymax), 1: (0, profileymax), 2: (0, profileymax)},
                                   capturebin = capturebin,
                                   xlabel = f'{chrom} [kbp]',
                                   annotations = annotations,
                                   xticknum = 10)

            fig.savefig(f'../TriC/{cell}_{capture}_profiles.pdf')

            plt.close(fig)

        # plot profile overlay
        for cell in cells:
            for origin in ['2way', '3plus', 'all']:
                profile = profiles[cell][origin]
                profile.loc[capturebin] = profile.loc[capturebin - 1]
                profiles[cell][origin] = profile

            tmpprofiles = {'2way': profiles[cell]['2way'], '3plus': profiles[cell]['3plus']}
            fig, axs = plotProfileOverlay(tmpprofiles,
                                          n_bins,
                                          chrom,
                                          (leftBound, rightBound),
                                          ('steelblue', 'gold'),
                                          (0, profileymax),
                                          capturebin = capturebin,
                                          xlabel = f'{chrom} [kbp]',
                                          annotations = annotations,
                                          xticknum = 10)

            fig.savefig(f'../TriC/{capture}_{cell}_profilesOverlay.pdf')
            plt.close(fig)