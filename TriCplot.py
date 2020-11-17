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
import logging
import argparse as ap
mpl.rcParams['pdf.fonttype'] = 42


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
    tab = pd.read_csv(annotation, sep = '\t')
    for i, locus in tab.iterrows():
        start, end = locus['start'] / scale, locus['end'] / scale
        x1 = xmin if start < xmin else (start - xmin) * increment
        x2 = xmax if end > xmax else (end - xmin) * increment
        x.append((x1 + x2)/2)

    ax.plot(x, [0.5] * len(x), '|')


def addAnnotationLine2D(ax, annotation, increment, scale, xmin, xmax, alternating = False):
    tab = pd.read_csv(annotation, sep = '\t')
    for i, locus in tab.iterrows():
        start, end = locus['start'] / scale, locus['end'] / scale
        x1 = xmin if start < xmin else (start - xmin) * increment
        x2 = xmax if end > xmax else (end - xmin) * increment

        if alternating:
            if i % 2 == 0:
                ax.add_line(Line2D([x1, x2], [0.425, 0.425], lw=5, solid_capstyle = 'butt'))

            else:
                ax.add_line(Line2D([x1, x2], [0.575, 0.575], lw=5, solid_capstyle = 'butt'))

        else:
            ax.add_line(Line2D([x1, x2], [0.5, 0.5], lw=5, solid_capstyle = 'butt'))

        if locus['name']:
            if i % 2 == 0:
                ax.text((x2 - x1) / 2 + x1, 0.25 if alternating else 0.325, locus['name'], ha='center', va='top', fontsize = 10)

            else:
                ax.text((x2 - x1) / 2 + x1, 0.75 if alternating else 0.675, locus['name'], ha='center', va='bottom', fontsize = 10)


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


def plotAnnotation(ax,
                   track,
                   ylabel,
                   alternating,
                   plottype,
                   ylim,
                   number_of_bins,
                   start,
                   end,
                   chrom,
                   scaling = 1000,
                   xticks = False):
    increment = number_of_bins / (end - start)
    ax.set_ylim(bottom=ylim[0], top=ylim[1])
    ax.set_yticks([])
    ax.set_ylabel(ylabel)

    if plottype == 'Line2D':
        addAnnotationLine2D(ax, track, increment, scaling, start, end, alternating)

    elif plottype == 'Marker':
        addAnnotationMarker(ax, track, increment, scaling, start, end)

    elif plottype == 'bigwig':
        addBigWigTrack(ax,
                       track,
                       chrom,
                       start * scaling,
                       end * scaling,
                       0,
                       number_of_bins)

    else:
        raise Exception("plottype not supported")

    if not xticks:
        for loc in ['left', 'top', 'right', 'bottom']:
            ax.spines[loc].set_visible(False)

        ax.set_xticks([])

    else:
        for loc in ['left', 'top', 'right']:
            ax.spines[loc].set_visible(False)

        ax.set_xlabel(chrom)

    ax.set_xlim((0, number_of_bins))

    return ax


def plotMatrix(ax,
               mat,
               cmap,
               xrange,
               scaling = 1000,
               capturebins = None,
               highlightbins = None,
               xlabel = None,
               xticknum = 0,
               cbarwidth = 5,
               vmin = 0,
               vmax = 50,
               mirror_horizontal = False,
               subplot_label = None):
    '''
    plotting function for triC results

    :param ax:                  plt.Axes object to generate the plot in
    :param mat:                 TriC matrix as generated by the CCseq pipeline
    :param cmap:                colormap to use for plotting
    :param xrange:              tuple of integer coordinates in bp of the genomic region plotted in the matrix
    :param scale:               determines the scale used (either 'none', 'kbp' or 'Mbp')
    :param capturebins:         bins containing the capture probes (optional)
    :param highlightbins:       list of tuples of startbin, endbin and highlightcolor for bins that should be highlighted
                                if endbin == startbin only the startbin will be highlighted
    :param xlabel:              xaxis label
                                the annotation given at index -1 will be the bottom most annotation
    :param xticknum:            number of xticks to plot. xticklabels will be interpolated with np.linspace
    :param cbarwidth:           width of the colorbar
    :param vmin:                minimum value of the colormap
    :param vmax:                maximum value of the colormap
    :param scaling:             divisor for scaling the xrange
    :param mirror_horizontal:   indicates if generated matrix plot should be mirrored at a horizontal line

    :return:                    plt.Axes
    '''

    xrange = (xrange[0]/scaling, xrange[1]/scaling)

    N = mat.shape[0]
    # Get the lower triangle of the matrix.
    C = np.triu(mat)
    # Mask the upper triangle
    C = np.ma.masked_array(C, C == 0)

    # Transformation matrix for rotating the heatmap.
    A = np.array([(y, x) for x in range(N, -1, -1) for y in range(N + 1)])
    t = np.array([[1, 0.5], [-1, 0.5]]) if not mirror_horizontal else np.array([[-1, 0.5], [1, 0.5]])
    A = np.dot(A, t)

    # Plot the heatmap triangle.
    X = A[:, 1].reshape(N + 1, N + 1)
    Y = A[:, 0].reshape(N + 1, N + 1)

    mesh = ax.pcolormesh(X, Y, np.flipud(C), cmap=cmap, vmin = vmin, vmax = vmax, zorder = 2)

    # mark bin containing capturesite in grey
    if capturebins:
        greymap = clr.LinearSegmentedColormap.from_list('greymap', ['Grey', 'Grey'], N = 256)
        capturemat = np.zeros(shape = mat.shape)
        for capturebin in capturebins:
            capturemat[capturebin, :] = 1
            capturemat[:, capturebin] = 1

        capturemat = np.triu(capturemat)
        capturemat = np.ma.masked_array(capturemat, capturemat == 0)
        capturemesh = ax.pcolormesh(X, Y, np.flipud(capturemat), cmap = greymap, vmin=0, vmax=1, zorder = 3)

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
            hlmesh = ax.pcolormesh(X, Y, np.flipud(hlmat), cmap=hlmap, vmin=0, vmax=1, alpha = 0.2, zorder = 1)

    # draw outlines of triangle plot
    vertices = np.array([[0, 0], [N/2, N], [N, 0]]) if not mirror_horizontal else np.array([[0, 0], [N/2, -N], [N, 0]])
    triangle = patches.Polygon(vertices, fill = False, edgecolor = 'black')
    ax.add_patch(triangle)

    ax.set_xlim(left = 0, right = N)

    if not mirror_horizontal:
        ax.set_ylim(bottom = 0, top = N)

    else:
        ax.set_ylim(bottom = -N, top = 0)

    if not mirror_horizontal:
        for loc in ['left', 'right', 'top']:
            ax.spines[loc].set_visible(False)

    else:
        for loc in ['left', 'right', 'bottom']:
            ax.spines[loc].set_visible(False)

        ax.xaxis.tick_top()

    ax.set_yticks([])
    ax.set_xticks([])

    if not xticknum:
        ax.set_xticks([])

    else:
        ax.set_xticks(np.linspace(0, N, xticknum))
        ax.set_xticklabels(['{val:,}'.format(val = i) for i in np.linspace(xrange[0], xrange[1], xticknum, dtype = int)])

    if xlabel:
        ax.set_xlabel(xlabel)

    # plot colorbar
    rect = patches.Rectangle((N - cbarwidth, N/2), cbarwidth, N/2, fill = False, edgecolor = 'black') \
           if not mirror_horizontal else \
           patches.Rectangle((N - cbarwidth, -N), cbarwidth, N / 2, fill=False, edgecolor='black')

    ax.add_patch(rect)


    cbarY = np.tile(np.linspace(N/2, N, cmap.N).reshape(-1, 1), 2) \
            if not mirror_horizontal else \
            np.tile(np.linspace(-N, -N/2, cmap.N).reshape(-1, 1), 2)
    cbarX = np.tile(np.array([N - cbarwidth, N]), (cbarY.shape[0], 1))
    cbarmesh = ax.pcolormesh(cbarX, cbarY, np.linspace(0, 1, cmap.N).reshape(-1, 1), cmap = cmap, vmin = 0, vmax = 1)

    ys = np.linspace(N/2, N, 5) if not mirror_horizontal else np.linspace(-N, -N/2, 5)
    for y, cmapval in zip(ys, np.linspace(vmin, vmax, 5)):
        ax.add_line(Line2D([N - cbarwidth - 1, N - cbarwidth], [y, y], color = 'black', lw = mpl.rcParams['patch.linewidth']))
        ax.text(N - cbarwidth - 1.5, y, '{:.01f}'.format(cmapval), ha = 'right', va = 'center')

    ax.text(N + 1, 3*N/4 if not mirror_horizontal else -3*N/4, 'RPM', ha = 'left', va = 'center', rotation = 90)

    if subplot_label:
        ax.text(0, N if not mirror_horizontal else -N, subplot_label,
                ha = 'left',
                va = 'top' if not mirror_horizontal else 'bottom', fs = 20)

    return ax


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


def plotProfileOverlay(ax,
                       profiledict,
                       number_of_bins,
                       xrange,
                       colors,
                       yrange = None,
                       scaling = 1000,
                       capturebins = None,
                       ylabel = None,
                       xlabel = None,
                       xticknum = 0):
    xrange = (xrange[0]/scaling, xrange[1]/scaling)

    if ylabel:
        ax.set_ylabel(ylabel)

    ax.set_xticks([])
    ax.set_xlim((0, number_of_bins))

    if not xticknum:
        ax.set_xticks([])

    else:
        ax.set_xticks(np.linspace(0, number_of_bins, xticknum))
        ax.set_xticklabels(['{val:,}'.format(val = i) for i in np.linspace(xrange[0], xrange[1], xticknum, dtype = int)])

    if yrange:
        ax.set_ylim(yrange)
        ax.set_yticks(yrange)

    if xlabel:
        ax.set_xlabel(xlabel)

    for loc in ['left', 'top', 'right']:
        ax.spines[loc].set_visible(False)

    for color, (profilename, profile) in zip(colors, profiledict.items()):
        ax.plot(np.arange(number_of_bins) + 0.5, profile, color = color, label = profilename)

        if capturebins:
            for capturebin in capturebins:
                ax.bar(capturebin + 0.5, axs[0].get_ylim()[1], align = 'center', width = 0.75, color = 'black')

    ax.legend(loc = 'upper right')
    return ax


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


def same_length(*args):
    l = len(args[0])
    for arg in args[1:]:
        if len(arg) != l:
            return False

    return True


def load_profiles(treatment_profile, control_profile, treatment_label, control_label, leftBound, rightBound, capturebins):
    profiles = {}
    for k, file in zip([treatment_label, control_label], [treatment_profile, control_profile]):
        profiletab = loadProfileTab(file, interval=(leftBound, rightBound))
        meanprofile = profiletab.loc[:, ~profiletab.columns.isin(['chr', 'start', 'end'])].mean(axis=1)
        totalnorm = 100000 / meanprofile.sum()
        #binnorm = 1000 / (meanprofile * totalnorm).max()
        for capturebin in capturebins:
            meanprofile.loc[capturebin] = 0  # setting capture site counts to 0

        profiles[k] = meanprofile * totalnorm  # * binnorm

    return profiles


wyorb = clr.LinearSegmentedColormap.from_list('wyorb', ['White', 'Yellow', 'Orange', 'Red', 'Black'], N=256)
gyorb = clr.LinearSegmentedColormap.from_list('gyorb', ['lightgrey', 'Yellow', 'Orange', 'Red', 'Black'], N=256)
bwr = plt.get_cmap('bwr')
#chrom, leftBound, rightBound, binsize = 'chr12', 114435000, 114669000, 1000

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('--treatment', '-t', nargs = '+',
                    help = 'one or more TriC interaction files for treatment condition'
                            'if more than one are given a mean is computed over all samples')
parser.add_argument('--treatment_label', '-tl',
                    help = 'name of the treatment condition to use for plot naming etc.'
                           'if not given the script tries to infer a name from the filename')
parser.add_argument('--control', '-c', nargs = '+',
                    help = 'one or more TriC interaction files for control condition'
                           'if more than one are given a mean is computed over all samples')
parser.add_argument('--control_label', '-control',
                    help = 'name of the control condition to use for plot naming etc.'
                           'if not given the script tries to infer a name from the filename')
parser.add_argument('--region', '-r', required = True,
                    help = 'region to consider for plotting in the format of chr:start-stop')
parser.add_argument('--capture_bins', required = True,
                    help = 'oligo file used for processing with CCseqBasic denoting the bins in the matrix and profile to mask')
parser.add_argument('--binsize', '-bs', default = 1000, type = int,
                    help = 'size of the bins the region is subdivided into')
parser.add_argument('--annotation', '-a', nargs = '*',
                    help = 'space-separated list of filenames that should be added to the figure'
                           'currently only bed and bigwig files are supported. Bed files can be in BED3 or BED4 format'
                           'which means chr, start, end or chr, start, end, name columns present with header.'
                           'if name is present the plotted annotations will be named')
parser.add_argument('--annotation_drawstyle', '-ad', nargs = '*',
                    help = 'space-separated list of drawstyles for each annotation file'
                           'in the same order as the files (has to be one of Line2D, Marker or bigwig'
                           'Line2D is usually use for bedfiles with small numbers of regions to display'
                           'e.g. genes in the locus, Marker is used for bedfiles with a large number of regions'
                           'e.g. NlaIII restriction sites and bigwig is reserved for bigwig files)')
parser.add_argument('--annotation_labels', '-al', nargs = '*',
                    help = 'names of the annotations used to name the generated plots'
                           'in the same order as the files given')
parser.add_argument('--annotation_yMin', nargs = '*',
                    help = 'minimum value of y axis of the annotation plots')
parser.add_argument('--annotation_yMax', nargs = '*',
                    help = 'maximum value of y axis of the annotation plots')
parser.add_argument('--alternating', nargs = '*',
                    help = 'space-separated list of 0 and 1 indicating if the the given annotation'
                           'should be drawn in an alternating fashing. Most suitable for dense annotations'
                           'to avoid label overlap')
parser.add_argument('--highlight_annotation', type = int,
                    help = 'the 1-based --annotation list index of the annotation to use for highlighting')
parser.add_argument('--highlight_feature', nargs = '*',
                    help = 'space-separated list of features in the --highlight_annotation to use for highlighting')
parser.add_argument('--treatment_3plus',
                    help = 'tab-separated table containing counts of reads with 3 or more valid restriction fragments in the treatment condition')
parser.add_argument('--control_3plus',
                    help = 'tab-separated table containing counts of reads with 3 or more valid restriction fragments in the control condition')
parser.add_argument('--profile_yMax',
                    help = 'maximum value of the y axis of the profile plot')
parser.add_argument('--figwidth', default = 6, type = float,
                    help = 'width of the generated figure. Height is computed accordingly.')
parser.add_argument('--outputFilePrefix', '-o', required = True,
                    help = 'prefix to use for output files')
args = parser.parse_args()

chrom, bounds = ':'.split(args.region)
leftBound, rightBound = [int(i) for i in '-'.split(bounds)]
n_bins = (rightBound - leftBound) // args.binsize

annotations = []
number_of_annotation_axes = 0
if args.annotations:
    assert same_length(args.annotation,
                       args.annotation_drawstyle,
                       args.annotation_labels,
                       args.annotation_yMin,
                       args.annotation_yMax,
                       args.alternating), \
           'If annotations are given, all arguments must have the same number of parameters'

    annotations = zip(args.annotation,
                      args.annotation_labels,
                      args.alternating,
                      args.annotation_drawstyle,
                      zip(args.annotation_yMin, args.annotation_yMax),
                      [False] * (len(args.annotations) - 1) + [True])       # xtick visibility

    number_of_annotation_axes = len(annotations)

treatments = [np.loadtxt(file) for file in args.treatment]
treatment_avg = computeAverageMatrix(treatments)

controls = [np.loadtxt(file) for file in args.control]
control_avg = computeAverageMatrix(controls)

# setting up figure layout
plt.rcParams['figure.constrained_layout.use'] = True
fig1 = plt.figure(dpi = 300)
matrix_subplot_height = args.figwidth / 2
annotation_height = 0.75
profile_height = 2.25
hspace = 0.2

profile_args = [args.treatment_3plus, args.control_3plus, args.profile_yMax]
if any(profile_args):
    assert all(profile_args), 'all profile arguments have to be set if one is used'
    number_of_annotation_axes += 1
    fig1_height = 2 * matrix_subplot_height + \
                  annotation_height * (number_of_annotation_axes - 1) + \
                  profile_height + \
                  hspace * (number_of_annotation_axes + 1)
    number_of_axes_fig1 = number_of_annotation_axes + 2
    height_ratios_fig1 = [matrix_subplot_height / fig1_height, profile_height / fig1_height] + \
                         [annotation_height / fig1_height] * (number_of_annotation_axes - 1) + \
                         [matrix_subplot_height / fig1_height]

    number_of_axes_fig2 = number_of_annotation_axes + 1
    fig2_height = fig1_height - matrix_subplot_height
    height_ratios_fig2 = [matrix_subplot_height / fig2_height, profile_height / fig2_height] + \
                         [annotation_height / fig2_height] * (number_of_annotation_axes - 1)

else:
    fig1_height = 2 * matrix_subplot_height + \
                  annotation_height * (number_of_annotation_axes) + \
                  hspace * (number_of_annotation_axes + 1)
    number_of_axes_fig1 = number_of_annotation_axes + 2
    height_ratios_fig1 = [matrix_subplot_height / fig1_height] + \
                         [annotation_height / fig1_height] * (number_of_annotation_axes) + \
                         [matrix_subplot_height / fig1_height]

    number_of_axes_fig2 = number_of_annotation_axes + 1
    fig2_height = fig1_height - matrix_subplot_height
    height_ratios_fig2 = [matrix_subplot_height / fig2_height] + \
                         [annotation_height / fig2_height] * (number_of_annotation_axes)

# compare figure
fig1.set_figwidth(args.figwidth)
fig1.set_figheight(fig1_height)
gridspec1 = gs.GridSpec(number_of_axes_fig1,
                        1,
                        height_ratios = height_ratios_fig1,
                        figure = fig1)
treatment_ax = fig1.add_subplot(gridspec1[0])
control_ax = fig1.add_subplot(gridspec1[-1])
profile_ax1 = fig1.add_subplot(gridspec1[1]) if any(profile_args) else None
annotation_axs1 = [fig1.add_subplot(gridspec1[i + 2]) for i in range(len(annotations))] if annotations else []

# diff figure
fig2 = plt.figure(dpi = 300)
fig2.set_figwidth(args.figwidth)
fig2.set_figheight(fig1_height)
gridspec2 = gs.GridSpec(number_of_axes_fig2,
                        1,
                        height_ratios = height_ratios_fig2,
                        figure = fig2)
diff_ax = fig2.add_subplot(gridspec2[0])
profile_ax2 = fig2.add_subplot(gridspec2[1]) if any(profile_args) else None
annotation_axs2 = [fig2.add_subplot(gridspec2[i + 2]) for i in range(len(annotations))] if annotations else []

oligotab = pd.read_csv(args.capture_bins,
                      sep = '\t',
                      header = None,
                      names = ['name', 'chrom', 'start', 'end'],
                      usecols = [0, 1, 2, 3])
capturebins = [getBinIndex(r['start'], leftBound, rightBound, args.binsize) for i, r in oligotab.iterrows()]

highlightbins = []
if args.highlight_annotation:
    highlightbins = getHighlightBinArgumentFromAnnotation(annotations[args.highlight_annotation - 1],
                                                          args.highlight_features,
                                                          leftBound,
                                                          rightBound,
                                                          args.binsize,
                                                          'cyan')

treatment_ax = plotMatrix(treatment_ax,
                          treatment_avg,
                          gyorb,
                          (leftBound, rightBound),
                          capturebins=capturebins,
                          highlightbins=highlightbins,
                          xticknum = 10,
                          vmax = 50,
                          subplot_label = args.treatment_label)

control_ax = plotMatrix(control_ax,
                        control_avg,
                        gyorb,
                        (leftBound, rightBound),
                        capturebins = capturebins,
                        highlightbins = highlightbins,
                        xticknum = 10,
                        vmax = 50,
                        mirror_horizontal = True,
                        subplot_label = args.control_label)

diff_ax = plotMatrix(diff_ax,
                     makeDiffMatrix(treatment_avg, control_avg),
                     bwr,
                     (leftBound, rightBound),
                     capturebins=capturebins,
                     highlightbins=highlightbins,
                     xticknum=10,
                     vmax=50,
                     mirror_horizontal=True,
                     subplot_label='-'.join((args.treatment_label, args.control_label)))

if any(profile_args):
    profiles = load_profiles(args.treatment_3plus, args.control_3plus,
                             args.treatment_label, args.control_label,
                             leftBound, rightBound,
                             capturebins)

    for profile_ax in [profile_ax1, profile_ax2]:
        profile_ax = plotProfileOverlay(profile_ax,
                                        profiles,
                                        n_bins,
                                        (leftBound, rightBound),
                                        yrange = args.profile_yMax,
                                        capturebins = capturebins,
                                        xticknum = 10,
                                        colors = ('steelblue', 'gold'))

if annotations:
    for annotation_axs in [annotation_axs1, annotation_axs2]:
        for ax, (track, label, alternating, plottype, ylim, visible_xticks) in zip(annotation_axs, annotations):
            ax = plotAnnotation(ax,
                                track,
                                label,
                                alternating,
                                plottype,
                                ylim,
                                n_bins,
                                leftBound,
                                rightBound,
                                chrom,
                                xticks = visible_xticks)


fig1.savefig(args.outputFilePrefix + '_compare.pdf')
fig2.savefig(args.outputFilePrefix + '_difference.pdf')