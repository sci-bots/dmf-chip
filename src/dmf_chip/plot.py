# coding: utf-8
from __future__ import absolute_import, unicode_literals, print_function
from collections import OrderedDict

import matplotlib as mpl
import matplotlib.collections
import matplotlib.patches
import matplotlib.pyplot as plt
import pandas as pd
import svg_model

from .core import (DEFAULT_DISTANCE_THRESHOLD, ELECTRODES_XPATH,
                   _extract_electrode_channels, get_all_intersections,
                   get_segments)

__all__ = ['draw', 'draw_w_segments']


def draw(svg_source, ax=None, labels=True):
    '''
    Draw the specified device, along with rays casted normal to the electrode
    line segments that intersect with a line segment of a neighbouring
    electrode.


    Parameters
    ----------
    svg_source : str or file-like or pandas.DataFrame
        File path, URI, or file-like object for SVG device file.

        If specified as ``pandas.DataFrame``, assume argument is in format
        returned by :func:`svg_model.svg_shapes_to_df`.
    ax : matplotlib.axes._subplots.AxesSubplot, optional
        Axis to draw on.
    labels : bool, optional
        Draw channel labels (default: ``True``).

    Returns
    -------
    `dict`
        Result `dict` includes::
        - ``axis``: axis to which the device was drawn
          (`matplotlib.axes._subplots.AxesSubplot`)
        - ``df_shapes``: table of electrode shape vertices (`pandas.DataFrame`)
        - ``electrode_channels``: mapping from channel number to electrode ID
          (`pandas.Series`).
        - ``channel_patches``: mapping from channel number to corresponding
          `matplotlib` electrode `Patch`.  May be used, e.g., to set color and
          alpha (`pandas.Series`).
    '''
    if not isinstance(svg_source, pd.DataFrame):
        df_shapes = svg_model.svg_shapes_to_df(svg_source,
                                               xpath=ELECTRODES_XPATH)
    else:
        df_shapes = svg_source
    electrode_channels = _extract_electrode_channels(df_shapes)

    # Compute center `(x, y)` for each electrode.
    electrode_centers = df_shapes.groupby('id')[['x', 'y']].mean()
    # Index by **channel number** instead of **electrode id**.
    electrode_centers.index = electrode_channels.reindex(electrode_centers
                                                         .index)

    patches = OrderedDict(sorted([(id_, mpl.patches
                                   .Polygon(df_shape_i[['x', 'y']].values,
                                            closed=False, label=id_))
                                  for id_, df_shape_i in
                                  df_shapes.groupby('id')]))
    channel_patches = pd.Series(patches.values(), index=electrode_channels
                                .reindex(patches.keys()))

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.set_aspect(True)

    # For colors, see: https://gist.github.com/cfobel/fd939073cf13a309d7a9
    for patch_i in patches.values():
        # Light blue
        patch_i.set_facecolor('#88bde6')
        # Medium grey
        patch_i.set_edgecolor('#4d4d4d')
        ax.add_patch(patch_i)

    ax.set_xlim(df_shapes.x.min(), df_shapes.x.max())
    ax.set_ylim(df_shapes.y.max(), df_shapes.y.min())

    if labels:
        for channel_i, center_i in electrode_centers.iterrows():
            ax.text(center_i.x, center_i.y, str(channel_i),
                    horizontalalignment='center', verticalalignment='center',
                    color='white', fontsize=10, bbox={'facecolor': 'black',
                                                      'alpha': 0.2, 'pad': 5})

    return {'axis': ax, 'electrode_channels': electrode_channels,
            'df_shapes': df_shapes, 'channel_patches': channel_patches}


def draw_w_segments(svg_source, ax=None,
                    distance_threshold=DEFAULT_DISTANCE_THRESHOLD):
    '''
    Draw the specified device, along with rays casted normal to the electrode
    line segments that intersect with a line segment of a neighbouring
    electrode.


    Parameters
    ----------
    svg_source : str or file-like or pandas.DataFrame
        File path, URI, or file-like object for SVG device file.

        If specified as ``pandas.DataFrame``, assume argument is in format
        returned by :func:`svg_model.svg_shapes_to_df`.
    ax : matplotlib.axes._subplots.AxesSubplot, optional
        Axis to draw on.
    distance_threshold : pint.quantity.Quantity, optional
        Maximum gap between electrodes to still be considered neighbours
        (default: ``%(distance_threshold)s``).

    Returns
    -------
    `dict`
        Return value from `draw()` with the following additional key::
        - ``df_intersections``: Indexed by `id, vertex_i, id_neighbour,
          vertex_i_neighbour`, where `id` is the corresponding electrode
          identifier (e.g., `electrode000`), `vertex_i` is the starting vertex
          of the line segment, `id_neighbour` is the identifier of the
          neighbouring electrode, and `vertex_i` is the starting vertex of the
          line segment of the neighbouring electrode (`pandas.DataFrame`).
    ''' % {'distance_threshold': DEFAULT_DISTANCE_THRESHOLD}
    result = draw(svg_source, ax=ax)
    df_shapes = result['df_shapes']
    ax = result['axis']

    df_intersections = \
        get_all_intersections(df_shapes, distance_threshold=distance_threshold)
    df_segments = get_segments(df_shapes,
                               distance_threshold=distance_threshold)

    for idx_i, segment_i in (df_intersections.reset_index([2, 3])
                             .join(df_segments, lsuffix='_neighbour')
                             .iterrows()):
        p = segment_i[['x_mid', 'y_mid']].values
        r = segment_i[['x_normal', 'y_normal']].values
        ax.plot(*zip(p, p + r), color='white')

    result['df_intersections'] = df_intersections

    return result
