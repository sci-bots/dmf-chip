import logging

import networkx as nx
import numpy as np
import pandas as pd
import six

from .load import draw
from .lib.tsp_local.base import TSP
from .lib.tsp_local.kopt import KOpt


logger = logging.getLogger(__name__)


__all__ = ['nearest_neighbor_tsp', 'compute_tour', 'draw_tour']


def nearest_neighbor_tsp(shortest_paths, starting_point=0):
    """
    Nearest neighbor TSP algorithm

    Adapted from (MIT License, Copyright (c) 2018, BraveDistribution):

    https://github.com/BraveDistribution/pytsp/blob/52144ef6283196cba01a2036987e9d7b982c6f99/pytsp/nearest_neighbor_tsp.py

    Args:
        shortest_paths: 2d numpy array
        starting_point: index of starting node

    Returns:
        tour approximated by nearest neighbor tsp algorithm
    Examples:
        >>> import numpy as np
        >>> shortest_paths = np.array([[  0, 300, 250, 190, 230],
        >>>                            [300,   0, 230, 330, 150],
        >>>                            [250, 230,   0, 240, 120],
        >>>                            [190, 330, 240,   0, 220],
        >>>                            [230, 150, 120, 220,   0]])
        >>> nearest_neighbor_tsp(shortest_paths)
    """
    number_of_nodes = len(shortest_paths)
    unvisited_nodes = list(range(number_of_nodes))
    unvisited_nodes.remove(starting_point)
    visited_nodes = [starting_point]

    while number_of_nodes > len(visited_nodes):
        neighbor_distances = pd.Series(shortest_paths[visited_nodes[-1]])
        neighbor_distances = neighbor_distances[(neighbor_distances > 0) &
                                                (neighbor_distances.index
                                                 .isin(set(unvisited_nodes)))]
        next_node = neighbor_distances.idxmin()
        visited_nodes.append(next_node)
        unvisited_nodes.remove(next_node)
    return visited_nodes


def compute_tour(chip_info, start_id=None, n=5):
    '''Compute low-cost tour to visit all chip electrodes at least once.

    Parameters
    ----------
    chip_info : dict
        See `dmf_chip.load()` return type.
    start_id : str or int, optional
        Electrode id (`str`) or channel number (`int`) to start/finish tour.
    n : int, optional
        Size of random tours tournament used to select best tour.

    Returns
    -------
    list[str]
        List of electrode IDs corresponding to
    '''
    electrode_channels = {e['id']: e['channels'][0]
                          for e in chip_info['electrodes']}
    channel_electrodes = {v: k for k, v in electrode_channels.items()}

    if isinstance(start_id, six.string_types):
        # Electrode id was specified.  Convert to channel number.
        start_id = electrode_channels[start_id]

    G = nx.Graph((electrode_channels[c['source']['id']],
                  electrode_channels[c['target']['id']])
                 for c in chip_info['connections'])
    nodes = list(G.nodes)

    # Compute shortest paths between all nodes.
    shortest_paths = np.array(nx.floyd_warshall_numpy(G))

    tours = []

    for r in sorted(np.random.choice(nodes, n, replace=False)):
        TSP.edges = {}  # Global cost matrix
        TSP.ratio = 10.  # Global ratio
        TSP.routes = {}  # Global routes costs

        start_i = nodes.index(r)

        tour = nearest_neighbor_tsp(np.array(shortest_paths),
                                    starting_point=start_i)

        TSP.setEdges(shortest_paths)
        lk = KOpt(tour)
        tour, cost = lk.optimise()
        tours.append((cost, tour))
        logger.info('start: %s, cost: %s' % (r, cost))

    tour_channels = [nodes[i] for i in sorted(tours)[0][1]]
    tour_channels = np.roll(tour_channels, -tour_channels.index(start_id))
    return [channel_electrodes[c] for c in tour_channels]


def draw_tour(chip_info, tour_ids, ax=None):
    '''Draw tour through electrodes specified by tour electrode ids.

    Parameters
    ----------
    chip_info : dict
        See `dmf_chip.load()` return type.
    tour_ids : list[str]
        List of electrode ids corresponding to
    ax : matplotlib axis, optional
        Axis to draw on.

    Returns
    -------
    dict
        Dictionary with the following keys::

         - ``axis`` (matplotlib axis)
         - ``patches`` (dict) : electrode patches, indexed by electrode id.
    '''
    result = draw(chip_info, ax=ax, labels=False)
    axis = result['axis']

    # For colors, see: https://gist.github.com/cfobel/fd939073cf13a309d7a9
    dark_green = '#059748'
    dark_orange = '#df5c24'

    df_centers = pd.DataFrame((e['pole_of_accessibility']
                               for e in chip_info['electrodes']),
                              index=(e['id'] for e in chip_info['electrodes']))
    df_source = df_centers.loc[tour_ids]
    df_target = df_centers.loc[tour_ids[1:] + [tour_ids[0]]]

    # # Draw route over electrodes as sequence of arrows.
    # # See: https://stackoverflow.com/a/7543518/345236
    q = axis.quiver(df_source.x.values, df_source.y.values,
                    df_target.x.values - df_source.x.values,
                    df_target.y.values - df_source.y.values,
                    scale=1.01, units='xy', angles='xy',
                    color=dark_green, alpha=.5)
    # # Ensure route is drawn on top layer of plot.
    q.set_zorder(30)

    # # Draw tour start as circle.
    x = df_source.x
    y = df_source.y
    s = axis.scatter(x[:1], y[:1], marker='o', s=15 ** 2,
                     edgecolor=dark_orange,
                     linewidth=2, facecolor='none', label='Route start')
    s.set_zorder(20)

    # # Draw tour end as square.
    x = df_target.x
    y = df_target.y
    s = axis.scatter(x[-1:], y[-1:], marker='s', s=10 ** 2, color=dark_orange,
                     linewidth=2, facecolor='none', label='Route end')
    s.set_zorder(20)
    return result
