# coding: utf-8
from __future__ import absolute_import, print_function
import hashlib
import itertools as it
import json

from click import echo, style as st_
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from svg_model.connections import extract_connections
import click
import json_tricks
import networkx as nx

from ..core import _resolve_source, chip_info, ELECTRODES_XPATH


@click.command()
@click.argument('chip_file')  #, help='Digital microfluidics chip file (`*.svg`).')
@click.option('--output', default='-', type=click.File('w'))
def info(chip_file, output):
    '''Display information about digital microfluidics chip design.'''
    chip_info_ = chip_info(chip_file)

    with open(chip_file, 'rb') as input_:
        hash_ = hashlib.sha256(input_.read())

    # json_tricks.dump(chip_info_, output, indent=4)

    summary = (('File', click.format_filename(chip_file)),
               ('SHA256', hash_.hexdigest()),
               ('Number of electrodes',
                str(chip_info_['electrode_channels'].index
                    .drop_duplicates().shape[0])),
               ('Number of channels',
                str(chip_info_['channel_electrodes'].index
                    .drop_duplicates().shape[0])))

    label_format = '%%-%ds' % (max(map(len, zip(*summary)[0])) + 2)

    for label, value in summary:
        echo(st_(label_format % (label + ':'), fg='magenta'), file=output,
             nl=False)
        echo(st_(value, fg='white'), file=output, nl=True)

    df_shapes = _resolve_source(chip_file, xpath=ELECTRODES_XPATH)
    electrode_polygons = {id_: Polygon(df_i[['x', 'y']].values)
                          for id_, df_i in df_shapes.groupby('id')}

    def find_shape(x, y):
        point = Point(x, y)
        for id_, polygon_i in electrode_polygons.items():
            if polygon_i.contains(point):
                return id_
        else:
            return None

    df_connections = extract_connections(chip_file, find_shape)

    g = nx.Graph(df_connections[['source', 'target']])
    g.add_nodes_from(df_shapes.id.values)
    subgraphs = list(nx.connected_component_subgraphs(g))

    def window(seq, n):
        '''
        Returns
        -------
        iter
            Sliding window (of width n) over data from the iterable::

                s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...
        '''
        it_ = iter(seq)
        result = tuple(it.islice(it_, n))
        if len(result) <= n:
            yield result
        for elem in it_:
            result = result[1:] + (elem,)
            yield result

    if len(subgraphs) > 1:
        # No connection exists between at least two electrodes on the chip.
        echo(st_('No connection exists between the following sets of '
                 'electrodes:\n', fg='red'))
        # Display list of electrodes in each subgraph (max 4 electrodes per
        # line).
        for s in subgraphs:
            echo(st_(' - %d electrode%s:\n' % (len(s), 's' if len(s) > 1
                                               else ''), fg='blue'))
            node_list = ('\n'.join(' ' * 7 + line
                                   for line in
                                   (',\n'.join(', '.join(w)
                                               for w in window(sorted(s.nodes),
                                                               4)))
                                   .splitlines()))
            echo(st_(node_list, fg='white'))
            echo('')
    else:
        echo(st_('At least one path exists between all electrodes.',
                 fg='blue'))


if __name__ == '__main__':
    info()
