# coding: utf-8
from __future__ import absolute_import, print_function
import itertools as it
import logging

from click import echo, style as st_
import click
import networkx as nx

from ..load import load

logger = logging.getLogger(__name__)


def info(chip_file, output):
    '''Display information about digital microfluidics chip design.'''
    chip_info = load(chip_file)

    channels_used = set(c for e in chip_info['electrodes']
                        if 'channels' in e for c in e['channels'])

    summary = (('File', click.format_filename(chip_file)),
               ('SHA256', chip_info['__metadata__']['sha256']),
               ('Pixels per inch', chip_info['__metadata__']['ppi']),
               ('Number of electrodes', len(chip_info['electrodes'])),
               ('Number of channels used', len(channels_used)))

    label_format = '%%-%ds' % (max(map(len, zip(*summary)[0])) + 2)

    for label, value in summary:
        echo(st_(label_format % (label + ':'), fg='magenta'), file=output,
             nl=False)
        echo(st_(str(value), fg='white'), file=output, nl=True)

    g = nx.Graph([(c['source']['id'], c['target']['id'])
                  for c in chip_info['connections']])
    g.add_nodes_from(e['id'] for e in chip_info['electrodes'])
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
