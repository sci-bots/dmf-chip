# coding: utf-8
from __future__ import absolute_import, print_function
import functools as ft
import hashlib
import io
import logging
import sys

from click import style as st_
import click
import click_log
import lxml.etree
import numpy as np
import pandas as pd
import six

from .bin.info import info as info_
from .core import ureg
from .edit import (DEFAULT_DISTANCE_THRESHOLD, write_connections_layer,
                   write_test_route)
from .load import load
from .tour import compute_tour

logger = logging.getLogger(__name__)


def handle_stdin(ctx, param, value):
    if value == '-':
        value = io.BytesIO(sys.stdin.read())
    return value


def mutex(*keys):
    def _wrapped(ctx, param, value):
        if value is not None and any(ctx.params.get(k) is not None
                                    for k in keys):
            raise click.BadParameter('at most one of the following args may be'
                                     ' specified: %s' %
                                     ', '.join('`%s`' % k for k in keys))
        return value
    return _wrapped


@click.group()
@click_log.simple_verbosity_option()
def dmf_chip():
    pass


# @click.group(invoke_without_command=True, help='Display information about '
             # 'digital microfluidics chip design.')
@dmf_chip.command(help='Display information about digital microfluidics chip '
                  'design.')
# @click.pass_context
@click.argument('chip_file', callback=handle_stdin)
@click.option('--output', default='-', type=click.File('w'))
# def info(ctx, chip_file, output):
def info(chip_file, output):
    info_(chip_file, output)


@click.group(help='Neighbour actions')
def neighbour():
    pass


@neighbour.command(help='Detect adjacent electrodes and add "Connections" '
                   'layer  to output SVG file containing lines; each '
                   'corresponding to a cached connection between two '
                   'electrodes.')
@click.argument('chip_file', callback=handle_stdin)
@click.argument('output', default='-', type=click.File('w'))
@click.option('--distance-threshold', help='Default: %s' %
              DEFAULT_DISTANCE_THRESHOLD)
def detect(chip_file, output, distance_threshold):
    if distance_threshold is None:
        distance_threshold = DEFAULT_DISTANCE_THRESHOLD
    else:
        distance_threshold = ureg.parse_expression(distance_threshold)
    doc = write_connections_layer(chip_file, distance_threshold)
    doc.write(output)
    logger.info(st_('Wrote connections layer to: ', fg='magenta') +
                st_('`%s`' % output.name, fg='white'))


@click.group(help='Route actions')
def route():
    pass


@route.command(help='Find low-cost test tour that passes through each '
               'electrode **at least** once and add test route to '
               '`<dmf:`ChipDesign><dmf:TestRoutes>`.')
@click.argument('chip_file', callback=handle_stdin)
@click.argument('output', default='-', type=click.File('w'))
@click.option('--start-id', help='Start electrode id',
              callback=mutex('start_id', 'start_channel'))
@click.option('--start-channel', help='Start channel number', type=int,
              callback=mutex('start_id', 'start_channel'))
@click.option('--seed', help='Random seed', type=int)
@click.option('--route-id', help='Test route id', default='default')
def compute(chip_file, output, start_id, start_channel, seed, route_id):
    doc = lxml.etree.parse(chip_file)
    root = doc.getroot()
    NSMAP = {k: v for k, v in root.nsmap.items() if k}
    _xpath = ft.partial(root.xpath, namespaces=NSMAP)
    routes = _xpath('//dmf:ChipDesign/dmf:TestRoutes/dmf:TestRoute[@id!=""]')
    if route_id in (r.attrib['id'] for r in routes if 'id' in r.attrib):
        raise click.BadParameter('a test route already exists with the id '
                                 '`%s`' % route_id, param_hint='route_id')

    if start_id is None and start_channel is None:
        raise click.BadParameter('one of the following args must be '
                                 'specified: %s' %
                                 ', '.join('`%s`' % arg
                                           for arg in ('start-id',
                                                       'start-channel')))
    elif start_channel is not None:
        start_id = start_channel

    chip_info = load(chip_file)

    if isinstance(start_channel, six.string_types):
        seed = np.fromstring(hashlib.sha256(start_channel).digest(),
                            dtype='uint32')[0]
    else:
        seed = start_channel

    with pd.util.testing.RNGContext(seed=seed):
        tour_ids = list(compute_tour(chip_info, start_id))

    electrode_channels = {e['id']: e['channels'][0]
                          for e in chip_info['electrodes']}
    tour_channels = [electrode_channels[i] for i in tour_ids]
    doc = write_test_route(chip_file, tour_channels, id_=route_id)
    doc.write(output)
    logger.info(st_('Wrote test route to: ', fg='magenta') +
                st_('`%s`' % output.name, fg='white'))


@route.command(help='Display routes', name='list')
@click.argument('chip_file', callback=handle_stdin)
def route__list(chip_file):
    chip_info = load(chip_file)

    for route_i in chip_info['__metadata__']['test-routes']:
        print(route_i)


@route.command(help='Add a test route', name='add')
@click.argument('chip_file', callback=handle_stdin)
@click.argument('output', default='-', type=click.File('w'))
@click.option('--id', required=True, help='Route id')
@click.option('-w', '--waypoint', multiple=True, help='Channel number. '
              'Multiple channel numbers may be specified, e.g., `-w 10 -w 23 '
              '...`', type=click.INT)
def route__add(chip_file, output, id, waypoint):
    doc = lxml.etree.parse(chip_file)
    root = doc.getroot()
    NSMAP = {k: v for k, v in root.nsmap.items() if k}
    _xpath = ft.partial(root.xpath, namespaces=NSMAP)
    routes = _xpath('//dmf:ChipDesign/dmf:TestRoutes/dmf:TestRoute[@id!=""]')
    if id in (r.attrib['id'] for r in routes if 'id' in r.attrib):
        raise click.BadParameter('a test route already exists with the id '
                                 '`%s`' % id, param_hint='id')
    doc = write_test_route(chip_file, list(waypoint), id_=id)
    doc.write(output)
    logger.info(st_('Wrote test route to: ', fg='magenta') +
                st_('`%s`' % output.name, fg='white'))


@route.command(help='Remove routes', name='remove')
@click.argument('chip_file', callback=handle_stdin)
@click.argument('output', default='-', type=click.File('w'))
@click.option('--id', required=True, multiple=True, help='id of route to '
              'remove. Multiple ids may be specified, e.g.,   `--id foo --id '
              'bar ...`')
def route__remove(chip_file, output, id):
    doc = lxml.etree.parse(chip_file)
    root = doc.getroot()
    NSMAP = {k: v for k, v in root.nsmap.items() if k}
    _xpath = ft.partial(root.xpath, namespaces=NSMAP)
    routes = _xpath('//dmf:ChipDesign/dmf:TestRoutes/dmf:TestRoute[@id!=""]')

    routes_by_id = {r.attrib['id']: r for r in routes if 'id' in r.attrib}

    for id_i in id:
        if id_i not in routes_by_id:
            raise click.BadParameter('no test route exists with the id `%s`' %
                                     id_i, param_hint='id')

    for id_i in id:
        route_i = routes_by_id[id_i]
        route_i.getparent().remove(route_i)
        logger.info(st_('Removed test route: ', fg='magenta') +
                    st_('`%s`' % id_i, fg='white'))
    doc.write(output)
    logger.info(st_('Wrote modified chip file to: ', fg='magenta') +
                st_('`%s`' % output.name, fg='white'))


dmf_chip.add_command(neighbour)
dmf_chip.add_command(route)


def main():
    for logger_name in (None, 'dmf_chip', 'TSP'):
        logger_ = logging.getLogger(logger_name)
        click_log.basic_config(logger_)
    click_log.basic_config(logger)
    dmf_chip()


if __name__ == '__main__':
    main()
