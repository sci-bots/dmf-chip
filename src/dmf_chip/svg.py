'''
Utility functions, etc. for working with SVG documents.
'''
import logging
import re

import numpy as np
import svg_model as sm

__all__ = ['shape_points', 'get_CTM', 'get_current_transformation_matrix']

CRE_TRANSFORM = re.compile(r'(?P<operation>(skew[XY]|scale|translate|rotate))'
                           r'\((?P<args>[^\)]+)\)')


def get_transform(op):
    '''
    See SVG `transform`_ documentation for matrix definitions.

    _transform:: https://www.w3.org/TR/SVG11/coords.html#TransformMatrixDefined
    '''
    # Start with identity matrix (no transform).
    T = np.eye(3)
    args = map(float, re.split(r'\s*[,\s]\s*', op['args']))

    if op['operation'] == 'matrix':
        T[0] = args[:3]
        T[1] = args[3:]
    elif op['operation'] == 'translate':
        if len(args) == 1:
            args.append(args[0])
        T[0, 2] = args[0]
        T[1, 2] = args[1]
    elif op['operation'] == 'scale':
        if len(args) == 1:
            args.append(args[0])
        T[0, 0] = args[0]
        T[1, 1] = args[1]
    elif op['operation'] == 'rotate':
        angle = (args[0] / 180.) * np.pi
        T = np.array([[np.cos(angle), -np.sin(angle), 0],
                      [np.sin(angle), np.cos(angle), 0],
                      [0, 0, 1]])
        if len(args) == 3:
            # Rotation point was specified; `(cx, cy)`.
            rotate_angle, cx, cy = args
            # Translate, perform rotation, and translate back.
            C = np.eye(3)
            C[:2, 2] = cx, cy
            C_reverse = C.copy()
            C_reverse[:2, 2] = -cx, -cy
            T = reduce(np.matmul, (C, T, C_reverse))
    elif op['operation'] == 'skewX':
        T[0, 1] = np.tan(args[0] / 180 * np.pi)
    elif op['operation'] == 'skewY':
        T[1, 0] = np.tan(args[0] / 180 * np.pi)

    return T


def get_current_transformation_matrix(element):
    transforms = []
    parents = [element.getparent()]

    while True:
        parent = parents[-1].getparent()
        if parent is None:
            break
        parents.append(parent)

    for i, parent in enumerate(parents[::-1]):
        if 'transform' in parent.attrib:
            transforms_i = [get_transform(match.groupdict())
                            for match in CRE_TRANSFORM
                            .finditer(parent.attrib['transform'])]
            transforms.extend(transforms_i)
    if transforms:
        return reduce(np.matmul, transforms)
    else:
        return np.eye(3, dtype=float)


def get_CTM(*args, **kwargs):
    '''
    Alias for :func:`get_current_transformation_matrix()`.
    '''
    return get_current_transformation_matrix(*args, **kwargs)


def shape_points(svg_element):
    '''
    Parameters
    ----------
    svg_element : lxml.etree.Element
        Either a ``<svg:path>`` or ``<svg:polygon>`` element.

    Returns
    -------
    list
        List of coordinates of points found in SVG shape element.

        Each point is represented by a dictionary with keys ``x`` and ``y``.
    '''
    if svg_element.tag.endswith('/svg}path'):
        # Decode `svg:path` vertices from [`"d"`][1] attribute.
        #
        # [1]: https://developer.mozilla.org/en-US/docs/Web/SVG/Attribute/d
        points = sm.shape_path_points(svg_element.attrib['d'])
        # Convert dictionary points to lists.
        points = [[p['x'], p['y']] for p in points]
    elif svg_element.tag.endswith('/svg}polygon'):
        # Decode `svg:polygon` vertices from [`"points"`][2] attribute.
        #
        # [2]: https://developer.mozilla.org/en-US/docs/Web/SVG/Attribute/points
        points = [map(float, v.split(','))
                  for i, v in enumerate(svg_element.attrib['points'].strip()
                                        .split(' '))]
    elif svg_element.tag.endswith('/svg}line'):
        # Decode `svg:polygon` vertices from [`"points"`][2] attribute.
        #
        # [2]: https://developer.mozilla.org/en-US/docs/Web/SVG/Attribute/points
        points = [map(float, p)
                  for p in [[svg_element.attrib['x1'],
                             svg_element.attrib['y1']],
                            [svg_element.attrib['x2'],
                             svg_element.attrib['y2']]]]
    else:
        raise NotImplementedError('Unsupported SVG tag: `%s`', svg_element.tag)
    points = np.asarray(points)

    T = get_CTM(svg_element)
    if not np.equal(T, np.eye(3)).all():
        if logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
            logging.debug('Applying transformation matrix: `%s`',
                          T[:2].ravel())
        padding = np.expand_dims(np.array([1] * len(points)), axis=1)
        points_padded = np.append(points, padding, axis=1)
        points = np.matmul(T, points_padded.T)[:2].T
    return points.tolist()