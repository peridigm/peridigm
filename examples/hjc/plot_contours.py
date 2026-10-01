#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Render HJC impact fields on a central section or the target's upper surface.

Requires NumPy, netCDF4 and Matplotlib. Input is serial/merged Exodus output.
Cell values are not smoothed. Mesh-node displacements are volume averages of
adjacent material-point displacements, for visualization only. No points are
deleted based on damage. The disk retains the y>=0 half; the sphere is whole.
With --top-view, show upper-surface target damage and omit the sphere.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from netCDF4 import Dataset, chartostring
import numpy as np


def section_faces(points, cells, cut=True):
    """Exterior and section faces, with interpolation weights for cut edges."""
    if cells.shape[1] == 8:
        local = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4],
                 [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]]
    elif cells.shape[1] == 4:
        local = [[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]]
    else:
        raise ValueError('expected the example HEX8/TETRA mesh')
    edges = sorted({tuple(sorted((face[i-1], face[i])))
                    for face in local for i in range(len(face))})
    all_faces = cells[:, local].reshape(-1, len(local[0]))
    _, first, count = np.unique(np.sort(all_faces, axis=1), axis=0,
                                return_index=True, return_counts=True)
    if not cut:
        indices = first[count == 1]
        exterior = all_faces[indices]
        return exterior, exterior, np.zeros((*exterior.shape, 1)), indices//len(local)
    faces, owners = [], []

    def crossing(a, b):
        fraction = -points[a, 1]/(points[b, 1]-points[a, 1])
        return (a, b, fraction)

    for index in first[count == 1]:
        face, clipped = all_faces[index], []
        for a, b in zip(np.roll(face, 1), face):
            inside_a, inside_b = points[a, 1] >= 0, points[b, 1] >= 0
            if inside_a != inside_b:
                clipped.append(crossing(a, b))
            if inside_b:
                clipped.append((b, b, 0.))
        if len(clipped) >= 3:
            faces.append(clipped)
            owners.append(index//len(local))
    for owner, cell in enumerate(cells):
        if not points[cell, 1].min() < 0 < points[cell, 1].max():
            continue
        cut = [(a, a, 0.) for a in cell if points[a, 1] == 0]
        for a, b in edges:
            a, b = cell[a], cell[b]
            if points[a, 1]*points[b, 1] < 0:
                cut.append(crossing(a, b))
        xyz = np.array([(1-f)*points[a]+f*points[b] for a, b, f in cut])
        center = xyz.mean(axis=0)
        order = np.argsort(np.arctan2(xyz[:, 2]-center[2], xyz[:, 0]-center[0]))
        faces.append([cut[i] for i in order])
        owners.append(owner)
    # Pad polygons with their final vertex for a compact vectorized mapping.
    width = max(map(len, faces))
    mapping = np.array([face+[face[-1]]*(width-len(face)) for face in faces])
    return mapping[:, :, 0].astype(int), mapping[:, :, 1].astype(int), \
        mapping[:, :, 2, None], np.array(owners)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('file', type=Path)
    parser.add_argument('--mesh', type=Path,
                        default=Path(__file__).resolve().parent.parent/'disk_impact/disk_impact.g')
    parser.add_argument('--times', nargs='+', type=float, default=[10., 25., 50.],
                        help='saved times in microseconds')
    parser.add_argument('--overview', action='store_true', help='show the entire half disk')
    parser.add_argument('--top-view', action='store_true',
                        help='target upper-surface damage, with a central close-up')
    parser.add_argument('--output', type=Path, default=Path('impact_contours.png'))
    args = parser.parse_args()
    with Dataset(args.mesh) as mesh, Dataset(args.file) as data:
        points = np.column_stack([np.asarray(mesh['coord'+axis][:]) for axis in 'xyz'])*1.e3
        element_names = list(chartostring(data['name_elem_var'][:]))
        node_names = list(chartostring(data['name_nod_var'][:]))
        times = np.asarray(data['time_whole'][:])*1.e6
        frames = [int(np.argmin(abs(times-time))) for time in args.times]
        if any(abs(times[frame]-time) > 1.e-5 for frame, time in zip(frames, args.times)):
            parser.error('--times must select saved frames')

        def element(name, block, frame):
            index = element_names.index(name)+1
            return np.asarray(data[f'vals_elem_var{index}eb{block}'][frame])

        blocks, offset = [], 0
        for block in ((1,) if args.top_view else (1, 2)):
            cells = np.asarray(mesh[f'connect{block}'][:])-1
            particles = np.asarray(data[f'connect{block}'][:]).ravel()-1
            # Element_Id survives MPI partitioning; never assume output order.
            ids = element('Element_Id', block, 0).astype(int)
            order = np.argsort(ids)
            if not np.array_equal(ids[order], np.arange(offset+1, offset+len(cells)+1)):
                raise ValueError('output Element_Id does not match the example mesh')
            offset += len(cells)
            volume = element('Volume', block, 0)[order]
            if not np.isfinite(volume).all() or np.any(volume <= 0):
                raise ValueError('invalid material-point volumes')
            weight = np.repeat(volume, cells.shape[1])
            nodal_weight = np.bincount(cells.ravel(), weights=weight, minlength=len(points))
            faces = section_faces(points, cells, cut=block == 1 and not args.top_view)
            if args.top_view:
                top = np.all(np.isclose(points[faces[0], 2], points[cells, 2].max()), axis=1)
                if not np.any(top):
                    raise ValueError('no planar upper surface in the target mesh')
                faces = tuple(values[top] for values in faces)
                print(f'{len(faces[3])} target upper-surface cells selected')
            blocks.append((block, cells, particles[order], order, weight, nodal_weight, faces))

        # Orthographic projection; all lengths and displacements use the same scale.
        azimuth, elevation = np.deg2rad([-70., 12.])
        right = np.array([-np.sin(azimuth), np.cos(azimuth), 0.])
        eye = np.array([np.cos(azimuth)*np.cos(elevation),
                        np.sin(azimuth)*np.cos(elevation), np.sin(elevation)])
        projection = np.stack((right, np.cross(eye, right), eye), axis=1)
        if args.top_view:
            projection = np.eye(3)
        # Match the white background and blue/gray/red fields in doc/slide-image-*.jpg.
        cmap = plt.get_cmap('coolwarm')
        fields = [('HJC_Damage', 1., Normalize(0., 1.), 'Damage', [0., .5, 1.]),
                  ('Von_Mises_Stress', 1.e-6, Normalize(0., 600.),
                   'Equivalent stress / MPa', [0., 300., 600.])]
        if args.top_view:
            fields = [fields[0], fields[0]]
        plt.rcParams.update({'font.size': 11, 'text.color': 'black',
                             'axes.labelcolor': 'black', 'font.family': 'DejaVu Sans'})
        figure, axes = plt.subplots(2, len(frames),
                                   figsize=(4.2*len(frames), 8. if args.top_view else 5.5),
                                   squeeze=False)
        figure.subplots_adjust(left=.018, right=.91, bottom=.065, top=.90,
                               wspace=.04, hspace=.15)

        for column, frame in enumerate(frames):
            displacement = np.column_stack([
                np.asarray(data[f'vals_nod_var{node_names.index("Displacement"+axis)+1}'][frame])
                for axis in 'XYZ'])*1.e3
            if not np.isfinite(displacement).all():
                raise ValueError('nonfinite displacement')
            geometry = []
            for block, cells, particles, order, weight, nodal_weight, faces in blocks:
                numerator = np.column_stack([
                    np.bincount(cells.ravel(), minlength=len(points), weights=weight*
                                np.repeat(displacement[particles, axis], cells.shape[1]))
                    for axis in range(3)])
                moved = points+numerator/np.maximum(nodal_weight[:, None], 1.e-30)
                a, b, fraction, owners = faces
                vertices = (1-fraction)*moved[a]+fraction*moved[b]
                geometry.append((vertices, owners, order, block))
            for row, (name, scale, norm, label, ticks) in enumerate(fields):
                polygons, colors = [], []
                for vertices, owners, order, block in geometry:
                    polygons.extend(vertices@projection)
                    if block == 1:
                        values = element(name, block, frame)[order]*scale
                        if not np.isfinite(values).all() or values.min() < norm.vmin-1.e-9 \
                                or values.max() > norm.vmax+1.e-9:
                            raise ValueError(f'{name} lies outside the fixed color scale')
                        colors.extend(cmap(norm(values[owners])))
                    else:
                        normal = np.cross(vertices[:, 1]-vertices[:, 0],
                                          vertices[:, 2]-vertices[:, 0])
                        length = np.linalg.norm(normal, axis=1)
                        brightness = .65+.28*abs(normal@np.array([-.4, -.5, .768]))/np.maximum(length, 1.e-30)
                        # Use the distinct gold sphere from the upstream impact gallery.
                        colors.extend(np.column_stack((brightness, brightness*.79,
                                                       brightness*.14, np.ones(len(vertices)))))
                depth = np.array([polygon[:, 2].mean() for polygon in polygons])
                order = np.argsort(depth)
                axis = axes[row, column]
                axis.add_collection(PolyCollection([polygons[i][:, :2] for i in order],
                                    facecolors=np.asarray(colors)[order], edgecolors='none',
                                    antialiased=False))
                limits = ((-42., 42.), (-18., 31.)) if args.overview \
                    else ((-18., 18.), (-7., 14.))
                if args.top_view:
                    extent = 40. if row == 0 else 15.
                    limits = ((-extent, extent), (-extent, extent))
                    if row == 0:
                        axis.add_patch(Rectangle((-15., -15.), 30., 30., fill=False,
                                                 edgecolor='#C9CED6', linewidth=.6))
                axis.set(xlim=limits[0], ylim=limits[1], aspect='equal')
                axis.set_axis_off()
                if row == 0:
                    axis.set_title(f'{times[frame]:g} μs', fontsize=12, pad=10)
            print(f'{times[frame]:g} us: fields rendered')
        for row, (_, _, norm, label, ticks) in enumerate(fields):
            position = axes[row, -1].get_position()
            label = ['Target damage', 'Impact region'][row] if args.top_view else label
            figure.text(.026, position.y1+.018, label, size=11)
            if args.top_view and row == 1:
                continue
            bar_position = (.935, .22, .012, .52) if args.top_view \
                else (.935, position.y0+.04, .012, position.height-.07)
            bar_axis = figure.add_axes(bar_position)
            bar = figure.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=bar_axis,
                                 ticks=ticks)
            bar.outline.set_visible(False)
            bar.ax.tick_params(size=0, labelsize=9, pad=6)
        bars = [(1, -14., -17., 5.), (0, -37., -37., 10.)] if args.top_view \
            else [(1, -38., -12., 5.)] if args.overview else [(1, -16., -5.5, 5.)]
        for row, xbar, ybar, length in bars:
            axes[row, 0].plot([xbar, xbar+length], [ybar, ybar], color='black',
                              lw=1.5, clip_on=False)
            axes[row, 0].text(xbar+length/2., ybar-1.5 if args.top_view or args.overview
                              else ybar-1.1, f'{length:g} mm', ha='center', va='top', size=9)
        figure.savefig(args.output, dpi=240, facecolor='white',
                       metadata={'Software': 'HJC impact visualization'})
        plt.close(figure)


if __name__ == '__main__':
    main()
