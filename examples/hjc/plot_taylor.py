#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Check and plot serial/merged Exodus output from generate_taylor.py.

Requires NumPy, netCDF4 and Matplotlib. Histories use all cylinder points;
contours show the reference y>=0 half at actual displaced coordinates.
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset, chartostring
import numpy as np


def read(path):
    case = json.loads((path.parent/'case.json').read_text())
    with Dataset(path) as data:
        enames = list(chartostring(data['name_elem_var'][:]))
        nnames = list(chartostring(data['name_nod_var'][:]))
        time = np.asarray(data['time_whole'][:])

        def element(name, block=1):
            return np.asarray(data[f'vals_elem_var{enames.index(name)+1}eb{block}'][:])

        def node(name):
            return np.asarray(data[f'vals_nod_var{nnames.index(name)+1}'][:])

        cylinder = np.asarray(data['connect1'][:]).ravel()-1
        wall = np.asarray(data['connect2'][:]).ravel()-1
        reference = np.column_stack([np.asarray(data['coord'+axis][:]) for axis in 'xyz'])
        displacement = np.stack([node('Displacement'+axis) for axis in 'XYZ'], axis=-1)
        velocity = np.stack([node('Velocity'+axis) for axis in 'XYZ'], axis=-1)
        contact = np.stack([node('Contact_Force'+axis) for axis in 'XYZ'], axis=-1)
        volume = element('Volume')[0]
        damage = element('HJC_Damage')
        q = element('Von_Mises_Stress')
        pressure = element('HJC_Pressure')
        ep = element('Equivalent_Plastic_Strain')
        vp = element('HJC_Plastic_Volume')
        gradient = np.stack([element('Deformation_Gradient'+i+j)
                             for i in 'XYZ' for j in 'XYZ'], axis=-1)
        jacobian = np.linalg.det(gradient.reshape(len(time), -1, 3, 3))
        acoustic = element('HJC_Acoustic_Modulus')
        for values in (time, reference, displacement, velocity, contact, volume,
                       damage, q, pressure, ep, vp, jacobian, acoustic):
            if not np.isfinite(values).all():
                raise ValueError(f'{path}: nonfinite output')
        if len(time) < 2 or abs(time[-1]-case['end_time']) > case['dt']:
            raise ValueError(f'{path}: incomplete history')
        if np.any(volume <= 0) or jacobian.min() <= 0:
            raise ValueError(f'{path}: invalid volume or deformation determinant')
        if damage.min() < 0 or damage.max() > 1:
            raise ValueError(f'{path}: damage outside [0,1]')
        for values in (damage, ep, vp):
            if values.min() < -1.e-12 or np.diff(values, axis=0).min() < -1.e-12:
                raise ValueError(f'{path}: irreversible history decreased')
        # Verlet output includes the final half-kick velocity on constrained
        # points. Their prescribed displacement, not that diagnostic velocity,
        # defines the fixed wall; only cylinder velocities enter the histories.
        if np.abs(displacement[:, wall]).max() > 1.e-14:
            raise ValueError(f'{path}: wall moved')
        positions = reference[cylinder][None, :, :] + displacement[:, cylinder]
        surface_gap = positions[:, :, 2].min(axis=1)-case['spacing']/2
        if surface_gap.min() < -.1*case['spacing']:
            raise ValueError(f'{path}: penetration exceeds 10% of spacing')
        mass = volume*case['density']
        momentum = (velocity[:, cylinder]*mass[None, :, None]).sum(axis=1)
        force = contact[:, cylinder].sum(axis=1)  # Includes point volume.
        impulse = np.vstack([np.zeros(3), np.cumsum(
            .5*(force[1:]+force[:-1])*np.diff(time)[:, None], axis=0)])
        residual = momentum-momentum[0]-impulse
        initial = np.linalg.norm(momentum[0])
        curves = np.column_stack((time, force[:, 2], momentum[:, 2]/mass.sum(),
                                  damage@volume/volume.sum(), surface_gap))
        summary = dict(
            case=case, minimum_J=float(jacobian.min()), maximum_damage=float(damage.max()),
            final_mean_damage=float(curves[-1, 3]), final_mean_axial_velocity_m_s=float(curves[-1, 2]),
            maximum_equivalent_stress_Pa=float(q.max()), peak_contact_force_N=float(force[:, 2].max()),
            maximum_wall_penetration_m=float(max(0., -surface_gap.min())),
            contact_force_balance_N=float(np.linalg.norm(contact.sum(axis=1), axis=1).max()),
            saved_frame_impulse_relative_residual=float(np.linalg.norm(residual, axis=1).max()/initial),
            acoustic_step_ratio=float(case['dt']*np.sqrt(acoustic.max()/case['density'])/case['spacing']),
            final_mean_plastic_strain=float(ep[-1]@volume/volume.sum()),
            final_mean_plastic_volume=float(vp[-1]@volume/volume.sum()))
    return curves, summary, (time, reference[cylinder], positions, damage, q)


def read_curves(path, saved, summary):
    """The small history file resolves contact transients at every timestep."""
    history = path.parent/'taylor_history.h'
    if not history.exists():
        return saved  # Allows inspection of field-only pilot runs.
    with Dataset(history) as data:
        names = list(chartostring(data['name_glo_var'][:]))
        values = np.asarray(data['vals_glo_var'][:])
        time = np.asarray(data['time_whole'][:])
        case = summary['case']
        count, dx = case['cylinder_points'], case['spacing']
        contact = np.stack([values[:, names.index('Cylinder_Contact_Density'+axis)]
                            for axis in 'XYZ'], axis=-1)*dx**3
        velocity = np.stack([values[:, names.index('Cylinder_Velocity_Sum'+axis)]
                             for axis in 'XYZ'], axis=-1)/count
        damage = values[:, names.index('Cylinder_Damage_Sum')]/count
        if damage.min() < 0 or damage.max() > 1 or np.diff(damage).min() < -1.e-12:
            raise ValueError(f'{history}: invalid damage history')
        impulse = np.vstack([np.zeros(3), np.cumsum(
            .5*(contact[1:]+contact[:-1])*np.diff(time)[:, None], axis=0)])
        mass = count*dx**3*case['density']
        residual = mass*(velocity-velocity[0])-impulse
        summary['step_impulse_relative_residual'] = float(
            np.linalg.norm(residual, axis=1).max()/(mass*case['speed']))
        summary['peak_contact_force_N'] = float(contact[:, 2].max())
        if not np.isfinite(values).all() or summary['step_impulse_relative_residual'] > 1.e-7:
            raise ValueError(f'{history}: nonfinite history or impulse imbalance')
        if abs(time[-1]-case['end_time']) > case['dt']:
            raise ValueError(f'{history}: incomplete history')
        return np.column_stack((time, contact[:, 2], velocity[:, 2], damage,
                                np.interp(time, saved[:, 0], saved[:, 4])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('files', type=Path, nargs='+')
    parser.add_argument('--labels', nargs='+')
    parser.add_argument('--output', type=Path, default=Path('taylor'))
    args = parser.parse_args()
    if args.labels and len(args.labels) != len(args.files):
        parser.error('one label is required per file')
    labels = args.labels or [path.parent.name for path in args.files]
    histories = []
    for path in args.files:
        curves, summary, frames = read(path)
        histories.append((read_curves(path, curves, summary), summary, frames))
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.edgecolor': '#A6ADB5',
                         'grid.color': '#E4E9ED', 'grid.linewidth': .6})
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.4), layout='constrained')
    colors = ['#3B4CC0', '#B40426', '#D2A329']
    summaries = []
    for index, ((curves, summary, _), label) in enumerate(zip(histories, labels)):
        summaries.append(summary)
        np.savetxt(args.output.with_name(args.output.name+f'_{index+1}.csv'), curves,
                   delimiter=',', header='time_s,contact_force_z_N,mean_velocity_z_m_s,mean_damage,surface_gap_m',
                   comments='')
        for column, ax in enumerate(axes, 1):
            ax.plot(curves[:, 0]*1.e6, curves[:, column]*(.001 if column == 1 else 1.),
                    lw=1.3, color=colors[index % len(colors)], label=label)
    for ax, label in zip(axes, ['Contact force / kN', 'Axial velocity / m/s', 'Mean damage']):
        ax.set(xlabel=r'Time / $\mu s$', ylabel=label)
        ax.grid(axis='y')
        ax.set_axisbelow(True)
    if len(histories) > 1:
        axes[0].legend(frameon=False)
    fig.savefig(args.output.with_name(args.output.name+'_response.png'), dpi=220)
    plt.close(fig)
    for index in range(1, len(histories)):
        base, other = histories[0][0], histories[index][0]
        # Compare on the coarser grid; do not invent finely sampled force
        # transients by interpolating a sparse history.
        sample = base if len(base) <= len(other) else other
        start, end = max(base[0, 0], other[0, 0]), min(base[-1, 0], other[-1, 0])
        sample = sample[(sample[:, 0] >= start) & (sample[:, 0] <= end)]
        error = {}
        for column, name in [(1, 'contact_force'), (2, 'velocity'), (3, 'mean_damage')]:
            first = np.interp(sample[:, 0], base[:, 0], base[:, column])
            second = np.interp(sample[:, 0], other[:, 0], other[:, column])
            scale = np.max(np.abs(second))
            error[name+'_peak_normalized_rmse'] = float(
                np.sqrt(np.mean((first-second)**2))/max(scale, 1.e-30))
        error['samples'] = len(sample)
        summaries[index]['comparison_to_first'] = error
    args.output.with_suffix('.json').write_text(json.dumps(summaries, indent=2)+'\n')
    time, reference, positions, damage, q = histories[0][2]
    targets = np.array([20., 40., 60.])*1.e-6
    frames = [int(np.argmin(np.abs(time-t))) for t in targets]
    if any(abs(time[f]-t) > 1.e-12 for f, t in zip(frames, targets)):
        parser.error('contours require saved frames at 20, 40 and 60 microseconds')
    upper = max(50., np.ceil(max(q[f].max() for f in frames)/5.e7)*50.)
    keep = reference[:, 1] >= 0
    fig = plt.figure(figsize=(10, 7), layout='constrained')
    axes = np.array([[fig.add_subplot(2, 3, 3*i+j+1, projection='3d')
                      for j in range(3)] for i in range(2)])
    dx = summaries[0]['case']['spacing']
    for j, frame in enumerate(frames):
        xyz = positions[frame, keep]*1.e3
        for i, (values, limit) in enumerate(((damage[frame, keep], 1.), (q[frame, keep]/1.e6, upper))):
            ax = axes[i, j]
            artist = ax.scatter(*xyz.T, c=values, s=5*(dx/.0005)**2,
                                cmap='coolwarm', vmin=0., vmax=limit,
                                linewidths=0, depthshade=False, rasterized=True)
            wx, wy = np.meshgrid([-7., 7.], [0., 6.])
            ax.plot_surface(wx, wy, np.zeros_like(wx), color='.7', alpha=.6, shade=False)
            ax.set(xlim=(-7, 7), ylim=(0, 7), zlim=(0, 22))
            ax.set_box_aspect((14, 7, 22))
            ax.set_proj_type('ortho')
            ax.view_init(elev=15, azim=-65)
            ax.set_axis_off()
            if i == 0:
                ax.set_title(r'$%g\,\mu s$' % (time[frame]*1.e6))
            if j == 2:
                fig.colorbar(artist, ax=axes[i, :], shrink=.65, ticks=[0., limit/2, limit],
                             label='Damage' if i == 0 else 'q / MPa')
    fig.savefig(args.output.with_name(args.output.name+'_contours.png'), dpi=220)
    plt.close(fig)


if __name__ == '__main__':
    main()
