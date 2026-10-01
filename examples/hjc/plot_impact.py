#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Plot the HJC disk-impact example and check its saved particle histories.

Requires numpy, netCDF4 and matplotlib. Pass one or more serial/epu-merged
Exodus files. Densities are those of disk_impact.yaml (2700 and 7700 kg/m^3).
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset, chartostring
import numpy as np


def read_history(path):
    with Dataset(path) as data:
        element_names = list(chartostring(data['name_elem_var'][:]))
        node_names = list(chartostring(data['name_nod_var'][:]))
        times = np.asarray(data['time_whole'][:])

        def element(name, block=1):
            return np.asarray(data[f'vals_elem_var{element_names.index(name)+1}eb{block}'][:])

        def node(name):
            return np.asarray(data[f'vals_nod_var{node_names.index(name)+1}'][:])

        disk = np.asarray(data['connect1'][:]).ravel()-1
        ball = np.asarray(data['connect2'][:]).ravel()-1
        disk_volume, ball_volume = element('Volume')[0], element('Volume',2)[0]
        velocity = np.stack([node('Velocity'+axis) for axis in 'XYZ'],axis=-1)
        contact = np.stack([node('Contact_Force'+axis) for axis in 'XYZ'],axis=-1)
        damage = element('HJC_Damage')
        gradient = np.stack([element('Deformation_Gradient'+i+j)
                             for i in 'XYZ' for j in 'XYZ'],axis=-1)
        jacobian = np.linalg.det(gradient.reshape(len(times),-1,3,3))
        pressure = element('HJC_Pressure')
        epsp = element('Equivalent_Plastic_Strain')
        plastic_volume = element('HJC_Plastic_Volume')
        for values in (velocity,contact,damage,jacobian,pressure,epsp,plastic_volume):
            if not np.isfinite(values).all():
                raise ValueError(f'{path}: nonfinite solver output')
        if jacobian.min() <= 0 or damage.min() < 0 or damage.max() > 1:
            raise ValueError(f'{path}: invalid deformation or damage bounds')
        if np.diff(damage,axis=0).min() < -1.e-12:
            raise ValueError(f'{path}: damage decreased')
        for values in (epsp,plastic_volume):
            if values.min() < -1.e-12 or np.diff(values,axis=0).min() < -1.e-12:
                raise ValueError(f'{path}: irreversible plastic history decreased')

        momentum = (velocity[:,disk]*disk_volume[None,:,None]*2700.).sum(axis=1)
        momentum += (velocity[:,ball]*ball_volume[None,:,None]*7700.).sum(axis=1)
        scale = np.linalg.norm(momentum[0])
        if scale <= 0:
            raise ValueError(f'{path}: expected nonzero initial sphere momentum')
        drift = np.linalg.norm(momentum-momentum[0],axis=1).max()/scale
        if drift > 1.e-6:
            raise ValueError(f'{path}: free-impact momentum balance failed ({drift:g})')
        ball_speed = (velocity[:,ball,2]*ball_volume).sum(axis=1)/ball_volume.sum()
        force = contact[:,ball,2].sum(axis=1)  # Contact_Force already includes volume.
        mean_damage = (damage*disk_volume).sum(axis=1)/disk_volume.sum()
        curves = np.column_stack((times,force,ball_speed,mean_damage))
        summary = {
            'file': str(path), 'last_time_s': float(times[-1]),
            'disk_points': len(disk), 'ball_points': len(ball),
            'minimum_disk_J': float(jacobian.min()),
            'maximum_disk_damage': float(damage.max()),
            'final_mean_disk_damage': float(mean_damage[-1]),
            'maximum_disk_pressure_Pa': float(pressure.max()),
            'maximum_disk_plastic_strain': float(epsp.max()),
            'maximum_disk_plastic_volume': float(plastic_volume.max()),
            'peak_contact_force_N': float(force.max()),
            'final_ball_velocity_m_s': float(ball_speed[-1]),
            'relative_momentum_drift': float(drift),
            'contact_force_balance_N': float(np.linalg.norm(contact.sum(axis=1),axis=1).max()),
        }
    return curves, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('files',type=Path,nargs='+')
    parser.add_argument('--labels',nargs='+')
    parser.add_argument('--output',type=Path,default=Path('impact_response'))
    args = parser.parse_args()
    if args.labels and len(args.labels) != len(args.files):
        parser.error('provide one label per input file')
    labels = args.labels or [path.stem for path in args.files]
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,
        'axes.spines.right':False,'axes.edgecolor':'#A6ADB5',
        'axes.labelcolor':'#34404D','xtick.color':'#59636F','ytick.color':'#59636F',
        'grid.color':'#E4E9ED','grid.linewidth':.6})
    figure,axes = plt.subplots(1,3,figsize=(11.5,3.4),layout='constrained')
    colors = ['#3B4CC0','#B40426','#D2A329','#333333']
    summaries = []
    for index,(path,label) in enumerate(zip(args.files,labels)):
        curves,summary = read_history(path)
        summaries.append(summary)
        np.savetxt(args.output.with_name(args.output.name+f'_{index+1}.csv'),curves,
                   delimiter=',',header='time_s,contact_force_z_N,ball_velocity_z_m_s,mean_disk_damage',comments='')
        for column,axis in enumerate(axes,1):
            scale = .001 if column == 1 else 1.
            axis.plot(curves[:,0]*1.e6,curves[:,column]*scale,
                      color=colors[index%len(colors)],lw=1.3,label=label)
    for axis,label in zip(axes,['Contact force z / kN','Sphere velocity z / m/s','Mean disk damage']):
        axis.set(xlabel='Time / microseconds',ylabel=label)
        axis.grid(axis='y')
        axis.set_axisbelow(True)
    axes[0].legend(frameon=False,fontsize=9)
    figure.savefig(args.output.with_suffix('.png'),dpi=220,bbox_inches='tight')
    plt.close(figure)
    args.output.with_suffix('.json').write_text(json.dumps(summaries,indent=2)+'\n',encoding='utf-8')


if __name__ == '__main__':
    main()
