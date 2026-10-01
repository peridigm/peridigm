#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Generate a 3-D HJC Taylor cylinder and fixed contact wall (SI units).

Uses the Python standard library. Run Peridigm taylor.xml from the generated
directory. Geometry is fixed; --spacing changes the Cartesian discretization.
"""
import argparse
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spacing', type=float, default=.0005)
    parser.add_argument('--dt', type=float, default=1.e-8)
    parser.add_argument('--end-time', type=float, default=6.e-5)
    parser.add_argument('--wall-factor', type=float, default=2.)
    parser.add_argument('--output', type=Path, default=Path('taylor'))
    args = parser.parse_args()
    if not all(math.isfinite(v) and v > 0 for v in
               (args.spacing, args.dt, args.end_time, args.wall_factor)):
        parser.error('spacing, dt, end-time and wall-factor must be finite and positive')
    n = round(.01 / args.spacing)
    if n < 4:
        parser.error('at least four subdivisions across the diameter are required')
    dx, gap, speed = .01 / n, .0005, 30.
    horizon = 3.01 * dx
    bulk, shear, rho = 4.e7 / .00124, 2.417e10, 2700.
    longitudinal = bulk + 4 * shear / 3
    # The native pair force visits both endpoints. For aligned equal-volume
    # points, acceleration/overlap = 18*K*dx^3/(pi*rho*horizon^5).
    spring = args.wall_factor * math.pi * 3.01**5 * longitudinal / 18
    if args.dt * math.sqrt(max(1., args.wall_factor) * longitudinal / rho) / dx > .15:
        parser.error('dt exceeds the example acoustic/contact estimate; reduce --dt')
    output_steps = round(1.e-6 / args.dt)
    if output_steps < 1 or abs(output_steps * args.dt - 1.e-6) > 1.e-12:
        parser.error('dt must divide the 1 microsecond output interval')
    args.output.mkdir(parents=True, exist_ok=False)
    mesh = []
    xy = [(i + .5) * dx - .005 for i in range(n)]
    for x in xy:
        for y in xy:
            if x*x + y*y < .005**2:
                for k in range(2*n):
                    mesh.append((x, y, gap + (k+.5)*dx, 1, dx**3))
    count = len(mesh)
    # A fixed, single-layer wall needs no correspondence shape tensor. Half
    # spacing in x/y reduces gaps in the native point contact surface; its
    # dx-thick integration cells have volume dx^3/4 and top surface z=0.
    pad = math.ceil(.003 / dx)
    wall_xy = [(i*.5+.5)*dx-.005 for i in range(-2*pad, 2*(n+pad))]
    for x in wall_xy:
        for y in wall_xy:
            mesh.append((x, y, -.5*dx, 2, dx**3/4))
    (args.output/'mesh.txt').write_text(''.join(
        f'{x:.16g} {y:.16g} {z:.16g} {block} {volume:.16g}\n'
        for x, y, z, block, volume in mesh), encoding='ascii')
    for name, ids in [('cylinder', range(1, count+1)),
                      ('wall', range(count+1, len(mesh)+1))]:
        (args.output/f'{name}_nodes.txt').write_text(
            ''.join(f'{i}\n' for i in ids), encoding='ascii')
    root = ET.Element('ParameterList')

    def group(parent, name):
        return ET.SubElement(parent, 'ParameterList', name=name)

    def value(parent, name, item):
        kind = ('bool' if isinstance(item, bool) else 'int' if isinstance(item, int)
                else 'double' if isinstance(item, float) else 'string')
        ET.SubElement(parent, 'Parameter', name=name, type=kind,
                      value=str(item).lower() if kind == 'bool' else str(item))

    disc = group(root, 'Discretization')
    for key, item in [('Type', 'Text File'), ('Input Mesh File', 'mesh.txt'),
                      ('Omit Bonds Between Blocks', 'All')]:
        value(disc, key, item)
    materials = group(root, 'Materials')
    concrete = group(materials, 'Concrete')
    for key, item in {
        'Material Model': 'HJC Correspondence', 'Density': rho,
        'Shear Modulus': shear, 'Hourglass Coefficient': .05,
        'A': .29, 'B': 2.06, 'C': .0013, 'N': .866,
        'Compressive Strength': 1.19e8, 'Tensile Strength': 8.2e6,
        'Reference Strain Rate': 1., 'EFMIN': .01, 'SFMAX': 5.,
        'Crushing Pressure': 4.e7, 'Crushing Volumetric Strain': .00124,
        'Locking Pressure': 1.2e9, 'Locking Plastic Volumetric Strain': .011,
        'D1': .04, 'D2': 1., 'K1': 1.287e10, 'K2': 1.631e10,
        'K3': 6.495e10}.items():
        value(concrete, key, item)
    wall = group(materials, 'Wall')
    for key, item in {'Material Model': 'Elastic', 'Density': rho,
                      'Bulk Modulus': bulk, 'Shear Modulus': shear}.items():
        value(wall, key, item)
    blocks = group(root, 'Blocks')
    for name, block_id in [('Concrete', 1), ('Wall', 2)]:
        block = group(blocks, name)
        value(block, 'Block Names', f'block_{block_id}')
        value(block, 'Material', name)
        value(block, 'Horizon', horizon)
    contact = group(root, 'Contact')
    value(contact, 'Search Radius', 3*dx)
    value(contact, 'Search Frequency', max(1, round(1.e-7/args.dt)))
    model = group(group(contact, 'Models'), 'Wall Contact')
    for key, item in {'Contact Model': 'Short Range Force', 'Contact Radius': dx,
                      'Spring Constant': spring, 'Friction Coefficient': 0.}.items():
        value(model, key, item)
    pair = group(group(contact, 'Interactions'), 'Cylinder Wall')
    for key, item in [('First Block', 'block_1'), ('Second Block', 'block_2'),
                      ('Contact Model', 'Wall Contact')]:
        value(pair, key, item)
    bc = group(root, 'Boundary Conditions')
    value(bc, 'Node Set Cylinder', 'cylinder_nodes.txt')
    value(bc, 'Node Set Wall', 'wall_nodes.txt')
    initial = group(bc, 'Initial Axial Velocity')
    for key, item in [('Type', 'Initial Velocity'), ('Node Set', 'Node Set Cylinder'),
                      ('Coordinate', 'z'), ('Value', str(-speed))]:
        value(initial, key, item)
    for axis in 'xyz':
        fixed = group(bc, 'Fixed Wall '+axis)
        for key, item in [('Type', 'Prescribed Displacement'), ('Node Set', 'Node Set Wall'),
                          ('Coordinate', axis), ('Value', '0.0')]:
            value(fixed, key, item)
    # Keep every-step global histories separate from the larger field output.
    # Force density is an existing field; multiply its block sum by dx^3 when
    # plotting. This avoids a derived Contact_Force compute-order dependency.
    computes = group(root, 'Compute Class Parameters')
    for variable, label in [('Contact_Force_Density', 'Cylinder_Contact_Density'),
                            ('Velocity', 'Cylinder_Velocity_Sum'),
                            ('HJC_Damage', 'Cylinder_Damage_Sum')]:
        compute = group(computes, label)
        for key, item in [('Compute Class', 'Block_Data'), ('Calculation Type', 'Sum'),
                          ('Block', 'block_1'), ('Variable', variable), ('Output Label', label)]:
            value(compute, key, item)
    solver = group(root, 'Solver')
    value(solver, 'Initial Time', 0.)
    value(solver, 'Final Time', args.end_time)
    value(group(solver, 'Verlet'), 'Fixed dt', args.dt)
    output = group(root, 'Output')
    value(output, 'Output File Type', 'ExodusII')
    value(output, 'Output Filename', 'taylor')
    value(output, 'Output Frequency', output_steps)
    fields = group(output, 'Output Variables')
    for name in ['Displacement', 'Velocity', 'Volume', 'Element_Id', 'Contact_Force',
                 'HJC_Pressure', 'Von_Mises_Stress', 'HJC_Damage', 'HJC_Plastic_Volume',
                 'Equivalent_Plastic_Strain', 'Deformation_Gradient', 'HJC_Acoustic_Modulus']:
        value(fields, name, True)
    histories = group(root, 'Output Histories')
    value(histories, 'Output File Type', 'ExodusII')
    value(histories, 'Output Filename', 'taylor_history')
    value(histories, 'Output Frequency', 1)
    fields = group(histories, 'Output Variables')
    for label in ['Cylinder_Contact_Density', 'Cylinder_Velocity_Sum', 'Cylinder_Damage_Sum']:
        value(fields, label, True)
    ET.indent(root)
    ET.ElementTree(root).write(args.output/'taylor.xml', encoding='utf-8', xml_declaration=True)
    case = dict(spacing=dx, cylinder_points=count, wall_points=len(mesh)-count,
                radius=.005, length=.02, gap=gap, speed=speed, density=rho,
                dt=args.dt, end_time=args.end_time, horizon=horizon,
                wall_factor=args.wall_factor, spring_constant=spring)
    (args.output/'case.json').write_text(json.dumps(case, indent=2)+'\n', encoding='utf-8')
    print(args.output/'taylor.xml')


if __name__ == '__main__':
    main()
