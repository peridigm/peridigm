#!/usr/bin/env python3
"""Generate a 3-D HJC compression/shear patch test (SI units, Python stdlib)."""
import argparse
from pathlib import Path
import xml.etree.ElementTree as ET

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--cells',type=int,default=8)
parser.add_argument('--dt',type=float,default=1.e-7)
parser.add_argument('--horizon-ratio',type=float,default=3.01)
parser.add_argument('--hourglass',type=float,default=.05)
parser.add_argument('--output',type=Path,default=Path('run'))
args=parser.parse_args()
if args.cells<3 or args.dt<=0 or args.horizon_ratio<1.01 or args.hourglass<0:
    parser.error('require cells>=3, dt>0, horizon-ratio>=1.01 and hourglass>=0')
args.output.mkdir(parents=True,exist_ok=False)
dx=.02/args.cells
mesh=[]
for k in range(args.cells):
    for j in range(args.cells):
        for i in range(args.cells):
            mesh.append(f'{(i+.5)*dx:.16g} {(j+.5)*dx:.16g} {(k+.5)*dx:.16g} 1 {dx**3:.16g}\n')
(args.output/'mesh.txt').write_text(''.join(mesh),encoding='ascii')
(args.output/'all_nodes.txt').write_text('\n'.join(str(i+1) for i in range(len(mesh)))+'\n',encoding='ascii')
root=ET.Element('ParameterList')
def group(parent,name):return ET.SubElement(parent,'ParameterList',name=name)
def value(parent,name,item):
    kind='bool' if isinstance(item,bool) else 'int' if isinstance(item,int) else 'double' if isinstance(item,float) else 'string'
    ET.SubElement(parent,'Parameter',name=name,type=kind,value=str(item).lower() if kind=='bool' else str(item))
disc=group(root,'Discretization')
value(disc,'Type','Text File');value(disc,'Input Mesh File','mesh.txt')
mat=group(group(root,'Materials'),'Concrete')
value(mat,'Material Model','HJC Correspondence')
for key,item in {
 'Density':2700.,'Shear Modulus':2.417e10,'Hourglass Coefficient':args.hourglass,
 'A':.29,'B':2.06,'C':.0013,'N':.866,'Compressive Strength':1.19e8,
 'Tensile Strength':8.2e6,'Reference Strain Rate':1.,'EFMIN':.01,'SFMAX':5.,
 'Crushing Pressure':4.e7,'Crushing Volumetric Strain':.00124,
 'Locking Pressure':1.2e9,'Locking Plastic Volumetric Strain':.011,
 'D1':.04,'D2':1.,'K1':1.287e10,'K2':1.631e10,'K3':6.495e10}.items():value(mat,key,item)
block=group(group(root,'Blocks'),'Cube')
value(block,'Block Names','block_1');value(block,'Material','Concrete')
value(block,'Horizon',args.horizon_ratio*dx)
bc=group(root,'Boundary Conditions');value(bc,'Node Set All','all_nodes.txt')
for axis,expression in [('x','(-0.015*x+0.04*y)'),('y','(-0.015*y+0.02*z)'),('z','(-0.015*z)')]:
    one=group(bc,'Affine '+axis)
    value(one,'Type','Prescribed Displacement');value(one,'Node Set','Node Set All')
    value(one,'Coordinate',axis);value(one,'Value',expression+'*sin(31415.92653589793*t)^2')
solver=group(root,'Solver');value(solver,'Initial Time',0.);value(solver,'Final Time',1.e-4)
value(group(solver,'Verlet'),'Fixed dt',args.dt)
out=group(root,'Output');value(out,'Output File Type','ExodusII')
value(out,'Output Filename','hjc');value(out,'Output Frequency',max(1,round(2.e-6/args.dt)))
fields=group(out,'Output Variables')
for name in ['Displacement','Cauchy_Stress','HJC_Pressure','Von_Mises_Stress',
             'HJC_Damage','HJC_Plastic_Volume','Equivalent_Plastic_Strain',
             'HJC_Density_Measure','HJC_Acoustic_Modulus','Element_Id']:
    value(fields,name,True)
ET.indent(root)
ET.ElementTree(root).write(args.output/'hjc.xml',encoding='utf-8',xml_declaration=True)
print(args.output/'hjc.xml')
