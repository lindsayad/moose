#!/usr/bin/env python3

import mms
import sympy

u = '0'
v = '0'
vel = u + '* e_i + ' + v + ' * e_k'

# Prescribed value: outlet
# Zero gradient: left right bottom
p = 'cos(y * pi / 2) * cos(x * pi)'

f_u, e_u = mms.evaluate('grad(p).dot(e_i)', u, variable='u', vel=vel, p=p, scalars=['mu', 'rho'], transformation='cylindrical', coordinate_names=('x', 'phi', 'y'))
f_v, e_v = mms.evaluate('grad(p).dot(e_k)', v, variable='v', vel=vel, p=p, scalars=['mu', 'rho'], transformation='cylindrical', coordinate_names=('x', 'phi', 'y'))
f_p, e_p = mms.evaluate('div(vel*rho)', p, variable='p', vel=vel, scalars=['rho'], transformation='cylindrical', coordinate_names=('x', 'phi', 'y'))

rho = sympy.Symbol('rho')

mms.print_hit(e_u, 'exact_u')
mms.print_hit(f_u, 'forcing_u', mu='${mu}', rho='${rho}')

mms.print_hit(e_v, 'exact_v')
mms.print_hit(f_v, 'forcing_v', mu='${mu}', rho='${rho}')

mms.print_hit(e_p, 'exact_p')
mms.print_hit(f_p, 'forcing_p', rho='${rho}')
