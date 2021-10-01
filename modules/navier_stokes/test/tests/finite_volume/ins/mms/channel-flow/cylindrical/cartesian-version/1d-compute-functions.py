#!/usr/bin/env python3

import mms
import sympy

u = 'cos(x)'
p = 'sin(x)'
vel = u + '* e_i'

f_u, e_u = mms.evaluate('div(vel*rho*u) - div(mu * grad(u)) + grad(p).dot(e_i)', u, variable='u', vel=vel, p=p, scalars=['mu', 'rho'])
f_p, e_p = mms.evaluate('div(vel*rho)', p, variable='p', vel=vel, scalars=['rho'])
flux_ux,_ = mms.evaluate('(vel*rho*u - mu * grad(u)).dot(e_i)', u, variable='u', vel=vel, scalars=['mu', 'rho'])
flux_px,_ = mms.evaluate('(vel*rho).dot(e_i)', p, variable='p', vel=vel, scalars=['rho'])

mms.print_hit(e_u, 'exact_u')
mms.print_hit(f_u, 'forcing_u', mu='${mu}', rho='${rho}')
mms.print_hit(flux_ux, 'flux_ux', mu='${mu}', rho='${rho}')

mms.print_hit(e_p, 'exact_p')
mms.print_hit(f_p, 'forcing_p', rho='${rho}')
mms.print_hit(flux_px, 'flux_px', rho='${rho}')
