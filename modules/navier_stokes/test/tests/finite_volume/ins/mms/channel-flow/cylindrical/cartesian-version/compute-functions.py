#!/usr/bin/env python3

import mms
import sympy

u = 'sin(pi*x)*cos(pi*y)'
v = 'cos(1.3*x)*cos(pi*y)'

vel = u + '* e_i + ' + v + ' * e_j'

p = 'sin(1.5*x)*cos(1.6*y)'

f_u, e_u = mms.evaluate('div(vel*rho*u) - div(mu * grad(u)) + grad(p).dot(e_i)', u, variable='u', vel=vel, p=p, scalars=['mu', 'rho'])
f_v, e_v = mms.evaluate('div(vel*rho*v) - div(mu * grad(v)) + grad(p).dot(e_j)', v, variable='v', vel=vel, p=p, scalars=['mu', 'rho'])
f_p, e_p = mms.evaluate('div(vel*rho)', p, variable='p', vel=vel, scalars=['rho'])

#
# Boundary fluxes
#

# u
advective_flux_u_left,_ = mms.evaluate('(vel*rho*u).dot(-e_i)', u, variable='u', vel=vel, scalars=['mu', 'rho'])
advective_flux_u_right,_ = mms.evaluate('(vel*rho*u).dot(e_i)', u, variable='u', vel=vel, scalars=['mu', 'rho'])
advective_flux_u_top,_ = mms.evaluate('(vel*rho*u).dot(e_j)', u, variable='u', vel=vel, scalars=['mu', 'rho'])
advective_flux_u_bottom,_ = mms.evaluate('(vel*rho*u).dot(-e_j)', u, variable='u', vel=vel, scalars=['mu', 'rho'])
diffusive_flux_u_left,_ = mms.evaluate('(- mu * grad(u)).dot(-e_i)', u, variable='u', vel=vel, scalars=['mu', 'rho'])
diffusive_flux_u_right,_ = mms.evaluate('(- mu * grad(u)).dot(e_i)', u, variable='u', vel=vel, scalars=['mu', 'rho'])
diffusive_flux_u_top,_ = mms.evaluate('(- mu * grad(u)).dot(e_j)', u, variable='u', vel=vel, scalars=['mu', 'rho'])
diffusive_flux_u_bottom,_ = mms.evaluate('(- mu * grad(u)).dot(-e_j)', u, variable='u', vel=vel, scalars=['mu', 'rho'])

# v
advective_flux_v_left,_ = mms.evaluate('(vel*rho*v).dot(-e_i)', v, variable='v', vel=vel, scalars=['mu', 'rho'])
advective_flux_v_right,_ = mms.evaluate('(vel*rho*v).dot(e_i)', v, variable='v', vel=vel, scalars=['mu', 'rho'])
advective_flux_v_top,_ = mms.evaluate('(vel*rho*v).dot(e_j)', v, variable='v', vel=vel, scalars=['mu', 'rho'])
advective_flux_v_bottom,_ = mms.evaluate('(vel*rho*v).dot(-e_j)', v, variable='v', vel=vel, scalars=['mu', 'rho'])
diffusive_flux_v_left,_ = mms.evaluate('(- mu * grad(v)).dot(-e_i)', v, variable='v', vel=vel, scalars=['mu', 'rho'])
diffusive_flux_v_right,_ = mms.evaluate('(- mu * grad(v)).dot(e_i)', v, variable='v', vel=vel, scalars=['mu', 'rho'])
diffusive_flux_v_top,_ = mms.evaluate('(- mu * grad(v)).dot(e_j)', v, variable='v', vel=vel, scalars=['mu', 'rho'])
diffusive_flux_v_bottom,_ = mms.evaluate('(- mu * grad(v)).dot(-e_j)', v, variable='v', vel=vel, scalars=['mu', 'rho'])

# p
flux_p_left,_ = mms.evaluate('(vel*rho).dot(-e_i)', p, variable='p', vel=vel, scalars=['rho'])
flux_p_right,_ = mms.evaluate('(vel*rho).dot(e_i)', p, variable='p', vel=vel, scalars=['rho'])
flux_p_top,_ = mms.evaluate('(vel*rho).dot(e_j)', p, variable='p', vel=vel, scalars=['rho'])
flux_p_bottom,_ = mms.evaluate('(vel*rho).dot(-e_j)', p, variable='p', vel=vel, scalars=['rho'])

mms.print_hit(e_u, 'exact_u')
mms.print_hit(f_u, 'forcing_u', mu='${mu}', rho='${rho}')
mms.print_hit(advective_flux_u_left, 'advective_flux_u_left', mu='${mu}', rho='${rho}')
mms.print_hit(advective_flux_u_right, 'advective_flux_u_right', mu='${mu}', rho='${rho}')
mms.print_hit(advective_flux_u_top, 'advective_flux_u_top', mu='${mu}', rho='${rho}')
mms.print_hit(advective_flux_u_bottom, 'advective_flux_u_bottom', mu='${mu}', rho='${rho}')
mms.print_hit(diffusive_flux_u_left, 'diffusive_flux_u_left', mu='${mu}', rho='${rho}')
mms.print_hit(diffusive_flux_u_right, 'diffusive_flux_u_right', mu='${mu}', rho='${rho}')
mms.print_hit(diffusive_flux_u_top, 'diffusive_flux_u_top', mu='${mu}', rho='${rho}')
mms.print_hit(diffusive_flux_u_bottom, 'diffusive_flux_u_bottom', mu='${mu}', rho='${rho}')

mms.print_hit(e_v, 'exact_v')
mms.print_hit(f_v, 'forcing_v', mu='${mu}', rho='${rho}')
mms.print_hit(advective_flux_v_left, 'advective_flux_v_left', mu='${mu}', rho='${rho}')
mms.print_hit(advective_flux_v_right, 'advective_flux_v_right', mu='${mu}', rho='${rho}')
mms.print_hit(advective_flux_v_top, 'advective_flux_v_top', mu='${mu}', rho='${rho}')
mms.print_hit(advective_flux_v_bottom, 'advective_flux_v_bottom', mu='${mu}', rho='${rho}')
mms.print_hit(diffusive_flux_v_left, 'diffusive_flux_v_left', mu='${mu}', rho='${rho}')
mms.print_hit(diffusive_flux_v_right, 'diffusive_flux_v_right', mu='${mu}', rho='${rho}')
mms.print_hit(diffusive_flux_v_top, 'diffusive_flux_v_top', mu='${mu}', rho='${rho}')
mms.print_hit(diffusive_flux_v_bottom, 'diffusive_flux_v_bottom', mu='${mu}', rho='${rho}')

mms.print_hit(e_p, 'exact_p')
mms.print_hit(f_p, 'forcing_p', rho='${rho}')
mms.print_hit(flux_p_left, 'flux_p_left', mu='${mu}', rho='${rho}')
mms.print_hit(flux_p_right, 'flux_p_right', mu='${mu}', rho='${rho}')
mms.print_hit(flux_p_top, 'flux_p_top', mu='${mu}', rho='${rho}')
mms.print_hit(flux_p_bottom, 'flux_p_bottom', mu='${mu}', rho='${rho}')
