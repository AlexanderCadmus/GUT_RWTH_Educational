### v2 save has sidebar


# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 10:58:38 2022

@author: Oxos
"""
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
import gempy as gp
import matplotlib.image as mpimg
import matplotlib as mpl
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
st.set_page_config(layout='wide')
###Introduction###
st.title('Das Kastner Modell')
st.markdown(r'''
Originalautor: Alexander Cadmus

https://github.com/AlexanderCadmus
''')
with st.sidebar:
    st.header('''Tunnel Eigenschaften''')
    d = st.slider('Tunnel Durchmesser (d) [m]', min_value=0.5, max_value=50.0, value=10.0, step=0.5, format=None, key=1, help=None, on_change=None, args=None, kwargs=None, disabled=False)
    rd = np.arange(d,d*5,1)
    h = st.slider('Tunnel Tiefe (h) [m]', min_value=d, max_value=500.0, value=60.0, step=1.0, format=None, key=2, help=None, on_change=None, args=None, kwargs=None, disabled=False)
    st.header('''__________________________________________________''')
    st.header('''Fels Eigenschaften''')
    v=st.slider('Querkontraktionszahl/Poissonzahl (v) [-]',value=0.2,min_value=0.0, max_value=0.5,step=0.01)
    K_null= v/(1-v)
    st.markdown('''Mit den Querkontraktionszahl kann die Erdruhedruckbeiwert oder koeffizient horizonal spannung gerechnet.''')
    st.latex(r'''K_0 = \left(\frac{v}{1-v}\right)''')
    st.markdown(np.round(K_null,3))
    gamma=st.slider('Effektiv Wichte Fels (γ) [kN/m3] ',value=28.0,min_value=0.0, max_value=50.0,step=0.5)
    st.header('''__________________________________________________''')
    st.header('''Block Eigenschaften''')
    theta=st.slider('Winkel von Vertical (θ) [°]',value=0,min_value=0, max_value=360,step=1)
    rm=st.slider('Entfernung von Tunnel Mitte [m]', value=rd.min()+10, min_value=rd.min(), max_value=rd.max(), step=1.0)
    rm_segment_temp = np.where(rd == rm)
    rm_segment = int(rm_segment_temp[0])
        
###########################################################################
###################-----Part 1-----########################################
###########################################################################
###########################################################################
###################-----Define Sampling Procedure-----##################### 
###########################################################################

# st.markdown('''___________________________________________''')
# st.header('''Part 1''')
# st.markdown('''___________________________________________''')
# st.markdown('''Notes''')


# with st.expander(r'''See Code: XXX'''):
#     metroplis_code='''Pv = gamma*h # Vertical stress
# Ph = K_null*Pv # Horizontal stress'''
#     st.code(metroplis_code, language='python')

Pv = gamma*h # Vertical stress
Ph = K_null*Pv # Horizontal stress
rd = np.arange(d,d*5,1) # distance (radius) form tunnel in [m]

# st.markdown('''___________________________________________''')
# st.header('''Part 2''')
# st.markdown('''___________________________________________''')
# st.markdown('''Notes''')

# with st.expander(r'''See Code: XXX'''):
#     metroplis_code='''def Kastner_rad(ro, r, vo, Pv, theta):
#         alpha = ro/r
#         Knull = vo/(1-vo)
#         # Radial Pressure
#         Pr = (Pv/2)*((1+Knull)*(1-alpha**2)+(1-Knull)*(1+3*alpha**4-4*alpha**2)*np.cos(2*np.radians(theta)))
#         return Pr'''
#     st.code(metroplis_code, language='python')
    

def Kastner_rad(ro, r, vo, Pv, theta):
    alpha = ro/r
    Knull = vo/(1-vo)
    # Radial Pressure
    Pr = (Pv/2)*((1+Knull)*(1-alpha**2)+(1-Knull)*(1+3*alpha**4-4*alpha**2)*np.cos(2*np.radians(theta)))
    return Pr

# st.markdown('''___________________________________________''')
# st.header('''Part 3''')
# st.markdown('''___________________________________________''')
# st.markdown('''Notes''')

# with st.expander(r'''See Code: XXX'''):
#     metroplis_code='''def Kastner_tang(ro, r, vo, Pv, theta):
#         alpha = ro/r
#         Knull = vo/(1-vo)
#         # Tangential Pressure
#         Pt = (Pv/2)*((1+Knull)*(1+alpha**2)-(1-Knull)*(1+3*alpha**4)*np.cos(2*np.radians(theta)))
#         return Pt'''
#     st.code(metroplis_code, language='python')


def Kastner_tang(ro, r, vo, Pv, theta):
    alpha = ro/r
    Knull = vo/(1-vo)
    # Tangential Pressure
    Pt = (Pv/2)*((1+Knull)*(1+alpha**2)-(1-Knull)*(1+3*alpha**4)*np.cos(2*np.radians(theta)))
    return Pt

# st.markdown('''___________________________________________''')
# st.header('''Part 4''')
# st.markdown('''___________________________________________''')
# st.markdown('''Notes''')

# with st.expander(r'''See Code: XXX'''):
#     metroplis_code='''def Kastner_schub(ro, r, vo, Pv, theta):
#         alpha = ro/r
#         Knull = vo/(1-vo)
#         # Schubspannungen
#         Trt = (-Pv/2)*(1-Knull)*(1-(3*alpha**4)+(2*alpha**2))*np.sin(2*np.radians(theta))
#         return Trt'''
#     st.code(metroplis_code, language='python')


def Kastner_schub(ro, r, vo, Pv, theta):
    alpha = ro/r
    Knull = vo/(1-vo)
    # Schubspannungen
    Trt = (-Pv/2)*(1-Knull)*(1-(3*alpha**4)+(2*alpha**2))*np.sin(2*np.radians(theta))
    return Trt

def setup_axes(fig, rect, rotation, axisScale, axisLimits, doShift):
    tr_rot = Affine2D().scale(axisScale[0], axisScale[1]).rotate_deg(rotation)

    # This seems to do nothing
    if doShift:
        tr_trn = Affine2D().translate(-90,-5)
    else:
        tr_trn = Affine2D().translate(0,0)

    tr = tr_rot + tr_trn

    grid_helper = floating_axes.GridHelperCurveLinear(tr, extremes=axisLimits)

    ax = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
    fig.add_subplot(ax)
    aux_ax = ax.get_aux_axes(tr)

    return ax, aux_ax

## This is a stable version
ro = d
vo = v
Krad=[]
Ktang=[]
Kschub=[]
for j in range(len(rd)):
    Krad.append(Kastner_rad(ro, rd[j], vo, Pv, theta))
    Ktang.append(Kastner_tang(ro, rd[j], vo, Pv, theta))
    Kschub.append(Kastner_schub(ro, rd[j], vo, Pv, theta))

#    print('Krad[-1] = ', Krad[-1])
#    print('Ktang[-1] = ', Ktang[-1])
#    print('Kschub[-1] = ', Kschub[-1])
fig, ax = plt.subplots(1, 3, figsize=(16, 5))

# Plot the tunnel face with investigation direction
point=(0,0)
dist = np.round(rd[rm_segment],0)
x, y = point
endy = (y + rd[rm_segment]*np.sin(np.radians(90-theta)))
endx = (x + rd[rm_segment]*np.cos(np.radians(90-theta)))
tunnel = plt.Circle((0,0),d, fill=True, edgecolor='k', facecolor='white')
rockmass = plt.Rectangle(((0-max(rd))*1.05,(0-max(rd))*1.05), (max(rd)*3+d), 
                         (max(rd)*3+d), fill=True, hatch='...',facecolor='grey')
floating_block_display_dims_factor = 1
floating_block_display = plt.Rectangle(((rd[rm_segment]*np.cos(np.radians(90-theta))),
                                        (rd[rm_segment]*np.sin(np.radians(90-theta)))), (floating_block_display_dims_factor), 
                                       (floating_block_display_dims_factor), fill=False, edgecolor='black',
                                       rotation_point='center', linewidth=10, angle=-theta)
ax[0].set_aspect(1)
ax[0].add_artist(rockmass)
ax[0].add_artist(tunnel)
ax[0].add_artist(floating_block_display)
ax[0].scatter(0,0,c='k')
ax[0].plot([x, endx], [y, endy],linestyle='dashed', c='black')
ax[0].set_xlim((-max(rd)*1.05), (max(rd))*1.05)
ax[0].set_ylim((-max(rd)*1.05), (max(rd))*1.05)
ax[0].set_xlabel('Breite [m]')
ax[0].set_ylabel('Hohe [m]')
ax[0].set_title('Tunnel Querschnitt')
fig.tight_layout()
ax[1].axhline(Pv, label=('Pv = '+str(Pv)), c='r', linestyle='--')
ax[1].axhline(Ph, label=('Ph = '+str(round(Ph,1))), c='b', linestyle='--')
ax[1].axvline(rd[rm_segment], label=('Segment Entfernung = '+str(dist)), c='y', linestyle='--')
ax[1].plot(rd, Krad, label='Radialspannung', c='g')
ax[1].plot(rd, Ktang, label='Tangentialspannung', c='m')
ax[1].plot(rd, Kschub, label='Schubspannung', c='k')
#ax[1].set_yticks(np.round(np.linspace(round(min(Krad+Ktang+Kschub),-1), round(max(Krad+Ktang+Kschub),-2), 10),-2))
ax[1].set_yticks(np.arange(-4000,20000,2000))
ax[1].set_ylim(-4000,20000)
ax[1].set_xlim(d, max(rd))
ax[1].set_xticks(np.arange(d, round(max(rd),-1), round(((max(rd)-d)/10),0)))
ax[1].set_xlabel('Entfenung von Tunnel-mitte [m]')
ax[1].set_ylabel('Spanung [kN/m$^2$]')
ax[1].set_title('Kastner Spanungen')
ax[1].grid()
ax[1].legend()
fig.tight_layout()
floating_block_dims_factor = 1
floating_block = plt.Rectangle((0-floating_block_dims_factor/2,0-floating_block_dims_factor/2), 
                               (floating_block_dims_factor), (floating_block_dims_factor), fill=True, hatch='.',
                               facecolor='grey',rotation_point='center', angle=(90-theta))
floating_block_frame = plt.Rectangle((0-floating_block_dims_factor/2,0-floating_block_dims_factor/2), 
                               (floating_block_dims_factor), (floating_block_dims_factor), fill=False,
                               edgecolor='black',rotation_point='center', angle=(90-theta), linewidth=2)
point_origin=(0,0)
def TwoD_euclidean_dist(startx,starty,endx,endy):
    euc_dist = (((endx-startx)**2)+(endy-starty)**(2))**(0.5)
    return euc_dist
# Radialspangung Pfile
x_arrow_r, y_arrow_r = point_origin
starty_arrow_r = y_arrow_r + floating_block_dims_factor/2 * np.sin(np.radians(90-theta))
startx_arrow_r = x_arrow_r + floating_block_dims_factor/2 * np.cos(np.radians(90-theta))
endy_arrow_r = starty_arrow_r + Krad[rm_segment]/4000 * np.sin(np.radians(90-theta))
endx_arrow_r = startx_arrow_r + Krad[rm_segment]/4000 * np.cos(np.radians(90-theta))
arrow_length_r = TwoD_euclidean_dist(startx_arrow_r,starty_arrow_r,endx_arrow_r,endy_arrow_r)
# Tangentialspangung Pfile
x_arrow_t, y_arrow_t = point_origin
starty_arrow_t = y_arrow_t + floating_block_dims_factor/2 * np.sin(np.radians(0-theta))
startx_arrow_t = x_arrow_t + floating_block_dims_factor/2 * np.cos(np.radians(0-theta))
endy_arrow_t = starty_arrow_t + Ktang[rm_segment]/4000 * np.sin(np.radians(0-theta))
endx_arrow_t = startx_arrow_t + Ktang[rm_segment]/4000 * np.cos(np.radians(0-theta))
arrow_length_t = TwoD_euclidean_dist(startx_arrow_t,starty_arrow_t,endx_arrow_t,endy_arrow_t)
# Schubspangung Pfile
x_arrow_s, y_arrow_s = point_origin
starty_arrow_s = y_arrow_s + (floating_block_dims_factor*np.sin(np.radians(45))) * np.sin(np.radians(135-theta))
startx_arrow_s = x_arrow_s + (floating_block_dims_factor*np.sin(np.radians(45))) * np.cos(np.radians(135-theta))
endy_arrow_s = starty_arrow_s + Kschub[rm_segment]/4000 * np.sin(np.radians(0-theta))
endx_arrow_s = startx_arrow_s + Kschub[rm_segment]/4000 * np.cos(np.radians(0-theta))
arrow_length_s = TwoD_euclidean_dist(startx_arrow_s,starty_arrow_s,endx_arrow_s,endy_arrow_s)
# Plotting rockmass segment
ax[2].set_aspect(1)
ax[2].add_artist(floating_block_frame)
ax[2].add_artist(floating_block)
ax[2].arrow(endx_arrow_r, endy_arrow_r, startx_arrow_r-endx_arrow_r, starty_arrow_r-endy_arrow_r,head_width=arrow_length_r/3 , head_length=arrow_length_r/2,
            length_includes_head=True, color='g', width=0.1, label='Radialspannung')
ax[2].arrow(endx_arrow_t, endy_arrow_t, startx_arrow_t-endx_arrow_t, starty_arrow_t-endy_arrow_t,head_width=arrow_length_t/3 , head_length=arrow_length_t/2,
            length_includes_head=True, color='m', width=0.1, label='Tangentialspannung')
ax[2].arrow(endx_arrow_s, endy_arrow_s, startx_arrow_s-endx_arrow_s, starty_arrow_s-endy_arrow_s,head_width=arrow_length_s/3 , head_length=arrow_length_s/2,
            length_includes_head=True, color='k', width=0.1, label='Schubspannung')
ax[2].set_xlim((-3), (3))
ax[2].set_ylim((-3), (3))
ax[2].set_xlabel('Breite [m]')
ax[2].set_ylabel('Hohe [m]')
ax[2].set_title('Tunnelwand-Segment '+str(dist)+' [m] entfernd von Tunnel-mitte')
ax[2].legend()
st.pyplot(fig)