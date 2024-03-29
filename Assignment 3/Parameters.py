'''
    Parameter list and init functions for the Animator
    For Assignment 3, TFY4345
    Created by: Haakon Svane
    Date: 13. November, 2019
'''

import matplotlib
matplotlib.use('TKAgg')                 # Default backend (support for both windows and MacOSX) to ensure optimal speed.
from matplotlib.pyplot import rcParams

FRAME_RATE = 60                         # Frames per second
T_MAX = 90                              # Total time of calculation

# Time-interval of each plot-segment (length of time-dependent plots)

FIGURE_RATIO = (18, 10)                 # Format: ([WIDTH], [HEIGHT])
FIGURE_SIZE_FACTOR = 6                  # Any number between 1-10 to set the scaling of the plot

# Parameters for the matplotlib plots to secure consistency on multiple devices
rcParams['figure.dpi'] = 100            # Sets the DPI of the figure. Should be left unchanged for best performance
rcParams["figure.figsize"] = [i/10*FIGURE_SIZE_FACTOR for i in FIGURE_RATIO]
#rcParams["toolbar"] = "None"
rcParams['axes.titlesize'] = "medium"
rcParams['axes.formatter.limits'] = -2, 2
rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False
rcParams['axes.spines.left'] = False
rcParams['axes.spines.bottom'] = False


g_acc = 9.80665                         # (m/s^2)
init_theta = 0.2                        # (rad)
init_theta_vel = 0                      # (rad/s)


# Declaration of initial values and physical constants
init_rod_length = g_acc                 # (m)
init_bob_mass = 1                       # (kg)
init_dt = 0.030                         # (s)
pivot_coords = (0, 0)                   # ([x]m, [y]m)

init_Om_d = 2/3
init_F_d = 0

init_Dtheta = 0.001


# Relative positioning of subplots and sliders
p_top = 0.94
p_bot = 0.08
p_left = 0.09
p_right = 0.82
h_space = 0.4
w_space = 0.25

slider_pad = 0.02
slider_w = 0.019
button_w = 3*slider_w+2*slider_pad+0.04

# Some default values related to the menu system
init_speed_index = 4                    # Default value for the animation speed setting. Index 4 refers to 'No Animation'
init_dampening_index = 3                # Default value for the dampening factor





