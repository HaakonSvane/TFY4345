'''
    A set of setup-functions used by the Animator class
    For Assignment 3, TFY4345
    Created by: Haakon Svane
    Date: 13. November, 2019
'''

from matplotlib.pyplot import axes, subplots, subplot, subplots_adjust, figure
from matplotlib.widgets import Circle, Line2D, Button, RadioButtons, Slider
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import matplotlib.gridspec as gridspec

import numpy as np

from ColorSchemer import *
from Parameters import *

class AnimatedActors:
    def __init__(self, parent, plot_color="default"):
        self.rod = None
        self.circle = None
        self.line_angle = None
        self.line_phase = None
        self.line_strob = None
        self.line_vel = None
        self.line_angle_diff = None
        self.line_vel_diff = None
        self.color = 0 if plot_color == "default" else plot_color
        self.parent = parent
        self.z_offset = 1 if not self.color == 0 else 0

    def init_pendulum_artists(self):
        self.circle = Circle((0, 0), -np.exp(-self.parent.bob_mass / 10) + 2.1, animated=True,
                             color=Color.BOB if not self.color else self.color, zorder=10 - self.z_offset)
        self.rod = Line2D([pivot_coords[0], self.circle.get_center()[0]],
                          [pivot_coords[1], self.circle.get_center()[1]], lw=3,
                          animated=True, color=Color.ROD if not self.color else self.color, zorder=5 - self.z_offset)

    def init_angle_artists(self, plot_frame):
        self.line_angle, = plot_frame.plot([], [], animated=True, linestyle='None',
                                           color=Color.PLOTLINE if not self.color else self.color,
                                           markersize=1 if self.z_offset == 0 else 1,
                                           zorder=5 - self.z_offset, marker=".")

    def init_phase_artists(self, plot_frame):
        self.line_phase, = plot_frame.plot([], [], animated=True, linestyle="None",
                                           color=Color.PLOTLINE if not self.color else self.color,
                                           markersize=1 if self.z_offset == 0 else 1,
                                           zorder=5 - self.z_offset, marker=".")

    def init_strob_artists(self, plot_frame):
        self.line_strob, = plot_frame.plot([], [], animated=True, marker="o", linestyle="None",
                                           color="darkred" if not self.color else "yellow",
                                           markersize=6 if self.z_offset == 0 else 8,
                                           zorder=5 - self.z_offset)

    def init_vel_artists(self, plot_frame):
        self.line_vel, = plot_frame.plot([], [], animated=True,
                                         color=Color.PLOTLINE if not self.color else self.color,
                                         lw=2 if self.z_offset == 0 else 3,
                                         zorder=5 - self.z_offset)

    def init_angle_diff_artists(self, plot_frame):
        self.line_angle_diff, = plot_frame.plot([], [], animated=True,
                                                color=Color.PLOTLINE if not self.color else self.color,
                                                lw=1)

    def init_vel_diff_artists(self, plot_frame):
        self.line_vel_diff, = plot_frame.plot([], [], animated=True,
                                              color=Color.PLOTLINE if not self.color else self.color,
                                              lw=1)

    def attach_pendulum_artists(self, plot_frame):
        plot_frame.add_patch(self.circle)
        plot_frame.add_line(self.rod)


    def get_artists(self):
        artists =  (self.line_angle, self.line_phase,
                    self.circle,     self.rod,
                    self.line_vel,   self.line_angle_diff,
                    self.line_strob, self.line_vel_diff)

        return [i for i in artists if i]


class AxisSetup:
    def __init__(self, parent_obj):
        self.parent = parent_obj
        self.fig = None

        self.pendulum_plot = None
        self.angle_plot = None
        self.phase_plot = None
        self.vel_plot = None
        self.angle_plot_diff = None
        self.vel_plot_diff = None

        self.ax_rod_length = None
        self.ax_bob_mass = None
        self.ax_dt = None
        self.ax_dampening = None
        self.ax_second_calc = None
        self.ax_anim_params = None
        self.ax_apply_run = None
        self.ax_reset = None

        self.slider_rod_length = None
        self.slider_bob_mass = None
        self.slider_dt = None
        self.radio_anim_params = None
        self.radio_damping = None
        self.button_apply_run = None
        self.button_reset_apply = None

        # Setting the text color to the current color-scheme
        rcParams['text.color'] = Color.TEXT
        rcParams['axes.labelcolor'] = Color.TEXT
        rcParams['xtick.color'] = Color.TEXT
        rcParams['ytick.color'] = Color.TEXT

    def plot_setup(self, pre_fig=None):

        # Setting the figure and unpacking the different subplots to variables

        if not self.parent.show_second:
            self.fig, ((self.pendulum_plot, self.angle_plot),
                       (self.phase_plot,    self.vel_plot)) = subplots(2, 2,
                                                                       gridspec_kw={'width_ratios': [1, 2.5]})
        else:
            self.fig = figure()
            self.gs = gridspec.GridSpec(ncols=5*4, nrows=3, figure=self.fig)
            self.pendulum_plot = self.fig.add_subplot(self.gs[0, 0:2*4])
            self.angle_plot = self.fig.add_subplot(self.gs[0, 2*4+1:])
            self.phase_plot = self.fig.add_subplot(self.gs[1, 0:2*4])
            self.vel_plot = self.fig.add_subplot(self.gs[1, 2 * 4 + 1:])
            self.angle_plot_diff = self.fig.add_subplot(self.gs[2, 1:10])
            self.vel_plot_diff = self.fig.add_subplot(self.gs[2, 11:])

        subplots_adjust(left=p_left, right=p_right, top=p_top, bottom=p_bot, hspace=h_space, wspace=w_space)
        self.fig.set_facecolor(Color.BACKGROUND)
        self.pendulum_plot.set(aspect='equal')
        self.phase_plot.set(aspect='equal')

        # Setting axes for the sliders and buttons
        self.ax_rod_length = axes([0.815 + slider_pad, 0.73, slider_w-0.005, p_top - p_bot - (0.75 - 0.10)],
                                  facecolor=Color.SLIDER_INACTIVE)
        self.ax_bob_mass = axes([0.815 + slider_pad * 2 + slider_w-0.005, 0.73, slider_w-0.005, p_top - p_bot - (0.75 - 0.10)],
                                facecolor=Color.SLIDER_INACTIVE)
        self.ax_dt = axes([0.815 + slider_pad * 3 + 2 * (slider_w-0.005), 0.73, slider_w-0.005, p_top - p_bot - (0.75 - 0.10)],
                          facecolor=Color.SLIDER_INACTIVE)

        self.ax_Fd = axes([0.815 + slider_pad * 4 + 3 * (slider_w-0.005), 0.73, slider_w-0.005, p_top - p_bot - (0.75 - 0.10)],
                          facecolor=Color.SLIDER_INACTIVE)

        self.ax_Om_d = axes([0.815 + slider_pad * 5 + 4 * (slider_w-0.005), 0.73, slider_w-0.005, p_top - p_bot - (0.75 - 0.10)],
                          facecolor=Color.SLIDER_INACTIVE)

        self.ax_dampening = axes([0.82 + slider_pad, p_bot + 0.315 + button_w, button_w, button_w],
                                 facecolor=Color.RADIO_BACKGROUND,
                                 title="Dampening type")
        self.ax_second_calc = axes([0.82 + slider_pad, p_bot + 0.34, button_w, button_w - 0.065],
                                   facecolor=Color.RADIO_BACKGROUND,
                                   title="Second trajectory")

        self.ax_Dtheta = axes([0.835 + slider_pad, p_bot + 0.30, button_w-0.055, slider_w],
                                   facecolor=Color.SLIDER_INACTIVE)

        self.ax_anim_params = axes([0.82 + slider_pad, p_bot + 0.105, button_w, button_w],
                                   facecolor=Color.RADIO_BACKGROUND, title="Animation speed")

        self.ax_apply_run = axes([0.82 + slider_pad, p_bot + 0.05, button_w, slider_w * 2])
        self.ax_reset = axes([0.82 + slider_pad, p_bot, button_w, slider_w * 2])

        self.slider_rod_length = Slider(self.ax_rod_length, "l", 5, 15, color=Color.SLIDER_ACTIVE,
                                        valinit=self.parent.rod_length,
                                        valstep=0.1, orientation="vertical")

        self.slider_bob_mass = Slider(self.ax_bob_mass, "m", 0.1, 10.0, color=Color.SLIDER_ACTIVE,
                                      valinit=self.parent.bob_mass,
                                      valstep=0.1, orientation="vertical")

        self.slider_dt = Slider(self.ax_dt, "dt", 0.001, 0.3, color=Color.SLIDER_ACTIVE, valinit=self.parent.dt,
                                valstep=0.001, orientation="vertical", valfmt="%1.3f")

        self.slider_Fd = Slider(self.ax_Fd, "$F_D$", 0.0, 2, color=Color.SLIDER_ACTIVE, valinit=self.parent.F_d,
                                valstep=0.1, orientation="vertical", valfmt="%1.1f")

        self.slider_Om_d = Slider(self.ax_Om_d, "$\Omega_D$", 0.0, 2*np.pi, color=Color.SLIDER_ACTIVE, valinit=self.parent.Om_d,
                                valstep=0.01, orientation="vertical", valfmt="%1.2f")

        self.radio_damping = RadioButtons(self.ax_dampening, ["underdamped", "overdamped", "critical damping", "none"],
                                          active=self.parent.dampening_index,
                                          activecolor=Color.RADIO_SELECTED)
        for i in self.radio_damping.labels: i.set_fontsize("small")

        self.radio_second_calc = RadioButtons(self.ax_second_calc, ["No", "Yes"],
                                              active=self.parent.show_second,
                                              activecolor=Color.RADIO_SELECTED)
        for i in self.radio_second_calc.labels: i.set_fontsize("small")

        self.slider_Dtheta = Slider(self.ax_Dtheta, "$\Delta\\theta$", 0.0, 0.01, color=Color.SLIDER_ACTIVE, valinit=self.parent.init_Dtheta,
                                valstep=0.0001, orientation="horizontal", valfmt="%1.4f")

        self.radio_anim_params = RadioButtons(self.ax_anim_params, ["$1/2$x", "1x", "4x", "10x", "No animation"],
                                              active=self.parent.speed_index,
                                              activecolor=Color.RADIO_SELECTED)
        for i in self.radio_anim_params.labels: i.set_fontsize("small")

        self.button_apply_run = Button(self.ax_apply_run, "Apply and run..", color=Color.BUTTON_OFF_HOVER,
                                       hovercolor=Color.BUTTON_ON_HOVER)
        self.button_reset_apply = Button(self.ax_reset, "Reset and run..", color=Color.BUTTON_OFF_HOVER,
                                         hovercolor=Color.BUTTON_ON_HOVER)

    def pendulum_plot_init(self):
        self.pendulum_plot.set_facecolor(Color.PLOT_BACKGROUND)
        self.pendulum_plot.set_xlabel('$x$')
        self.pendulum_plot.set_ylabel('$y$')
        self.pendulum_plot.set_title('Pendulum display')

        fac = np.ceil(self.parent.rod_length) * (1 + 0.2)
        self.pendulum_plot.set_xlim(-fac, fac)
        self.pendulum_plot.set_ylim(-fac, fac)
        self.pendulum_plot.grid(b=True, which='both')
        self.pendulum_plot.xaxis.set_major_locator(MultipleLocator(2 * fac / 4))
        self.pendulum_plot.yaxis.set_major_locator(MultipleLocator(2 * fac / 4))
        self.pendulum_plot.xaxis.set_minor_locator(AutoMinorLocator(2))
        self.pendulum_plot.yaxis.set_minor_locator(AutoMinorLocator(2))

    def angle_plot_init(self):
        self.angle_plot.set_facecolor(Color.PLOT_BACKGROUND)
        self.angle_plot.set_xlabel('$t$')
        self.angle_plot.set_ylabel('$\Theta$ (rad)')
        self.angle_plot.set_title('Angle $\Theta$ (rad) over time $t$')

        fac = np.max(np.abs(self.parent.solutions[0][1])) * (1 + 0.1) if np.max(
            np.abs(self.parent.solutions[0][1])) < np.pi else 2 / 3 * np.pi
        self.angle_plot.set_xlim(0, T_MAX)
        self.angle_plot.set_ylim(-fac, fac)
        self.angle_plot.grid(b=True, which='both')
        self.angle_plot.xaxis.set_minor_locator(AutoMinorLocator(2))

    def phase_plot_init(self):
        self.phase_plot.set_facecolor(Color.PLOT_BACKGROUND)
        self.phase_plot.set_xlabel('$\Theta$ (rad)')
        self.phase_plot.set_ylabel('$\omega$')
        self.phase_plot.set_title('Stroboscobic/phase space plot')

        fac_1 = np.max(np.abs(self.parent.solutions[0][2]))*(1+0.1) if np.max(np.abs(self.parent.solutions[0][2])) < 5 else 5
        fac_2 = np.max(np.abs(self.parent.solutions[0][1])) * (1 + 0.1) if np.max(np.abs(self.parent.solutions[0][1]))*(1+0.1) < np.pi else np.pi
        fac = fac_1 if fac_1 > fac_2 else fac_2
        self.phase_plot.set_xlim(-fac, fac)
        self.phase_plot.set_ylim(-fac, fac)
        self.phase_plot.grid(b=True, which='both')
        self.phase_plot.xaxis.set_major_locator(MultipleLocator(2 * fac / 4))
        self.phase_plot.yaxis.set_major_locator(MultipleLocator(2 * fac / 4))
        self.phase_plot.xaxis.set_minor_locator(AutoMinorLocator(2))
        self.phase_plot.yaxis.set_minor_locator(AutoMinorLocator(2))

    def vel_plot_init(self):
        self.vel_plot.set_facecolor(Color.PLOT_BACKGROUND)
        self.vel_plot.set_xlabel('$t$')
        self.vel_plot.set_ylabel('$\omega$')
        self.vel_plot.set_title('Angular velocity $\omega$ over time $t$')

        if self.parent.show_second:
            fac = (np.max(np.abs(np.concatenate((self.parent.solutions[0][2], self.parent.solutions[1][2]), axis=0))) * (1 + 0.1)
                   if np.max(np.abs(np.concatenate((self.parent.solutions[0][2], self.parent.solutions[1][2]), axis=0))) < 100
                   else 100)
        else:
            fac = (np.max(np.abs(self.parent.solutions[0][2]))*(1+0.05) if np.max(np.abs(self.parent.solutions[0][2])) < 100 else 100)
        self.vel_plot.set_xlim(0, T_MAX)
        self.vel_plot.set_ylim(-fac, fac)
        self.vel_plot.grid(b=True, which='both')
        self.vel_plot.xaxis.set_minor_locator(AutoMinorLocator(2))

    def angle_diff_plot_init(self):
        self.angle_plot_diff.set_facecolor(Color.PLOT_BACKGROUND)
        self.angle_plot_diff.set_xlabel('$t$')
        self.angle_plot_diff.set_ylabel('$\Delta \Theta$')
        self.angle_plot_diff.set_title('Angular difference $\Delta \Theta$')

        fac1 = np.max(self.parent.angle_diffs) * (1 + 0.25) if np.max(self.parent.angle_diffs) < 100 else 100
        fac2 = np.min(self.parent.angle_diffs) * (1 + 0.25)

        self.angle_plot_diff.set_yscale("log")
        self.angle_plot_diff.set_xlim(0, T_MAX)
        self.angle_plot_diff.set_ylim(fac2, fac1+10)
        self.angle_plot_diff.grid(b=True, which='both')
        #self.angle_plot_diff.xaxis.set_minor_locator(AutoMinorLocator(2))

    def vel_diff_plot_init(self):
        self.vel_plot_diff.set_facecolor(Color.PLOT_BACKGROUND)
        self.vel_plot_diff.set_xlabel('$t$')
        self.vel_plot_diff.set_ylabel('$\Delta \omega$', rotation='horizontal')
        self.vel_plot_diff.set_title('Angular velocity difference $\Delta \omega$')
        self.vel_plot_diff.yaxis.set_label_coords(-0.03, 1.02)

        fac = np.max(np.abs(self.parent.vel_diffs)) * (1 + 0.15) if np.abs(np.max(self.parent.vel_diffs)) < 100 else 100


        self.vel_plot_diff.set_xlim(0, T_MAX)
        self.vel_plot_diff.set_ylim(-fac, fac)
        self.vel_plot_diff.grid(b=True, which='both')
        self.vel_plot_diff.xaxis.set_minor_locator(AutoMinorLocator(2))