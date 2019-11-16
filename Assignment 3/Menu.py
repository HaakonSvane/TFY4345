'''
    Animation plotter with interactive parameter buttons
    For Assignment 3, TFY4345
    Created by: Haakon Svane
    Date: 13. November, 2019
'''

from matplotlib.animation import FuncAnimation
from matplotlib import pyplot as plt

from AnimatorHelper import *
from Solver import *


# TODO: Remove instance-variables that are not instance-dependent
# TODO: Fix frame stuff

class Menu:
    def __init__(self, color_sch="default"):
        set_color_scheme(color_sch)
        self.color_scheme = color_sch

        self.ani = None

        self.bob_mass = init_bob_mass
        self.rod_length = init_rod_length
        self.dt = init_dt
        self.Om_d = init_Om_d
        self.F_d = init_F_d
        self.show_second = False
        self.dampening_index = init_dampening_index
        self.speed_index = init_speed_index
        self.init_theta = init_theta
        self.init_Dtheta = init_Dtheta
        self.initialized = False

        self.s = AxisSetup(self)
        self.s.plot_setup()

        self.aa_first = AnimatedActors(self)
        self.aa_second = AnimatedActors(self, "darkgreen")
        self.acts = (self.aa_first, self.aa_second)
        self.artists = None

        self.solve_first = Solver(method="runge_kutta")
        self.solve_second = Solver(method="runge_kutta")
        self.solve_dummy = Solver(method="dummy")

        self.solution_first = None
        self.solution_second = None
        self.angle_diffs = None
        self.vel_diffs = None

        self.plot_indexes = [0, 0]

    def _ani_init(self):
        if not self.initialized:
            self.initialized = True
            self.solution_first = self.solve_first.solve(self.init_theta, self.rod_length,
                                                         self.bob_mass, self.dt,
                                                         0.5,
                                                         self.F_d, self.Om_d)
            if self.show_second:
                self.solution_second = self.solve_second.solve(self.init_theta+self.init_Dtheta,self.rod_length,
                                                               self.bob_mass, self.dt,
                                                               0.5,
                                                               self.F_d, self.Om_d)

            self.solutions = [self.solution_first, self.solution_second]
            self.s.pendulum_plot_init()
            self.s.angle_plot_init()
            self.s.phase_plot_init()
            self.s.vel_plot_init()
            if self.show_second:
                self.angle_diffs, self.vel_diffs = self.solve_second.calc_differences(self.solution_first[1],
                                                                                         self.solution_first[2])
                self.s.angle_diff_plot_init()
                self.s.vel_diff_plot_init()

            self.aa_first.init_pendulum_artists()
            self.aa_first.attach_pendulum_artists(self.s.pendulum_plot)
            self.aa_first.init_angle_artists(self.s.angle_plot)
            self.aa_first.init_phase_artists(self.s.phase_plot)
            self.aa_first.init_strob_artists(self.s.phase_plot)
            self.aa_first.init_vel_artists(self.s.vel_plot)

            self.aa_second.init_pendulum_artists()
            self.aa_second.attach_pendulum_artists(self.s.pendulum_plot)
            if not self.show_second:
                self.aa_second.rod.set_visible(False)
                self.aa_second.circle.set_visible(False)
            self.aa_second.init_angle_artists(self.s.angle_plot)
            self.aa_second.init_phase_artists(self.s.phase_plot)
            self.aa_second.init_strob_artists(self.s.phase_plot)
            self.aa_second.init_vel_artists(self.s.vel_plot)
            if self.show_second:
                self.aa_second.init_angle_diff_artists(self.s.angle_plot_diff)
                self.aa_second.init_vel_diff_artists(self.s.vel_plot_diff)

            self.artists = (*self.aa_first.get_artists(), *self.aa_second.get_artists())

        return self.artists

    def _ani_update(self, anim_t):
        for i, j in enumerate([n for n in self.solutions if n]):
            # Checks whether the animation time has surpassed the data time or not.
            # Correctly adjusts for frame-skips
            if anim_t >= self.solutions[i][0][self.plot_indexes[i]]:
                self.plot_indexes[i] += int(
                    (anim_t - self.solutions[i][0][self.plot_indexes[i]]) // (self.dt if i == 0 else (1/FRAME_RATE)))
                if self.plot_indexes[i] >= self.solutions[i][0].size or self.speed_index == 4:
                    self.plot_indexes[i] = self.solutions[i][0].size - 1

            self.acts[i].line_angle.set_data(j[0][:self.plot_indexes[i]], j[1][:self.plot_indexes[i]])

            self.acts[i].line_phase.set_data(j[1][:self.plot_indexes[i]], j[2][:self.plot_indexes[i]])
            a = self.solve_first.strob_times if i == 0 else self.solve_second.strob_times
            ind = a[np.where(a <= self.plot_indexes[i])[0]]
            self.acts[i].line_strob.set_data(j[1][ind], j[2][ind])

            self.acts[i].circle.center = (self.rod_length * np.sin(j[1][self.plot_indexes[i]]),
                                          -self.rod_length * np.cos(j[1][self.plot_indexes[i]]))
            self.acts[i].rod.set_data([pivot_coords[0], self.acts[i].circle.get_center()[0]],
                                      [pivot_coords[1], self.acts[i].circle.get_center()[1]])
            self.acts[i].line_vel.set_data(j[0][:self.plot_indexes[i]], j[2][:self.plot_indexes[i]])

            if self.show_second and i == 0:
                self.acts[1].line_angle_diff.set_data(j[0][:self.plot_indexes[i]], self.angle_diffs[:self.plot_indexes[i]])
                self.acts[1].line_vel_diff.set_data(j[0][:self.plot_indexes[i]], self.vel_diffs[:self.plot_indexes[i]])


        return self.artists

    def _button_reset(self, event):
        self.bob_mass = init_bob_mass
        self.rod_length = init_rod_length
        self.dt = init_dt
        self.F_d = init_F_d
        self.Om_d = init_Om_d
        self.dampening_index = init_dampening_index
        self.speed_index = init_speed_index
        self.init_Dtheta = init_Dtheta
        self.show_second = False
        self.reset()

    def _button_apply(self, event):
        self.bob_mass = self.s.slider_bob_mass.val
        self.rod_length = self.s.slider_rod_length.val
        self.dt = self.s.slider_dt.val
        self.F_d = self.s.slider_Fd.val
        self.Om_d = self.s.slider_Om_d.val
        self.init_Dtheta = self.s.slider_Dtheta.val
        if not self.dampening_index == -1:
            self.dampening_index = [i.get_text() for i in self.s.radio_damping.labels].index(
                self.s.radio_damping.value_selected)
        self.speed_index = [i.get_text() for i in self.s.radio_anim_params.labels].index(
            self.s.radio_anim_params.value_selected)
        self.show_second = [i.get_text() for i in self.s.radio_second_calc.labels].index(
            self.s.radio_second_calc.value_selected)
        self.reset()

    def _speed_from_index(self, index):
        if index == 0:
            return 2
        elif index == 1:
            return 1
        elif index == 2:
            return 1 / 4
        elif index == 3:
            return 1 / 10
        elif index == 4:
            return 1 / (T_MAX * FRAME_RATE) + 1  # If the user wants no animation, create only two frames
        else:
            print("Invalid index %s given. Setting to default (1).." % index)
            return 1

    def start_ani(self):
        self.ani = FuncAnimation(self.s.fig, self._ani_update,
                                 frames=np.linspace(0, T_MAX,
                                                    T_MAX * FRAME_RATE * self._speed_from_index(self.speed_index)),
                                 init_func=self._ani_init, blit=True, interval=(1 / FRAME_RATE) * 1e3, repeat=False)
        self.s.button_reset_apply.on_clicked(self._button_reset)
        self.s.button_apply_run.on_clicked(self._button_apply)
        plt.show()

    def reset(self):
        self.ani.event_source.stop()
        plt.close(self.s.fig)

        del self.s
        self.s = None

        del self.aa_first
        self.aa_first = None
        del self.aa_second
        self.aa_second = None

        del self.ani
        self.ani = None
        self.initialized = False
        self.solution_first = None
        self.solution_second = None
        if self.show_second:
            rcParams["figure.figsize"][1] *= 1.4
        else:
            rcParams["figure.figsize"] = [i / 10 * FIGURE_SIZE_FACTOR for i in FIGURE_RATIO]
        self.plot_indexes = [0, 0]

        self.s = AxisSetup(self)
        self.s.plot_setup()

        self.aa_first = AnimatedActors(self)
        self.aa_second = AnimatedActors(self, "darkgreen")
        self.acts = (self.aa_first, self.aa_second)
        self.artists = None

        self.start_ani()
