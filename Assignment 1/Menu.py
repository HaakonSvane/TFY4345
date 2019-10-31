'''
    Animation plotter with interactive parameter buttons
    For Assignment 1, TFY4345
    Created by: Haakon Svane
    Date: 9. September, 2019
'''

from matplotlib.animation import FuncAnimation
from matplotlib import pyplot as plt

from AnimatorHelper import *
from Solver import *


# TODO: Remove instance-variables that are not instance-dependent

class Menu:
    def __init__(self, color_sch="default"):
        set_color_scheme(color_sch)
        self.color_scheme = color_sch

        self.ani = None

        self.bob_mass = init_bob_mass
        self.rod_length = init_rod_length
        self.dt = init_dt
        self.show_analytical = False
        self.method_index = init_method_index
        self.speed_index = init_speed_index
        self.initialized = False

        self.s = AxisSetup(self)
        self.s.plot_setup()

        self.aa_calc = AnimatedActors(self)
        self.aa_anal = AnimatedActors(self, "darkgreen")
        self.acts = (self.aa_calc, self.aa_anal)
        self.artists = None

        self.sol_euler = Solver("euler")
        self.sol_euler_cromer = Solver("euler_cromer")
        self.sol_runge_kutta = Solver("runge_kutta")
        self.sol_analytical = Solver("analytical")
        self.sol_dummy = Solver("dummy")

        self.solution = None
        self.solution_analytical = None
        self.angle_diffs = None
        self.energy_diffs = None

        self.plot_indexes = [0, 0]

    def _ani_init(self):
        if not self.initialized:
            self.initialized = True
            if self.method_index == -1:
                self.solution = self.sol_dummy.solve(self.rod_length, self.bob_mass, self.dt)
            elif self.method_index == 0:
                self.solution = self.sol_euler.solve(self.rod_length, self.bob_mass, self.dt)
            elif self.method_index == 1:
                self.solution = self.sol_euler_cromer.solve(self.rod_length, self.bob_mass, self.dt)
            elif self.method_index == 2:
                self.solution = self.sol_runge_kutta.solve(self.rod_length, self.bob_mass, self.dt)
            else:
                print("Error in Animator: method_index '%s' does not point to any solver method. "
                      "Reverting to default: 'euler'" % self.method_index)
                self.solution = self.sol_euler.solve(self.rod_length, self.bob_mass, self.dt)
            if self.show_analytical:
                self.solution_analytical = self.sol_analytical.solve(self.rod_length, self.bob_mass, self.dt)

            self.solutions = [self.solution, self.solution_analytical]
            self.s.pendulum_plot_init()
            self.s.angle_plot_init()
            self.s.phase_plot_init()
            self.s.energy_plot_init()
            if self.show_analytical:
                self.angle_diffs, self.energy_diffs = self.sol_analytical.calc_differences(self.solution[1],
                                                                                           self.solution[3])
                self.s.angle_diff_plot_init()
                self.s.energy_diff_plot_init()

            self.aa_calc.init_pendulum_artists()
            self.aa_calc.attach_pendulum_artists(self.s.pendulum_plot)
            self.aa_calc.init_angle_artists(self.s.angle_plot)
            self.aa_calc.init_phase_artists(self.s.phase_plot)
            self.aa_calc.init_energy_artists(self.s.energy_plot)

            self.aa_anal.init_pendulum_artists()
            self.aa_anal.attach_pendulum_artists(self.s.pendulum_plot)
            if not self.show_analytical:
                self.aa_anal.rod.set_visible(False)
                self.aa_anal.circle.set_visible(False)
            self.aa_anal.init_angle_artists(self.s.angle_plot)
            self.aa_anal.init_phase_artists(self.s.phase_plot)
            self.aa_anal.init_energy_artists(self.s.energy_plot)
            if self.show_analytical:
                self.aa_anal.init_angle_diff_artists(self.s.angle_plot_diff)
                self.aa_anal.init_energy_diff_artists(self.s.energy_plot_diff)

            self.artists = (*self.aa_calc.get_artists(), *self.aa_anal.get_artists())

        return self.artists

    def _ani_update(self, anim_t):
        for i, j in enumerate([n for n in self.solutions if n]):

            # Checks whether the animation time has surpassed the data time or not.
            # Correctly adjusts for frame-skips
            if anim_t >= self.solutions[i][0][self.plot_indexes[i]]:
                self.plot_indexes[i] += int(
                    (anim_t - self.solutions[i][0][self.plot_indexes[i]]) //(self.dt if i == 0 else (1/FRAME_RATE)))
                if self.plot_indexes[i] >= self.solutions[i][0].size or self.speed_index == 4:
                    self.plot_indexes[i] = self.solutions[i][0].size - 1

            self.acts[i].line_angle.set_data(j[0][:self.plot_indexes[i]], j[1][:self.plot_indexes[i]])
            self.acts[i].line_phase.set_data(j[1][:self.plot_indexes[i]], j[2][:self.plot_indexes[i]])
            self.acts[i].circle.center = (self.rod_length * np.sin(j[1][self.plot_indexes[i]]),
                                          -self.rod_length * np.cos(j[1][self.plot_indexes[i]]))
            self.acts[i].rod.set_data([pivot_coords[0], self.acts[i].circle.get_center()[0]],
                                      [pivot_coords[1], self.acts[i].circle.get_center()[1]])
            self.acts[i].line_energy.set_data(j[0][:self.plot_indexes[i]], j[3][:self.plot_indexes[i]])

            if self.show_analytical and i == 0:
                self.acts[1].line_angle_diff.set_data(j[0][:self.plot_indexes[i]], self.angle_diffs[:self.plot_indexes[i]])
                self.acts[1].line_energy_diff.set_data(j[0][:self.plot_indexes[i]], self.energy_diffs[:self.plot_indexes[i]])


        return self.artists

    def _button_reset(self, event):
        self.bob_mass = init_bob_mass
        self.rod_length = init_rod_length
        self.dt = init_dt
        self.method_index = init_method_index
        self.speed_index = init_speed_index
        self.run = 0
        self.show_analytical = False
        self.reset()

    def _button_apply(self, event):
        self.bob_mass = self.s.slider_bob_mass.val
        self.rod_length = self.s.slider_rod_length.val
        self.dt = self.s.slider_dt.val
        if not self.method_index == -1:
            self.method_index = [i.get_text() for i in self.s.radio_calc_method.labels].index(
                self.s.radio_calc_method.value_selected)
        self.speed_index = [i.get_text() for i in self.s.radio_anim_params.labels].index(
            self.s.radio_anim_params.value_selected)
        self.run = 0
        self.show_analytical = [i.get_text() for i in self.s.radio_show_calc.labels].index(
            self.s.radio_show_calc.value_selected)
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

        del self.aa_calc
        self.aa_calc = None
        del self.aa_anal
        self.aa_anal = None

        del self.ani
        self.ani = None
        self.initialized = False
        self.solution = None
        self.solution_analytical = None
        if self.show_analytical:
            rcParams["figure.figsize"][1] *= 1.4
        else:
            rcParams["figure.figsize"] = [i / 10 * FIGURE_SIZE_FACTOR for i in FIGURE_RATIO]
        self.plot_indexes = [0, 0]

        self.s = AxisSetup(self)
        self.s.plot_setup()

        self.aa_calc = AnimatedActors(self)
        self.aa_anal = AnimatedActors(self, "darkgreen")
        self.acts = (self.aa_calc, self.aa_anal)
        self.artists = None

        self.start_ani()
