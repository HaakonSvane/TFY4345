'''
    Solver and timer class to time and numerically solve the differential equations
    For Assignment 1, TFY4345
    Created by: Haakon Svane, Birk Karlsen-Bæck
    Date: 17. September, 2019
'''

import time

import numpy as np

from Parameters import *


# Class to ease time-measurement of the different solvers
class Timer:
    def __init__(self):
        self.init_abs_time = 0
        self.init_proc_time = 0
        self.dt = 0

    def start(self):
        self.init_abs_time = time.perf_counter()
        self.init_proc_time = time.process_time()

    def reset(self):
        self.init_abs_time = 0
        self.init_proc_time = 0

    def stop_and_get(self):
        print(
            "\t\t\tPerformance ('absolute') time elpapsed:\t\t%s (s)\n\t\t\tProcessor time elapsed:\t\t\t\t\t\t%s (s)" %
            (round(time.perf_counter() - self.init_abs_time, 5), round(time.process_time() - self.init_proc_time, 5)))
        self.reset()


# Class containing the different solvers.
class Solver:
    def __init__(self, method="euler"):
        self.function = None
        self.method = method

        if self.method == "euler":
            self.function = self._euler
        elif self.method == "euler_cromer":
            self.function = self._euler_cromer
        elif self.method == "runge_kutta":
            self.function = self._runge_kutta
        elif self.method == "analytical":
            self.function = self._analytical_solution
        elif self.method == "dummy":
            self.function = self._dummy_func
        else:
            print("Invalid solver '%s' given. Reverting to default: 'euler'" % self.method)
            self.function = self._euler

        self.timed = True
        self.t = Timer()

        self.mass = 0
        self.length = 0
        self.dt = 0

    def solve(self, rod_len, bob_mass, dt):
        if self.timed:
            self.mass = bob_mass
            self.length = rod_len
            self.dt = dt

            self.t.start()
            (t_array, y_array, y_vel_array, energy_array) = self.function()
            print("\n---------Stats from using solver '%s'---------"
                  "(mass=%s, length=%s, dt=%s, total time steps=%s)---------" % (
                  self.method, self.mass, self.length, self.dt, round(T_MAX / self.dt)))
            self.t.stop_and_get()
            print("-" * 106 + "")
            return t_array, y_array, y_vel_array, energy_array

    # Backend solvers.
    # Note that the array_definitions has been moved to these solvers so that time-measurement is more accurate.
    def _euler(self):
        t_arr = np.arange(0, T_MAX, self.dt)
        y_arr = np.empty(t_arr.size)
        y_vel_arr = np.empty(t_arr.size)
        energy_arr = np.empty(t_arr.size)

        # Set the initial conditions
        y_arr[0] = init_theta
        y_vel_arr[0] = init_theta_vel

        # Calculate the differential equation using the Euler Method.
        for i in range(1, t_arr.size):
            y_arr[i] = y_arr[i - 1] + y_vel_arr[i - 1] * self.dt
            y_vel_arr[i] = y_vel_arr[i - 1] - (g_acc / self.length) * y_arr[i - 1] * self.dt

        formula_energy = lambda y, y_vel: 1 / 2 * self.mass * self.length ** 2 * y_vel ** 2 + self.mass * g_acc * self.length * (1 - np.cos(y))
        f_energy = np.vectorize(formula_energy, otypes=[np.float])
        energy_arr = f_energy(y_arr, y_vel_arr)

        return t_arr, y_arr, y_vel_arr, energy_arr

    def _euler_cromer(self):
        t_arr = np.arange(0, T_MAX, self.dt)
        y_arr = np.empty(t_arr.size)
        y_vel_arr = np.empty(t_arr.size)
        energy_arr = np.empty(t_arr.size)

        # Set the initial conditions
        y_arr[0] = init_theta
        y_vel_arr[0] = init_theta_vel

        # Calculate the differential equation using the Euler-Cromer Method
        for i in range(1, t_arr.size):
            y_vel_arr[i] = y_vel_arr[i - 1] - (g_acc / self.length) * y_arr[i - 1] * self.dt
            y_arr[i] = y_arr[i - 1] + y_vel_arr[i] * self.dt

        formula_energy = lambda y, y_vel: 1 / 2 * self.mass * self.length ** 2 * y_vel ** 2 + self.mass * g_acc * self.length * (
                    1 - np.cos(y))
        f_energy = np.vectorize(formula_energy, otypes=[np.float])
        energy_arr = f_energy(y_arr, y_vel_arr)

        return t_arr, y_arr, y_vel_arr, energy_arr

    def _runge_kutta(self):
        t_arr = np.arange(0, T_MAX, self.dt)
        y_arr = np.empty(t_arr.size)
        y_vel_arr = np.empty(t_arr.size)
        energy_arr = np.empty(t_arr.size)

        # Set the initial conditions
        y_arr[0] = init_theta
        y_vel_arr[0] = init_theta_vel

        for i in range(1, t_arr.size):
            k1y = (y_vel_arr[i - 1]) * self.dt
            k1v = - (g_acc / self.length) * y_arr[i - 1] * self.dt

            k2y = (y_vel_arr[i - 1] + k1v * (1 /2)) * self.dt
            k2v = - (g_acc / self.length) * (y_arr[i - 1] + k1y * (1 / 2)) * self.dt

            k3y = (y_vel_arr[i - 1] + k2v * (1 / 2)) * self.dt
            k3v = - (g_acc / self.length) * (y_arr[i - 1] + k2y * (1 / 2)) * self.dt

            k4y = (y_vel_arr[i - 1] + k3v) * self.dt
            k4v = - (g_acc / self.length) * (y_arr[i - 1] + k3y) * self.dt

            y_arr[i] = y_arr[i - 1] + (1 / 6) * (k1y + 2 * k2y + 2 * k3y + k4y)
            y_vel_arr[i] = y_vel_arr[i - 1] + (1 / 6) * (k1v + 2 * k2v + 2 * k3v + k4v)

        formula_energy = lambda y,y_vel: 1 / 2 * self.mass * self.length ** 2 * y_vel ** 2 + self.mass * g_acc * self.length * (1 - np.cos(y))
        f_energy = np.vectorize(formula_energy, otypes=[np.float])
        energy_arr = f_energy(y_arr, y_vel_arr)

        return t_arr, y_arr, y_vel_arr, energy_arr

    def _analytical_solution(self):
        formula_y = lambda t: init_theta*np.cos(np.sqrt(g_acc/self.length)*t)
        formula_y_der = lambda t: -init_theta*np.sqrt(g_acc/self.length)*np.sin(np.sqrt(g_acc/self.length)*t)
        formula_energy = lambda y, y_vel: 1/2*self.mass*self.length**2*y_vel**2 + self.mass*g_acc*self.length*(1-np.cos(y))

        f = np.vectorize(formula_y, otypes=[np.float])
        f_der = np.vectorize(formula_y_der, otypes=[np.float])
        f_energy = np.vectorize(formula_energy, otypes=[np.float])

        t_arr = np.arange(0, T_MAX, 1/FRAME_RATE)
        y_arr = f(t_arr)
        y_vel_arr = f_der(t_arr)
        energy_arr = f_energy(y_arr, y_vel_arr)

        return t_arr, y_arr, y_vel_arr, energy_arr
        

    # Dummy function for debugging and development
    def _dummy_func(self):
        formula_energy = lambda y,y_vel: 1 / 2 * self.mass * self.length ** 2 * y_vel ** 2 + \
                                         self.mass * g_acc * self.length * (1 - np.cos(y))
        f_energy = np.vectorize(formula_energy, otypes=[np.float])
        t_arr = np.arange(0, T_MAX, self.dt)
        y_arr = 3*np.pi/5 * np.cos(t_arr) * np.exp(-t_arr / 20)
        y_vel_arr = -3*np.pi/100 * np.exp(-t_arr / 20) * (20 * np.sin(t_arr) + np.cos(t_arr))
        energy_arr = f_energy(y_arr, y_vel_arr)

        return t_arr, y_arr, y_vel_arr, energy_arr

    def calc_differences(self, calced_angles, calced_energies):
        if not (calced_angles.size == calced_energies.size):
            print("Error: The two arrays are not of the same size. Returning nothing.")
            return

        t_arr = np.linspace(0, T_MAX, calced_angles.size)
        formula_y = lambda t: init_theta * np.cos(np.sqrt(g_acc / self.length) * t)
        formula_y_der = lambda t: -init_theta * np.sqrt(g_acc / self.length) * np.sin(np.sqrt(g_acc / self.length) * t)
        formula_energy = lambda y, y_vel: 1 / 2 * self.mass * self.length ** 2 * y_vel ** 2 + self.mass * g_acc * self.length * (1 - np.cos(y))

        f = np.vectorize(formula_y, otypes=[np.float])
        f_der = np.vectorize(formula_y_der, otypes=[np.float])
        f_energy = np.vectorize(formula_energy, otypes=[np.float])

        y_arr = f(t_arr)
        y_vel_arr = f_der(t_arr)
        energy_arr = f_energy(y_arr, y_vel_arr)

        return (calced_angles-y_arr, calced_energies-energy_arr)
