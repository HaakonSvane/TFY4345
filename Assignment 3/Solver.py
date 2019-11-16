'''
    Solver and timer class to time and numerically solve the differential equations
    For Assignment 3, TFY4345
    Created by: Haakon Svane, Birk Karlsen-Bæck
    Date: 13. November, 2019
'''

import time
import math
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
    def __init__(self, method="runge_kutta"):
        self.function = None

        self.method = method
        self.dampening = None

        if self.method == "runge_kutta":
            self.function = self._runge_kutta
        elif self.method == "dummy":
            self.function = self._dummy_func
        else:
            print("Invalid solver '%s' given. Reverting to default: 'runge_kutta'" % self.method)
            self.function = self._runge_kutta

        self.t = Timer()

        self.mass = 0
        self.length = 0
        self.dt = 0
        self.q = 0
        self.omega_d = 0
        self.Fd = 0
        self.q_scale = 0.90

        self.t_arr = None
        self.y_arr = None
        self.y_vel_arr = None
        self.energy_arr = None
        self.strob_times = []

    def solve(self, init_theta, rod_len, bob_mass, dt, dampening, F_d, Om_d):
        self.init_theta = init_theta
        self.length = rod_len
        self.mass = bob_mass
        self.dt = dt
        self.dampening = dampening
        self.Fd = F_d
        self.omega_d = Om_d


        if self.dampening == "underdamped":
            self.q = 2 * np.sqrt(g_acc / self.length * (1 - self.q_scale))
        elif self.dampening == "overdamped":
            self.q = 2 * np.sqrt(g_acc / self.length * (1 + self.q_scale))
        elif self.dampening == "critical damping":
            self.q = 2 * np.sqrt(g_acc / self.length)
        elif self.dampening == "none":
            self.q = 0
        else:
            self.q = self.dampening

        self.t.start()
        print("\n---------Stats from using solver '%s' (%s)---------"
              "(mass=%.1f, length=%.1f, dt=%.3f, Nt=%s, Fd=%.1f, Od=%.1f, q=%.3f, O=%.3f)---------" % (
                  self.method,self.dampening, self.mass, self.length, self.dt, round(T_MAX / self.dt),
                  self.Fd, self.omega_d, self.q, np.sqrt(g_acc/self.length)))
        self.function()
        self.t.stop_and_get()
        print("-" * 145 + "")

        self.strob_times = []
        if self.omega_d:
            N = int(self.omega_d*T_MAX/(2*np.pi))
            for n in range(N+1):
                self.strob_times.append(self._find_nearest_time_value(2 * np.pi * n / self.omega_d))
            self.strob_times = np.array(self.strob_times)
        else:
            self.strob_times=np.arange(self.t_arr.size)

        return self.t_arr, self.y_arr, self.y_vel_arr

    def _calc_energy(self, y_array, y_vel_array):
        formula_energy = lambda y, y_vel: 1 / 2 * self.mass * self.length ** 2 * y_vel ** 2 + \
                                          self.mass * g_acc * self.length * (1 - np.cos(y))
        return formula_energy(y_array, y_vel_array)

    # Backend solvers.
    # Note that the array_definitions has been moved to these solvers so that time-measurement is more accurate.
    def _runge_kutta(self):
        self.t_arr = np.arange(0, T_MAX, self.dt)
        self.y_arr = np.empty(self.t_arr.size)
        self.y_vel_arr = np.empty(self.t_arr.size)

        # Set the initial conditions
        self.y_arr[0] = self.init_theta
        self.y_vel_arr[0] = init_theta_vel

        self.y_arr[0] = self.init_theta
        self.y_vel_arr[0] = init_theta_vel

        for i in range(1, self.t_arr.size):
            k1y = (self.y_vel_arr[i - 1]) * self.dt
            k1v = self.dt * (-(g_acc / self.length) * (np.sin(self.y_arr[i - 1])) - self.q * (
                self.y_vel_arr[i - 1]) + self.Fd * np.sin(
                self.omega_d * self.t_arr[i - 1]))

            k2y = (self.y_vel_arr[i - 1] + k1v * (1 / 2)) * self.dt
            k2v = self.dt * (-(g_acc / self.length) * np.sin(self.y_arr[i - 1] + (1 / 2) * k1y) - self.q * (
                    self.y_vel_arr[i - 1] + (1 / 2) * k1v) + self.Fd * np.sin(
                self.omega_d * (self.t_arr[i - 1] + self.dt / 2)))

            k3y = (self.y_vel_arr[i - 1] + k2v * (1 / 2)) * self.dt
            k3v = self.dt * (-(g_acc / self.length) * np.sin(self.y_arr[i - 1] + (1 / 2) * k2y) - self.q * (
                    self.y_vel_arr[i - 1] + (1 / 2) * k2v) + self.Fd * np.sin(
                self.omega_d * (self.t_arr[i - 1] + self.dt / 2)))

            k4y = (self.y_vel_arr[i - 1] + k3v) * self.dt
            k4v = self.dt * (-(g_acc / self.length) * np.sin(self.y_arr[i - 1] + k3y) - self.q * (
                    self.y_vel_arr[i - 1] + k3v) + self.Fd * np.sin(self.omega_d * (self.t_arr[i - 1] + self.dt)))

            self.y_arr[i] = self.y_arr[i - 1] + (1 / 6) * (k1y + 2 * k2y + 2 * k3y + k4y)
            self.y_vel_arr[i] = self.y_vel_arr[i - 1] + (1 / 6) * (k1v + 2 * k2v + 2 * k3v + k4v)


            if self.y_arr[i] > np.pi:
                self.y_arr[i] -= 2*np.pi
            elif self.y_arr[i] < -np.pi:
                self.y_arr[i] += 2*np.pi



    def _find_nearest_time_value(self, value):
        idx = np.searchsorted(self.t_arr, value, side="left")
        if idx > 0 and (idx == self.t_arr.size or math.fabs(value - self.t_arr[idx - 1]) < math.fabs(value - self.t_arr[idx])):
            return idx - 1
        else:
            return idx

    # Dummy function for debugging and development
    def _dummy_func(self):
        self.t_arr = np.arange(0, T_MAX, self.dt)
        self.y_arr = 3 * np.pi / 5 * np.cos(self.t_arr) * np.exp(-self.t_arr / 20)
        self.y_vel_arr = -3 * np.pi / 100 * np.exp(-self.t_arr / 20) * (20 * np.sin(self.t_arr) + np.cos(self.t_arr))
        self.energy_arr = self._calc_energy(self.y_arr, self.y_vel_arr)

    def calc_differences(self, calced_theta, calced_theta_vel):
        if not (calced_theta.size == calced_theta_vel.size):
            print("Error: The two arrays are not of the same size. Returning nothing.")
            return
        r = np.where(calced_theta-self.y_arr>=np.pi, 2*np.pi-(calced_theta-self.y_arr), calced_theta-self.y_arr)
        r = np.where(r<-np.pi, 2*np.pi+r, r)
        return np.abs(r), calced_theta_vel - self.y_vel_arr
