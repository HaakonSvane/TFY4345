'''
    Solver and timer class to time and numerically solve the differential equations
    For Assignment 2, TFY4345
    Created by: Haakon Svane, Birk Karlsen-Bæck
    Date: 30. October, 2019
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
    def __init__(self, method="runge_kutta"):
        self.function = None

        self.method = method
        self.dampening = None

        if self.method == "runge_kutta":
            self.function = self._runge_kutta
        elif self.method == "analytical":
            self.function = self._analytical_solution
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
        self.omega_d = 1
        self.Fd = 0
        self.q_scale = 0.90

    def solve(self, rod_len, bob_mass, dt, dampening, F_d, Om_d):
        self.mass = bob_mass
        self.length = rod_len
        self.dt = dt
        self.Fd = F_d
        self.omega_d = Om_d
        self.dampening = dampening

        if self.dampening == "underdamped":
            self.q = 2 * np.sqrt(g_acc / self.length * (1 - self.q_scale))
        elif self.dampening == "overdamped":
            self.q = 2 * np.sqrt(g_acc / self.length * (1 + self.q_scale))
        elif self.dampening == "critical damping":
            self.q = 2 * np.sqrt(g_acc / self.length)
        elif self.dampening == "none":
            self.q = 0
        else:
            print("DAMPENING %s DOES NOT EXIST!" % dampening)

        self.t.start()
        print("\n---------Stats from using solver '%s' (%s)---------"
              "(mass=%.1f, length=%.1f, dt=%.3f, Nt=%s, Fd=%.1f, Od=%.1f, q=%.3f, O=%.3f)---------" % (
                  self.method,self.dampening, self.mass, self.length, self.dt, round(T_MAX / self.dt),
                  self.Fd, self.omega_d, self.q, np.sqrt(g_acc/self.length)))
        self.function()
        self.t.stop_and_get()
        print("-" * 145 + "")
        return self.t_arr, self.y_arr, self.y_vel_arr, self.energy_arr

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
        energy_arr = np.empty(self.t_arr.size)

        # Set the initial conditions
        self.y_arr[0] = init_theta
        self.y_vel_arr[0] = init_theta_vel

        for i in range(1, self.t_arr.size):
            k1y = (self.y_vel_arr[i - 1]) * self.dt
            k1v = self.dt * (-(g_acc / self.length) * (self.y_arr[i - 1]) - self.q * (
            self.y_vel_arr[i - 1]) + self.Fd * np.sin(
                self.omega_d * self.t_arr[i - 1]))

            k2y = (self.y_vel_arr[i - 1] + k1v * (1 / 2)) * self.dt
            k2v = self.dt * (-(g_acc / self.length) * (self.y_arr[i - 1] + (1 / 2) * k1y) - self.q * (
                    self.y_vel_arr[i - 1] + (1 / 2) * k1v) + self.Fd * np.sin(
                self.omega_d * (self.t_arr[i - 1] + self.dt / 2)))

            k3y = (self.y_vel_arr[i - 1] + k2v * (1 / 2)) * self.dt
            k3v = self.dt * (-(g_acc / self.length) * (self.y_arr[i - 1] + (1 / 2) * k2y) - self.q * (
                    self.y_vel_arr[i - 1] + (1 / 2) * k2v) + self.Fd * np.sin(
                self.omega_d * (self.t_arr[i - 1] + self.dt / 2)))

            k4y = (self.y_vel_arr[i - 1] + k3v) * self.dt
            k4v = self.dt * (-(g_acc / self.length) * (self.y_arr[i - 1] + k3y) - self.q * (
                    self.y_vel_arr[i - 1] + k3v) + self.Fd * np.sin(self.omega_d * (self.t_arr[i - 1] + self.dt)))

            self.y_arr[i] = self.y_arr[i - 1] + (1 / 6) * (k1y + 2 * k2y + 2 * k3y + k4y)
            self.y_vel_arr[i] = self.y_vel_arr[i - 1] + (1 / 6) * (k1v + 2 * k2v + 2 * k3v + k4v)

        self.energy_arr = self._calc_energy(self.y_arr, self.y_vel_arr)

        return self.t_arr, self.y_arr, self.y_vel_arr, self.energy_arr

    def _analytical_solution(self):

        omega = g_acc / self.length
        if self.dampening == "critical damping":
            A1 = init_theta + (16 * self.Fd * self.omega_d * self.q) / ((self.q ** 2 + 4 * self.omega_d ** 2) ** 2)
            A2 = init_theta_vel + (1 / 2) * self.q * A1 + (4 * self.Fd * self.omega_d) / (
                        4 * self.omega_d ** 2 + self.q ** 2) - (
                         8 * self.Fd * self.omega_d * self.q ** 2) / ((self.q ** 2 + 4 * self.omega_d ** 2) ** 2)

            formula_y = lambda t: A1 * np.exp(-(self.q * t) / 2) + A2 * t * np.exp(-(self.q * t) / 2) - (
                        (2 * self.Fd) / ((self.q ** 2 + 4 * self.omega_d ** 2) ** 2)) * (
                                          (
                                                      self.q ** 3 * t - 2 * self.q ** 2 + 4 * self.q * t * self.omega_d ** 2 + 8 * self.omega_d ** 2) * np.sin(
                                      self.omega_d * t)
                                          - 2 * self.omega_d * (
                                                      self.q ** 2 * t - 4 * self.q + 4 * t * self.omega_d ** 2) * np.cos(
                                      self.omega_d * t)) + (
                                          (2 * self.Fd * t) / (self.q ** 2 + 4 * self.omega_d ** 2)) * (
                                              self.q * np.sin(self.omega_d * t) - 2 * self.omega_d * np.cos(
                                          self.omega_d * t))

            formula_y_der = lambda t: -(1 / 2) * self.q * A1 * np.exp(-(self.q * t) / 2) - (
                        1 / 2) * self.q * A2 * t * np.exp(-(self.q * t) / 2) + A2 * np.exp(
                -(self.q * t) / 2) - ((2 * self.Fd) / ((4 * self.omega_d ** 2 + self.q ** 2) ** 2)) * (
                                                  self.omega_d * self.q *
                                                  (
                                                              4 * self.omega_d ** 2 * t + self.q ** 2 * t - 4 * self.q) * np.cos(
                                              self.omega_d * t) + (8 * self.omega_d ** 4 * t
                                                                   + 2 * self.omega_d ** 2 * self.q * (
                                                                               self.q * t - 2) + self.q ** 3) * np.sin(
                                              self.omega_d * t)) + (
                                              (2 * self.Fd) / (4 * self.omega_d ** 2 + self.q ** 2)) * (
                                                  (2 * self.omega_d ** 2 * t + self.q) * np.sin(self.omega_d * t)
                                                  + self.omega_d * (self.q * t - 2) * np.cos(self.omega_d * t))

        elif self.dampening == "overdamped":
            phi = np.sqrt(self.q ** 2 - 4 * omega)
            k1 = (1 / 2) * (phi - self.q)
            k2 = (1 / 2) * (phi + self.q)
            A_ = (4 * self.Fd * self.omega_d) / (phi * ((self.q - phi) ** 2 + 4 * self.omega_d ** 2))
            B_ = (4 * self.Fd * self.omega_d) / (phi * ((self.q + phi) ** 2 + 4 * self.omega_d ** 2))
            C_ = (2 * self.Fd * self.omega_d * (self.q - phi)) / (
                        phi * ((self.q - phi) ** 2 + 4 * self.omega_d ** 2))
            D_ = (2 * self.Fd * self.omega_d * (self.q + phi)) / (
                        phi * ((self.q + phi) ** 2 + 4 * self.omega_d ** 2))
            A1 = (A_ * k2 - B_ * k2 - C_ + D_ + k2 * init_theta + init_theta_vel) / (k1 + k2)
            A2 = (A_ * k1 - B_ * k1 + C_ - D_ + k1 * init_theta - init_theta_vel) / (k1 + k2)

            formula_y = lambda t: A1 * np.exp((phi - self.q) * t / 2) + A2 * np.exp(-(phi + self.q) * t / 2) + (
                    (2 * self.Fd) / (phi * ((self.q - phi) ** 2 + 4 * self.omega_d ** 2))) * (
                                          (self.q - phi) * np.sin(self.omega_d * t) - 2 * self.omega_d * np.cos(
                                      self.omega_d * t)) - (
                                          (2 * self.Fd) / (phi * ((self.q + phi) ** 2 + 4 * self.omega_d ** 2))) * (
                                          (self.q + phi) * np.sin(self.omega_d * t) - 2 * self.omega_d * np.cos(
                                      self.omega_d * t))

            formula_y_der = lambda t: (1 / 2) * (phi - self.q) * A1 * np.exp((phi - self.q) * t / 2) - (1 / 2) * (
                        phi + self.q) * A2 * np.exp(
                -(phi + self.q) * t / 2) + (
                                              (2 * self.Fd * self.omega_d) / (
                                                  phi * ((self.q - phi) ** 2 + 4 * self.omega_d ** 2))) * (
                                              (self.q - phi) * np.cos(self.omega_d * t) + 2 * self.omega_d * np.sin(
                                          self.omega_d * t)) - (
                                              (2 * self.Fd * self.omega_d) / (
                                                  phi * ((self.q + phi) ** 2 + 4 * self.omega_d ** 2))) * (
                                              (self.q + phi) * np.cos(self.omega_d * t) + 2 * self.omega_d * np.sin(
                                          self.omega_d * t))

        elif self.dampening == "underdamped" or self.dampening == "none":
            phi = np.sqrt(omega - (1 / 4) * self.q ** 2)
            K1 = self.q ** 2 + 4 * (phi - self.omega_d) ** 2
            K2 = self.q ** 2 + 4 * (phi + self.omega_d) ** 2
            A1 = init_theta + ((self.Fd * self.q) / phi) * (1 / K1 - 1 / K2)
            A2 = (init_theta_vel / phi + (1 / 2) * ((self.q * A1) / phi) +
                  ((2 * self.Fd) / (phi ** 2)) * ((phi - self.omega_d) ** 2 / K1 - (phi + self.omega_d) ** 2 / K2) -
                  ((2 * self.Fd) / phi) * ((phi - self.omega_d) / K1 - (phi + self.omega_d) / K2))

            formula_y = lambda t: A1 * np.exp(-(self.q * t) / 2) * np.cos(phi * t) + A2 * np.exp(
                -(self.q * t) / 2) * np.sin(phi * t) - (
                                          self.Fd / phi) * np.cos(phi * t) * (
                                          (self.q * np.cos((phi - self.omega_d) * t) + 2 * (
                                                      phi - self.omega_d) * np.sin((phi - self.omega_d) * t)) / K1
                                          - (self.q * np.cos((self.omega_d + phi) * t) + 2 * (
                                              self.omega_d + phi) * np.sin((self.omega_d + phi) * t)) / K2) + (
                                          self.Fd / phi) * np.sin(phi * t) * (
                                          (2 * (phi - self.omega_d) * np.cos(
                                              (phi - self.omega_d) * t) - self.q * np.sin(
                                              (phi - self.omega_d) * t)) / K1
                                          + (self.q * np.sin((self.omega_d + phi) * t) - 2 * (
                                              self.omega_d + phi) * np.cos((self.omega_d + phi) * t)) / K2)

            formula_y_der = lambda t: (- (1 / 2) * self.q * A1 * np.exp(-(self.q * t) / 2) * np.cos(phi * t) -
                                       A1 * np.exp(-(self.q * t) / 2) * phi * np.sin(phi * t) -
                                       (1 / 2) * self.q * A2 * np.exp(-(self.q * t) / 2) * np.sin(phi * t) +
                                       A2 * np.exp(-(self.q * t) / 2) * phi * np.cos(phi * t) -
                                       (self.Fd / phi) * np.cos(phi * t) *
                                       ((phi - self.omega_d) * (
                                                   2 * (phi - self.omega_d) * np.cos((phi - self.omega_d) * t) -
                                                   self.q * np.sin((phi - self.omega_d) * t)) / K1 -
                                        (phi + self.omega_d) * (
                                                    2 * (phi + self.omega_d) * np.cos((phi + self.omega_d) * t) -
                                                    self.q * np.sin((phi + self.omega_d) * t)) / K2) +
                                       self.Fd * np.sin(phi * t) *
                                       ((self.q * np.cos((phi - self.omega_d) * t) + 2 * (
                                                   phi - self.omega_d) * np.sin((phi - self.omega_d) * t)) / K1
                                        - (self.q * np.cos((self.omega_d + phi) * t) + 2 * (
                                                           self.omega_d + phi) * np.sin(
                                                   (self.omega_d + phi) * t)) / K2) +
                                       self.Fd * np.cos(phi * t) *
                                       ((2 * (phi - self.omega_d) * np.cos(
                                           (phi - self.omega_d) * t) - self.q * np.sin(
                                           (phi - self.omega_d) * t)) / K1
                                        + (self.q * np.sin((self.omega_d + phi) * t) - 2 * (
                                                           self.omega_d + phi) * np.cos(
                                                   (self.omega_d + phi) * t)) / K2) +
                                       (self.Fd / phi) * np.sin(phi * t) *
                                       ((phi - self.omega_d) * (-self.q * np.cos((phi - self.omega_d) * t) -
                                                                2 * (phi - self.omega_d) * np.sin(
                                                   (phi - self.omega_d) * t)) / K1 +
                                        (phi + self.omega_d) * (self.q * np.cos((phi + self.omega_d) * t) +
                                                                2 * (phi + self.omega_d) * np.sin(
                                                           (phi + self.omega_d) * t)) / K2))


        self.t_arr = np.arange(0, T_MAX, self.dt)
        self.y_arr = formula_y(self.t_arr)
        self.y_vel_arr = formula_y_der(self.t_arr)
        self.energy_arr = self._calc_energy(self.y_arr, self.y_vel_arr)

        return self.t_arr, self.y_arr, self.y_vel_arr, self.energy_arr

    # Dummy function for debugging and development
    def _dummy_func(self):
        self.t_arr = np.arange(0, T_MAX, self.dt)
        self.y_arr = 3 * np.pi / 5 * np.cos(self.t_arr) * np.exp(-self.t_arr / 20)
        self.y_vel_arr = -3 * np.pi / 100 * np.exp(-self.t_arr / 20) * (20 * np.sin(self.t_arr) + np.cos(self.t_arr))
        self.energy_arr = self._calc_energy(self.y_arr, self.y_vel_arr)
        return self.t_arr, self.y_arr, self.y_vel_arr, self.energy_arr

    def calc_differences(self, calced_angles, calced_energies):
        if not (calced_angles.size == calced_energies.size):
            print("Error: The two arrays are not of the same size. Returning nothing.")
            return
        return calced_angles - self.y_arr, calced_energies - self.energy_arr
