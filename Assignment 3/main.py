'''
    Main module of the project. Creates a menu instance and runs it.
    For Assignment 3, TFY4345
    Created by: Haakon Svane
    Date: 13. November, 2019

    Please refer to the README.txt should you have any intentions of changing the parameters or color_schemes
'''


from Menu import *

def main():
    m = Menu(color_sch='default')
    m.start_ani()

    # Code for the poincaré section. Uncomment and set the parameter T_MAX to a high value (100000) for best result
    # s = Solver("runge_kutta")
    # t_arr, y_arr, y_vel_arr = s.solve(init_theta, init_rod_length, init_bob_mass, init_dt, 0.5, 1.2, 2/3)
    # ind = s.strob_times
    # plt.plot(y_arr[ind], y_vel_arr[ind], marker=".", linestyle="None", markersize=1, color="black")
    # plt.title("Poincaré section of the chaotic system over $t = 100000$ s")
    # plt.ylabel("$\omega (t)$")
    # plt.xlabel("$\\theta (t)$")
    # plt.show()
main()