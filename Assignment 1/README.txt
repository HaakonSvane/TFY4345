This package contains our solution for Assignment 1 in TFY4345.
In order to avoid crashes and (unexpected) bugs, matplotlib should be updated to its newest version.
The code has currently been tested (to a varying degree) in Python 3.6.* and 3.7.* running on MacOS Mojave (10.14.*).

Bugs:
    Depending on the platform from which you are executing the code, strange stuff may occur. There are also some known
    bugs that I did not care to dwell deep enough into. Some of the known bugs are:
        *   When running a new animation from the menu (hitting "Apply and run.." or "Reset and run.."), some machines
            do not properly close the previous figure window. This probably has to do with threads (callback
            functions) not exiting in time.

        *   On certain interactive elements of the menu, hovering over them or changing its state causes the plots in the
            axis to disappear for a brief amount of time. This is as far as I know, a "feature" of matplotlib.

        *   When recalculating, an error message may appear in the console reading
            "invalid command name [NUMBER]_on_timer". This is due to callback timer not properly disconnecting on time.

Modules:
    The package contains several modules aimed to aid the readability and development of the source code.
    The purpose of the module are as follows:
        'Menu.py':              For handling the animations and button-inputs
        'AnimatorHelper.py':    Module with a helper-class to improve readability of 'Animator.py'. Contains
                                parameters for the plots.
        'ColorSchemer':         Contains classes for creating and getting color-schemes
        'main.py':              Main module. Gathers the different modules and sets up the workflow of the project.
        'Parameters':           Contains parameters for matplotlib, 'Animator.py', 'Solver.py' and general constants
                                used in the project.
        'Solver.py':            Contains the source code for the different solvers used in the project
                                (Euler, Euler-Cromer and Runge-Kutta (4th order)).

Color-schemes:
    At instantiation of an Animator object, the optional argument 'color_sch' can be set to swap between color-schemes.
    Supported color-schemes are as follows:

        'default':              Standard white/gray scheme, if no scheme is passed as parameter for the animator,
                                it defaults to this.

        'synthwave_sunrise':    Only for the brave. Readability = 0, Visual stimuli = 9.5/10.

        'municipal_melancholy': As boring as it can be in terms of color usage.

        'lucrative_luxury':'    A range of flat colors for those who are tired of 'default',
                                but not brave enough for 'synthwave_sunrise' . It just works!

Parameters:
    The 'Parameters.py' module contains physical constants and initial values to be used for the first
    run of the program. When soft-resetting the program (pressing 'Reset and run' instead of 'Apply and run'), the
    program also reverts back to these values.
    Perhaps the most interesting parameter in this module is the one named 'T_MAX' which allows you to adjust the total
    time interval that the calculations are done over. Increasing this value will give you a better understanding of
    the divergent/convergent properties of the different solvers at the cost of computational time and effort.

Animation speed and framerate:
    The framerate is a parameter found in the 'Parameters.py' module. This is decoupled from the animation speed
    to ensure a steady stream of processed data. Since the timestep 'dt' is changeable while the framerate is constant,
    different values for 'dt' will affect the smoothness of the animations appearing on screen. This is due to there
    being less points for the plots to respond to per animation frame. In other words, blame the math, not the
    programmer.

