'''
    Main module of the project. Creates a menu instance and runs it.
    For Assignment 1, TFY4345
    Created by: Haakon Svane
    Date: 16. September, 2019

    Please refer to the README.txt should you have any intentions of changing parameters or color_schemes
'''


from Menu import *

def main():
    m = Menu(color_sch='default')
    m.start_ani()

main()