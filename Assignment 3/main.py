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

main()