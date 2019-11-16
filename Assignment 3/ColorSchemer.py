'''
    Handling the colorschemes of the interactive figure
    For Assignment 3, TFY4345
    Created by: Haakon Svane
    Date: 12. September, 2019
'''

import sys

# Scheme class for storing the color-values of each scheme. Tuple-lengths must equal number of members in 'Color' class
class __schemes:
    '''
    Order for tuple-indexes are the same as for the Color class. That is:
        0: Text color
        1: Figure background color
        2: Plot background color
        3: Plot-line color
        4: Rod color
        5: Bob color
        6: Slider, selected portion
        7: Slider, non-selected portion
        8: Buttons, non-hover
        9: Buttons, on hover
        10: Radio button, background
        11: Radio button, selected button
    '''

    default = (
        "black",
        "white",
        "white",
        "black",
        "darkgray",
        "black",
        "darkgray",
        "lightgray",
        "darkgray",
        "lightgray",
        "lightgray",
        "darkgray"
    )

    synthwave_sunrise = (
        "#ffd319",
        "#8c1eff",
        "#f222ff",
        "#ffd319",
        "#ff901f",
        "#ffd319",
        "#ff2975",
        "#ff901f",
        "#ff2975",
        "#f222ff",
        "#ff2975",
        "#ff901f"
    )

    municipal_melancholy = (
        "black",
        "peachpuff",
        "white",
        "black",
        "darkgray",
        "black",
        "darkorange",
        "bisque",
        "bisque",
        "antiquewhite",
        "bisque",
        "darkorange"
    )

    lucrative_luxury = (
        "#ffe5a9",
        "#525266",
        "#cbbeb5",
        "#ff6666",
        "#525266",
        "#423f3b",
        "#ffe5a9",
        "#423f3b",
        "#423f3b",
        "#cbbeb5",
        "#423f3b",
        "#ffe5a9",
    )

# Class container for storing the active colors
class Color:
    TEXT = 0
    BACKGROUND = 0
    PLOT_BACKGROUND = 0
    PLOTLINE = 0
    ROD = 0
    BOB = 0
    SLIDER_ACTIVE = 0
    SLIDER_INACTIVE = 0
    BUTTON_OFF_HOVER = 0
    BUTTON_ON_HOVER = 0
    RADIO_BACKGROUND = 0
    RADIO_SELECTED = 0

# Function for setting the color scheme. Takes a string (scheme-name) as input and updates the members of 'Color'.
def set_color_scheme(scheme):
    try:
        getattr(__schemes, scheme)
    except:
        print("ERROR: Could not find color scheme '%s', reverting to default." %scheme)
        scheme = "default"
    j = 0
    for i in vars(Color):
        if not i.startswith("__"):
            try:
                setattr(Color, i, getattr(__schemes, scheme)[j])
            except: sys.exit("SYSTEM EXIT:\nMismatch between color scheme '%s' and 'Color' class members.\nCheck length of '%s'-tuple and 'Color' members."%(scheme, scheme))
            j+=1

