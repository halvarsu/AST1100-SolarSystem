from atmosphereAnalyser import *
import sys
import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-dp", "--dont_plot", 
                        action = 'store_true', default = False)
    parser.add_argument("-s", "--save_fig", 
                        action = 'store_true', default = False)
    parser.add_argument('--val', help="plot normalized value of least "\
                        "square function", action = 'store_true', 
                        default = False)
    parser.add_argument('-n', '--noise', help="plot noise", 
                        action = 'store_true', default = False)
    parser.add_argument('-o', '--one_gas', help="only calc one gas", 
                        default = '')
    parser.add_argument('-i', '--iteration', help="choose iteration",
                        type=int, default = 3)
    parser.add_argument('-v', '--verbose', help="print more info",
                        action = 'store_true', default =False)
    return parser.parse_args()

args = get_args()
find_spectral_lines(args)
