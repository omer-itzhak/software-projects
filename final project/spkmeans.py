import sys
from enum import Enum
# how to connect when on a different file?
class goal(Enum):
    wam = 1
    ddg = 2
    lnorm = 3
    jacobi = 4

def args_parsing():
    argv = sys.argv
    if len(argv) != 4:
        invalid_input()
    k = argv[1]
    if not k.isnumeric():
        invalid_input()
    k = int(k)
    if k == 0:
        k = eigengap_heuristic()
    g = argv[2]
    try:
        #i want the enum to be defined on a different file
        g = goal[g]
    except:
        invalid_input()
    input_file_name = argv[3]
    
    

    


def eigengap_heuristic():
    return


def invalid_input():
    print("Invalid Input!")
    exit(1)
