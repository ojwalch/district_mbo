import sys
import run_state
import config
import os.path

def print_syntax():
    print('Syntax:')
    print('python demo_gif.py [state]')
    print('python demo_gif.py [state] [filename.shp (in data/)]')
    print('python demo_gif.py [state] [filename.shp (in data/)] [# nearest neighbors]')
    print('python demo_gif.py [state] [filename.shp (in data/)] [# nearest neighbors] [centroid coefficient]')
    print('python demo_gif.py [state] [filename.shp (in data/)] [# nearest neighbors] [centroid coefficient] [temperature]\n')
    print('python demo_gif.py [state] [filename.shp (in data/)] [# nearest neighbors] [centroid coefficient] [temperature] [annealing rate]\n')
    
    print('Only state is required. When other parameters are not specified, the following are used:')
    print('- Initialization is random')
    print('- k (# of nearest neighbors) is ' + str(k))
    print('- alpha (centroid coefficient) is ' + str(ms_param) + '\n')
    print('- temperature is ' + str(temp) + '  and annealing rate is ' + str(annealing))



# Set default parameters
k =  150
max_iter = 700
init_state = config.INIT_RAND
alt_map = ''
ms_param = 2
stopping_criteria = 0.00
lb_frac = 0.985
temp =  0.05
annealing = 0.99
verbose = 0
driving_distance = 0
mode = config.MODE_LOG

if(len(sys.argv) < 2):
    print_syntax()
else:
    state = str(sys.argv[1]).upper()
    
    if(state in config.state_dict):
        
        if(len(sys.argv) > 2):
            file_path = sys.argv[2]
            if(os.path.exists(file_path)):
                if(file_path.endswith('.shp')):
                    init_state = config.INIT_ALT
                    alt_map = file_path
                else:
                    print('Provided file ' + file_path + ' does not have extension .shp. Please provide a shapefile.')

            else:
                print('Could not find file at path: ' + file_path)
                print('Please confirm path.')

        if(len(sys.argv) > 3):
            try:
                temp_k = int(sys.argv[3])
                if temp_k <= 0:
                    print('Positive integer value required for k')
                    print_syntax()
                else:
                    k = temp_k
            except ValueError:
                print('Positive integer value required for k')
                print_syntax()

        if(len(sys.argv) > 4):
            try:
                temp_ms = float(sys.argv[4])
                if temp_ms < 0:
                    print('Nonnegative value required for alpha')
                    print_syntax()
                else:
                    ms_param = temp_ms
            except ValueError:
                print('Nonnegative value required for alpha')
                print_syntax()

        if(len(sys.argv) > 5):
            try:
                temp_temp = float(sys.argv[5])
                if temp_temp < 0:
                    print('Nonnegative value required for temperature')
                    print_syntax()
                else:
                    temp = temp_temp
            except ValueError:
                print('Nonnegative value required for temperature')
                print_syntax()


        if(len(sys.argv) > 6):
            try:
                temp_annealing = float(sys.argv[6])
                if temp_annealing < 0:
                    print('Nonnegative value required for annealing')
                    print_syntax()
                else:
                    annealing = temp_annealing
            except ValueError:
                print('Nonnegative value required for annealing')
                print_syntax()

        run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,alt_map)
        print('GIF saved to gifs/' + state + '_output.gif')
    else:
        print('Invalid state! Please enter a two letter state abbreviation.\n\n')
        print_syntax()
    
    


