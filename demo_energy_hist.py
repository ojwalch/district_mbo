import sys
import run_state
import config
import os.path
import numpy as np
import matplotlib.pyplot as plt

def print_syntax():
    print('Syntax:')
    print('python demo_energy_hist.py [state] [filename.shp (in data/)]')
    print('python demo_energy_hist.py [state] [filename.shp (in data/)] [# nearest neighbors]')
    print('python demo_energy_hist.py [state] [filename.shp (in data/)] [# nearest neighbors] [centroid coefficient]')
    print('python demo_energy_hist.py [state] [filename.shp (in data/)] [# nearest neighbors] [centroid coefficient] [num runs]')

    print('Only state and new districting file are required. When other parameters are not specified, the following are used:')
    print('- k (# of nearest neighbors) is ' + str(k))
    print('- alpha (centroid coefficient) is ' + str(ms_param) + '\n')
    print('- number of runs ' + str(num_runs) + '\n')


# Set default parameters
k =  150
max_iter = 500
init_state = config.INIT_RAND
alt_map = ''
ms_param = 2
stopping_criteria = 0.00
lb_frac = 0.985
temp =  0.05
annealing = 0.99
verbose = 0
driving_distance = 0
mode = config.MODE_BEGIN_END
num_runs = 30

if(len(sys.argv) < 3):
    print_syntax()
else:
    state = str(sys.argv[1]).upper()
    
    if(state in config.state_dict):
        
        if(len(sys.argv) > 2):
            file_path = sys.argv[2]
            if(os.path.exists('data/' + file_path)):
                if(file_path.endswith('.shp')):
                    init_state = config.INIT_ALT
                    alt_map = file_path
                else:
                    print('Provided file ' + file_path + ' does not have extension .shp. Please provide a shapefile.')
        
            else:
                print('Could not find file ' + file_path + ' in data/ folder.')
                print('Please copy to that folder.')
                sys.exit()


            if(len(sys.argv) > 3):
                try:
                    temp_k = int(sys.argv[3])
                    if temp_k <= 0:
                        print('Positive integer value required for k')
                        print_syntax()
                        sys.exit()

                    else:
                        k = temp_k
                except ValueError:
                    print('Positive integer value required for k')
                    print_syntax()
                    sys.exit()


            if(len(sys.argv) > 4):
                try:
                    temp_ms = float(sys.argv[4])
                    if temp_ms < 0:
                        print('Nonnegative value required for alpha')
                        print_syntax()
                        sys.exit()

                    else:
                        ms_param = temp_ms
                except ValueError:
                    print('Nonnegative value required for alpha')
                    print_syntax()
                    sys.exit()


            if(len(sys.argv) > 5):
                try:
                    temp_num_runs = int(sys.argv[5])
                    if temp_num_runs <= 0:
                        print('Positive value required for number of runs')
                        print_syntax()
                        sys.exit()

                    else:
                        num_runs = temp_num_runs
                except ValueError:
                    print('Positive value required for number of runs')
                    print_syntax()
                    sys.exit()



            start_energies_alt = []
            final_energies = []
            
            # Run provided map.
            init_state = config.INIT_ALT
            num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,alt_map)
            energy = np.fromfile("energy",dtype=np.float64)
            start_energies_alt.append(energy[0])

            # Run from random num_runs times
            init_state = config.INIT_RAND

            for i in range(0,num_runs):
                num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,"")
                energy = np.fromfile("energy",dtype=np.float64)
                final_energy = energy[num_iter - 1]
                
                for j in range(0,num_iter - 1):
                    if(energy[j] > 0 and energy[j + 1] == 0):
                        final_energy = energy[j]
                
                final_energies.append(final_energy)
                plt.close("all")


            figure_energy = plt.figure()
            axis_energy = figure_energy.add_subplot(111)


            bins = np.linspace(0.02, 0.05, 100)
            axis_energy.hist(final_energies, bins, facecolor='blue',alpha=0.5,label='Results from flow')
            axis_energy.hist(start_energies_alt, bins, facecolor='yellow', alpha=0.5,label='Provided mapping')
            axis_energy.legend(loc='upper right')

            plt.xlabel('Energy', fontsize=18)
            plt.ylabel('Count', fontsize=18)
            axis_energy.spines['right'].set_visible(False)
            axis_energy.spines['top'].set_visible(False)
            axis_energy.xaxis.set_ticks_position('bottom')
            axis_energy.yaxis.set_ticks_position('left')

            plt.savefig(state + '_energy_hist.png')

            print('Energy histogram saved to ' + state + '_energy_hist.png')
    else:
        print('Invalid state! Please enter a two letter state abbreviation.\n\n')
        print_syntax()
        sys.exit()




