import sys
import run_state
import config
import os.path
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont
import math


def print_syntax():
    print('Syntax:')
    print('python demo_nearby_maps.py [state]')

    print('python demo_nearby_maps.py [state] [filename.shp (in data/)]')
    print('python demo_nearby_maps.py [state] [filename.shp (in data/)] [# nearest neighbors]')
    print('python demo_nearby_maps.py [state] [filename.shp (in data/)] [# nearest neighbors] [centroid coefficient]')
    print('python demo_nearby_maps.py [state] [filename.shp (in data/)] [# nearest neighbors] [centroid coefficient] [temperature]')
    print('python demo_nearby_maps.py [state] [filename.shp (in data/)] [# nearest neighbors] [centroid coefficient] [temperature] [annealing]')
    print('python demo_nearby_maps.py [state] [filename.shp (in data/)] [# nearest neighbors] [centroid coefficient] [temperature] [annealing] [max iter]')

    print('Only state file is required. When other parameters are not specified, the following are used:')
    print('- Initialized from current mapping')
    print('- k (# of nearest neighbors) is ' + str(k))
    print('- alpha (centroid coefficient) is ' + str(ms_param) + ' (high to prevent splitting in the short number of iterations)')
    print('- temperature is ' + str(temp) + '  and annealing rate is ' + str(annealing))
    print('- max number for iterations is ' + str(max_iter))

# Set default parameters
k =  150
init_state = config.INIT_CURRENT
ms_param = 2.5
stopping_criteria = 0.00
lb_frac = 0.985
temp =  0.01
annealing = 0.5
verbose = 0
driving_distance = 0
mode = config.MODE_BEGIN_END
max_iter = 3


if(len(sys.argv) < 2):
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
                temp_temp = float(sys.argv[5])
                if temp_temp < 0:
                    print('Nonnegative value required for temperature')
                    print_syntax()
                    sys.exit()
                
                else:
                    temp = temp_temp
            except ValueError:

                print('Nonnegative value required for temperature')
                print_syntax()
                sys.exit()
            
            
            
        if(len(sys.argv) > 6):
            try:
                temp_annealing = float(sys.argv[6])
                if temp_annealing < 0:
                    print('Nonnegative value required for annealing')
                    print_syntax()
                    sys.exit()
                
                else:
                    annealing = temp_annealing
            except ValueError:
                print('Nonnegative value required for annealing')
                print_syntax()
                sys.exit()

        if(len(sys.argv) > 7):
            try:
                temp_iter = int(sys.argv[7])
                if temp_iter < 0:
                    print('Positive integer value required for number of iterations')
                    print_syntax()
                    sys.exit()
                
                else:
                    max_iter = temp_iter
            except ValueError:
                print('Positive integer value required for number of iterations')
                print_syntax()
                sys.exit()



        num_runs = 16

        for i in range(0,num_runs):
            if i == 0:
                p = run_state.read(state,'')

            num_iter = run_state.run_with_data(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,'',p)
            
            img = Image.open('output/flow_' + state + '_' +  str(num_iter - 1).rjust(3, '0') + ".png")
            
            
            if i == 0:
                img = Image.open('output/flow_' + state + '_' +  str(0).rjust(3, '0') + ".png")
                
                w, h = img.size
                
                x_pad = 20
                y_pad = 10
                
                new_im = Image.new('RGBA', (int(4*w) + (5)*x_pad, int(4*h) + 5*y_pad))
                draw = ImageDraw.Draw(new_im)
                    
            x0 = int(x_pad) + int((i % 4)*(w + x_pad))
            y0 = int(y_pad) + int(math.floor(i/4)*(h + y_pad))
        
            new_im.paste(img,(x0,y0))
            plt.close("all")

        new_im.save(state + '_nearby.png')
        print('Nearby maps image saved to ' + state + '_nearby.png')

    else:
        print('Invalid state! Please enter a two letter state abbreviation.\n\n')
        print_syntax()
        sys.exit()
