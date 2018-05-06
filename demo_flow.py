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
    print('python demo_flow.py [state]')
    
    print('python demo_flow.py [state] [filename.shp (in data/)]')
    print('python demo_flow.py [state] [filename.shp (in data/)] [# nearest neighbors]')
    print('python demo_flow.py [state] [filename.shp (in data/)] [# nearest neighbors] [centroid coefficient]')
    
    print('Only state file is required. When other parameters are not specified, the following are used:')
    print('- Initialized from current mapping')
    print('- k (# of nearest neighbors) is ' + str(k))
    print('- alpha (centroid coefficient) is ' + str(ms_param) + ' (high to prevent splitting in the short number of iterations)')

# Set default parameters
k =  150
max_iter = 8
init_state = config.INIT_CURRENT
ms_param = 1
stopping_criteria = 0.00
lb_frac = 0.985
temp =  0
annealing = 0.5
verbose = 0
driving_distance = 0
mode = config.MODE_ALL
txt_buff = 60
alt_map = ''

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

        p = run_state.read(state,alt_map)
        num_iter = run_state.run_with_data(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,alt_map,p)

        for i in range(0,min(num_iter,max_iter)):
            
            if i == 0:
                img = Image.open('output/flow_' + state + '_' +  str(0).rjust(3, '0') + ".png")

                w, h = img.size
                
                x_pad = 20
                y_pad = 90
                
                new_im = Image.new('RGBA', (int(4*w) + (5)*x_pad, int(2*h) + 3*y_pad))
                draw = ImageDraw.Draw(new_im)
                font = ImageFont.truetype("/Library/Fonts/Arial Unicode.ttf",70)
            else:
                img = Image.open('output/flow_' + state + '_' +  str(i).rjust(3, '0') + ".png")
            
            x0 = int(x_pad) + int((i % 4)*(w + x_pad))
            y0 = int(y_pad) + int(math.floor(i/4)*(h + y_pad))

            new_im.paste(img,(x0,y0))
            draw.text((x0 + w/2 - txt_buff, y0 - y_pad),'i = ' + str(i),font=font,fill='#000000')
            
            plt.close("all")

        new_im.save(state + '_flow.png')
        print('Flow figure saved to ' + state + '_flow.png')

    else:
        print('Invalid state! Please enter a two letter state abbreviation.\n\n')
        print_syntax()
        sys.exit()

