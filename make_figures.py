import run_state
import sys
from PIL import Image, ImageDraw, ImageFont
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import config
from subprocess import call
import math


# Run MD with and without driving distance as weight
def make_driving_dist():
    
    state = 'MD'
    k =  150
    max_iter = 300
    init_state = config.INIT_RAND
    ms_param = 1
    stopping_criteria = 0.00
    lb_frac = 0.980
    temp =  .1
    annealing = 0.985
    verbose = 0
    driving_distance = 0
    mode = config.MODE_BEGIN_END

    # Run without driving distance
    num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,'')

    img = Image.open('output/flow_' + state + '_' +  str(num_iter - 1).rjust(3, '0') + '.png')
    w, h = img.size

    x_pad = 200
    y_pad = 30

    new_im = Image.new('RGBA', (int(2*w) + 2*x_pad, int(h) + y_pad))
    draw = ImageDraw.Draw(new_im)

    font = ImageFont.truetype('/Library/Fonts/Arial Unicode.ttf',75)
    new_im.paste(img,(int(x_pad*0.5),int(y_pad)))
    draw.text((0,int(2*y_pad)),'A',font=font,fill='#000000')

    # Run with driving distance
    driving_distance = 1
    num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,'')
    
    img = Image.open('output/flow_' + state + '_' + str(num_iter - 1).rjust(3, '0') + '.png')
    w, h = img.size

    new_im.paste(img,(w + int(x_pad*1.5),y_pad))
    draw.text((w + int(x_pad),int(2*y_pad)),'B',font=font,fill='#000000')

    new_im.save('figures/figure_driving_dist.png')


# Run VA at 5 different choices of nearest neighbors
def make_knn():
    
    state = 'VA'
    ms_param =  1
    max_iter = 300
    init_state = config.INIT_CURRENT
    stopping_criteria = 0.00
    lb_frac = 0.980
    temp =  0
    annealing = 0.985
    verbose = 0
    driving_distance = 0
    mode = config.MODE_BEGIN_END
    
    k_range = [5, 10, 50, 100, 200, 500]
    
    for i in range(0,len(k_range)):
        
        k = k_range[i]
        num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,'')
        
        img = Image.open('output/flow_' + state + '_' +  str(num_iter - 1).rjust(3, '0') + '.png')
        
        if i == 0:
            w, h = img.size
            
            x_pad = 50
            y_pad = 30
            
            new_im = Image.new('RGBA', (int(3*w) + (4)*x_pad, int(2*h) + 3*y_pad))
            draw = ImageDraw.Draw(new_im)
        
        if(i < 3):
            x0 = int(x_pad*0.5) + i*(w + x_pad)
            y0 = y_pad
        else:
            x0 = int(x_pad*0.5) + (i - 3)*(w + x_pad)
            y0 = 2*y_pad + int(h)
    
        new_im.paste(img,(x0,y0))
        
        title_string =  'k = ' + str(k)
        font=ImageFont.truetype('/Library/Fonts/Arial Unicode.ttf',40)
        
        draw.text((int(x0 + (w + x_pad)/2), y0),title_string,font=font,fill='#000000')

    new_im.save('figures/figure_knn.png')


# Run VA at 5 different choices of centroid distance parameter
def make_ms():
    
    state = 'VA'
    k =  150
    max_iter = 300
    init_state = config.INIT_CURRENT
    stopping_criteria = 0.00
    lb_frac = 0.980
    temp =  0
    annealing = 0.985
    verbose = 0
    driving_distance = 0
    mode = config.MODE_BEGIN_END
    
    ms_range = [0, 0.01, 0.1, 0.5, 1, 10]
    
    for i in range(0,len(ms_range)):
        
        ms_param = ms_range[i]
        num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,'')
    
        img = Image.open('output/flow_' + state + '_' +  str(num_iter - 1).rjust(3, '0') + '.png')

        if i == 0:
            w, h = img.size
        
            x_pad = 50
            y_pad = 30
        
            new_im = Image.new('RGBA', (int(3*w) + (4)*x_pad, int(2*h) + 3*y_pad))
            draw = ImageDraw.Draw(new_im)

        if(i < 3):
            x0 = int(x_pad*0.5) + i*(w + x_pad)
            y0 = y_pad
        else:
            x0 = int(x_pad*0.5) + (i - 3)*(w + x_pad)
            y0 = 2*y_pad + int(h)

        new_im.paste(img,(x0,y0))

        title_string = u'\u03B1' + ' = ' + str(ms_param) #
        font=ImageFont.truetype('/Library/Fonts/Arial Unicode.ttf',40)

        draw.text((int(x0 + (w + x_pad)/2), y0),title_string,font=font,fill='#000000')


    new_im.save('figures/figure_ms.png')


# Show how energy changes with iteration number for 4 states, initializing from random
def make_energy_traj():
    
    states = ['MI','PA','NC','VA']
    labels = ['A','B','C','D']

    for state in states:
        k =  150
        max_iter = 100
        init_state = config.INIT_RAND
        stopping_criteria = 0.00
        lb_frac = 0.980
        temp =  0.1
        annealing = 0.95
        verbose = 0
        driving_distance = 0
        mode = config.MODE_NONE
        ms_param = 1
        
        num_runs = 20
        
        roi = 100
        figure_energy = plt.figure()
        axis_energy = figure_energy.add_subplot(111)

        colors = iter(cm.rainbow(np.linspace(0, 1, num_runs)))

        for i in range(0,num_runs):
            num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,'')
            energy = np.fromfile('energy',dtype=np.float64)

            axis_energy.plot(energy[0:roi], 'k--',color=next(colors))

        axis_energy.spines['right'].set_visible(False)
        axis_energy.spines['top'].set_visible(False)
        axis_energy.xaxis.set_ticks_position('bottom')
        axis_energy.yaxis.set_ticks_position('left')

        plt.xlabel('Iteration #', fontsize=16)
        plt.ylabel('Energy', fontsize=16)

        plt.savefig('figures/' + state + '_figure_energy.png')

    for i in range(0,len(states)):

        img = Image.open('figures/' + states[i] + '_figure_energy.png')

        w, h = img.size
        x_pad = 50
        y_pad = 30
        
        if i == 0:
            new_im = Image.new('RGBA', (int(2*w) + 3*x_pad, int(2*h) + 3*y_pad))
            draw = ImageDraw.Draw(new_im)
            font=ImageFont.truetype('/Library/Fonts/Arial Unicode.ttf',85)

        x_ind = i % 2
        y_ind = math.floor(i/2)
        new_im.paste(img,(int(x_pad + x_ind*(w + x_pad)),int(y_pad + y_ind*(h + y_pad))))
        draw.text((int(x_pad/2 + x_ind*(w + x_pad)), int(y_pad + y_ind*(h + y_pad/2))),labels[i],font=font,fill='#000000')



    new_im.save('figures/figure_combined_energy.png')



# Compare energies from PA districting (2016, and Remedial) with outputs from the MBO algorithm
def make_energy_hist():
    
    state = 'PA'
    k =  150
    max_iter = 100
    init_state = config.INIT_RAND
    stopping_criteria = 0.00
    lb_frac = 0.980
    temp =  0
    annealing = 0.985
    verbose = 0
    driving_distance = 0
    mode = config.MODE_NONE
    ms_param = 2
    
    final_energies = []
    start_energies_current = []
    start_energies_alt = []
    
    # Run current configuration
    init_state = config.INIT_CURRENT
    num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,'')
    energy = np.fromfile('energy',dtype=np.float64)
    start_energies_current.append(energy[0])

    # Run remedial
    init_state = config.INIT_ALT
    num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,'Remedial Plan Shapefile.shp')
    energy = np.fromfile('energy',dtype=np.float64)
    start_energies_alt.append(energy[0])

    # Run from random num_runs times
    num_runs = 30
    init_state = config.INIT_RAND
    for i in range(0,num_runs):
        num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,'')
        energy = np.fromfile('energy',dtype=np.float64)
        final_energy = energy[num_iter - 1]
        
        for j in range(0,num_iter - 1):
            if(energy[j] > 0 and energy[j + 1] == 0):
                final_energy = energy[j]

        final_energies.append(final_energy)
        
    figure_energy = plt.figure()
    axis_energy = figure_energy.add_subplot(111)


    bins = np.linspace(0.02, 0.05, 100)
    axis_energy.hist(final_energies, bins, facecolor='blue',alpha=0.5,label='Results from MBO')
    axis_energy.hist(start_energies_current, bins, facecolor='red', alpha=0.5,label='2016 Districting')
    axis_energy.hist(start_energies_alt, bins, facecolor='green', alpha=0.5,label='2018 Remedial Plan')
    axis_energy.legend(loc='upper right')

    plt.xlabel('Energy', fontsize=18)
    plt.ylabel('Count', fontsize=18)
    axis_energy.spines['right'].set_visible(False)
    axis_energy.spines['top'].set_visible(False)
    axis_energy.xaxis.set_ticks_position('bottom')
    axis_energy.yaxis.set_ticks_position('left')

    plt.savefig('figures/figure_energy_hist.png')


# Make the PA figure, varying parameters
def make_hist_grid():
    state = 'PA'
    k_range = [50, 100, 150, 200, 250]
    ms_range = [0, 1, 2]

    for k_ind in range(0,len(k_range)):
        for ms_ind in range(0,len(ms_range)):
            
            k = k_range[k_ind]
            ms_param = ms_range[ms_ind]
            
            max_iter = 500
            init_state = config.INIT_RAND
            stopping_criteria = 0.00
            lb_frac = 0.980
            temp =  0.1
            annealing = 0.985
            verbose = 0
            driving_distance = 0
            mode = config.MODE_NONE
            
            final_energies = []
            start_energies_current = []
            start_energies_alt = []
            
            # Run current configuration
            init_state = config.INIT_CURRENT
            num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,'')
            energy = np.fromfile('energy',dtype=np.float64)
            start_energies_current.append(energy[0])
            
            # Run remedial
            init_state = config.INIT_ALT
            num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,'Remedial Plan Shapefile.shp')
            energy = np.fromfile('energy',dtype=np.float64)
            start_energies_alt.append(energy[0])
            
            # Run (initializing with random) num_runs times
            num_runs = 20
            init_state = config.INIT_RAND
            for i in range(0,num_runs):
                
                if i == 0:
                    p = run_state.read(state,'')
            
                num_iter = run_state.run_with_data(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,'',p)
                energy = np.fromfile('energy',dtype=np.float64)
                final_energy = energy[num_iter - 1]
                
                for j in range(0,num_iter - 1):
                    if(energy[j] > 0 and energy[j + 1] == 0):
                        final_energy = energy[j]
            
                final_energies.append(final_energy)

            figure_energy = plt.figure()
            axis_energy = figure_energy.add_subplot(111)
            
            bins = np.linspace(0.01, 0.05, 100)

            hist_temp,be = np.histogram(final_energies,bins=bins)
            
            weights = np.ones_like(final_energies)/float(max(hist_temp))

            axis_energy.hist(final_energies, bins, weights=weights, facecolor='blue',alpha=0.5,label='Results from MBO')
            axis_energy.hist(start_energies_current, bins, facecolor='red', alpha=0.5,label='2016 Districting')
            axis_energy.hist(start_energies_alt, bins, facecolor='green', alpha=0.5,label='2018 Remedial Plan')
            
            plt.xlabel('Energy', fontsize=18)
            plt.ylabel('Count', fontsize=18)
            axis_energy.spines['right'].set_visible(False)
            axis_energy.spines['top'].set_visible(False)
            axis_energy.xaxis.set_ticks_position('bottom')
            axis_energy.yaxis.set_ticks_position('left')
            axis_energy.set_xlim([0.01,0.05])
            axis_energy.set_ylim([0,1])
            save_name = 'output/energy_grid' + str(k) + '_' + str(ms_param) +'.png'
            plt.savefig(save_name)
            plt.close('all')


            img = Image.open(save_name)
    
            if k_ind == 0 and ms_ind == 0:
                w, h = img.size
            
                x_pad = 50
                y_pad = 30
            
                num_rows = len(ms_range)
                num_cols = len(k_range)
                new_im = Image.new('RGBA', (int(num_cols*w) + (num_cols + 1)*x_pad, int(num_rows*h) + (num_rows + 1)*y_pad))
                draw = ImageDraw.Draw(new_im)
    
            x0 = int(x_pad) + int(k_ind*(w + x_pad))
            y0 = int(y_pad) + int(ms_ind*(h + y_pad))

            new_im.paste(img,(x0,y0))

    new_im.save('figures/figure_energy_grid.png')


# Make a GIF showing flow in PA from Remedial plan.
def make_pa_gif():
    state = 'PA'
    k =  150
    max_iter = 300
    init_state = config.INIT_ALT
    ms_param = 1
    stopping_criteria = 0.00
    lb_frac = 0.980
    temp =  0
    annealing = 0.985
    verbose = 0
    driving_distance = 0
    mode = config.MODE_LOG
    alt_map = 'Remedial Plan Shapefile.shp'
    run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,alt_map)
    call('mv gifs/PA_output.gif figures/figure_PA_output.gif', shell=True)


# Show 2016 and Remedial plans before and after MBO.
def make_pa_before_and_after():
    state = 'PA'
    k =  150
    max_iter = 300
    stopping_criteria = 0.00
    lb_frac = 0.980
    temp =  0
    annealing = 0.985
    verbose = 0
    driving_distance = 0
    mode = config.MODE_BEGIN_END
    ms_param = 1
    initial_state_range = [config.INIT_CURRENT, config.INIT_ALT]

    labels = ['A','B','C','D']
    
    for i in range(0,len(initial_state_range)):
        
        init_state = initial_state_range[i]
        
        alt_map = ''
        if(init_state == config.INIT_ALT):
            alt_map = 'Remedial Plan Shapefile.shp'
        num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,alt_map)
        
        
        begin_img = Image.open('output/flow_' + state + '_' +  str(0).rjust(3, '0') + '.png')
        end_img = Image.open('output/flow_' + state + '_' +  str(num_iter - 1).rjust(3, '0') + '.png')
        
        w, h = begin_img.size
        x_pad = 50
        y_pad = 30

        if i == 0:
            new_im = Image.new('RGBA', (int(2*w) + 3*x_pad, int(2*h) + 3*y_pad))
            draw = ImageDraw.Draw(new_im)
            font = ImageFont.truetype('/Library/Fonts/Arial Unicode.ttf',85)

        new_im.paste(begin_img,(x_pad,y_pad + i*(h + y_pad)))
        new_im.paste(end_img,(2*x_pad + w,y_pad + i*(h + y_pad)))
        draw.text((0, int(y_pad + i*(h + y_pad))),labels[i*2],font=font,fill='#000000')
        draw.text((int(w + x_pad), int(y_pad + i*(h + y_pad))),labels[i*2 + 1],font=font,fill='#000000')

    new_im.save('figures/figure_pa.png')

# Generate maps by perturbing slightly and running for a small number of iterations.
def make_nearby_maps():
    state = 'OH'
    k =  5
    max_iter = 5
    init_state = config.INIT_CURRENT
    ms_param = 5
    stopping_criteria = 0.00
    lb_frac = 0.980
    temp =  0.05
    annealing = 0.5
    verbose = 0
    driving_distance = 0
    mode = config.MODE_BEGIN_END

    num_runs = 16

    for i in range(0,num_runs):
        if i == 0:
            p = run_state.read(state,'')
       
        num_iter = run_state.run_with_data(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,'',p)

        img = Image.open('output/flow_' + state + '_' +  str(num_iter - 1).rjust(3, '0') + '.png')


        if i == 0:
            img = Image.open('output/flow_' + state + '_' +  str(0).rjust(3, '0') + '.png')

            w, h = img.size
            
            x_pad = 20
            y_pad = 10
            
            new_im = Image.new('RGBA', (int(4*w) + (5)*x_pad, int(4*h) + 5*y_pad))
            draw = ImageDraw.Draw(new_im)

        x0 = int(x_pad) + int((i % 4)*(w + x_pad))
        y0 = int(y_pad) + int(math.floor(i/4)*(h + y_pad))
        

        new_im.paste(img,(x0,y0))

        plt.close('all')

        
    new_im.save('figures/figure_nearby_maps.png')


# Show flow steps
def make_flow_figure():
    state = 'NC'
    k =  150
    max_iter = 8
    init_state = config.INIT_CURRENT
    ms_param = 1
    stopping_criteria = 0.00
    lb_frac = 0.980
    temp =  0
    annealing = 0.5
    verbose = 0
    driving_distance = 0
    mode = config.MODE_ALL
    txt_buff = 60

    p = run_state.read(state,'')
    num_iter = run_state.run_with_data(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,'',p)
    
    for i in range(0,min(num_iter,max_iter)):
        
        if i == 0:
            img = Image.open('output/flow_' + state + '_' +  str(0).rjust(3, '0') + '.png')
            
            w, h = img.size
            
            x_pad = 20
            y_pad = 90
            
            new_im = Image.new('RGBA', (int(4*w) + (5)*x_pad, int(2*h) + 3*y_pad))
            draw = ImageDraw.Draw(new_im)
            font = ImageFont.truetype('/Library/Fonts/Arial Unicode.ttf',70)

        
        else:
            img = Image.open('output/flow_' + state + '_' +  str(i).rjust(3, '0') + '.png')
        
     
        x0 = int(x_pad) + int((i % 4)*(w + x_pad))
        y0 = int(y_pad) + int(math.floor(i/4)*(h + y_pad))
        
        new_im.paste(img,(x0,y0))
        draw.text((x0 + w/2 - txt_buff, y0 - y_pad),'i = ' + str(i),font=font,fill='#000000')

        plt.close('all')

    new_im.save('figures/figure_flow.png')


## This is where the functions are called. Uncomment to run. ##

#make_pa_gif()                  # Make flow gif
#make_pa_before_and_after()     # Make GIF showing differences between start and flow for 2016 and Remedial maps
#make_energy_traj()             # Show energy trajectories as iteration changes
#make_driving_dist()            # Show how driving time can be used as an alternative distance measure
#make_ms()                      # Show effects of alpha (centroid distance)
#make_knn()                     # Show effects of changing k (number of nearest neighbors)
#make_energy_hist()             # Make energy histogram
#make_nearby_maps()             # Generate grid of nearby mappings
#make_hist_grid()               # Make grid showing how alpha and k do not change relative preference for 2016 and Remedial mappings
make_flow_figure()              # Make figure showing flow
