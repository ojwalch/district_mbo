import run_state
import sys
from PIL import Image, ImageDraw, ImageFont
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import config
from subprocess import call


def make_driving_dist():
    
    # Run MD with and without driving distance as weight
    state = 'MD'
    k =  150
    max_iter = 300
    init_state = config.INIT_RAND
    ms_param = 1
    stopping_criteria = 0.00
    lb_frac = 0.985
    temp =  .1
    annealing = 0.985
    verbose = 1
    driving_distance = 0
    mode = config.MODE_BEGIN_END

    # Run without driving distance
    num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,"")

    img = Image.open('output/flow_' + state + '_' +  str(num_iter - 1).rjust(3, '0') + ".png")
    w, h = img.size

    x_pad = 100
    y_pad = 30

    new_im = Image.new('RGBA', (int(2*w) + 2*x_pad, int(h) + y_pad))
    new_im.paste(img,(int(x_pad*0.5),y_pad))

    # Run with driving distance
    driving_distance = 1
    num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,"")
    
    img = Image.open('output/flow_' + state + '_' + str(num_iter - 1).rjust(3, '0') + ".png")
    w, h = img.size

    new_im.paste(img,(w + int(x_pad*1.5),y_pad))

    new_im.save('figure_driving_dist.png')



def make_ms():
    
    # Run VA at 5 different choices of MS Parameter
    state = 'VA'
    k =  150
    max_iter = 300
    init_state = config.INIT_RAND
    stopping_criteria = 0.00
    lb_frac = 0.985
    temp =  0.1
    annealing = 0.985
    verbose = 0
    driving_distance = 0
    mode = config.MODE_BEGIN_END
    
    ms_range = [0, 0.01, 0.1, 0.5, 1, 10]
    
    for i in range(0,len(ms_range)):
        
        ms_param = ms_range[i]
        num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,"")
    
        img = Image.open('output/flow_' + state + '_' +  str(num_iter - 1).rjust(3, '0') + ".png")

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

        title_string = u"\u03B1" + " = " + str(ms_param) #
        font=ImageFont.truetype("/Library/Fonts/Arial Unicode.ttf",40)

        draw.text((int(x0 + (w + x_pad)/2), y0),title_string,font=font,fill='#000000')


    new_im.save('figure_ms.png')



def make_energy_traj():
    
    # Show how energy changes with iteration number, initializing from current configuration
    
    state = 'MI'
    k =  150
    max_iter = 100
    init_state = config.INIT_RAND
    stopping_criteria = 0.00
    lb_frac = 0.985
    temp =  0.1
    annealing = 0.95
    verbose = 1
    driving_distance = 0
    mode = config.MODE_NONE
    ms_param = 1
    
    num_runs = 20
    
    roi = 100
    figure_energy = plt.figure()
    axis_energy = figure_energy.add_subplot(111)

    colors = iter(cm.rainbow(np.linspace(0, 1, num_runs)))

    for i in range(0,num_runs):
        num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,"")
        energy = np.fromfile("energy",dtype=np.float64)

        axis_energy.plot(energy[0:roi], 'k--',color=next(colors))

    figure_energy.suptitle('Energy changes with iteration count', fontsize=20)
    plt.xlabel('Iteration number', fontsize=18)
    plt.ylabel('Energy', fontsize=16)

    plt.savefig('figure_energy.png')



def make_energy_hist():
    
    # Compare energies from PA districting (2016, and Remedial) and
    state = 'PA'
    k =  150
    max_iter = 100
    init_state = config.INIT_RAND
    stopping_criteria = 0.00
    lb_frac = 0.985
    temp =  0
    annealing = 0.985
    verbose = 1
    driving_distance = 0
    mode = config.MODE_NONE
    ms_param = 2
    
    final_energies = []
    start_energies_current = []
    start_energies_alt = []
    
    # Run current configuration
    init_state = config.INIT_CURRENT
    num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,"")
    energy = np.fromfile("energy",dtype=np.float64)
    start_energies_current.append(energy[0])

    # Run remedial
    init_state = config.INIT_ALT
    num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,"Remedial Plan Shapefile.shp")
    energy = np.fromfile("energy",dtype=np.float64)
    start_energies_alt.append(energy[0])

    # Run frmo random num_runs times
    num_runs = 30
    init_state = config.INIT_RAND
    for i in range(0,num_runs):
        num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,"")
        energy = np.fromfile("energy",dtype=np.float64)
        final_energy = energy[num_iter - 1]
        
        for j in range(0,num_iter - 1):
            if(energy[j] > 0 and energy[j + 1] == 0):
                final_energy = energy[j]

        final_energies.append(final_energy)
        
    figure_energy = plt.figure()
    axis_energy = figure_energy.add_subplot(111)


    bins = np.linspace(0.02, 0.05, 100)
    axis_energy.hist(final_energies, bins, facecolor='blue', alpha=0.5,label='Results from MBO')
    axis_energy.hist(start_energies_current, bins, facecolor='red', alpha=0.5,label='2016 Districting')
    axis_energy.hist(start_energies_alt, bins, facecolor='green', alpha=0.5,label='2018 Remedial Plan')
    axis_energy.legend(loc='upper right')


    plt.xlabel('Energy', fontsize=18)
    plt.ylabel('Count', fontsize=18)
    axis_energy.spines['right'].set_visible(False)
    axis_energy.spines['top'].set_visible(False)
    axis_energy.xaxis.set_ticks_position('bottom')
    axis_energy.yaxis.set_ticks_position('left')

    plt.savefig('figure_energy_hist.png')


def make_pa_gif():
    # Make a GIF showing flow in PA from Remedial plan.
    
    state = 'PA'
    k =  150
    max_iter = 300
    init_state = config.INIT_ALT
    ms_param = 1
    stopping_criteria = 0.00
    lb_frac = 0.985
    temp =  0
    annealing = 0.985
    verbose = 0
    driving_distance = 0
    mode = config.MODE_LOG
    alt_map = "Remedial Plan Shapefile.shp"
    run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,alt_map)
    call("mv PA_output.gif figure_PA_output.gif", shell=True)


def make_pa_before_and_after():
    
    # Show 2016 and Remedial plans before and after MBO.

    state = 'PA'
    k =  150
    max_iter = 300
    stopping_criteria = 0.00
    lb_frac = 0.985
    temp =  0
    annealing = 0.985
    verbose = 1
    driving_distance = 0
    mode = config.MODE_BEGIN_END
    ms_param = 1
    initial_state_range = [config.INIT_CURRENT, config.INIT_ALT]

    for i in range(0,len(initial_state_range)):
        
        init_state = initial_state_range[i]
        
        alt_map = ""
        if(init_state == config.INIT_ALT):
            alt_map = "Remedial Plan Shapefile.shp"
        num_iter = run_state.run(state,k,max_iter,init_state,ms_param,stopping_criteria,lb_frac,temp,annealing,verbose,driving_distance,mode,alt_map)
        
        
        begin_img = Image.open('output/flow_' + state + '_' +  str(0).rjust(3, '0') + ".png")
        end_img = Image.open('output/flow_' + state + '_' +  str(num_iter - 1).rjust(3, '0') + ".png")
        
        w, h = begin_img.size
        x_pad = 50
        y_pad = 30

        if i == 0:
            new_im = Image.new('RGBA', (int(2*w) + 3*x_pad, int(2*h) + 3*y_pad))

        new_im.paste(begin_img,(x_pad,y_pad + i*(h + y_pad)))
        new_im.paste(end_img,(2*x_pad + w,y_pad + i*(h + y_pad)))
        
            

    new_im.save('figure_pa.png')

make_pa_gif()
make_pa_before_and_after()
make_energy_traj()     # energy trajectories
make_driving_dist()    # driving distance
make_ms()              # mumford-shah
make_energy_hist()      # histograms
