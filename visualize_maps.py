import matplotlib.pyplot as plt
import shapefile
from shapely.geometry import shape
import numpy as np
import imageio
import config
import math

# Plot mappings
def make_pics(state,shape_data,num_lab,num_iter,mode):
    
    # Get [num_districts] unique colors
    cm = plt.get_cmap('Paired')
    colors = [cm(1.*i/num_lab) for i in range(num_lab)]
    
    # Sorts colors by geographic location
    geosorted_labels  = np.fromfile("geosorted_labels",dtype=np.uint32)
    filenames = []
    
    shape_dict = {}     # GEOID -> shape data

    for tract in shape_data.shapeRecords():
        aff_geo_id = int(tract.record[3][config.id_0:]) # GEOID
        shape_dict[aff_geo_id] = tract
    
    if mode == config.MODE_BEGIN_END: # First and last only
        picture_range = [0, num_iter-1]

    if mode == config.MODE_LOG: # Logarithm
        st = 10
        base = 1.5;
    
        picture_range = range(0,int(math.log(num_iter,base)))
        for i in range(0,len(picture_range)):
            picture_range[i] = int(pow(base,picture_range[i])) - 1

    if mode == config.MODE_ALL: # All
        picture_range = range(0,num_iter)
            
    
    for frame in picture_range:
        name = config.temp_folder + state + config.unit_data_suffix + '_step_' + str(frame)
 
        # Read label file & reshape into a matrix
        labels  = np.fromfile(name,dtype=np.uint32)
        labels = np.reshape(labels,(-1,num_lab+1))
        
        # Configure plot
        fig = plt.figure(frameon=False,facecolor="white")
        ax = fig.add_axes([0, 0, 1, 1])
        ax.axis('off')
        
        # For each tract, get GEOID and check dictionary for shape
        for row in labels:
            
            id = row[0]
            if id in shape_dict:
                
                # Get tract information & labeling
                tract = shape_dict[id]
                tractLabels = row[1:]
                
                # Get index of max population to assign color
                idx = np.argmax(tractLabels, axis=None)
                
                # Resort according to latitude
                idx, = np.where(geosorted_labels == idx)
                idx = idx[0]
                
                # Plot
                x = [i[0] for i in tract.shape.points[:]]
                y = [i[1] for i in tract.shape.points[:]]
                plt.fill(x,y,facecolor=colors[idx],edgecolor=colors[idx])
    
        # Configure axis
        plt.rcParams['axes.facecolor']='white'
        plt.rcParams['savefig.facecolor']='white'
        
        fig.patch.set_visible(False)
        
        ax.set_aspect('equal', 'box')
        ax.axis('off')
        
        # Write to file
        name = "output/flow_" + state + "_" +  str(frame).rjust(3, '0') + ".png"
        with open(name, 'w') as outfile:
            fig.canvas.print_png(outfile)
        filenames.append(name)

    # Save GIF
    images = []

    for filename in filenames:
        images.append(imageio.imread(filename))

    # Reverse for looping
    filenames.reverse()
    for filename in filenames:
        images.append(imageio.imread(filename))

    imageio.mimsave(state + '_output.gif', images,'GIF',duration=5.0/len(filenames))



