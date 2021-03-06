usage = '''

This script can be used to run a simulation of the GLIF cell model using the Allen SDK 

To run a simulation with model ID 473875489 with 120pA input run:

    python run_glif.py 473875489 120
    
To run all the models used in this repo type:

    python run_glif.py -all

'''
import allensdk.core.json_utilities as json_utilities
from allensdk.model.glif.glif_neuron import GlifNeuron
import sys

from pyneuroml import pynml

def run_one_cell(glif_dir, curr_pA, dt=5e-7, show_plot=True):

    # initialize the neuron
    neuron_config = json_utilities.read('%s/neuron_config.json'%glif_dir)
    neuron = GlifNeuron.from_dict(neuron_config)

    # important! set the neuron's dt value for your stimulus in seconds
    neuron.dt = dt

    # make a short square pulse. stimulus units should be in Amps.
    stimulus = [ 0.0 ] * int(0.1/neuron.dt) + [ curr_pA * 1e-12 ] * int(1/neuron.dt) + [ 0.0 ] * int(0.1/neuron.dt)
    times = [ i*neuron.dt for i in range(len(stimulus)) ]

    # simulate the neuron
    output = neuron.run(stimulus)

    print "run.."

    voltage = output['voltage']
    threshold = output['threshold']
    spike_times = output['interpolated_spike_times']

    info = "Model %s; %spA stimulation; dt=%ss; %i spikes"%(glif_dir,curr_pA,dt,len(spike_times))
    print(info)
    
    v_file = open('%s/original.v.dat'%glif_dir,'w')
    th_file = open('%s/original.thresh.dat'%glif_dir,'w')
    
    for i in range(len(times)):
        t = times[i]
        v = voltage[i]
        th = threshold[i]
        v_file.write('%s\t%s\n'%(t,v))
        th_file.write('%s\t%s\n'%(t,th))
    v_file.close()
    th_file.close()
    
    pynml.generate_plot([times],
                            [voltage], 
                            "Membrane potential; %s"%info, 
                            colors = ['k'], 
                            xaxis = "Time (s)", 
                            yaxis = "Voltage (V)", 
                            grid = True,
                            show_plot_already=False,
                            save_figure_to='%s/MembranePotential_%ipA.png'%(glif_dir,curr_pA))
                
    pynml.generate_plot([times],
                            [threshold], 
                            "Threshold; %s"%info, 
                            colors = ['r'], 
                            xaxis = "Time (s)", 
                            yaxis = "Voltage (V)", 
                            grid = True,
                            show_plot_already=show_plot,
                            save_figure_to='%s/Threshold_%ipA.png'%(glif_dir,curr_pA))



if __name__ == '__main__':
    
    if '-all' in sys.argv:
        run_one_cell('473875489', 120, dt=5e-7, show_plot=False)
        run_one_cell('480629471', 50, dt=5e-7, show_plot=False)
        run_one_cell('480629475', 50, dt=5e-7, show_plot=False)
        run_one_cell('480633674', 120, dt=5e-7, show_plot=False)
        run_one_cell('486557295', 160, dt=5e-7, show_plot=False)
        run_one_cell('472451425', 180, dt=5e-7, show_plot=False)
        run_one_cell('472308324', 150, dt=5e-7, show_plot=False)
        run_one_cell('472455459', 120, dt=5e-7, show_plot=False)
        
        
        exit()
        
    if len(sys.argv)==3:

        glif_dir = sys.argv[1]
        curr_pA = float(sys.argv[2])
        if len(sys.argv)==4:
            dt=float(sys.argv[3])
        else:
            dt=5e-6
        show_plot = '-nogui' not in sys.argv
        run_one_cell(glif_dir, curr_pA, dt=dt, show_plot=show_plot)
    
    else:
        print(usage)