
import os
import os.path
import json
import sys

import neuroml
from pyneuroml import pynml

sys.path.append("../data")
import data_helper as DH

cell_types = ['Izh', 'HH']

for cell_type in cell_types:
    
    for dataset_id in DH.CURRENT_DATASETS:

        target_sweep_numbers = DH.DATASET_TARGET_SWEEPS[dataset_id]

        net_id = "network_%s_%s"%(dataset_id, cell_type)
        net = neuroml.Network(id=net_id, type="networkWithTemperature", temperature=DH.SIMULATION_TEMPERATURE)

    
        net_doc = neuroml.NeuroMLDocument(id=net.id)
        net_doc.networks.append(net)
        
        net_doc.includes.append(neuroml.IncludeType('%s.cell.nml'%cell_type))

        number_cells = len(target_sweep_numbers)
        pop = neuroml.Population(id="Pop0",
                            component="RS",
                            size=number_cells,
                            type="populationList")
        net.populations.append(pop)
        for i in range(number_cells):
            location = neuroml.Location(x=100*i,y=0,z=0)
            pop.instances.append(neuroml.Instance(id=i,location=location))

        print target_sweep_numbers
        f = "../data/%s_analysis.json"%dataset_id
        with open(f, "r") as json_file:
            data = json.load(json_file) 

        id = data['data_set_id']
        sweeps = data['sweeps']

        print("Looking at data analysis in %s (dataset: %s)"%(f,id))

        index = 0
        for s in target_sweep_numbers:
            current = float(sweeps['%i'%s]["sweep_metadata"]["aibs_stimulus_amplitude_pa"])
            print("Sweep %s (%s pA)"%(s, current))
            
            stim_amp = "%s pA"%current
            input_id = ("input_%i"%s)
            pg = neuroml.PulseGenerator(id=input_id,
                                        delay="270ms",
                                        duration="1000ms",
                                        amplitude=stim_amp)
            net_doc.pulse_generators.append(pg)

            input_list = neuroml.InputList(id=input_id,
                                     component=pg.id,
                                     populations=pop.id)
            input = neuroml.Input(id='0', 
                                  target="../%s[%i]"%(pop.id, index), 
                                  destination="synapses")
            index+=1
            input_list.input.append(input)
            net.input_lists.append(input_list)

        net_file_name = 'prototypes/RS/%s.net.nml'%net_id
        pynml.write_neuroml2_file(net_doc, net_file_name)