from allensdk.api.queries.biophysical_api import BiophysicalApi

import sys
import os
import json

import tables

import pickle

import pprint
pp = pprint.PrettyPrinter(indent=4)


def download():

    with open('/Users/katedoxey/Desktop/research/AIBS/parse_models/all_active_models.txt', 'rb') as fp:
        all_active_model_ids = pickle.load(fp)
    with open('/Users/katedoxey/Desktop/research/AIBS/parse_models/perisomatic_models.txt', 'rb') as fp:
        perisomatic_model_ids = pickle.load(fp)

    # perisomatic_model_ids.remove(486508789)
    # perisomatic_model_ids.remove(483108490)
    # perisomatic_model_ids.remove(486508702)
    # perisomatic_model_ids.remove(487245118)
    # perisomatic_model_ids.remove(487246046)
    # perisomatic_model_ids.remove(480630344)
    # perisomatic_model_ids.remove(489932551)
    # perisomatic_model_ids.remove(482520370)
    # perisomatic_model_ids.remove(489932682)
    # perisomatic_model_ids.remove(487245719)
    # perisomatic_model_ids.remove(486509958)
    # perisomatic_model_ids.remove(473871592)
    # perisomatic_model_ids.remove(488083972)
    # perisomatic_model_ids.remove(482657528)

    neuronal_model_ids = perisomatic_model_ids

    print("---- Downloading %i cell models..."%len(neuronal_model_ids))

    for neuronal_model_id in neuronal_model_ids:
        if os.path.isdir('%i'%neuronal_model_id):
            continue
        else:
            print("---- Downloading cell model: %s..."%neuronal_model_id)
            try:

                bp = BiophysicalApi('http://api.brain-map.org')
                bp.cache_stimulus = False # change to False to not download the large stimulus NWB file
                working_directory='%i'%neuronal_model_id
                bp.cache_data(neuronal_model_id, working_directory=working_directory)
                print("---- Saved model into %s, included NWB file: %s"%(working_directory,bp.cache_stimulus))

                with open(working_directory+'/manifest.json', "r") as json_file:
                    manifest_info = json.load(json_file)

                metadata={}
                # exp_id = int(manifest_info["biophys"][0]["model_file"][1][:9])
                model_file = manifest_info["biophys"][0]["model_file"][1]
                if model_file[0] == 'f':
                    exp_id = int(model_file[4:13])
                else:
                    exp_id = int(model_file[:9])
                metadata['exp_id'] = exp_id

                metadata['URL'] = 'http://celltypes.brain-map.org/mouse/experiment/electrophysiology/%s'%exp_id

                for m in manifest_info['manifest']:
                    if m['key']=="output_path":
                        nwb_file = working_directory+'/'+m['spec']

                if os.path.isfile(nwb_file):
                    print("---- Extracting metadate from NWB file: %s"%(nwb_file))
                    h5file=tables.open_file(nwb_file,mode='r')
                    metadata['AIBS:aibs_dendrite_type'] = str(h5file.root.general.aibs_dendrite_type.read())
                    metadata['AIBS:aibs_cre_line'] = str(h5file.root.general.aibs_cre_line.read())
                    metadata['AIBS:aibs_specimen_id'] = str(h5file.root.general.aibs_specimen_id.read())
                    metadata['AIBS:aibs_specimen_name'] = str(h5file.root.general.aibs_specimen_name.read())
                    metadata['AIBS:intracellular_ephys:Electrode 1:location'] = str(h5file.root.general.intracellular_ephys._v_children['Electrode 1'].location.read())
                    metadata['AIBS:session_id'] = str(h5file.root.general.session_id.read())
                    metadata['AIBS:subject:age'] = str(h5file.root.general.subject.age.read())
                    metadata['AIBS:subject:description'] = str(h5file.root.general.subject.description.read())
                    metadata['AIBS:subject:genotype'] = str(h5file.root.general.subject.genotype.read())
                    metadata['AIBS:subject:sex'] = str(h5file.root.general.subject.sex.read())
                    metadata['AIBS:subject:species'] = str(h5file.root.general.subject.species.read())
                else:
                    print("---- Can't find NWB file: %s!"%(nwb_file))

                print('    Metadata:')
                pp.pprint(metadata)
                with open(working_directory+'/metadata.json', 'w') as f:
                    json.dump(metadata, f, indent=4)



            except IndexError:
                print("Problem!")



if __name__ == '__main__':


    test = '-test' in sys.argv

    if test:
        '''
        dataset_ids = []
        ct = CellTypesApi()

        all_cells = ct.list_cells(require_morphology=True, require_reconstruction=True, reporter_status=None)

        count = 0
        for cell in all_cells:
            id = cell['id']
            print("===============================\nCell %i/%i (%s):"%(count, len(all_cells), id))
            dataset_ids.append(id)
            pp.pprint(cell)
            count+=1


        print("Found %i datasets: %s"%(len(dataset_ids),dataset_ids))
        '''
        all_biophys_models = []


        count = 0

        dataset_ids = [471141261,464198958,325941643,479704527]

        for dataset_id in dataset_ids:

            url = 'http://celltypes.brain-map.org/api/v2/data/Specimen/query.xml?criteria=model::Specimen,rma::criteria,%5Bid$eq' + \
                  str(dataset_id) + \
                  '%5D,specimen_types%5Bid$eq305008011%5D,rma::include,donor(transgenic_lines%5Btransgenic_line_type_code$eq%27D%27%5D),structure,cell_soma_locations,ephys_features,ephys_sweeps,ephys_result(well_known_files),neuronal_models(neuronal_model_template,neuronal_model_runs(well_known_files%5Bwell_known_file_type_id$eq481007198%5D))'

            import urllib2
            xml = urllib2.urlopen(url).read()

            #print xml

            import xml.etree.ElementTree as ET

            root = ET.fromstring(xml)

            print("\n\n Got XML: %s for dataset %s (%i/%i)"%(root.tag, dataset_id,count,len(dataset_ids)))
            for model in root.findall('./specimens/specimen/neuronal-models/neuronal-model'):
                for child in model:
                    print(child.tag, child.attrib, child.text)
                id = model.find('./id').text
                name = model.find('./name').text
                print("===========\n  Model: %s, %s"%(id,name))
                if 'Biophysical - perisomatic' in name:
                    print("Has biophysical/perisomatic model!")
                    all_biophys_models.append(id)
            count+=1

        print("Found %i biophysical models: %s"%(len(all_biophys_models),all_biophys_models))

    else:

        download()
