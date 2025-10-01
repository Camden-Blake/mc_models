import os
import sys
sys.path.append("..")

import openmc
from beavrs_settings import set_settings, get_beavrs_model

beavrs_model = get_beavrs_model()

## Swap TSLs for OTF data
tsl_change_map = {
    "c_H_in_H2O" : "tsl-HinH2O_OTF"
}

for mat_name, mat_data in beavrs_model.mats.items():
    if mat_data._sab:  # Check if _sab is not empty
        for i, item in enumerate(mat_data._sab):
            if item[0] in tsl_change_map:
                print(f"{mat_name} : Changing {item[0]} to {tsl_change_map[item[0]]}")
                mat_data._sab[i] = (tsl_change_map[item[0]], item[1])

beavrs_model.write_openmc_geometry()
beavrs_model.write_openmc_materials()
beavrs_model.write_openmc_settings()

model = openmc.Model.from_xml(geometry="geometry.xml", materials="materials.xml", settings="settings.xml")
os.remove('geometry.xml')
os.remove('materials.xml')
os.remove('settings.xml')

model.settings.temperature = {'method': 'interpolation'}

set_settings(model)
model.export_to_model_xml()