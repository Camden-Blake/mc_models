import os
import sys
sys.path.append("..")

import openmc
from beavrs_settings import set_settings, get_beavrs_model

beavrs_model = get_beavrs_model()

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

