import mcnpy
import openmc

geometry_file = "geometry.xml"
materials_file = "materials.xml"
settings_file = "settings.xml"

mcnp_deck = mcnpy.translate_mcnp_openmc.translate_file(geometry_file, materials_file, settings_file)

mcnp_deck.write("beavrs.mcnp")