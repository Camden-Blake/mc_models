import mcnpy
import mcnpy.translate_mcnp_openmc
import openmc

mcnp_file = "rcf.mcnp"
mcnp_deck = mcnpy.Deck().read(mcnp_file)

geo, mat = mcnpy.translate_mcnp_openmc.mcnp_to_openmc(mcnp_deck)

model = openmc.Model(geo, mat)

model.export_to_model_xml()