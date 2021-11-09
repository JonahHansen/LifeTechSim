import matplotlib.pyplot as plt
from cycler import cycler
from pokepalette.pokemon_colours import pokemon_colours_dict as poke


def pokemon_colours(name):
    colours = poke[name].split(",")
    pokemon_cycler = cycler('color',colours)
    plt.rc('axes', prop_cycle=pokemon_cycler)
    return

#Pokemon colours!
