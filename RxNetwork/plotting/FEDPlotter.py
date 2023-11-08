from network.Network import Network
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
from data.Calculator import Reaction

class FED_plotter:

    def __init__(self, reaction, FED_energy) -> None:
        self.reaction = reaction
        self.FED_energy = FED_energy
        # self.network = Network(reaction)
        self.reaction_data = {"NRR": {"N2":(0,0), # tuple represents (index, pathway)
                                      "N2*": (1,0),
                                      "N2H*": (2,0),
                                      "NNH2*": (3,0),
                                      "N*": (4,0),
                                      "NH*": (5,0),
                                      "NH2*": (6,0),
                                      "NH3*": (7,0),
                                      "2NH3": (8,0),
                                      "NHNH*": (3,1),
                                      "NHNH2*": (4,1),
                                      "NH2NH2*": (5,1)},
                                "HER": {"2H+":(0,0),
                                        "H*": (1,0),
                                        "H2": (2,0)},
                                }
    
    def plot(self, surface:str, bias:str, color="#f00000", linewidth=1, graph_objects=None, label_str=None):
        '''
        graph_objects is a tuple of (fig, ax) that can be passed through to add multiple plots to the same axis
        '''
        reaction = Reaction(self.reaction)
        initial_state, final_state = reaction.terminal_to_states()
        custom_line = [Line2D([0], [0], color=color, lw=4)]
        if label_str == None:
            custom_label = f"{surface.split('_')[0]} ({surface.split('_')[1]}) {bias}"
        elif label_str != None:
            custom_label = label_str
        state_width = 1
        connector_width = 1/2
        ticks = []
        tick_labels = []
        for i in range(self.reaction_data[self.reaction][final_state][0]+1):
            ticks.append(i*(state_width+connector_width)+state_width/2)
            state = list(self.reaction_data[self.reaction].keys())[list(self.reaction_data[self.reaction].values()).index((i,0))]
            tick_labels.append(f"{state}")


        if graph_objects == None:
            fig, ax = plt.subplots(dpi=300, figsize=(10,5))
        elif graph_objects != None: #ability to pass through an axis for adding multiple plots
            fig, ax = graph_objects 

        states = self.get_states()
        for state in states:
            state_index = self.reaction_data[self.reaction][state][0]
            pathway_index = self.reaction_data[self.reaction][state][1]
            if state_index != 0:
                # this loop will go back in the pathway until it finds a state that is in the FED data.
                for i in range(1, state_index+1):
                    previous_state_tuple = (state_index-i, pathway_index)
                    # the line below is a way of finding the keys associeted with a value, where the value is previous_state_tuple
                    previous_state = list(self.reaction_data[self.reaction].keys())[list(self.reaction_data[self.reaction].values()).index(previous_state_tuple)]
                    if previous_state in list(self.FED_energy.keys()):
                        previous_energy = self.FED_energy[previous_state]
                        previous_state_tuple = (state_index-i, pathway_index)
                        state_subtractor = i
                        break
                    else:
                        pass
                # previous_state = list(self.reaction_data[self.reaction].keys())[list(self.reaction_data[self.reaction].values()).index(previous_state_tuple)]
                # previous_energy = self.FED_energy[previous_state]
                state_x_coords = [state_index*(state_width+connector_width), state_index*(state_width+connector_width)+state_width]
                state_y_coords = [self.FED_energy[state], self.FED_energy[state]]
                connector_x_coords = [(state_index - state_subtractor + 1)*(state_width) + (state_index-state_subtractor)*connector_width, 
                                      state_index*(state_width+connector_width)]
                connector_y_coords = [previous_energy, self.FED_energy[state]]
                ax.hlines(state_y_coords[0], state_x_coords[0], state_x_coords[1], color=color, linewidth=linewidth)
                ax.plot(connector_x_coords, connector_y_coords, color=color, linewidth=linewidth, linestyle="--")
            elif state_index == 0:
                state_x_coords = [state_index*(state_width)+state_index, state_index*(state_width+connector_width)+state_width]
                state_y_coords = [self.FED_energy[state], self.FED_energy[state]]
                ax.hlines(state_y_coords[0], state_x_coords[0], state_x_coords[1], color=color, linewidth=linewidth, label=custom_label)
        
        # ax.legend(custom_line, [f"{surface.split('_')[0]} {surface.split('_')[1]} {bias}"])
        h, l = ax.get_legend_handles_labels()
        ax.legend(h, l)
        ax.set_xticks(ticks)
        ax.set_xticklabels(tick_labels)
        return fig, ax

    def get_states(self) -> list:
        return [i for i in self.FED_energy.keys()]



    