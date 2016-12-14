# Copyright 2016 Autodesk Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

__all__ = ["EnergyModel", "energy_model", "BOLTZMANN_CONSTANT", "convert_temperature_K_to_C"]

import math

DEFAULT_TEMPERATURE_IN_KELVIN = 37.0 + 273.15
# 37 degrees C, in Kelvin.

BOLTZMANN_CONSTANT = 0.0019872041  
# kcal/(mol*K)

def convert_temperature_K_to_C( temperature_in_K ):
    return temperature_in_K - 273.15


def str_by_twos( iterable ):
    """Iterate over a string by consecutive pairs. Used for stack energy
calculations and maybe should be local there, but there may be other areas of
use. Looked for an equivalent function in itertools and didn't find anything
obvious."""
    iterable = iter(iterable)
    cur_item  = iterable.next()
    next_item = iterable.next()
    yield cur_item + next_item
    for item in iterable:
        cur_item = next_item
        next_item = item
        yield cur_item + next_item
    


class EnergyModel(object):
    def __init__(self):
        # TODO: (JMS 4/8/16) Add docstring.
        
        # These energies are from the following papers:
        # For Watson-Crick pairs, Table 1 and 2 in:
        # John SantaLucia Jr., "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics", Proceedings of the National Academy of Sciences 95.4 (1998): 1460-1465.
        # For GT wobble pairs, Table 5 in:
        # Hatim T. Allawi and John SantaLucia Jr., "Thermodynamics and NMR of Internal G-T Mismatches in DNA", Biochemistry 36 (1997): 10581-10594.

        # We use the following convention for storage, which is similar to the one used by Nupack (J. N. Zadeh, C. D. Steenberg, J. S. Bois, B. R. Wolfe, M. B. Pierce, A. R. Khan, R. M. Dirks, N. A. Pierce. NUPACK: analysis and design of nucleic acid systems. J Comput Chem, 32:170-173, 2011.)
        #
        # For a given NN stack: 
        #   Stacking 5' W Y 3'
        #            3' X Z 5'
        #   WX = AT CG GC TA GT TG (row index)
        #   YZ = AT CG GC TA GT TG (column index)
        # 
        self.stack_dG_37 = [
            [-1.00, -1.44, -1.28, -0.88,  0.71,  0.07],
            [-1.45, -1.84, -2.17, -1.28, -0.47, -0.32],
            [-1.30, -2.24, -1.84, -1.44,  0.08, -0.59],
            [-0.58, -1.30, -1.45, -1.00,  0.43,  0.34],
            [ 0.34, -0.59, -0.32,  0.07,  0.74,  1.15],
            [ 0.43,  0.08, -0.47,  0.71,  0.52,  0.74]
            ]
        # units of kcal/mol

        self.stack_dH = [
            [-7.9, -8.4,  -7.8, -7.2,  1.0, -2.5],
            [-8.5, -8.0, -10.6, -7.8, -4.1, -2.8],
            [-8.2, -9.8,  -8.0, -8.4,  3.3, -4.4],
            [-7.2, -8.2,  -8.5, -7.9, -0.1, -1.3],
            [-1.3, -4.4,  -2.8, -2.5,  5.8,  4.1],
            [-0.1,  3.3,  -4.1,  1.0, -1.4,  5.8]
        ]
        # units of kcal/mol

        self.stack_dS = [
            [-22.2, -22.4, -21.0, -20.4,   0.9,  -8.3],
            [-22.7, -19.9, -27.2, -21.0, -11.7,  -8.0],
            [-22.2, -24.4, -19.9, -22.4,  10.4, -12.3],
            [-21.3, -22.2, -22.7, -22.2,  -1.7,  -5.3],
            [ -5.3, -12.3,  -8.0,  -8.3,  16.3,   9.5],
            [ -1.7,  10.4, -11.7,   0.9,  -6.2,  16.3]
        ]
        # units of cal / K mol      <--- IMPORTANT

        self.pair_types = {
            "AT": 0,
            "TA": 3,
            "CG": 1,
            "GC": 2,
            "GT": 4,
            "TG": 5,
            "at": 0,
            "ta": 3,
            "cg": 1,
            "gc": 2,
            "gt": 4,
            "tg": 5
        }
        
        self.temperature_in_K = DEFAULT_TEMPERATURE_IN_KELVIN
    # end: def __init__()


    def pair_type( self, base_1, base_2 ):
        pair = base_1 + base_2
        if pair not in self.pair_types:
            return -1
        else:
            return self.pair_types[pair]
        
    def stack_energy( self, sequence_1, sequence_2 ):
        """ Computes the energy of a single base pair stack, or a helical region consisting only of stacks.

        sequence_1:   5'-XY-3'     "XY"
        sequence_2:   3'-WZ-5'     "ZW"

        Computes the nearest neighbor energy of the above stack. If the sequences are more than 2 bases, e.g.:

        sequence_1:   5'-GGCCA-3'    "GGCCA"
        sequence_2:   3'-CCGGT-5'    "TGGCC" 

        it will compute the NN energies for each stack in turn and return the sum. In this case, it would compute
        the energy of the four stacks below and return the sum:

        5'-GG-3'    5'-GC-3'    5'-CC-3'    5'-CA-3'
        3'-CC-5'    3'-CG-5'    3'-GG-5'    3'-GT-5'


        Args:
            sequence_1 (str): base identifiers for the top strand, 5' to 3' order
            sequence_2 (str): base identifiers for the bottom strand, 5' to 3' order

        Returns:
            3-tuple containing the computed dG_37, dH, and dS.
        """

        # TODO (JMS 4/8/16): Need to remove the duplicate dG calculations;
        # currently I get different results for calculating at a specific
        # temperature (37 deg C) vs what the reference dG has. It doesn't appear
        # to be a rounding/floating point issue, but something is off. However,
        # it does seem to agree roughly with Nupack calculations using the same
        # energy model.


        dG_37 = 0.0
        dG_check = 0.0
        dH = 0.0
        dS = 0.0

        for top,bot in zip( str_by_twos( sequence_1), str_by_twos(reversed(sequence_2))):
            pair_1 = self.pair_type(top[0],bot[0])
            pair_2 = self.pair_type(top[1],bot[1])

            dH_pair = self.stack_dH[pair_1][pair_2]
            dS_pair = self.stack_dS[pair_1][pair_2]
            dG_37_pair = dH_pair - self.temperature_in_K * dS_pair / 1000.0   # the division by 1000.0 is the unit change back to kcal/mol
            dG_check_pair = self.stack_dG_37[pair_1][pair_2]

            dG_check += dG_check_pair
            dG_37 += dG_37_pair
            dH += dH_pair
            dS += dS_pair

        return (dG_37, dG_check, dH, dS)


    def melting_temperature( self, dH, dS, staple_conc = 100e-9, scaffold_conc = 10e-9, sodium_conc = 1.0, magnesium_conc = 20e-3):
        """
        
        Units: dH is in kcal/mol, dS is in cal/mol (note difference). All concentrations are in M.
        """
        
        # more details on this derivation will be added to the doc string later.
        # TODO JMS 4/1.

        eff_staple_conc_at_tm = staple_conc - .5 * scaffold_conc
        conc_derived_term =  BOLTZMANN_CONSTANT * math.log( eff_staple_conc_at_tm )
        denominator = dS/1000.0 + conc_derived_term  # important unit conversion for dS
        return dH / denominator


energy_model = EnergyModel()


# TODO: Convert the below code to a unit test on the energy computation? (JMS 11/21/16)
#
# if __name__ == "__main__":
#     a = EnergyModel()
#     seq_1 = "G"*3
#     seq_2 = "C"*3
#     print a.stack_energy(seq_1,seq_2)
#     _,_,dH,dS = a.stack_energy(seq_1, seq_2)
#     print a.melting_temperature(dH,dS)
#     print a.melting_temperature(dH,dS) - 273.15
