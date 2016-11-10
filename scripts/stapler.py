import os
import sys
import json
import logging
import argparse
import numpy as np

try:
    from nanodesign.converters.converter import Converter,ConverterFileFormats
except ImportError:
    import sys
    base_path = os.path.abspath( os.path.join( os.path.dirname(os.path.abspath( __file__)), '../'))
    sys.path.append(base_path)
    from nanodesign.converters.converter import Converter,ConverterFileFormats
    from nanodesign.converters.cadnano.reader import CadnanoReader
    from nanodesign.converters.cadnano.writer import CadnanoWriter
    from nanodesign.converters.cadnano.convert_design import CadnanoConvertDesign
    from nanodesign.data.parameters import DnaParameters
    sys.path = sys.path[:-1]

#reduce logging level TODO
#logging.disable(logging.ERROR) davep

class Stapler(object):
    """
        dna_structure         = dna_structure (dna_structure.py) from converter.cadnano_convert_design.create_structure
        path_list             = list of non-scaffold strands of strand domain lengths
        crossovers            = list of all crossovers, [ start helix_num, cadnano position start, start scaffold base_id, end helix_num, end cadnano position, end scaffold base_id, strand id, domain id, path_list id]
        crossovers_sorted     = list of all crossovers, [ start helix_num, cadnano position start, start scaffold base_id, end helix_num, end cadnano position, end scaffold base_id, strand id, domain id, path_list id], sorted by smaller scaffold base number in scaffold strand 0
        crossovers_joint      = list of all unique crossover locations, number of currently closed crossovers (1 : single crossover, 2 : double crossover)
        path_crossover_list   = list of paths of domains, [ crossover_joint index of crossover at 3' end of domain , [ crossovers entry for crossover at 3' at end of domain ]
        strand_index          = strand id of path i
        
        template              = template staple design, list of domain lengths 5' to 3'
        temperature           = monte carlo temperature ( in kT ? )
        energy                = list of paths: energies of paths of current state
        path_probabilities    = probablity to chose path i for step
        uncovered_domain_penalty = energy penalty for domains that are not covered by any template
        template_overhang_penalty_factor = energy penalty for template domains not associated to any path domain = length of template domain * template_overhang_penalty_factor
        temperature             = system temperature for monte carlo probability calculation
        step_probabilities    = probabily of each step type during monte carlo
        standard_domain_length = preffered length of domains during breaking that are not covered by any template
        
        template_types        = list of paths of template types placed
        template_positions    = list of paths of template positions
        template_offsets      = list of template offsets relative to template position
        break_positions       = list of paths of strand break positions
        single_strand_side_exception = special case handling, list of paths, 1 if only strand in a circular strand and template location is identical to single domain break location
            in the case of a circular strand, with a single break, a state with a break location i that is even (on a domain, not a crossover) and an identical template location i (on the broken segment),
            the state is not uniquely identifiable. The template can either be placed left of the break, or right of the break. In this case:
            if single_strand_side_exception = 1, place on the left of break
            if single_strand_side_exception = 0, place on the right of break
        
        crossover connectivity handling:
            when path i crossover j is changed, change crossovers_joint[path_crossover_list[i][j][0]]. If change towards 0, connectivity has been broken, change is forbidden
            double crossovers map to the same crossovers_joint field
    """

    def __init__(self,dna_structure,cadnano_design):
        self.dna_structure = dna_structure
        self.cadnano_design = cadnano_design
        self.path_list = []
        self.crossovers = []
        self.crossovers_sorted = []
        self.crossovers_joint = [0]
        self.path_crossover_list = []
        self.strand_index = []
        
        self.template = []
        self.energy = []
        self.path_probabilities = []
        self.uncovered_domain_penalty = 0.0
        self.template_overhang_penalty_factor = 0.0
        self.temperature = 1.0
        self.step_probabilities = []
        self.standard_domain_length = 1
        
        self.template_types = []
        self.template_positions = []
        self.template_offsets = []
        self.break_positions = []
        self.single_strand_side_exception = []
        self._logger = self._setup_logging()

        self._generate_paths()
        self._join_crossovers()

    def _setup_logging(self):
        """ Set up logging."""
        logger = logging.getLogger("stapler")
        logger.setLevel(logging.INFO)
    
        # create console handler and set format
        console_handler = logging.StreamHandler()
        formatter = logging.Formatter('[%(name)s] %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
        return logger
    
    def _generate_paths(self):
        """ create path_list, crossovers, crossovers_sorted, path_crossover_list, strand_index
        """
        #self._logger.setLevel(logging.DEBUG)
        self._logger.debug("=====================  generate paths =====================")
        path_number = 0

        for strand_nr, strand in enumerate(self.dna_structure.strands):
            if not strand.is_scaffold:
                path = []                      #list of domain lengths
                path_crossovers = []           #list of crossover placeholders
                self._logger.debug("")
                self._logger.debug("---------- strand %d ---------" % strand.id)
                self._logger.debug("Number of domains %d " % len(strand.domain_list))
                
                for domain_nr, domain in enumerate(strand.domain_list):
                    path.append(len(domain.base_list))      #add domain length
                    
                    # Collect crossover information, add to crossover list.
                    if strand.is_circular or domain_nr < ( len(strand.domain_list) - 1 ):
                        helix1 = domain.helix.lattice_num
                        helix2 = strand.domain_list[( domain_nr + 1 ) % len(strand.domain_list)].helix.lattice_num
                        base1 = domain.base_list[-1].p
                        base2 = strand.domain_list[( domain_nr + 1 ) % len(strand.domain_list)].base_list[0].p
                        """ TODO
                            scaffold strand is required to be strand 0
                        """
                        scaff_base1 = self.dna_structure.strands[0].get_base_index(domain.base_list[-1].across)
                        scaff_base2 = self.dna_structure.strands[0].get_base_index(strand.domain_list[( domain_nr + 1 ) % len(strand.domain_list)].base_list[0].across)
                       
                        self.crossovers.append([helix1, base1, scaff_base1, helix2, base2, scaff_base2, strand.id, domain.id, path_number, domain_nr])
                        path_crossovers.append([[],self.crossovers[-1]])  #add crossover placeholder for path crossover list, crossover information

                self.path_list.append(path)
                self.path_crossover_list.append(path_crossovers)
                path_number+=1
                self.strand_index.append(strand_nr)
    
        self.crossovers_sorted=sorted(self.crossovers, key=lambda val: (min(val[2], val[5])))       #sort crossovers according to smaller scaffold base number

    def _join_crossovers(self):
        """ create entries in path_crossover_list[path][domain][0] referincing unique crossover position in crossovers_joint, add +1 to crossovers_joint for each closed crossover in unique crossover position
        """
        counter = 0
        for index, crossover in enumerate(self.crossovers_sorted[:-1]):
            self.crossovers_joint[counter] += 1
            self.path_crossover_list[ crossover[8] ][ crossover[9] ][0] = counter

            if ( min(crossover[2], crossover[5]) + 1 != min(self.crossovers_sorted[index + 1][2], self.crossovers_sorted[index + 1][5]) ) or ( max(crossover[2], crossover[5]) - 1 != max(self.crossovers_sorted[index + 1][2], self.crossovers_sorted[index + 1][5] ) ):
                counter += 1
                self.crossovers_joint.append(0)
        #add last crossover of structure
        if self.crossovers_sorted:              #check if any crossovers exist at all
            self.crossovers_joint[counter] += 1
            self.path_crossover_list[ self.crossovers_sorted[-1][8] ][ self.crossovers_sorted[-1][9] ][0] = counter
    
    def initialize_system(self):
        """ initialize system: place 1 template in random location with random offset on each path
            add a single strand break in all circular strands to make them linear/physical

        self.template_positions             = domains on which template are projected onto, list of paths of positions
        self.template_offsets               = template domain that is attached at template_position, list of paths of offsets
        self.break_positions                = locations of strand breaks. locations are alternating: domain, crossover | first location is domain 0, list of paths of positions
        self.single_strand_side_exception   = track exception case of circular strand with single strand and identical break and template position
        
        template position can be position from ( 0 ) to ( nr_domains - 1 )
        template offset can be from ( 0 ) to ( len(template) - 1 )
        
        template position:      -0-
        template offset:        -2-
        
        template:       -0- -1- -2-
        strand:                 -0- -1- -2- -3- -4-
        
        template position:                      -4-
        template offset:                        -1-
        
        template:                           -0- -1- -2-
        strand:                 -0- -1- -2- -3- -4-
        
        circular strands: template wraps around:
        
        template position:                      -4-
        template offset:                        -1-
        
        template:               -2-         -0- -1-
        strand:                 -0- -1- -2- -3- -4-
        
        template position:          -2-
        template offset:            -2-
        
        template:               -1- -2-         -0-
        strand:                 -0- -1- -2- -3- -4-
        
        possible break positions:
        
        template:               -1- -2-         -0-
        strand:                 -0- -1- -2- -3- -4-
        break pos:               0 1 2 3 4 5 6 7 8 (9) <- static in linear strand, variable in circular strand

        TODO:   states are entropically weighted with their number of domains within the target unbroken path part, because of multiple representations for the same system state.
                Therefore overhanging states are aditionally penalized
        """
        
        for path_nr, path in enumerate(self.path_list):
            template_type = np.random.randint( 0 , len( self.template ) )
            self.template_types.append( [ template_type ] )   #type of template placed on current part of current path
            self.template_positions.append( [ np.random.randint( 0 , len( path ) ) ] )   #random placement of initial template on each path
            self.template_offsets.append([ np.random.randint( 0, len( self.template[ template_type ] ) ) ] ) #random offset of template to position
            
            if self.dna_structure.strands[self.strand_index[path_nr]].is_circular:  #add one break to circular strands
                break_pos_temp = np.random.randint( 0, 2 * len(path) )
                
                # If crossver break, check if double crossover, if not, reroll position.
                while ( ( break_pos_temp % 2 ) == 1 ) and ( self.crossovers_joint[ self.path_crossover_list[path_nr][ break_pos_temp // 2 ][0] ] < 2 ):
                    break_pos_temp = np.random.randint( 0, 2 * len(path) )
            
                # Add break, reduce crossover counter in crossovers_joint if crossover break.
                self.break_positions.append([break_pos_temp])
                if ( break_pos_temp % 2 ) == 1:
                    self.crossovers_joint[ self.path_crossover_list[path_nr][ break_pos_temp // 2 ][0] ] -= 1

                if self.break_positions[path_nr][0] / 2 == self.template_positions[path_nr][0] and self.break_positions[path_nr][0] % 2 == 0:
                    self.single_strand_side_exception.append( np.random.randint(0, 2) )
                else:
                    self.single_strand_side_exception.append(0)
        
            else:
                self.break_positions.append([ 2 * len(path) - 1 ])       #place a permanent break at final possible break position for linear strands
                self.single_strand_side_exception.append(0)

        # Calculate total energy of system.
        for path_nr in range(len(self.path_list)):
            self.energy.append( self._calculate_energy(self.template_positions,self.template_offsets,self.break_positions,self.single_strand_side_exception,self.template_types,path_nr) )
        
        domain_number = 0       #calculate total number of path domains
        for path in self.path_list:
            domain_number += len(path)

        for path in self.path_list: #set probabilities for each path to be chosen during monte-carlo steps. paths weighted with number of domains
            self.path_probabilities.append( float( len(path) ) / domain_number )
                
        print "Starting state energy",sum(self.energy),self.energy #TODO

    def _calculate_energy(self,template_positions,template_offsets,break_positions,single_strand_side_exception,template_types,path_nr):
        """ calculate energy of given path in current system state
            fully unbound template domains add ( template_overhang_penalty_factor * length ) to energy
            mismatches in lengths between domain and template domain add ( 1.0 * length ) to energy
        """
        
        energy=0
        
        path = self.path_list[path_nr]
           
        path_copy=path[:]    #make a copy of the path domain lengths
        
        for template_nr,template_pos in enumerate(template_positions[path_nr]): #project template domains onto current unbroken strand section domains
            
            template_off = template_offsets[path_nr][template_nr]
            
            template_type = template_types[path_nr][template_nr]
            
            """ start_pos is first domain of current unbroken strand section
                end_pos is last domain of current unbroken strand section
                start_pos_break is break index left of current unbroken strand section
                end_pos_break is break index right of current unbroken strand section
            """
            
            start_pos_break = break_positions[path_nr][ ( template_nr - 1 ) % len(break_positions[path_nr]) ]
            start_pos = ( start_pos_break + 1 ) // 2
                
            end_pos_break = break_positions[path_nr][template_nr]
            end_pos = end_pos_break // 2

            template_copy=self.template[template_type][:]   #make a copy of the template domain lengths

            """ current_pos is postion in the coordinate system of the template strand, some examples:
                linear strand, template_pos = 0, template_offset = 2 , no added breaks
                
                template position:      -0-
                template offset:        -2-
                current_pos:             2   3   4   5   6
                template:       -0- -1- -2- -3-
                strand:                 -0- -1- -2- -3- -4-
                
                linear strand, template_pos = 2, template_offset = 0, added break at break_pos 7
                template position:              -2-
                template offset:                -0-
                current_pos:            -2  -1   0   1
                template:                       -0- -1-   -2- -3-
                strand:                 -0- -1- -2- -3- | -4-
                
                circular strand, template_pos = 0, template_offset = 2, added break at break_pos 3
                
                template position:      -0-
                template offset:        -2-
                current_pos:             2   3    -1   0   1
                template:               -2- -3-       -0- -1-
                strand:                 -0- -1- | -2- -3- -4-
                
                circular strand, template_pos = 4, template_offset = 1, added break at break_pos 3
                
                template position:                        -4-
                template offset:                          -1-
                current_pos:             2   3    -1   0   1
                template:               -2- -3-       -0- -1-
                strand:                 -0- -1- | -2- -3- -4-
                
                circular strand, template_pos = 1, template_offset = 1, added break at break_pos 3
                
                template position:          -1-
                template offset:            -0-
                current_pos:            -1   0    -4  -3  -2
                template:                   -0-   -1- -2- -3-
                strand:                 -0- -1- | -2- -3- -4-
                
                circular strand, template_pos = 2, template_offset = 1, added break at break_pos 3
                
                template position:                -2-
                template offset:                  -1-
                current_pos:             4   5     1   2   3
                template:                   -0-   -1- -2- -3-
                strand:                 -0- -1- | -2- -3- -4-
                
                circular strand, template_pos = 1, template_offset = 0, added break at break_pos 2, exception = 1
                
                template position:          -1-
                template offset:            -0-
                current_pos:            -1   0  -5  -4  -3  -2
                template:                   -0- -1- -2- -3-
                strand:                 -0- -1 | 1- -2- -3- -4-
                
                circular strand, template_pos = 1, template_offset = 0, added break at break_pos 2, exception = 0
                
                template position:              -1-
                template offset:                -0-
                current_pos:             4   5   0   1   2   3
                template:                       -0- -1- -2- -3-
                strand:                 -0- -1 | 1- -2- -3- -4-
            """
            
            #shift positions if wrap around case
            if start_pos_break >= end_pos_break:
                if ( template_pos >= start_pos ) and ( single_strand_side_exception[path_nr] == 0 ):
                    end_pos += len(self.path_list[path_nr])
                else:
                    start_pos -= len(self.path_list[path_nr])
        
            for current_pos in xrange(start_pos, end_pos + 1):
                #current_pos in template coordinate system
                current_template_pos = current_pos - template_pos + template_off
                
                #shift positions back to real coordinates
                current_pos_shift = current_pos % len(self.path_list[path_nr])
            
                #template reaches current domain
                if ( current_template_pos >= 0 ) and ( current_template_pos < len(self.template[template_type]) ):
                    
                    path_copy[ current_pos_shift ] -= self.template[ template_type ][ current_template_pos ]
                    template_copy[ current_template_pos ] -= self.template[ template_type ][ current_template_pos ]
                #domain not covered, penalize
                else:
                    energy += self.uncovered_domain_penalty
    
            for domain in template_copy:
                energy += abs(domain) * self.template_overhang_penalty_factor

        for domain in path_copy:
            energy += abs(domain)

        return energy
            
    def _take_step(self):
        """ do a monte carlo type step
            change_type determines type of change in current step:
                0 = add break to path
                1 = remove break to path
                2 = move existing break
                3 = move existing template
                4 = change template type
        """
        
        # Determine path to change.
        current_path = np.random.choice(range(len(self.path_list)), None, p = self.path_probabilities)
        current_strand = self.strand_index[current_path]
        
        # Determine type of change.
        change_type = np.random.choice((0,1,2,3,4), None, p = self.step_probabilities)
        
        # Copy system state.
        crossovers_joint_temp = self.crossovers_joint[:]
        
        template_positions_temp = []
        for element in self.template_positions:
            template_positions_temp.append(element[:])
        
        template_offsets_temp = []
        for element in self.template_offsets:
            template_offsets_temp.append(element[:])
        
        template_types_temp = []
        for element in self.template_types:
            template_types_temp.append(element[:])
        
        break_positions_temp = []
        for element in self.break_positions:
            break_positions_temp.append(element[:])
        
        single_strand_side_exception_temp = self.single_strand_side_exception[:]
        
        #keep track if step was forbidden because not physical
        forbidden = False
       
        if change_type == 0:        #add new break to path
            
            # List of non-broken locations in path.
            break_positions_available = [pos for pos in range( 2 * len(self.path_list[current_path]) ) if pos not in break_positions_temp[current_path]]

            if not break_positions_available:   #no positions available
                print "current path is broken at every position" #TODO
                forbidden = True
            else:
                break_pos_new = np.random.choice( break_positions_available )          #roll new break position
                
                # Check if single crossover in unique crossover location.
                if ( ( break_pos_new % 2 ) == 1 ) and ( crossovers_joint_temp[ self.path_crossover_list[current_path][ break_pos_new // 2 ][0] ] < 2 ):
                    forbidden = True
                else:
                    insertion_index = -1

                    # Reduce crossover counter in crossovers_joint_temp if crossover break.
                    if ( break_pos_new % 2 ) == 1:
                        crossovers_joint_temp[ self.path_crossover_list[current_path][ break_pos_new // 2 ][0] ] -= 1

                    # Unbroken (single break) strand, insert break directly at index 0 of break_positions.
                    if len(break_positions_temp[current_path]) == 1:
                        break_positions_temp[current_path].insert(0, break_pos_new)
                        insertion_index = 0
                        break_pos_left = break_positions_temp[current_path][1]
                        break_pos_right = break_pos_left
                
                    # Find insertion index for new break_pos
                    else:
                        #find break positions larger than new break position
                        breaks_larger = [[index, break_pos] for index, break_pos in enumerate(break_positions_temp[current_path]) if break_pos > break_pos_new]
                    
                        if not breaks_larger:   #new break_pos is largest break_pos so far (must be circular strand), insert after previously largest break_pos
                        
                            #find break positions smaller than new break position
                            breaks_smaller = [[index, break_pos] for index, break_pos in enumerate(break_positions_temp[current_path]) if break_pos < break_pos_new]
                        
                            break_smaller_largest = max(breaks_smaller, key = lambda value: value[1])
                            break_positions_temp[current_path].insert(break_smaller_largest[0] + 1, break_pos_new)
                            insertion_index = break_smaller_largest[0] + 1
                            break_pos_left = break_smaller_largest[1]
                            break_pos_right = break_positions_temp[current_path][ ( break_smaller_largest[0] + 2 ) % len(break_positions_temp[current_path]) ]
                    
                        else:  #new break_pos is not largest break_pos so far, insert before smallest break_pos larger than new break_pos
                        
                            break_larger_smallest = min(breaks_larger, key = lambda value: value[1])
                            break_positions_temp[current_path].insert(break_larger_smallest[0], break_pos_new)
                            insertion_index = break_larger_smallest[0]
                            break_pos_left = break_positions_temp[current_path][ ( break_larger_smallest[0] - 1 ) % len(break_positions_temp[current_path]) ]
                            break_pos_right = break_larger_smallest[1]

                    #current unbroken path wraps around, map values to continous range
                    if break_pos_left >= break_pos_right:
                        break_pos_right += 2 * len(self.path_list[current_path])
                        if break_pos_new < break_pos_left:
                            break_pos_new += 2 * len(self.path_list[current_path])
                    
                    #possible domains for template_positions (modulo to shift domains back to real range if positions were shifted to continuous range after wrapping around)
                    domains_left = [domains % len(self.path_list[current_path]) for domains in range( ( break_pos_left + 1 ) // 2, ( break_pos_new // 2 ) + 1 )]
                    domains_right = [domains % len(self.path_list[current_path]) for domains in range( ( break_pos_new + 1 ) // 2, ( break_pos_right // 2 ) + 1 )]

                    template_pos_current = template_positions_temp[current_path][ insertion_index % len(template_positions_temp[current_path]) ]  #save current template_pos
                    template_off_current = template_offsets_temp[current_path][ insertion_index % len(template_positions_temp[current_path]) ]    #save current template_off
                    template_type_current = template_types_temp[current_path][ insertion_index % len(template_positions_temp[current_path]) ]    #save current template_type
                    
                    template_type_new = np.random.randint( 0, len(self.template) )      #roll new template type
                    template_off_new = np.random.randint( 0 , len(self.template[template_type_new]) )  #roll new template_off

                    # existing template_pos is left of new break_pos
                    if template_positions_temp[current_path][insertion_index % len(template_positions_temp[current_path])] in domains_left:
                        if single_strand_side_exception_temp[current_path] == 1:
                            """ this case should only occur if the following are true: circular strand, single break, break in domain, template domain identical to break domain, template to the left of break
                                
                                template is found in left side, but single_strand_exception = 1 means it is left of its break, = right side from new break
                                insert new template_pos left of current template
                                """
                            template_pos_new = np.random.choice( domains_left )             #roll new template_pos
                        
                            template_positions_temp[current_path].insert(insertion_index, template_pos_new )    #insert new template_pos left
                            template_offsets_temp[current_path].insert(insertion_index, template_off_new )      #insert new template_off left
                            template_types_temp[current_path].insert(insertion_index, template_type_new )       #insert new template type left
                    
                            template_positions_temp[current_path][ ( insertion_index + 1 ) % len(template_positions_temp[current_path]) ] = template_pos_current    #insert current template_pos right
                            template_offsets_temp[current_path][ ( insertion_index + 1 ) % len(template_positions_temp[current_path]) ] = template_off_current      #insert current template_off right
                            template_types_temp[current_path][ ( insertion_index + 1 ) % len(template_positions_temp[current_path]) ] = template_type_current      #insert current template_type right
                    
                            single_strand_side_exception_temp[current_path] = 0   #change single_strand_exception to 0
                        else:
                            """ template is found in left side, no single_strand exception, put new template in right side, insert right of current template
                                """
                            template_pos_new = np.random.choice( domains_right )            #roll new template_pos

                            template_positions_temp[current_path].insert(insertion_index, template_pos_current )    #insert current template_pos left
                            template_offsets_temp[current_path].insert(insertion_index, template_off_current )      #insert current template_off left
                            template_types_temp[current_path].insert(insertion_index, template_type_current )      #insert current template_off left
                        
                            template_positions_temp[current_path][ ( insertion_index + 1 ) % len(template_positions_temp[current_path]) ] = template_pos_new    #insert new template_pos right
                            template_offsets_temp[current_path][ ( insertion_index + 1 ) % len(template_positions_temp[current_path]) ] = template_off_new    #insert new template_off right
                            template_types_temp[current_path][ ( insertion_index + 1 ) % len(template_positions_temp[current_path]) ] = template_type_new    #insert new template_off right
                            #existing template_pos is right of new break_pos
                    else:
                        template_pos_new = np.random.choice( domains_left )             #roll new template_pos
                    
                        template_positions_temp[current_path].insert(insertion_index, template_pos_new )    #insert new template_pos left
                        template_offsets_temp[current_path].insert(insertion_index, template_off_new )      #insert new template_off lefft
                        template_types_temp[current_path].insert(insertion_index, template_type_new )      #insert new template_off lefft
                    
                        template_positions_temp[current_path][ ( insertion_index + 1 ) % len(template_positions_temp[current_path]) ] = template_pos_current    #insert current template_pos right
                        template_offsets_temp[current_path][ ( insertion_index + 1 ) % len(template_positions_temp[current_path]) ] = template_off_current      #insert current template_off right
                        template_types_temp[current_path][ ( insertion_index + 1 ) % len(template_positions_temp[current_path]) ] = template_type_current      #insert current template_off right

        if change_type == 1:    #remove break from path
            
            if len(break_positions_temp[current_path]) < 2:
                #print "current path has no free breaks left to remove" #TODO
                forbidden = True
            else:
                #roll break to remove
                break_index = np.random.randint( 0 , len(break_positions_temp[current_path]) - ( not self.dna_structure.strands[current_strand].is_circular )  )
            
                #roll which template (left or right of break) to delete
                template_left_or_right = np.random.choice((0,1))
                template_index_to_delete = ( break_index + template_left_or_right ) % len(template_positions_temp[current_path])
            
                #circular strand and single break will remain: check for single_strand_side_exception status
                if ( len(template_positions_temp[current_path]) == 2 ) and self.dna_structure.strands[current_strand].is_circular:
                    template_index_to_keep = 1 - template_index_to_delete
                    break_index_to_keep = 1 - break_index
                    #remaining break and template are on same domain, save on which side remaining template is relative to remaining break: 0 = right, 1 = left
                    if ( template_positions_temp[current_path][template_index_to_keep] == ( break_positions_temp[current_path][break_index_to_keep] // 2 ) ) and ( ( break_positions_temp[current_path][break_index_to_keep] % 2 ) == 0 ):
                        single_strand_side_exception_temp[current_path] = 1 - template_left_or_right
            
                #increase crossover counter in crossovers_joint_temp if crossover break
                if ( break_positions_temp[current_path][break_index] % 2 ) == 1:
                        crossovers_joint_temp[ self.path_crossover_list[current_path][ break_positions_temp[current_path][break_index] // 2 ][0] ] += 1
                
                #delete chosen break and template
                del break_positions_temp[current_path][break_index]
                del template_positions_temp[current_path][template_index_to_delete]
                del template_offsets_temp[current_path][template_index_to_delete]
                del template_types_temp[current_path][template_index_to_delete]

        if change_type == 2:    #move break
        
            if ( len(break_positions_temp[current_path]) - ( not self.dna_structure.strands[current_strand].is_circular ) ) == 0:
                #can't move, is linear strand with only ending break
                forbidden = True
            else:
                #roll break to move
                break_index = np.random.randint( 0 , len(break_positions_temp[current_path]) - ( not self.dna_structure.strands[current_strand].is_circular )  )
                break_pos = break_positions_temp[current_path][break_index]
                
                """ TODO
                    choose step direction and distance
                    currently included stepsize > 1 to prevent getting stuck. This may not be necessary
                """
                step = np.random.choice((-5,-4,-3,-2,-1,1,2,3,4,5),p=(0.025,0.025,0.05,0.15,0.25,0.25,0.15,0.05,0.025,0.025)) #roll step distance
                
                #break_pos left and right of current break_pos
                break_limit_left = break_positions_temp[current_path][ ( break_index - 1 ) % len(break_positions_temp[current_path]) ]
                break_limit_right = break_positions_temp[current_path][ ( break_index + 1 ) % len(break_positions_temp[current_path]) ]
                
                template_limit_left = 2 * template_positions_temp[current_path][ break_index ]
                template_limit_right = 2 * template_positions_temp[current_path][ ( break_index + 1 ) % len(break_positions_temp[current_path]) ]
                
                # if limiting breaks are on limiting domains, move limits closer (modulo wrapping omitted because done in next step anyway
                limit_left = template_limit_left + ( template_limit_left == break_limit_left )
                limit_right = template_limit_right - ( template_limit_right == break_limit_right )
                
                #shift limit by path length * 2 if wrap around case
                if template_limit_left >= template_limit_right:
                    if break_pos > limit_left  or ( single_strand_side_exception_temp[current_path] == 1 ):
                        limit_right += 2 * len(self.path_list[current_path])
                    else:
                        limit_left -= 2 * len(self.path_list[current_path])
            
                break_pos_new = break_pos + step

                #check if new break_pos is outside of limits, or if break is trapped between two templates on same domain
                if ( break_pos_new > limit_right ) or ( break_pos_new < limit_left ) or ( ( break_pos == template_limit_left ) and ( break_pos == template_limit_right ) and ( len(break_positions_temp[current_path]) > 1 ) ):
                    forbidden = True
                else:
                    #shift new break_pos back to real coordinates
                    break_pos_new = break_pos_new % ( 2 * len(self.path_list[current_path]) )
                
                    #check if crossover is broken, and if yes if its a single crossover
                    if ( break_pos_new % 2 ) == 1:  #two step check, because linear strands do not have crossover entry for last segment
                        if crossovers_joint_temp[ self.path_crossover_list[current_path][ break_pos_new // 2 ][0] ] == 1:
                            forbidden = True
                    if not forbidden:
                        #adjust unique crossover position counters
                        if ( break_pos % 2 ) == 1:
                            crossovers_joint_temp[ self.path_crossover_list[current_path][ break_pos // 2 ][0] ] += 1
                
                        if ( break_pos_new % 2 ) == 1:
                            crossovers_joint_temp[ self.path_crossover_list[current_path][ break_pos_new // 2 ][0] ] -= 1
                
                        #change break_pos
                        break_positions_temp[current_path][break_index] = break_pos_new

                        """" check if circular strand and only one break, control single_strand_side_exception status
                            if = 1, all moves will change state to 0 (even a full circle around a single-segment path)
                            if = 0, if new domain = template domain, if move to the left, set to 1, if move to the right, set to 0
                        """
                        if self.dna_structure.strands[current_strand].is_circular and ( len(break_positions_temp[current_path]) == 1 ):
                            if single_strand_side_exception_temp[current_path] == 1:
                                single_strand_side_exception_temp[current_path] = 0
                            elif ( template_positions_temp[current_path][0] == ( break_positions_temp[current_path][0] // 2 ) ) and ( ( break_positions_temp[current_path][0] % 2 ) == 0 ):
                                if step < 0:
                                    single_strand_side_exception_temp[current_path] = 1

        if change_type == 3:    #move template
            
            #roll template to move
            template_index = np.random.randint( 0 , len(template_positions_temp[current_path]) )
            template_pos = template_positions_temp[current_path][template_index]

            #break_pos left and right of current template_pos
            limit_right = break_positions_temp[current_path][ template_index ]
            limit_left = break_positions_temp[current_path][ ( template_index - 1 ) % len(break_positions_temp[current_path]) ]

            step_pos = 0
            step_off = 0
            while ( step_pos == 0 ) and ( step_off == 0 ):
                
                """ TODO
                    choose step direction and distance
                    currently included stepsize > 1 to prevent getting stuck. This may not be necessary
                """
                step_pos = np.random.choice((-5,-4,-3,-2,-1,0,1,2,3,4,5),p=(0.025,0.025,0.05,0.1,0.2,0.2,0.2,0.1,0.05,0.025,0.025)) #roll step distance template_pos
                step_off = np.random.choice((-5,-4,-3,-2,-1,0,1,2,3,4,5),p=(0.025,0.025,0.05,0.1,0.2,0.2,0.2,0.1,0.05,0.025,0.025)) #roll step distance template_off
        
            #shift limit by path length * 2 if wrap around case
            if limit_left >= limit_right:
                if template_pos * 2 >= limit_left  and ( single_strand_side_exception_temp[current_path] == 0 ):
                    limit_right += 2 * len(self.path_list[current_path])
                else:
                    limit_left -= 2 * len(self.path_list[current_path])
                            
            limit_left = ( limit_left + 1 ) // 2
            limit_right = limit_right // 2
            template_pos_new = template_pos + step_pos
            template_off_new = template_offsets_temp[current_path][template_index] + step_off

            #check if new position and new offset are in physical limits
            if ( template_pos_new > limit_right ) or ( template_pos_new < limit_left ) or ( template_off_new < 0 ) or ( template_off_new >= len(self.template[ template_types_temp[current_path][template_index] ]) ):
                forbidden = True
            else:
                #shift new template_pos back to real coordinates
                template_pos_new = template_pos_new % len(self.path_list[current_path])
                    
                #change template_pos, template_off
                template_positions_temp[current_path][template_index] = template_pos_new
                template_offsets_temp[current_path][template_index] = template_off_new
                        
                """" check if circular and only break, control single_strand_side_exception status
                    if = 1, all moves will change state to 0 (even a full circle around a single-segment path)
                    if = 0, if new domain = template domain, if move to the left, set to 1, if move to the right, set to 0
                """
                if self.dna_structure.strands[current_strand].is_circular and ( len(template_positions_temp[current_path]) == 1 ):
                    if single_strand_side_exception_temp[current_path] == 1:
                        single_strand_side_exception_temp[current_path] = 0
                    elif ( template_positions_temp[current_path][0] == ( break_positions_temp[current_path][0] // 2 ) ) and ( ( break_positions_temp[current_path][0] % 2 ) == 0 ):
                        if step_pos > 0:
                            single_strand_side_exception_temp[current_path] = 1

        if change_type == 4:    #change template type
        
            template_type_new = np.random.randint( 0, len(self.template) )  #roll new template type
            template_index = np.random.randint( 0 , len(template_positions_temp[current_path]) ) #roll template index
            
            #check if current template offset is compatible with new template type = (offset < length)
            if template_offsets_temp[current_path][template_index] < len(self.template[template_type_new]):
                template_types_temp[current_path][template_index] = template_type_new
            else:
                forbidden = True

        energy_new = self._calculate_energy(template_positions_temp,template_offsets_temp,break_positions_temp,single_strand_side_exception_temp,template_types_temp,current_path)
        energy_diff = energy_new - self.energy[current_path]
        
        p = np.exp(- energy_diff / self.temperature)
        
        random_draw = np.random.random()
        
        if p > random_draw:
            self.crossovers_joint = crossovers_joint_temp
            self.template_positions = template_positions_temp
            self.template_offsets = template_offsets_temp
            self.template_types = template_types_temp
            self.break_positions = break_positions_temp
            self.single_strand_side_exception = single_strand_side_exception_temp
            self.energy[current_path] = energy_new
            acceptance = ( True and not forbidden )
        else:
            acceptance = False

        return acceptance, change_type
    
    def generate(self,nr_steps,nr_steps_timescale,temperature_adjust_rate):
        #approximate average acceptance rate for each step type over last 10000 steps
        acceptance_rate = [0.0,0.0,0.0,0.0,0.0]
        #approximate average change in system state energy per step over last 10000 steps
        energy_rate = 0.0
        #count number of steps since system state energy didn't change
        counter = 0
        
        for step_nr in range(nr_steps):
            energy_previous = sum(self.energy)
            acceptance_current, type = self._take_step()
            energy_new = sum(self.energy)
            if energy_new == energy_previous:
                counter += 1
            else:
                counter = 0
            
            if counter == nr_steps_timescale:
                break
                    
            energy_rate = ( 1.0 - 1.0/nr_steps_timescale ) * energy_rate + ( 1.0/nr_steps_timescale ) * ( energy_new - energy_previous )
            acceptance_rate[type] = ( 1.0 - 1.0/nr_steps_timescale ) * acceptance_rate[type] + ( 1.0/nr_steps_timescale ) * acceptance_current
            if energy_rate > 0.0:
                self.temperature *= temperature_adjust_rate
            if step_nr % 1000 == 0:
                print "step",step_nr,"acceptance_rate",acceptance_rate,"step type",type,"energy",sum(self.energy),"temp",self.temperature,"e rate",energy_rate #TODO
        
        self._break_json()

        return

    def _break_json(self):
        """ write all breaks into cadnano_design
        """
        self._logger.setLevel(logging.DEBUG)
        self._logger.debug("=====================  break json =====================")
        
        for path_nr, path_breaks in enumerate(self.break_positions):
            self._logger.debug("")
            self._logger.debug("---------- path %d ----------" % path_nr)
            strand_nr = self.strand_index[ path_nr ]
            self._logger.debug("Strand number %d" % strand_nr)
            self._logger.debug("Number of path breaks %d" % len(path_breaks))
            is_circular = self.dna_structure.strands[strand_nr].is_circular
        
            # Break all breaks if circular strand, all breaks - 1 if linear strand.
            for break_nr,break_current in enumerate(path_breaks[0:len(path_breaks) - 1 + is_circular]): 
                domain = break_current // 2
                self._logger.debug("")
                self._logger.debug("----- Domain %d----- " % domain)
                self._logger.debug("break_current %d " % break_current)

                if break_current % 2 == 0:           #domain break
                    template_pos_left = self.template_positions[path_nr][break_nr]

                    # Shift coordinates if wrap around.
                    if ( template_pos_left > domain ) or ( ( template_pos_left == domain ) and ( len(path_breaks) == 1 ) and ( self.single_strand_side_exception[path_nr] == 0 ) ):
                        template_pos_left -= len(self.path_list[path_nr])

                    # Template index at break domain.
                    template_domain_at_break_left = domain - template_pos_left + self.template_offsets[path_nr][break_nr]
                    
                    # Check if template covers domain to be broken, from left.
                    if ( template_domain_at_break_left>= 0 ) and ( template_domain_at_break_left < len(self.template[ self.template_types[path_nr][break_nr] ]) ):
                        domain_length_left = self.template[ self.template_types[path_nr][break_nr] ][template_domain_at_break_left]
                    else:
                        # Set preferred domain length if not covered by template.
                        domain_length_left = self.standard_domain_length
                    
                    # Shift coordinates if wrap around.
                    template_pos_right = self.template_positions[path_nr][( break_nr + 1) % len(path_breaks)]
                    if ( template_pos_right < domain ) or ( self.single_strand_side_exception[path_nr] == 1 ):
                        template_pos_right += len(self.path_list[path_nr])

                    # Template index at break domain.
                    template_domain_at_break_right = domain - template_pos_right + self.template_offsets[path_nr][break_nr]
                    
                    # Check if template covers domain to be broken, from right.
                    if ( template_domain_at_break_right>= 0 ) and ( template_domain_at_break_right < len(self.template[ self.template_types[path_nr][break_nr] ]) ):
                        domain_length_right = self.template[ self.template_types[path_nr][break_nr] ][template_domain_at_break_right]
                    else:
                        # Set preferred domain length if not covered by template.
                        domain_length_right = self.standard_domain_length
                    
                    domain_length = self.path_list[path_nr][domain]
                    
                    # Both template domains don't overfill the domain.
                    # TODO Different behavior could be implemented as well.

                    if domain_length_left + domain_length_right <= domain_length:
                        break_base = domain_length_left - 1
                    else:
                        # shorter template domain is longer than half the domain, break in half
                        if min(domain_length_left,domain_length_right) * 2 > domain_length:
                            break_base = ( domain_length // 2 ) - 1
                        #left domain is shorter and shorter than half the domain, keep full
                        elif domain_length_left <= domain_length_right:
                            break_base = domain_length_left - 1
                        #right domain is shorter and shorter than half the domain, keep full
                        else:
                            break_base = domain_length - domain_length_right - 1
                    
                    # TODO this should not occur in standard cadnano designs, a one-base domain
                    if len(self.dna_structure.strands[strand_nr].domain_list[domain].base_list) < 2:
                        self._logger.warning("Domain %d number of domain bases < 2" % (domain))
                        #print "strand nr",strand_nr,"path nr",path_nr,"domain nr",domain,"has < 2 bases, cannot break" 
                    else:
                        helix = self.dna_structure.strands[strand_nr].domain_list[domain].helix.lattice_num
                        position1 = self.dna_structure.strands[strand_nr].domain_list[domain].base_list[ break_base ].p
                        position2 = self.dna_structure.strands[strand_nr].domain_list[domain].base_list[ break_base + 1 ].p
                        helix_index = self.dna_structure.structure_helices_map[ helix ].load_order
                        self.cadnano_design.helices[helix_index].staple_strands[position2].initial_strand=-1
                        self.cadnano_design.helices[helix_index].staple_strands[position2].initial_base=-1
                        self.cadnano_design.helices[helix_index].staple_strands[position1].final_strand=-1
                        self.cadnano_design.helices[helix_index].staple_strands[position1].final_base=-1

                # Crossover break.

                else:         
                    crossover = self.path_crossover_list[path_nr][domain][1]
                    helix_index_3 = self.dna_structure.structure_helices_map[ crossover[3] ].load_order
                    helix_index_0 = self.dna_structure.structure_helices_map[ crossover[0] ].load_order
                    self.cadnano_design.helices[ helix_index_3 ].staple_strands[crossover[4]].initial_strand=-1
                    self.cadnano_design.helices[ helix_index_3 ].staple_strands[crossover[4]].initial_base=-1
                    self.cadnano_design.helices[ helix_index_0 ].staple_strands[crossover[1]].final_strand=-1
                    self.cadnano_design.helices[ helix_index_0 ].staple_strands[crossover[1]].final_base=-1

                #__if break_current % 2 == 0:

            #__for break_nr,break_current in enumerate(path_breaks[0:len(path_breaks) - 1 + is_circular]):

        #__for path_nr, path_breaks in enumerate(self.break_positions)

    #__def _break_json

def main():
    
    converter = Converter()

    # TODO: I removed some abspath stuff here so that it just tries to expand the first argument, that way this can be called from different relative paths. Remove this comment or edit to match a consistent style for argument handling.
    if (len(sys.argv) != 2):
        sys.stderr.write("**** ERROR: Wrong number of arguments.\n") 
        sys.stderr.write("Usage: stapler.py <filename>\n")
        sys.stderr.write("Output: If <filename> is path/name.json, will be placed in path/name_recode.json")
        sys.exit(1)
    
    file_full_path_and_name = os.path.abspath( os.path.expanduser( sys.argv[1] ))
    file_name = os.path.basename( file_full_path_and_name )
    file_path = os.path.dirname ( file_full_path_and_name )

    output_file_name = file_name[:-5]+"_recode.json"
    output_file_full_path_and_name = os.path.join( file_path, output_file_name)
    
    cadnano_reader = CadnanoReader()
    converter.cadnano_design = cadnano_reader.read_json(file_full_path_and_name)

    dna_parameters = DnaParameters()
    converter.cadnano_convert_design = CadnanoConvertDesign(dna_parameters)
    converter.dna_structure = converter.cadnano_convert_design.create_structure(converter.cadnano_design)
    converter.dna_structure.get_domains()

    #TODO testing variables:

    stapler = Stapler(converter.dna_structure,converter.cadnano_design)
    stapler.template = [[7,7,14,7,7,7],[7,14,7,7,7],[14,7,7,7],[7,7,7,14,7,7],[7,7,7,14,7],[7,7,7,14],[7,7,14,7,7],[7,14,7,7],[14,7,7],[7,7,14,7],[7,7,14],[7,14,7]]
    #stapler.template = [[7,7,7,14,7,7,7]]
    #stapler.template = [[7,7,7,14,7],[7,7,14,7,7],[7,14,7,7,7]]
    #stapler.template = [[16,16,16],[8,16],[16,8]]
    #stapler.template = [[8,16,8]]
    
    stapler.uncovered_domain_penalty = 10.0
    stapler.template_overhang_penalty_factor = 0.5
    stapler.step_probabilities = (0.1,0.1,0.35,0.35,0.1)
    stapler.standard_domain_length = 5

    stapler.initialize_system()
    
    stapler.temperature = 10.0
    nr_steps = 10000000
    nr_steps_timescale = 10000
    stapler.generate(nr_steps, nr_steps_timescale, 0.999975)
    
    """ TODO
       structure with only two staples with two segments each, one double crossover, crashes:
       **** ERROR: Reached a visited base.Traceback (most recent call last):
       File "testScript3.py", line 357, in <module>
       main()
       File "testScript3.py", line 346, in main
       converter.dna_structure = converter.cadnano_convert_design.create_structure(converter.cadnano_design)
       File "/Users/jeanphilippe/Google Drive/tilibit other/Todo/nanodesign/nanodesign_transition/converters/cadnano/convert_design.py", line 170, in create_structure
       self._set_strands_colors(strands)
       File "/Users/jeanphilippe/Google Drive/tilibit other/Todo/nanodesign/nanodesign_transition/converters/cadnano/convert_design.py", line 1124, in _set_strands_colors
       for strand in strands:
       TypeError: 'NoneType' object is not iterable
    """
    
    converter.dna_structure = converter.cadnano_convert_design.create_structure(converter.cadnano_design)
    converter.dna_structure.get_domains()
    print "staple domains:"
    for strands in converter.dna_structure.strands[1:]:
        lengths = []
        for domain in strands.domain_list:
            lengths.append(len(domain.base_list))
        print lengths

    # Write a caDNAno JSON file.
    cadnano_writer = CadnanoWriter(converter.dna_structure)
    cadnano_writer.write(output_file_full_path_and_name)

if __name__ == '__main__':
    main()
