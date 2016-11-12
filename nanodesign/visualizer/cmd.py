#!/usr/bin/env python
""" This module is used to read, execute and generate commands for visualization. 

    Visualization commands are used to perform all visualization operations. Commands are generated
    by the popup menu and executed here. 

    Commands are space-separated name-value pairs with the format:

        <entity>  <name1>=<value1>  <name2>=<value2>  ... <nameN>=<valueN>

            where <entity>  = domain | graphics | helix | strand 

    Operations performed using the popup menu are written to the file 'vis.cmd'. These commands may be read in at the
    start of the visualization session or included on the command line.
"""
import logging
import os
import re
from .menu import VisMenuEntity
import model

class VisCommandEntity:
    """ Defines the types of commands recognized by the visualizer. """
    ATOMIC_STRUCTURE = 'atomic_structure' 
    COMMENT  = '#'
    DOMAIN   = 'domain'
    GRAPHICS = 'graphics'
    HELIX    = 'helix'
    MODEL    = 'model'
    STRAND   = 'strand'

class VisCommand(object):
    """ This class implements command processing. """
    NONE = 'None'
    ALL = 'All'
    COMMAND_DELIM = ';'
    SHOW_ATOMIC_STRUCTURE = 'atomic_structure name=%s  rep=%s  show=%s' 
    SHOW_HELIX = 'helix name=%s  rep=%s  show=%s' 
    SHOW_STRAND = 'strand name=%s  rep=%s  show=%s' 
    SHOW_MODEL_REP  = 'model rep=%s show=%s' 

    def __init__(self, model, cmd_file_name, commands):
        """ 
            Initialize VisCommand object.

            Arguments:
                model (VisModel): The visualization model object used to interface with the DNA design structure
                    and manage the visualization of all representions.
                cmd_file_name (String): The name of the file to read commands from.
                commands(string): The string of commands given on the command line.
        """
        self.model = model
        #self.update_menu = False
        self.update_menu = True
        self._logger = logging.getLogger(__name__)
        # Open command logging file.
        path = os.getcwd()
        self.file_name = os.path.join(path, "vis.cmd")
        self.file = open(self.file_name, 'w')
        self.cmd_file_name = cmd_file_name
        self.commands = commands
        # The dict used to map command entities to their parsing functions.
        self.parsing_map = {  
            VisCommandEntity.ATOMIC_STRUCTURE : self.proc_atomic_struct_cmd, 
            VisCommandEntity.COMMENT : self.proc_comment, 
            VisCommandEntity.GRAPHICS : self.proc_graphics_cmd, 
            VisCommandEntity.HELIX : self.proc_helix_cmd, 
            VisCommandEntity.MODEL : self.proc_model_cmd, 
            VisCommandEntity.STRAND : self.proc_strand_cmd 
        }

    def write_cmd(self, cmd):
        """ Write a command to a file. """
        self.file.write(cmd + "\n")

    def execute_file_cmds(self):
        """ Execute the commands in a file. """
        if not self.cmd_file_name:
            return
        self._logger.info("Reading commands from \'%s\' " % self.cmd_file_name)
        self.update_menu = True
        with open(self.cmd_file_name) as cmd_file:
            for line in cmd_file:
                self._logger.info("Command \'%s\'" % line)
                cmd = line.rstrip()
                self.execute_cmd(cmd)
                self.write_cmd(cmd)
            #__for line in cmd_file
        #__with open(self.cmd_file_name) as cmd_file

    def execute_cmds(self, commands=None):
        """ Execute the commands from a string. """
        if not (self.commands or commands):
            return
        if self.commands:
            self._logger.info("Execute commands from command line.")
            cmd_list = self.commands.split(VisCommand.COMMAND_DELIM)
        else:
            cmd_list = commands.split(VisCommand.COMMAND_DELIM)
        self.update_menu = True
        for cmd in cmd_list:
            self._logger.info("Command \'%s\'" % cmd)
            self.execute_cmd(cmd)
            self.write_cmd(cmd)

    def execute_cmd(self, cmd):
        """ Execute a command. """
        tokens = cmd.split()
        if not tokens:
            return
        entity = tokens[0]
        if entity in self.parsing_map:
            self.parsing_map[entity](cmd, tokens[1:])
        else:
            self._logger.error("Unknown entity \'%s\' " % entity) 

    def proc_comment(self, cmd, tokens):
        pass

    def proc_helix_cmd(self, cmd, tokens):
        """ Process a 'helix' command. """
        args = [ 'color', 'name', 'names', 'rep', 'show' ]
        name = None
        rep = None
        attributes = []
        names = []
        for token in tokens:
            arg_name,arg_value = token.split('=')
            if arg_name not in args:
                self._logger.error("Unknown helix argument \'%s\' " % arg_name) 
                return
            if arg_name == 'name':
                names.append(arg_value)
            elif arg_name == 'names':
                names = self.parse_helix_ids(arg_value)
            elif arg_name == 'rep':
               rep = arg_value
            elif arg_name == 'show':
               show = (arg_value == 'true')
               attributes.append(('show',show))
            elif arg_name == 'color':
               color = self.parse_color(arg_value)
               attributes.append(('color',color))
        #__for token in tokens
        for name in names:
            self.model.show_helix(name, rep, attributes)
        # Update the menu but after graphics is up and fully initialized (delayed=True).
        if (show and self.update_menu):
            delay = True
            self.model.menu.update_submenu_selection(VisMenuEntity.HELIX, name, rep, delay)

    def generate_helix_cmd(self, helix_name, helix_rep, show):
        """ Generate a helix command. """
        if helix_name == VisCommand.NONE: 
            cmd = VisCommand.SHOW_HELIX % ('All', helix_rep, 'false')
        elif helix_name == VisCommand.ALL: 
            cmd = VisCommand.SHOW_HELIX % (helix_name, helix_rep, 'true')
        else:
            if show: 
                show_str = 'true' 
            else:
                show_str = 'false' 
            cmd = VisCommand.SHOW_HELIX % (helix_name, helix_rep, show_str)
        self.write_cmd(cmd)
        self.execute_cmd(cmd)

    def generate_model_cmd(self, rep, show):
        """ Generate a model command. """
        SHOW_MODEL_REP  = 'model rep=%s show=%s' 
        if show: 
            show_str = 'true' 
        else:
            show_str = 'false' 
        rep = rep.replace(" ", "_")
        cmd = VisCommand.SHOW_MODEL_REP % (rep, show_str)
        self.write_cmd(cmd)
        self.execute_cmd(cmd)

    def generate_strand_cmd(self, strand_name, strand_rep, show):
        """ Generate a strand command. """
        if strand_name == VisCommand.NONE:
            cmd = VisCommand.SHOW_STRAND % ('All', strand_rep, 'false')
        elif strand_name == VisCommand.ALL:
            cmd = VisCommand.SHOW_STRAND % (strand_name, strand_rep, 'true')
        else:
            if show: 
                show_str = 'true' 
            else:
                show_str = 'false' 
            cmd = VisCommand.SHOW_STRAND % (strand_name, strand_rep, show_str)
        self.write_cmd(cmd)
        self.execute_cmd(cmd)

    def generate_atomic_struct_cmd(self, atom_struct_name, atom_struct_rep, show):
        """ Generate an atomic struct command. """
        if atom_struct_name == VisCommand.NONE:
            cmd = VisCommand.SHOW_ATOMIC_STRUCTURE % ('All', atom_struct_rep, 'false')
        elif atom_struct_name == VisCommand.ALL:
            cmd = VisCommand.SHOW_ATOMIC_STRUCTURE % (atom_struct_name, atom_struct_rep, 'true')
        else:
            if show: 
                show_str = 'true' 
            else:
                show_str = 'false' 
            cmd = VisCommand.SHOW_ATOMIC_STRUCTURE % (atom_struct_name, atom_struct_rep, show_str)
        self.write_cmd(cmd)
        self.execute_cmd(cmd)

    def proc_model_cmd(self, cmd, tokens):
        """ Process a 'model' command. """
        args = [ 'rep', 'show' ]
        bbox = None
        for token in tokens:
            arg_name,arg_value = token.split('=')
            if arg_name not in args:
                self._logger.error("Unknown model argument \'%s\' " % arg_name)
                return
            if arg_name == 'rep':
               rep = arg_value
            elif arg_name == 'show':
               show = (arg_value == 'true')
        #__for token in tokens
        rep = rep.replace("_", " ")
        if rep == model.VisModelRepType.BOUNDING_BOX:
            self.model.show_bounding_box(show)
        elif rep == model.VisModelRepType.GEOMETRY:
            self.model.show_structure_geometry(show)
        elif rep == model.VisModelRepType.HELIX_NUMBERS:
            self.model.show_helix_numbers(show)
        elif rep == model.VisModelRepType.HELIX_PROJECTION:
            self.model.show_helix_projections(show)
        else:
            self._logger.error("Unknown model rep \'%s\' " % rep)
            return

        if (show and self.update_menu):
            delay = True
            self.model.menu.update_selection(VisMenuEntity.MODEL, rep, delay)

    def proc_strand_cmd(self, cmd, tokens):
        """ Process a 'strand' command. """
        args = [ 'color', 'line_width', 'name', 'names', 'rep', 'show' ]
        name = None
        rep = None
        attributes = []
        names = []
        for token in tokens:
            arg_name,arg_value = token.split('=')
            if arg_name not in args:
                self._logger.error("Unknown strand argument \'%s\' " % arg_name)
                return
            if arg_name == 'name':
                names.append(arg_value)
            elif arg_name == 'names':
                names = self.parse_strand_ids(arg_value)
            elif arg_name == 'rep':
               rep = arg_value
            elif arg_name == 'show':
               show = (arg_value == 'true')
               attributes.append(('show',show))
            elif arg_name == 'color':
               color = self.parse_color(arg_value)
               attributes.append(('color',color))
            elif arg_name == 'line_width':
               line_width = float(arg_value)
               attributes.append(('line_width',line_width))
        #__for token in tokens
        for name in names:
            self.model.show_strand(name, rep, attributes)
        # Update the menu but after graphics is up and fully initialized (delayed=True).
        if (show and self.update_menu):
            delay = True
            self.model.menu.update_submenu_selection(VisMenuEntity.STRAND, name, rep, delay)

    def proc_atomic_struct_cmd(self, cmd, tokens):
        """ Process a 'atomic_structure' command. """
        args = [ 'name', 'rep', 'show' ]
        name = None
        rep = None
        show = True
        for token in tokens:
            arg_name,arg_value = token.split('=')
            if arg_name not in args:
                self._logger.error("Unknown atomic structure argument \'%s\' " % arg_name)
                return
            if arg_name == 'name':
               name = arg_value
            elif arg_name == 'rep':
               rep = arg_value
            elif arg_name == 'show':
               show = (arg_value == 'true')
        #__for token in tokens
        self.model.show_atomic_struct(name, rep, show)
        # Update the menu but after graphics is up and fully initialized (delayed=True).
        if (show and self.update_menu):
            delay = True
            self.model.menu.update_submenu_selection(VisMenuEntity.ATOMIC_STRUCTURE, name, rep, delay)
    #__def proc_atomic_struct_cmd(self, tokens)

    def generate_graphics_cmd(self, name, value):
        """ Generate a graphic command. """
        if name == "center":
            cmd = "graphics  center=(%g,%g,%g)" % (value[0],value[1],value[2])
        self.write_cmd(cmd)

    def proc_graphics_cmd(self, cmd, tokens):
        """ Process a 'graphics' command. """
        args = [ 'center' ]
        #print("proc_graphics_cmd tokens %s" % str(tokens))
        #s = ''.join(tokens)
        #print("proc_graphics_cmd s %s" % s)
        for token in tokens:
            arg_name,arg_value = token.split('=')
            if arg_name not in args:
                self._logger.error("Unknown graphics argument \'%s\' " % arg_name)
                return
            if arg_name == 'center':
                #print("proc_graphics_cmd value %s" % arg_value)
                values = arg_value[arg_value.find("(")+1:arg_value.find(")")].split(",")
                #print("proc_graphics_cmd values %s" % values)
                point = [float(x) for x in values]
                self.model.graphics.set_center(point)
        #__for token in tokens
    #__def proc_graphics_cmd(self, cmd, tokens)

    def parse_color(self, value):
        pattern = re.compile(r"[,()]")
        color_tokens = pattern.split(value)
        return [ float(color) for color in color_tokens if color != '']
    #__def parse_color(self, value)

    def parse_helix_ids(self, value):
        """ Parse a list of helix IDs. """
        pattern = re.compile(r"[,\[\]]")
        tokens = pattern.split(value)
        helix_ids = []
        for s in tokens:
            if "-" in s:
                ids = self.parse_id_range(s)
                for id in ids:
                    helix_ids.append(str(id))
            elif s:
                helix_ids.append(s)
        #__for s in tokens
        return helix_ids
    #__def parse_helix_ids

    def parse_strand_ids(self, value):
        """ Parse a list of strand IDs. """
        pattern = re.compile(r"[,\[\]]")
        tokens = pattern.split(value)
        strand_ids = []

        if (tokens[0] == "start_helices") or (tokens[0] == "in_helices"):
            helix_ids = set()
            for s in tokens[1:]:
                if "-" in s:
                    ids = self.parse_id_range(s)
                    for id in ids: 
                        helix_ids.add(id)
                elif s != '':
                    helix_ids.add(int(s))
            #__for s in tokens[1:]

            if (tokens[0] == "start_helices"):
                for strand in self.model.strands.values():
                    if strand.dna_strand.tour[0].h in helix_ids:
                        strand_ids.append(strand.name)
                #__for strand in self.model.strands

            elif (tokens[0] == "in_helices"):
                for strand in self.model.strands.values():
                    if not strand.dna_strand.is_scaffold:
                        for base in strand.dna_strand.tour:
                            if base.h in helix_ids:
                                strand_ids.append(strand.name)
                                break
                #__for strand in self.model.strands.values()

            #__if (tokens[0] == "start_helices")

        #__if tokens[0] == "start_helices or ..."

        else:
            for s in tokens:
                strand_ids.append(s)
        #__if tokens[0] == "start_helices"
        return strand_ids
    #__def parse_strand_ids

    def parse_id_range(self, value):
        """ Parse a range of integer IDs of the form '<start> - <end>'. """
        rtoks = value.split("-")
        start = int(rtoks[0])
        end = int(rtoks[1])+1
        ids = []
        for id in xrange(start,end):
            ids.append(id)
        return ids
    #__def parse_id_range

#__class VisCommand(object)

