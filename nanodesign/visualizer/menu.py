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

"""This module is used to create a GLUT popup menu. 

   The main menu contains items for different graphics representions of certain
   entities from a DNA structure design: virtual helices, strands, domains,
   etc. A submenu of entity names is defined for each representation.  The
   entity names are used to selectively visualize components of the DNA
   structure design. The submenus for virtual helix representions contain
   integer IDs for the helix 'number' obtained from the caDNAno design file.

   Most submenus contain 'All' and 'None' entries that are used to selected all
   or none of the entities from the submenu. When an entity from a submenu is
   selected its visibility is modified depending on its current selection
   state. If it is already selected then it is set to be hidden. The entity
   name in the submenu is updated with a '+' if it not selected.

   Menu selection generates a command that is then executed to perform the visualization operation.
"""
import inspect
import logging
from .helix import VisHelixRepType
from .strand import VisStrandRepType
from .atomic_struct import VisAtomicStructureRepType
import model

try:
    from OpenGL.GL import *
    from OpenGL.GLU import *
    from OpenGL.GLUT import *

except ImportError as e:
    print "Could not import PyOpenGL."
    raise e

class VisMenuItem:
    """ This class defines the two menu items that appear in all submenus. """
    ALL  = "All"
    NONE = "None"

class VisMenuEntity:
    """ This class defines the entities that can be visualized. """ 
    UNKNOWN = 'unknown'
    ATOMIC_STRUCTURE = 'atomic_structure'
    HELIX    = 'helix'
    MODEL    = 'model'
    STRAND   = 'strand'

class VisMenuItemSelect(object):
    """ This class is used to process the selection of submenu item for a given representation from submenus. 

        Attributes:
            selection (VisRepMenu): The menu selection object that this representation is part of. 
            rep (String): The representation this menu item selects. The string is from the values of representation 
                type classes (e.g. VisHelixRepType).

        When a submenu item it selected the callback method in this class is called with the item number selected.
        The item number and the representation is then used to process the selection of an entiy of that representation.
    """ 
    def __init__(self, selection, rep):
        self.selection = selection
        self.rep = rep

    def callback(self, item):
        """ This function is the callback given to a GLUT menu. 

            Arguments:
                item (int): The menu item number of the selection. 

            This function is executed when an item from a menu is selected. 
        """
        self.selection.select_rep(item, self.rep)
        return 0

class VisMenuMainEntry(object):
    """ This class is used to manage the main menu selections. 
    """
    def __init__(self, name, menu_id, callback):
        self.name = name
        self.menu_id = menu_id
        self.callback = callback
        self.selected = False

class VisSubMenu(object):
    """ This class is used to manage a submenu for a set of representations.

        Attributes:
            menu (VisMenu): The menu object used to manage the visualizer menus.
            command_generate (Method): The function used to generate commands for this representation.
            rep_names (List[String]): The list of representation names displayed in the submenu.
            selections (Dict[bool]): The dictionary of boolean values that determine if a representation is selected. 
                It is indexed by representation name.
            menuID (int): The menu ID created by GLUT when creating a menu using glutCreateMenu.  he menu ID is used to 
                set the current menu when updating submenu names when selected.

        The menu selections for a set of representation are used to selectively visualize representations without 
        entities (e.g. helix, strand, etc.). A submenu is created for representations that can be selected. A boolean 
        is stored for all the representations. If the boolean is True the the representation is seleced and is visible.
    """
    def __init__(self, menu, menu_name, rep, command_generate):
        """ Initialize a VisSubMenu object.

            Arguments:
                menu (VisMenu): The menu object used to manage the visualizer menus.
                menu_name (String): The name of the submenu. 
                rep (Class): The class that defines an entity's representation types (e.g. VisModelRepType).
                command_generate (Method): The function used to generate commands for this representation.

            The type names for a given representaion are obtained from its class attibutes. A dict 
            stores a boolean value for each representaion. A submenu is created with an item for each
            representaion. A VisMenuItemSelect object is used to handle the callbacks from the submenu.
        """
        self.menu = menu
        self.menu_name = menu_name
        self.command_generate = command_generate
        # Get the representation type names.
        attributes = inspect.getmembers(rep, lambda a:not(inspect.isroutine(a)))
        rep_types = [a for a in attributes if not(a[0].startswith('__') and a[0].endswith('__'))]
        # Create the selections dict, set callbacks and create submenus for each rep.
        self.selections = {}
        self.menuIDs = {}
        self.rep_names = []
        for rep_type in rep_types:
            rep = rep_type[1]
            if rep == VisMenuEntity.UNKNOWN:
                continue
            self.selections[rep] = False
            self.rep_names.append(rep)
        #__for rep_type in rep_types
        menu_item_select = VisMenuItemSelect(self, menu_name)
        self.menuID = self.menu._add_submenu_item(menu_item_select.callback, self.rep_names)

    def add_submenu(self):
        glutAddSubMenu(self.menu_name, self.menuID)

    def select_rep(self, item, rep):
        """ Select a representation.

            Arguments:
                rep (String): The name of the representation selected.
                item (int): The menu item number representation selected.

            Selecting a representation will cause it to be displayed or hidden from view. A command is generated
            and then executed to perform the visualization of the representation.
        """
        selected_rep = self.rep_names[item]
        show = self.menu.set_selections(self.rep_names, self.selections, item, rep, self.menuID)
        self.command_generate(selected_rep, show)

    def update_selection(self, name, rep):
        """ Update a selection for the given representation. This is used to update a menu from a command.

            Arguments:
                name (String): Not used. 
                rep (String): The name of the representation.
        """
        for item,rep_name in enumerate(self.rep_names):
            if rep_name == rep:
                show = self.menu.set_selections(self.rep_names, self.selections, item, rep, self.menuID)
                break
        #__for item,rep_name in enumerate(self.rep_names)


class VisRepMenu(object):
    """ This class is used to manage the submenus for a set of representations. 

        Attributes:
            menu (VisMenu): The menu object used to manage the visualizer menus. 
            names (List[String]): The list of entity names displayed for the representation's submenu. 
            command_generate (Method): The function used to generate commands for this representation.
            selections (Dict[Dict[bool]]): The dictionary of boolean values that determine if an entity of
                representation is selected. It is indexed by representation name and entity name. 
            menuIDs (Dict[int]): The dictionary of the menu IDs created by GLUT when creating a menu using glutCreateMenu.
                The menu IDs are used to set the current menu when updating submenu names when selected.

        The menu selections for a given representation are used to selectively visualize the individual
        entities (helix, strand, etc.) of a model. Submenus are create for each representation containing
        entity names (e.g. helix number or strand name) that can be selected. A boolean is stored for all the entities
        defined for each representation. If the boolean is True the the entity is seleced and is visible. 
    """
    def __init__(self, menu, rep, entity_names, command_generate):
        """ Initialize a VisRepMenu object.

            Arguments:
                menu (VisMenu): The menu object used to manage the visualizer menus. 
                entity_names (List[String]): The list of entity names displayed for the representation's submenu. 
                rep (Class): The class that defines an entity's representation types (e.g. VisHelixRepType).
                command_generate (Method): The function used to generate commands for this representation.

            The type names for a given representaion are obtained from its class attibutes. A dict of dicts
            stores a boolean value for each entity of each representaion. A submenu is created for each
            representaion and uses VisMenuItemSelect object to handle the callbacks from the submenu.
        """
        self.menu = menu
        self.entity_names = entity_names
        self.command_generate = command_generate
        # Get the representation type names. 
        attributes = inspect.getmembers(rep, lambda a:not(inspect.isroutine(a)))
        rep_types = [a for a in attributes if not(a[0].startswith('__') and a[0].endswith('__'))]
        # Create the selections dict, set callbacks and create submenus for each rep.
        self.selections = {}
        self.menuIDs = {}
        for rep_type in rep_types:
            rep = rep_type[1]
            if rep == VisMenuEntity.UNKNOWN:
                continue
            self.selections[rep] = {}
            # Set the selection flag for each entity name.
            for name in entity_names:
                self.selections[rep][name] = False
            menu_item_select = VisMenuItemSelect(self, rep)
            self.menuIDs[rep] = self.menu._add_submenu_item(menu_item_select.callback, entity_names)
        #__for rep_type in rep_types

    def add_submenu(self, rep, menu_name):
        """ Create a GLUT representation submenu. 

            Arguments:
                rep (String): The name of a representation, taken from a representation class such as VisHelixRepType.
                menu_name (String): The menu name that is displayed in the popup menu. 
        """
        glutAddSubMenu(menu_name, self.menuIDs[rep])

    def select_rep(self, item, rep):
        """ Select a representation. 

            Arguments:
                rep (String): The name of the representation selected.
                item (int): The menu item number representation selected. 

            Selecting a representation will cause it to be displayed or hidden from view. A command is generated
            and then executed to perform the visualization of the representation.
        """
        selected_name = self.entity_names[item]
        show = self.menu.set_selections(self.entity_names, self.selections[rep], item, rep, self.menuIDs[rep])
        self.command_generate(selected_name, rep, show)

    def update_selection(self, name, rep):
        """ Update a selection for the given representation and name. This is used to update a menu from a command.

            Arguments:
                name (String): The entity name.
                rep (String): The name of the representation.
        """
        for item,entity_name in enumerate(self.entity_names):
            if entity_name == name:
                show = self.menu.set_selections(self.entity_names, self.selections[rep], item, rep, self.menuIDs[rep])
                break
        #__for item,item_name in enumerate(self.entity_names)

class VisMenu(object):
    """ This class is used to manage the visualizer menus. 

        Attributes:
            atomic_struct_names (List[String]): The list of atomic structure names to be displayed for the atomic structure 
                representation submenus. 
            delayed_updates (List[Tuple]): A list of menu selection updates from commands that are applied after the 
                graphics loop has started. We neeed to do this because there seemed to be some menu state issues.
            helix_names (List[String]): The list of helix entity names to be displayed for the helix representation submenus. 
            item_count (int): The count of main menu items. This is used to index menu_items[].
            menu_items (List[int]): The list of main (non-submenu) menu items. 
            selections (Dict[VisRepMenu]): The dict of selection objects used to manage menu selections.
            strand_names (List[String]): The list of strand names to be displayed for the strand representation submenus. 
            updated (bool): If True then the menu has been updated by commands.

        The visualizer popup menu displays menus for different graphics representions of certain entities from a DNA 
        structure design: virtual helices, strands, domains, etc. A submenu of entity names is defined for each representation. 
        The entity names are used to selectively visualize components of the DNA structure design. Entity selections are
        managed using a the 'selections' attribute, a dictionary of VisRepMenu objects.
    """
    menu_items = {} 
    item_count = 1
    updated = False
    delayed_updates = []
    delayed_submenu_updates = []
    selected_symbol = " - on"

    def __init__(self, command, helix_names, strand_names, atomic_struct_names):
        self.command = command
        self.helix_names = helix_names 
        self.strand_names = strand_names
        self.atomic_struct_names = atomic_struct_names 
        self.selections = {}
        self._logger = logging.getLogger(__name__)
        self._create_menu_items()

    def _create_menu_items(self):
        """ Create items for the popup menu. 

            To create the menu items we first create items for the submenus for each representation and with the 
            associated entity names. The main menu items are then created using titles for each representation. 
        """

        # Create model sub-menu.
        model_select = VisSubMenu(self, "Model", model.VisModelRepType, self.command.generate_model_cmd)
        self.selections[VisMenuEntity.MODEL] = model_select

        # Create helix sub-menus.
        helix_select = VisRepMenu(self, VisHelixRepType, self.helix_names, self.command.generate_helix_cmd)
        self.selections[VisMenuEntity.HELIX] = helix_select

        # Create strand sub-menus.
        strand_select = VisRepMenu(self, VisStrandRepType, self.strand_names, self.command.generate_strand_cmd)
        self.selections[VisMenuEntity.STRAND] = strand_select

        # Create atomic stucture sub-menus.
        if len(self.atomic_struct_names) > 2:
            atomic_select = VisRepMenu(self, VisAtomicStructureRepType, self.atomic_struct_names, 
                self.command.generate_atomic_struct_cmd)
            self.selections[VisMenuEntity.ATOMIC_STRUCTURE] = atomic_select
        else:
            self.selections[VisMenuEntity.ATOMIC_STRUCTURE] = None 
            atomic_select = None 

        # Create main menu items. Most of the main menu items will be a submenu of entity names.
        self.menu = glutCreateMenu(self.main_callback)

        # Create the submenus for atomic structure representations. 
        if atomic_select: 
            atomic_select.add_submenu(VisAtomicStructureRepType.BACKBONE, "Atomic structure backbone") 
            atomic_select.add_submenu(VisAtomicStructureRepType.BONDS,    "Atomic structure bonds") 
            atomic_select.add_submenu(VisAtomicStructureRepType.CHECK,    "Atomic structure check") 

        # Create the submenu for the model menu. This only has model representations for its submenu.
        model_select.add_submenu()

        # Create the submenus for strand representations. 
        strand_select.add_submenu(VisStrandRepType.CONNECTORS,  "Strand connectors")
        strand_select.add_submenu(VisStrandRepType.DOMAINS,     "Strand domains")
        strand_select.add_submenu(VisStrandRepType.FRAMES,      "Strand frames")
        strand_select.add_submenu(VisStrandRepType.PATH,        "Strand path")
        strand_select.add_submenu(VisStrandRepType.TEMPERATURE, "Strand temperature")

        # Create the submenus for the virtual helix representations. 
        helix_select.add_submenu(VisHelixRepType.BASE_POSITIONS,    "Virtual helix base positions") 
        helix_select.add_submenu(VisHelixRepType.COORDINATE_FRAMES, "Virtual helix coordinate frames") 
        helix_select.add_submenu(VisHelixRepType.DESIGN_CROSSOVERS, "Virtual helix design crossovers")
        helix_select.add_submenu(VisHelixRepType.COORDINATES,       "Virtual helix DNA helix P coordinates")
        helix_select.add_submenu(VisHelixRepType.DOMAINS,           "Virtual helix domains")
        helix_select.add_submenu(VisHelixRepType.GEOMETRY,          "Virtual helix geometry")
        helix_select.add_submenu(VisHelixRepType.INSERTS_DELETES,   "Virtual helix inserts and deletes")
        helix_select.add_submenu(VisHelixRepType.MAXIMAL_CROSSOVERS,"Virtual helix maximal crossovers")
        helix_select.add_submenu(VisHelixRepType.PAIRED_GEOMETRY,   "Virtual helix paired geometry")
        helix_select.add_submenu(VisHelixRepType.STRANDS,           "Virtual helix strands")
        helix_select.add_submenu(VisHelixRepType.TEMPERATURE,       "Virtual helix temperature")

        # Add the quit menu item.
        self._add_menu_item("Quit", self.quit_callback) 
        glutAttachMenu(GLUT_RIGHT_BUTTON)

    @staticmethod
    def help():
        """ Display the visualizer menu help. """
        print(" ====================== Menu =========================")
        print(" Atomic structure backbone - Show the atomic structure backbone P atoms for the selected strand.")
        print(" Atomic structure bonds - Show the atomic structure atom bonds as lines and DNA base planes as filled polygons.")
        print(" Atomic structure check - Check the bond lengths between backbone P atoms. Bond lengths that are too short or long are highlighted.")

        print(" Model - Bounding box - Show the bounding box of the design.")
        print("         Geometry - Show the geometry of the design representing dsDNA as solid cylinders.")
        print("         Virtual helix numbers - Show the cadnano virtual helix numbering displayed on in a plane on the bounding box.")
        print("         Virtual helix projection - Show the virtual helix axes projected onto two planes on the bounding box.")
        print(" Strand domains - Show the domains defined for a strand.")
        print(" Strand frames - Show the coordinate frames for a strand. An arrow shows the direction of a strand within a virtual helix.")
        print(" Strand geometry - Show a strand as a series of straight lines. A sphere shows the start of the strand.")
        print(" Strand temperature - Show the melting temperature for domains defined for a strand.")

        print(" Virtual helix base positions - Show the location of base positions along a helix axis.")
        print(" Virtual helix coordinate frames - Show the virtual helix coordinate frames. An arrow on the helix axis shows its polarity.")
        print(" Virtual helix design crossovers - Show the virtual helix crossovers present in the design.")
        print(" Virtual helix geometry - Show the virtual helix as a transparent cylinder.")
        print(" Virtual helix DNA helix P coordinates - Show the approximate coordinates of P atoms of a DNA helix.")
        print(" Virtual helix domains - Show the domains defined for the virtual helix.")
        print(" Virtual helix inserts and deletes - Show the locations of base insertions and deletions for the virtual helix.")
        print(" Virtual helix strands - Show the strands that start in a virtual helix.")
        print(" Virtual helix temperature - Show the melting temperature for domains defined for the virtual helix.")
        print("\n")

    def update(self):
        """ Update the menu with selections given from commands. 
          
            This is needed because the menus must be defined and GLUT finished initializing
            before the menus can be updated with selections from commands.
        """
        if self.updated:
            return
        self.updated = True
        for update in self.delayed_submenu_updates:
            method = update[0]
            type = update[1]
            name = update[2]
            rep = update[3]
            method(type, name, rep)

        for update in self.delayed_updates:
            method = update[0]
            type = update[1]
            name = update[2]
            rep = update[3]
            method(type, rep)

    def _add_menu_item(self, name, callback):
        """ Add a main menu item that does not have a submenu. """
        glutAddMenuEntry(name, VisMenu.item_count)
        entry = VisMenuMainEntry(name, VisMenu.item_count, callback)
        VisMenu.menu_items[VisMenu.item_count] = entry
        VisMenu.item_count += 1 

    def _add_submenu_item(self, callback, name_list):
        """ Add a submenu item. 

            Arguments:
                callback (Function): The callback function executed when an item is seleced from the submenu.
                name_list(List[String]): The list of entity names for the submenu. 
        """
        sub_menu = glutCreateMenu(callback)
        n = 0
        for name in name_list:
            glutAddMenuEntry(str(name),  n)
            n += 1
        VisMenu.item_count += 1 
        return sub_menu 

    def update_selection(self, type, rep, delay=False):
        """ Update a selection for the given rep name. """
        if delay:
            self.delayed_updates.append((self.update_selection, type, None, rep))
            return
        self.selections[type].update_selection(type,rep)

    def update_submenu_selection(self, type, name, rep, delay=False):
        """ Update a selection for the given rep and entity name. 

            Arguments:
                type (String): The menu entity type taken from VisMenuEntity.
                name (String): The entity name to update.
                rep (String): The representation name.
                delay (bool): If False then delay the menu update until the graphics has been initialized.
        """
        if delay:
            self.delayed_submenu_updates.append((self.update_submenu_selection, type, name, rep))
            return
        self.selections[type].update_selection(name,rep)

    def set_selections(self, names, selections, item, rep, menu=-1):
        """ Set the selection of an entity submenu item.

            Arguments:
                names (List[String]): The list of submenu entity names.            
                selections (Dict[String]): A dict of boolean values for each entity name. 
                item (int): The number of the submenu item selected. 
                rep (String): The name of the representation. 
                menu (int): The submenu ID. If -1 then don't set the menu.

            Returns a boolean indicating the entity's visibility. 

            When an entity from a submenu is selected its visibility is modified depending on its current selection
            state. If it is already selected then it is set to be hidden. The entity name in the submenu is updated 
            with a '+' if it not selected. 
        """
        name = names[item]
        cmenu = glutGetMenu()
        if menu != -1:
            glutSetMenu(menu)
        # Select all of the submenu entities.
        if name == VisMenuItem.ALL: 
            for item,name in enumerate(names):
                if item < 2:
                    continue
                entry_name = name + VisMenu.selected_symbol
                selections[name] = True
                glutChangeToMenuEntry(item+1, entry_name, item)
        # Select none  of the submenu entities.
        elif name == VisMenuItem.NONE: 
            for item,name in enumerate(names):
                if item < 2:
                    continue
                selections[name] = False
                glutChangeToMenuEntry(item+1, name, item)
        # Select an entity from the submenu entities.
        else:
            if not selections[name]:
                entry_name = name + VisMenu.selected_symbol
                selections[name] = True
            else:
                entry_name = name
                selections[name] = False
            glutChangeToMenuEntry(item+1, entry_name, item)
        return selections[name]
 
    def main_callback(self, item):
        """ This is the callback function for menus that are not submenus. """
        entry = VisMenu.menu_items[item]
        VisMenu.menu_items[item].callback(entry)
        return 0

    def quit_callback(self, entry):
        """ This the callback function for the 'quit' menu selection. """
        sys.exit(0)


