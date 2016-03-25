__all__ = ["Domain"]

class Domain(object):
    """ This class stores information for a DNA domain. """
    def __init__(self, id, helix):
        self.id = id
        self.helix = helix
        self.strand = None
        self.base_list = []
        self.color = [0.5,0.5,0.5]
        self.connected_strand = -1
        self.connected_domain = -1

    def get_end_points(self):
        base_pos = self.helix.helix_axis_nodes
        #print"######## base_pos %d" % len(base_pos)
        base1 = self.base_list[0].p 
        base2 = self.base_list[-1].p
        point1 = base_pos[base1]
        point2 = base_pos[base2]
        return point1,point2
