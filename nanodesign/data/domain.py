__all__ = ["Domain"]

from .energymodel import energy_model, convert_temperature_K_to_C

class Domain(object):
    """ This class stores information for a DNA domain. """
    def __init__(self, id, helix, strand, bases):
        self.id = id
        self.helix = helix
        self.strand = strand
        self.base_list = bases
        self._color = None
        
        self.sequence = "".join([b.seq for b in self.base_list])

        self.connected_strand = -1
        self.connected_domain = -1


    @property
    def color( self ):
        """ The color property is related to what Cadnano thinks the color should be in
            the input file, or can be arbitrarily assigned. There is an open question as to
            whether this should be directly part of the Domain object or should be metadata
            connected to it."""
        if self._color is not None:
            return self._color
        elif self.strand is not None:
            return self.strand.color
        else:
            return [0.5,0.5,0.5]

    @color.setter
    def color( self, newval ):
        self._color = newval
        # should we also set the strand's color? Or can the domain have a different color from the strand?

    def get_end_points(self):
        """ Get the end coordinates of the domain. """
        base_pos = self.helix.helix_axis_coords
        base1 = self.base_list[0]
        base2 = self.base_list[-1]
        point1 = base1.coordinates
        point2 = base2.coordinates
        #print("[get_end_points] base 1 p %d  h %d  coord %g %g %g" % (base1.p, base1.h, point1[0], point1[1], point1[2]))
        #print("                 base 2 p %d  h %d  coord %g %g %g" % (base2.p, base2.h, point2[0], point2[1], point2[2]))
        return point1,point2

    def melting_temperature(self):
        # If we are not paired, return a nonphysical melting temperature.
        if self.connected_domain == -1: 
            return -500.0

        # If we have any "N" bases in our sequence, also return a nonphysical melting temperature:
        if "N" in self.sequence or 'n' in self.sequence:
            return -501.0

        def rev_complement( seq ):
            comp = {"G":"C","C":"G",
                    "A":"T","T":"A", 
                    "g":"c","c":"g",
                    "a":"t","t":"a",
                    "N":"N","n":"n"}
            rev = []
            for i in xrange(len(seq)):
                rev.append(comp[seq[len(seq)-i-1]])
            return "".join( rev )

        
        _,_,dH,dS = energy_model.stack_energy( self.sequence, rev_complement( self.sequence ))
        return convert_temperature_K_to_C( energy_model.melting_temperature( dH, dS ) )
        
