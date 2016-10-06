#####################
#CLASSES:
#####################

class Model(object):
    """ Represents a model """    
    def __init__(self,
                 domain = set(),
                 relation = dict(),
                 atoms = set(),
                 valuation = dict()
                 ): 
        self.domain = domain
        self.relation = relation
        self.atoms = atoms
        self.valuation = valuation

    def set_domain(self, domain):
        assert isinstance(domain, set)
        for codon in domain:
            assert isinstance(codon, str)
        self.domain = domain
    def get_domain(self):
        return self.domain

    def set_relation(self, relation):
        assert isinstance(relation, dict)
        self.relation = relation
    def get_relation(self):
        return self.relation

    def set_atoms(self, atoms):
        assert isinstance(atoms, set)
        for atom in atoms:
            assert isinstance(atom, str)
        self.atoms = atoms
    def get_atoms(self):
        return self.atoms
    
    def set_valuation(self, valuation):
        assert isinstance(valuation, dict)
        self.valuation = valuation
    def get_valuation(self):
        return self.valuation

class Node(object):
    """ Represents a node """    
    def __init__(self,
                 codon = None,
                 children = list(),
                 formula = None,
                 operator = None
                 ): 
        self.codon = codon
        self.children = children
        self.formula = formula
        self.operator = operator

    def set_codon(self, codon):
        assert isinstance(codon, str)
        self.codon = codon
    def get_codon(self):
        return self.codon

    def set_children(self, children):
        assert isinstance(children, list)
        self.children = children
    def get_children(self):
        return self.children

    def set_formula(self, formula):
        self.formula = formula
    def get_formula(self):
        return self.formula

    def set_operator(self, operator):
        self.operator = operator
    def get_operator(self):
        return self.operator

    def display(self):
        print 'Node: %s' % self
        print 'Codon: %s' % self.codon
        print 'Children: %s' % self.children
        print 'Formula: %s' % self.formula
        print 'Operator: %s' % self.operator
