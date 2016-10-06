from classes import Model
from classes import Node
import string
import itertools

#####################
#NOTES:
#####################

#Convention wrt syntax: Brackets around each negation, conjunction, diamond, etc., e.g., (E&(~E)), ((~E)&E), (~(E&E)), (~((~E)&E)), (dE), (~((dW)&E)))
#Convention wrt to naming functions: 'give_' prefix indicates that something is returned; 'build_' prefix indicates that class object is changed

#####################
#FUNCTIONS:
#####################

def give_modeltruth(model, formula, verbose): #Return truthvalue of formula in model
    for codon in model.get_domain():
        truthvalue = give_codontruth(model, codon, formula, verbose)
        if truthvalue != '1':
            print '%s is FALSE in the Simple Model' % formula
            print ''
            return truthvalue
    print '%s is TRUE in the Simple Model' % formula
    print ''
    return truthvalue
            
def give_codontruth(model, codon, formula, verbose): #Return truthvalue of formula at codon in model
    try:
        assert codon in model.get_domain()
    except AssertionError:
        print 'Error: Not a codon!'
        print ''
        return None
    try:
        assert isinstance(formula, str)
    except AssertionError:
        print 'Error: Not a formula!'
        print ''
        return None
    try:
        for symbol in formula:
            assert (symbol in model.get_atoms() or symbol in ['(', ')', '~', '&', 'v', '>', 'd', 'b'])
    except AssertionError:
        print 'Error: Formula not well-formed!'
        print ''
        return None
    try:
        root = Node()
        root.set_codon(codon)
        root.set_formula(formula)
        tree = build_tree(model, root, verbose)       
        build_evaluation(model, tree, verbose) 
        truthvalue = tree.get_formula()
        if truthvalue == '1':
            pass
            print '%s is TRUE at %s in the Simple Model' % (formula, codon)
            print ''
        else:
            print '%s is FALSE at %s in the Simple Model' % (formula, codon)
            print ''
        return truthvalue
    except IndexError:
        print 'Error: Formula not well-formed (brackets)!'
        print ''
        return None

#--------------------
#MODEL:
#--------------------

def give_simplemodel(): #Return the simple model
    simplemodel = Model()

    #Setup domain:
    domain = set()         
    for bases in itertools.product('ACGT', repeat=3):
        codon = ''.join(bases)
        domain.add(codon)
    simplemodel.set_domain(domain)

    #Setup relation:
    relation = dict()       
    alphabet = ['A','C','G','T']
    for codon in simplemodel.get_domain():
        position = 0
        substitutions = set()
        while position < len(codon):
            for letter in alphabet:
                substitution = list(codon)
                substitution[position] = letter
                substitution = ''.join(substitution)
                substitutions.add(substitution)
            position += 1
        substitutions.remove(codon)
        relation[codon] = list(substitutions)
    simplemodel.set_relation(relation)

    #Setup atoms:
    simplemodel.set_atoms(set(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']))

    #Setup valuation:
    simplemodel.set_valuation({'A': ['GCA', 'GCC', 'GCG', 'GCT'],
                               'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
                               'N': ['AAC', 'AAT'],
                               'D': ['GAC', 'GAT'],
                               'C': ['TGC', 'TGT'],
                               'Q': ['CAA', 'CAG'],
                               'E': ['GAA', 'GAG'],
                               'G': ['GGA', 'GGC', 'GGG', 'GGT'],
                               'H': ['CAC', 'CAT'],
                               'I': ['ATA', 'ATC', 'ATT'],
                               'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
                               'K': ['AAA', 'AAG'],
                               'M': ['ATG'],
                               'F': ['TTC', 'TTT'],
                               'P': ['CCA', 'CCC', 'CCG', 'CCT'],
                               'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
                               'T': ['ACA', 'ACC', 'ACG', 'ACT'],
                               'W': ['TGG'],
                               'Y': ['TAC', 'TAT'],
                               'V': ['GTA', 'GTC', 'GTG', 'GTT'],
                               '*': ['TAA', 'TAG', 'TGA']})
    return simplemodel

#--------------------
#SEMANTICS:
#--------------------

def build_tree(model, root, verbose): #Update root into tree (recursive)
    operator = give_operator(root.get_formula())
    
    if operator == 'atom':
        root.set_operator(operator)

    if operator == 'negation':
        child = Node()
        child.set_codon(root.get_codon())
        child.set_formula(root.get_formula()[2:len(root.get_formula())-1])
        root.set_children([child])
        root.set_operator(operator)

    if operator == 'conjunction':
        build_children(operator, root)

    if operator == 'disjunction':
        build_children(operator, root)

    if operator == 'implication':
        build_children(operator, root)
    
    if (operator == 'diamond' or operator == 'box'):
        substitutions = model.get_relation()[root.get_codon()]
        children = list()
        for substitution in substitutions:
            child = Node()
            child.set_codon(substitution)
            child.set_formula(root.get_formula()[2:len(root.get_formula())-1])
            children.append(child)
        root.set_children(children)
        root.set_operator(operator)         

    for child in root.get_children():
        build_tree(model, child, verbose)

    return root

def build_evaluation(model, node, verbose): #Update node with truthvalue (recursive)
    if verbose == 1:
        print '%s at %s (%s)' % (node.get_formula(), node.get_codon(), node.get_operator())
        
    if node.get_operator() == 'atom':
        if node.get_codon() in model.get_valuation()[node.get_formula()]:
            node.set_formula('1')
            node.set_operator('resolved')
        else:
            node.set_formula('0')
            node.set_operator('resolved')            

    if node.get_operator() == 'negation':
        for child in node.get_children():
            while child.get_operator() != 'resolved':
                build_evaluation(model, child, verbose)
            if child.get_formula() == '1':
                node.set_formula('0')
                node.set_operator('resolved')
            else:
                node.set_formula('1')
                node.set_operator('resolved')
         
    if (node.get_operator() == 'conjunction' or node.get_operator() == 'box'):
        counter = 0
        for child in node.get_children():
            while child.get_operator() != 'resolved':
                build_evaluation(model, child, verbose)
            if child.get_formula() == '1':
                counter += 1
        if counter == len(node.get_children()):
            node.set_formula('1')
            node.set_operator('resolved')
        else:
            node.set_formula('0')
            node.set_operator('resolved')

    if (node.get_operator() == 'disjunction' or node.get_operator() == 'diamond'):
        counter = 0
        for child in node.get_children():
            while child.get_operator() != 'resolved':
                build_evaluation(model, child, verbose)
            if child.get_formula() == '1':
                counter += 1
                break
        if counter > 0:
            node.set_formula('1')
            node.set_operator('resolved')
        else:
            node.set_formula('0')
            node.set_operator('resolved')

    if node.get_operator() == 'implication':
        antecedentchild = node.get_children()[0]
        while antecedentchild.get_operator() != 'resolved':
            build_evaluation(model, antecedentchild, verbose)
        if antecedentchild.get_formula() == '0':
            node.set_formula('1')
            node.set_operator('resolved')
        else:
            consequentchild = node.get_children()[1]
            while consequentchild.get_operator() != 'resolved':
                build_evaluation(model, consequentchild, verbose)
            if consequentchild.get_formula() == '1':
                node.set_formula('1')
                node.set_operator('resolved')
            else:
                node.set_formula('0')
                node.set_operator('resolved')

    if verbose == 1:
        print '%s at %s (%s)' % (node.get_formula(), node.get_codon(), node.get_operator())
        print '.'

    return node
   
#--------------------
#UTILITY:
#--------------------

def build_children(operator, node): #Update node with children for 2-ary operators
    if operator == 'conjunction':
        symbol = '&'
    if operator == 'disjunction':
        symbol = 'v'
    if operator == 'implication':
        symbol = '>'
    child1 = Node()
    child2 = Node()
    child1.set_codon(node.get_codon())
    child2.set_codon(node.get_codon())
    child1.set_formula(give_parts(symbol, node.get_formula())[0])
    child2.set_formula(give_parts(symbol, node.get_formula())[1])
    node.set_children([child1, child2])
    node.set_operator(operator)

def give_operator(formula): #Return type of formula (highest-ranking operator)
    if '(' not in formula:
        return 'atom'
    opened = give_leftmost('(', formula)
    if formula[opened+1] == '~':
        return 'negation'
    if formula[opened+1] == 'd':
        return 'diamond'
    if formula[opened+1] == 'b':
        return 'box'
    try:
        give_parts('&', formula)
        return 'conjunction'
    except:
        pass
    try:
        give_parts('v', formula)
        return 'disjunction'
    except:
        pass
    try:
        give_parts('>', formula)
        return 'implication'
    except:
        pass
    
def give_leftmost(symbol, formula): #Return leftmost symbol in formula
    assert isinstance(formula, str)
    assert isinstance(formula, str)
    position = 0
    while True:
            if formula[position] == symbol:
                return position
            else:
                position += 1

def give_parts(symbol, formula): #Return both conjuncts of a conjuction (resp. disjuncts of disjunction resp. antecedent/consequent of implication)
    position = give_leftmost('(', formula)
    counter = 0
    while True:
        if formula[position] == '(':
            counter += 1
        if formula[position] == ')':
            counter -= 1
        if formula[position] == symbol and counter == 1:
            return [formula[1:position], formula[position+1:len(formula)-1]]
        position += 1
        


