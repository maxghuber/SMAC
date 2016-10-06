from functions import give_codontruth
from functions import give_simplemodel

def SMAC():
    print '-' * 72
    print '|' + ' ' * 11 + 'SMAC (Simple Model Aminoacid Checker) version 1' + ' ' * 11 + '|'
    print '-' * 72
    print '| X := Atom | (~X) | (X&X) | (XvX) | (X>X) | (dX) | (bX)' + ' ' * 15 + '|'
    print '| Atoms: A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, *' + ' ' + '|'
    print '| Examples: E, (~W), (~(E&(~E))), (d((bD)>(d((dE)v(~K)))))' + ' ' * 13 + '|'
    print '-' * 72
    
    simplemodel = give_simplemodel()
    print 'Verbose(Y/N)'
    answer = raw_input()
    if answer == 'Y':
        verbose = 1
    else:
        verbose = 0
    while True:
        print 'Specifiy codon:'
        codon = raw_input()
        print 'Specifiy formula:'
        formula = raw_input()
        give_codontruth(simplemodel, codon, formula, verbose)
        print 'Again? (Y/N)'
        answer = raw_input()
        if answer != 'Y':
            raise SystemExit(0)

SMAC()
