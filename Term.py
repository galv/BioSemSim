'''
Created on Jun 25, 2013

@author: s-galvez
'''

import Ontology
import AnnotationSet
import math
import sets

class Term(Ontology.OboTerm):

    
    def __init__(self, id, name, ontol):
        super(Term, self).__init__(id, name, ontol)
        #self.originalAncestor = namespace #= Use the root ancestor of one of the three ontologies.
        self.ancestors = sets.Set() #Change to general ancestors, possibly.
        self.descendants= sets.Set()
        self.descendantsIC = 0
        self.annotationsIC = None
        #self.dag = {} #Use for Wang method
        self.annotation_count = 0
        self.annotatedObjects = {}
        self.wang_s = {} #Dict mapping a node's ancestors to their semantic contributions.
        self.wang_sv = 0
        self.ancestorGraph ="" #String
    #Can use to override OboTerm method. 
    #def __str__(self):
    def getIC(self, type):
        if type == "annotations":
            return self.annotationsIC
        else:
            return self.descendantsIC
        
    def __str__(self):
        return self.id
    