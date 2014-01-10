'''
Created on Jun 25, 2013

@author: s-galvez
'''
import DAG
#import Ontology
import math
import AnnotationSet

class ForwardClosure(DAG.Closure):
    '''
    classdocs
    '''
    
    def __init__(self, annotationSet, ontology):
        super(ForwardClosure, self).__init__()
        self.annotations = annotationSet.annotations
        self.totalNodes = {}
        self.totalUniqueAnnotations = {}
        self.ontology = ontology
        for namespace in ontology.getNamespaces():
            self.totalUniqueAnnotations[namespace] = 0
            self.totalNodes[namespace] = 0
        
        self.annotationSet = annotationSet
        for node in ontology.getNodes(): #Note that each node is a Term object.
            self.totalNodes[node.namespace] += 1
            self.totalUniqueAnnotations[node.namespace] += len(annotationSet.annotations[node.id])
    def afterNode(self,dag,node,path):
        node.descendants = self.closure[node]
        #if not node.is_obsolete == ['true']: #Filter out unannotated terms and obsolete terms
        node.annotations = self.annotationSet.getAnnotationsForTerm(node.id)
        for child in dag.getChildren(node):
            node.annotations |= child.annotations
        if len(node.annotations) == 0:
            pass
        else:
            node.annotationsIC = -math.log(1.0*len(node.annotations)/self.totalUniqueAnnotations[node.namespace])
        node.descendantsIC = -math.log(1.0*len(node.descendants)/self.totalNodes[node.namespace])
        
        for a in node.descendants:
            if a == node:
                a.wang_s[node] = 1
            else:
                max = 0
                for c in dag.getChildren(node):
                    if a in c.descendants:
                        #Choose edge weight based on 
                        #edge = dag.getEdge(a,c) #ontology instead of dag?
                        #dag.getNode(a)
                        tuple = dag.nodes[a]
                        #edge = aDict[c]
                        in_edges = tuple[0]
                        #out_edges = tuple[1]
                        #edge = dag.getEdge(node,a)
                        #print edge
                        
                        edge = dag.getEdge(node, c)
                        
                        #print in_edges
                        #edge = in_edges[node]                      
                        if edge == "is_a":
                            w = 0.8
                        elif edge == "part_of":
                            w = 0.6
                        else:
                            w = 0.1 #Not too sure what to do with "regulates" edges...
                        newval = w * a.wang_s[c]
                        if newval > max:
                            max = newval
                a.wang_s[node] = max
        
        
        
class ReverseClosure(DAG.Closure):
        
    def __init__(self, reversed = True):
        super(ReverseClosure, self).__init__()
        #self.annotationSet = annotationSet.annotations
    def afterNode(self,dag,node,path):
        node.ancestors = self.closure[node]
        
        node.ancestorGraph = []
        
        for ancestor in node.ancestors:
            parents = dag.getParents(ancestor)
            #print parents
            for parent in parents:
                edge = dag.getEdge(parent, ancestor)
                #node.ancestorGraph += "\n" + ancestor.id.replace("GO:", "") + " -> " + parent.id.replace("GO:", "")
                addendum = "\n" + ancestor.id.replace("GO:", "").replace("MP:", "") + " -> " + parent.id.replace("GO:", "").replace("MP:", "") + "[len = .3, weight = 2"
                
                if edge == "part_of":
                    addendum += ", color = deepskyblue2"
                elif edge == "negatively_regulates":
                    addendum += ", color = red"
                elif edge == "positively_regulates":
                    addendum += ", color = green"
                elif edge == "regulates":
                    addendum += ", color = purple"
                else:
                    addendum += ", color = black"
                addendum += "]"
                
                node.ancestorGraph.append(addendum)
        '''
        #print node
        #print "ancestors", node.ancestors
        for ancestor in node.ancestors:
            parents = dag.getParents(ancestor)
            for parent in parents:
                edge = dag.getEdge(parent, ancestor)
                print "Edge,", edge
                node.ancestorGraph += "\n"  + parent.id.replace("GO:", "").replace("MP:", "") + "->" + node.id.replace("GO:", "").replace("MP:", "") + " [minlen = .5"
                
                if edge == "part_of":
                    node.ancestorGraph += ", color = blue"
                elif edge == "negatively_regulates":
                    node.ancestorGraph += ", color = red"
                elif edge == "positively_regulates":
                    node.ancestorGraph += ", color = green"
                elif edge == "regulates":
                    node.ancestorGraph += ", color = purple"
                else:
                    node.ancestorGraph += ", color = black"
                node.ancestorGraph += "]"
        '''
        if not node.is_obsolete == ['true']:
            for a in node.ancestors:
                #print a.id
                node.wang_sv += node.wang_s[a]
        
class TinyClosure(ReverseClosure):
    def afterNode(self,dag,node,path):
        ancestors = self.closure[node]
        node.ancestorGraph = ""
        
        for ancestor in ancestors:
            parents = dag.getParents(ancestor)
            for parent in parents:
                edge = dag.getEdge(parent, ancestor)
                print "Edge,", edge
                node.ancestorGraph += "\n"  + parent.id.replace("GO:", "").replace("MP:", "") + "->" + node.id.replace("GO:", "").replace("MP:", "") + " [minlen = .5"
                
                if edge == "part_of":
                    node.ancestorGraph += ", color = blue"
                elif edge == "negatively_regulates":
                    node.ancestorGraph += ", color = red"
                elif edge == "positively_regulates":
                    node.ancestorGraph += ", color = green"
                elif edge == "regulates":
                    node.ancestorGraph += ", color = purple"
                else:
                    node.ancestorGraph += ", color = black"
                node.ancestorGraph += "]"
                
                