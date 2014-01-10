'''
Created on Jul 18, 2013

@author: s-galvez
'''

import ObjectSimilarity
import AnnotationSet
import MyClosure

class DAGPrinter(object):
    '''
    classdocs
    '''

    def __init__(self):
        pass
    def outputGraph(self, object1, object2,root, icType, ontology, annotationSet):
        #outputString = "digraph test {\n\tgraph [ratio=fill];\n\tnode [label=\"\N\", color=black, fillcolor=white, fontcolor=blue, fontsize=10, shape=box, style=filled];\n\tedge [fontsize=8];\n\tgraph [bb=\"0,0,162,446\"];\n\taccall [label=\"all\\nall\", fontname=Courier, pos=\"72,18\", width=\"0.75\", height=\"0.5\"];\n"
        if root == "F":
            rootNode = "0003674"
        elif root == "C":
            rootNode = "0005575"
        elif root == "P":
            rootNode = "0008150"
        else:
            rootNode = "0000001"
        #print root
        #print object1.object_ID, object2.object_ID #size = \"7.75,10.25\"
        #Deleted: size = \"7.75,10.25\",  root = " + rootNode + ",
        outputString = "digraph {\n\tgraph [size = \"7.75,10.25\", ranksep = 1];\n\tgraph [ label = \"" + object1.symbol + "//" + object1.object_ID + " (blue), " + object2.symbol + "//" + object2.object_ID + " (pink), " + "Intersection green\", "  +"fontsize = 50];"
        intersection = ObjectSimilarity.ObjectSimilarity.cachedIntersections[(object1.object_ID, object2.object_ID, root)]
        #for node in intersection:
        #    print node, AnnotationSet.AnnotationSet.TERMS[node].namespace
        alreadyAdded = {}
        i = 0
        while i < len(intersection.keys()):
            term = AnnotationSet.AnnotationSet.TERMS[intersection.keys()[i]]
            outputString += "\t" + term.id.replace("GO:", "").replace("MP:", "") + " [label=\"" + term.id +"\\n" + term.name + "\\n" + str(term.getIC(icType)) + "\"," + "fontname=Courier, " + "style=filled, " + "fillcolor=greenyellow, " + "fontsize = 50];\n"
            alreadyAdded[term.id] = i
            i +=1
        j = 0
        while j < len(object1.annotationDict[root].keys()):
            if not alreadyAdded.has_key(object1.annotationDict[root].keys()[j]):
                term = AnnotationSet.AnnotationSet.TERMS[object1.annotationDict[root].keys()[j]]
                outputString += "\t" + term.id.replace("GO:", "").replace("MP:", "") + " [label=\"" + term.id +"\\n" + term.name + "\\n" + str(term.getIC(icType)) + "\"," + "fontname=Courier, " + "style=filled, " + "fillcolor=cyan2, " + "fontsize = 50];\n"
                alreadyAdded[term.id] = i
                i += 1
                
            j += 1
        k = 0
        while k < len(object2.annotationDict[root].keys()):
            if not alreadyAdded.has_key(object2.annotationDict[root].keys()[k]):
                term = AnnotationSet.AnnotationSet.TERMS[object2.annotationDict[root].keys()[k]]
                outputString += "\t" + term.id.replace("GO:", "").replace("MP:", "") + " [label=\"" + term.id +"\\n" + term.name + "\\n" + str(term.getIC(icType)) + "\"," + "fontname=Courier, " + "style=filled, " + "fillcolor=lightpink, " + "fontsize = 50];\n"
                alreadyAdded[term.id] = i
                i +=1
            k += 1
        #Create ancestor nodes
        alreadyAddedNonAncestors = alreadyAdded
        for term_ID in alreadyAddedNonAncestors.keys():
            ancestors = AnnotationSet.AnnotationSet.TERMS[term_ID].ancestors
            for term in ancestors:
                outputString += "\t"  + term.id.replace("GO:", "").replace("MP:", "") + " [label=\"" + term.id +"\\n" + term.name + "\\n" + str(term.getIC(icType)) + "\"," + "fontname=Courier, " + "fontsize = 50];\n"
                alreadyAdded[term.id] = i
                i += 1
        #Create edges
        alreadyAddedEdges ={}        
        for term_id in alreadyAdded.keys():
            term = AnnotationSet.AnnotationSet.TERMS[term_id]
            edges = term.ancestorGraph
            for edge in edges:
                if not alreadyAddedEdges.has_key(edge):
                    outputString += edge
                    alreadyAddedEdges[edge] = 0
        
        outputString += "\n}"
        return outputString #Print to file later.