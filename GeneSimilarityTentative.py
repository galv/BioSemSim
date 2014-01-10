'''
Created on Jun 20, 2013

# 6-20-13: NEw idea: make all methods in GeneSimilarity static.

@author: s-galvez
'''

import math
import Term
import AnnotationSet

class GeneSimilarity(object):
    '''
    classdocs
    '''
    TOTALCCNODES = len(Term.Term("GO:0005575", "C").descendantIDs)
    TOTALMFNODES = len(Term.Term("GO:0003674", "F").descendantIDs)
    TOTALBPNODES = len(Term.Term("GO:0008150", "P").descendantIDs)
    ALLANNOTATIONS = AnnotationSet.AnnotationSet.ALLANNOTATIONS #Dictionary file

    def __init__(self):
        '''
        self.g1 = g1
        self.g2 = g2
        #self.annotationSet = annotationSet
        #self.g1Annotations
        #self.g2Annotations
        #self.fileGOID = fileGOID
        self.allAnnotations = allAnnotations
        '''
        self.intersection = {}
        #self.findIntersectionNodes()
        self.union = {}
        #self.findUnionNodes()

        #Intersection of the nodes of each gene.
        
    def findIntersectionNodes(self, g1, g2):
        #print self.g1.annotations
        #print self.g2.annotations
        #pseudo-code: intersection[ID] = term
        
        for annotation1 in self.g1.geneSpecificAnnotations:
            for annotation2 in self.g2.geneSpecificAnnotations:
                if annotation1.term.GO_ID == annotation2.term.GO_ID:
                    self.intersection[annotation1.term.GO_ID] = annotation1.term
                    #print self.findIC(annotation1.term)

        
    def findUnionNodes(self, g1, g2):
        for annotation in self.g1.geneSpecificAnnotations:
            self.union[annotation.term.GO_ID] = annotation.term
        for annotation in self.g2.geneSpecificAnnotations:
            GO_IDkeys = self.union.keys
            if GO_IDkeys.count(annotation.term.GO_ID) == 0:
                self.union.append(annotation.term)
                #print self.findIC(annotation.term)
                
    def findIC(self, term):
        numberAnnotations = 0
        for annotation in self.allAnnotations.annotations:
            if annotation.term.GO_ID == term.GO_ID:
                numberAnnotations = numberAnnotations + 1
        p = ((float)(numberAnnotations))/len(self.allAnnotations.annotations)
        #print numberAnnotations
        #print len(self.allAnnotations)
        return -math.log(p)
    
    def findchildrenbasedIC(self, term):
        totalChildren = len(term.descendantIDs)
        #print len(term.descendantIDs)
        totalNodes = 0
        if term.ontology == "C":
            totalNodes = GeneSimilarity.TOTALCCNODES
        if term.ontology == "F":
            totalNodes = GeneSimilarity.TOTALMFNODES
        if term.ontology == "P":
            totalNodes = GeneSimilarity.TOTALBPNODES
        p = (1.0*totalChildren)/totalNodes
        #print GeneSimilarity.TOTALCCNODES
        return -math.log(p)
        
    #Creates similarity measurement which equals magnitude of union of g1 and g2 divided by the magnitude of the intersection of g1 and g2. 
class JaccardSimilarity(GeneSimilarity):
    @staticmethod
    def findSimilarity(g1, g2):
        findIntersectionNodes(g1,g2)
        totalNodesNumber = len(g1.geneSpecificAnnotations) + len(g2.geneSpecificAnnotations) - len(self.intersection) #equals number of common modes. (intersection)
        print len(g1.geneSpecificAnnotations)
        print len(g2.geneSpecificAnnotations)
        print len(self.intersection)
        '''
        for annotation in self.g1.geneSpecificAnnotations:
            print annotation.term.GOID
        for annotation in self.g2.geneSpecificAnnotations:
            print annotation.term.GOID
        '''
        #print self.totalNodesNumber
        #print len(self.intersection)
        #print self.intersection
        #print self.totalNodesNumber
        return (float)(len(self.intersection))/self.totalNodesNumber
    
    #Creates similarity measurement like original method, except that method obtains and compares all algorithms  degree degrees away from the nodes genes are originally annotated to. 
    #All of these nodes are then analyzed.  
    def findExtendedJaccardSimilarity(self, degree):
        pass
class DiceSimilarity(GeneSimilarity):
    def findSimilarity(self):
        return 2.0*len(self.intersection)/(len(self.g1.geneSpecificAnnotations) + len(self.g2.geneSpecificAnnotations))

class GICSimilarity(GeneSimilarity):
    def findSimilarity(self):
        numerator = 0.0
        for node in self.intersection:
            numerator = numerator + self.findIC(node)
        denominator = 0.0
        for node in self.union:
            denominator = denominator + self.findIC(node)
        return numerator/denominator
'''
#This class is used to separate transitive closure-based code from methods not needing the transitive closure.    
class TransitiveClosureBasedGeneSimilarity(GeneSimilarity):
    def __init__(self):
'''     
class ResnikSimilarity(GeneSimilarity):

    #Note: this isn't a gene similarity measure so much as a term similarity measure. Just to note. Should use it on
    def findResnikMeasure(self):
        MICAic = 0 #Information content of most informative common ancestor.
        for g1Annotation in GeneSimilarity.g1.geneSpecificAnnotations:
            for g2Annotation in GeneSimilarity.g2.ancestorIDs:
                if g1Annotation.GO_ID == g2Annotation.GO_ID and self.findchildrenbasedIC(g1Annotation.term) > MICAic:
                    MICAic = self.findchildrenbasedIC(g1Annotation.term)
        return MICAic

    
    '''
    def findLinMeasure(self):
        
    def findJiangMeasure(self):
        
    def findRelMeasure(self):
    ''' 

