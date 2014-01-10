'''
Created on Jun 13, 2013

@author: s-galvez
'''

import math
import AnnotationSet
import sys
import sets
class ObjectSimilarity(object):
    '''
    classdocs
    '''
    cachedIntersections = {} #Dict mapping tuples of two genes' ids to intersection of those two genes.
    
    def __init__(self):
        pass
    def findIntersectionNodes(self, 
                              annotations1, 
                              annotations2,
                              g1,
                              g2,
                              root):
        #Finds intersection of two sets of nodes.
        if ObjectSimilarity.cachedIntersections.has_key((g1.object_ID, g2.object_ID, root)):
            pass
        else:
            intersection = {}
            if len(annotations1.keys()) < len(annotations2.keys()):
                for key1 in annotations1.keys():
                    if key1 in annotations2.keys():
                        intersection[key1] = annotations1[key1]
            else:
                for key2 in annotations2.keys():
                    if key2 in annotations1.keys():
                        intersection[key2] = annotations1[key2]
            
            ObjectSimilarity.cachedIntersections[(g1.object_ID, g2.object_ID, root)] = intersection
        return ObjectSimilarity.cachedIntersections[(g1.object_ID,g2.object_ID, root)]
        
    def findUnionNodes(self, annotations1, annotations2):
        #Finds union of two nodes.
        union = {}
        for key in annotations1.keys():
            union[key] = annotations1[key]
        for key in annotations2.keys():
            union[key] = annotations2[key]
        return union

    #Creates similarity measurement which equals magnitude of union of g1 and g2 divided by the magnitude of the intersection of g1 and g2. 
class JaccardSimilarity(ObjectSimilarity):
    def findSimilarity(self, g1, #First object
                        g2, #Second object
                         root): #"F", "P", or "C"
    #Returns intersection/union of g1 and g2. Between 0 and 1. Higher valeus indicate greater similarity.
        intersection = self.findIntersectionNodes(g1.annotationDict[root],g2.annotationDict[root],g1,g2, root)
        union = self.findUnionNodes(g1.annotationDict[root],g2.annotationDict[root])
        try:
            return (1.0*len(intersection))/len(union)
        except ZeroDivisionError:
            return 0
#Can make more efficient by finding union then applying equation to find length of intersection.
class JaccardSimilarityExtended(ObjectSimilarity):
    
    cached = {}# Key: (g, root) Value: Second dictionary. See findAncestors.
    cachedIntersections = {} #Dict mapping tuple of (object1 id, object2 id, root) to dictionary mapping term ids to terms that make up a genes ancestors. 
    
    def findIntersectionNodes(self,g1Ancestors,g2Ancestors,g1,g2, root):
        if JaccardSimilarityExtended.cachedIntersections.has_key((g1.object_ID, g2.object_ID, root)):
            pass
        else:
            intersection = {}
            if len(g1Ancestors.keys()) < len(g2Ancestors.keys()):
                for key1 in g1Ancestors.keys():
                    if key1 in g2Ancestors.keys():
                        intersection[key1] = g1Ancestors[key1]
            else:
                for key2 in g2Ancestors.keys():
                    if key2 in g1Ancestors.keys():
                        intersection[key2] = g1Ancestors[key2]
            
            JaccardSimilarityExtended.cachedIntersections[(g1.object_ID, g2.object_ID, root)] = intersection
        return JaccardSimilarityExtended.cachedIntersections[(g1.object_ID,g2.object_ID, root)]
    def findSimilarity(self, g1, g2, root):
    #Returns index between 0 and 1 (least to most similar based on intersection of each gene's set of terms' total set of ancestors.
        g1Ancestors = self.findAncestors(g1,root)
        g2Ancestors = self.findAncestors(g2,root)
        intersection = self.findIntersectionNodes(g1Ancestors,g2Ancestors,g1,g2, root)
        union = self.findUnionNodes(g1Ancestors, g2Ancestors)
        try:
            return 1.0*len(intersection)/len(union)
        except ZeroDivisionError:
            return 0 #If a gene has no annotations, assume it has the minimum similarity.
        #try:
        #    return 0.5*(((float)(len(intersection)))/len(g1Ancestors) + ((float)(len(intersection)))/len(g2Ancestors))
    def findAncestors(self,g, #Object of interest
                      root): #Ontology of the object of interest.
    #Return dict mapping term id to term object. Dict contains all ancestor terms of g's terms.
        gAncestors = {}
        if not JaccardSimilarityExtended.cached.has_key((g,root)):
            for term_ID in g.annotationDict[root].keys():
                term = AnnotationSet.AnnotationSet.TERMS[term_ID]
                for ancestor in term.ancestors:
                    gAncestors[ancestor.id] = ancestor #ancestor is an instance of Term.
            JaccardSimilarityExtended.cached[(g, root)] = gAncestors
        return JaccardSimilarityExtended.cached[(g,root)]
    
class GICSimilarityExtended(JaccardSimilarityExtended):
    def __init__(self, icType):
        self.icType = icType
    def findSimilarity(self, g1, g2, root):
        try:
            g1Ancestors = self.findAncestors(g1,root)
            g2Ancestors = self.findAncestors(g2,root)
            intersection = self.findIntersectionNodes(g1Ancestors,g2Ancestors,g1,g2, root)
            union = self.findUnionNodes(g1Ancestors, g2Ancestors)
            unionIC = 0
            for term_id in union.keys():
                unionIC += union[term_id].getIC(self.icType)
            intersectionIC = 0
            for term_id in intersection.keys():
                intersectionIC += AnnotationSet.AnnotationSet.TERMS[term_id].getIC(self.icType)
            return intersectionIC/unionIC
        except ZeroDivisionError:
            return 0 #See JaccardSimilarityExtended for reasoning behind this.
        
class DiceSimilarity(ObjectSimilarity):
    def findSimilarity(self, g1, #First object being compared 
                       g2, #Second object being compared
                        root):
    #Return a similarity value between 0 (least similar) and 1 (most similar) for two genes.
        intersection =  self.findIntersectionNodes(g1.annotationDict[root],g2.annotationDict[root],g1,g2, root)
        return 2.0*len(intersection)/(len(g1.annotationDict[root]) + len(g2.annotationDict[root]))

class GICSimilarity(ObjectSimilarity):
    def __init__(self, icType):
        self.icType = icType
    def findSimilarity(self, g1, g2, root):
    #Return similarity index between 0 (least similar) and 1 (most similar) for two genes based on their terms IC's.
        numerator = 0.0
        intersection = self.findIntersectionNodes(g1.annotationDict[root],g2.annotationDict[root],g1,g2, root)
        union = self.findUnionNodes(g1.annotationDict[root],g2.annotationDict[root])
        
        for key in intersection.keys():
            numerator = numerator + AnnotationSet.AnnotationSet.TERMS[intersection[key].term_ID].getIC(self.icType) #self.findChildrenBasedIC(AnnotationSet.AnnotationSet.TERMS[intersection[key].term_ID])
        denominator = 0.0
        for key in union.keys():
            denominator = denominator + AnnotationSet.AnnotationSet.TERMS[union[key].term_ID].getIC(self.icType)
        try:
            return numerator/denominator
        except ZeroDivisionError:
            return 0
class ResnikSimilarity(ObjectSimilarity):        
    def __init__(self, icType):
        self.icType = icType
    def termSimilarity(self,t1,t2):
    #Returns most informative common ancestor of the two nodes.
        ic_max = 0
        for term1 in t1.ancestors:
            if term1 in t2.ancestors and term1.getIC(self.icType) >= ic_max:
                ic_max = term1.getIC(self.icType)
                mica = term1
        
        return (ic_max, mica)
    
class ResnikSimilarityAvg(ResnikSimilarity):
#May be non-functionnal right now! Don't use without validating against GOSemSim or another package!
    def findSimilarity(self,g1,g2,root):
    #Returns average of all pair-based term values.
        avgR = 0.0
        count2 = len(g1.annotationDict[root])*len(g2.annotationDict[root])
        if count2 == 0:
            return 0 #Sometimes, genes with annotations have no annotations in a specific GO ontology. e.g., http://www.informatics.jax.org/marker/MGI:2387613 Say such have zero similarity
        tupleList = []
        for key1 in g1.annotationDict[root].keys():
            for key2 in g2.annotationDict[root].keys():
                tuple = self.termSimilarity(AnnotationSet.AnnotationSet.TERMS[g1.annotationDict[root][key1].term_ID], AnnotationSet.AnnotationSet.TERMS[g2.annotationDict[root][key2].term_ID])
                test = tuple[0]
                tupleList.append(tuple)
                avgR += test
                g2.resnikAvgMICAs.append(tuple)
        return avgR/count2

class ResnikSimilarityMax(ResnikSimilarity):
    def findSimilarity(self,g1,g2,root):
    #Returns maximum of all pair-based term values.
        tupleList = []
        max = 0.0
        for key1 in g1.annotationDict[root].keys():
            for key2 in g2.annotationDict[root].keys():
                tuple = self.termSimilarity(AnnotationSet.AnnotationSet.TERMS[g1.annotationDict[root][key1].term_ID], AnnotationSet.AnnotationSet.TERMS[g2.annotationDict[root][key2].term_ID])
                tupleList.append(tuple)
                test = tuple[0]
                if test > max:
                    max = test
                    maxTuple = tuple
        try:
            g2.resnikMaxMICA = maxTuple[1].id
        except UnboundLocalError:
            g2.resnikMaxMICA = ""
        return max
#BMA method is supposed to be the best REsnik method. This is not working right now. Must fix before using.
class ResnikSimilarityBMA(ResnikSimilarity):
    def findSimilarity(self,g1,g2,root):
        maxList = []
        len = 0
        if len(g1.annotationDict[root]) > len(g2.annotationDict[root]):
            count = len(g1.annotationDict[root]) 
        else:
            count = len(g2.annotationDict[root])
            
        for key1 in g1.annotationDict[root].keys():
            max = 0
            for key2 in g2.annotationDict[root].keys():
                test = self.termSimilarity(AnnotationSet.AnnotationSet.TERMS[g1.annotationDict[root][key1].term_ID], AnnotationSet.AnnotationSet.TERMS[g2.annotationDict[root][key2].term_ID])
                if test > max:
                    max = test
        return

class LinSimilarityAvg(ResnikSimilarityAvg):
    def termSimilarity(self,t1,t2):
        
        ic_t1 = t1.getIC(self.icType)
        ic_t2 = t2.getIC(self.icType)
        return 2*super(LinSimilarityAvg, self).termSimilarity(t1,t2)/(ic_t1 + ic_t2)
    
class LinSimilarityMax(ResnikSimilarityMax):
    def termSimilarity(self,t1,t2):
        
        ic_t1 = t1.id
        ic_t2 = t2.id
        return 2*super(LinSimilarityMax, self).termSimilarity(t1,t2)/(ic_t1 + ic_t2)
class RelSimilarityAvg(LinSimilarityAvg):
    def termSimilarity(self,t1,t2):
        #Need to find MICA!
        super(RelSimilarityAvg, self).termSimilarity(t1,t2)*(1-len())
    
class WangSimilarity(ObjectSimilarity):
    
    cache = {} #Dict mapping two terms' IDs to termSimilarity returned value.
    
    def findSimilarity(self, g1, g2, root):
        part1 = 0
        for term_ID1 in g1.annotationDict[root].keys():
            part1 += self.termTermsSimilarity(AnnotationSet.AnnotationSet.TERMS[term_ID1], g2.annotationDict[root])
        part2 = 0
        for term_ID2 in g2.annotationDict[root].keys():
            part2 += self.termTermsSimilarity(AnnotationSet.AnnotationSet.TERMS[term_ID2], g1.annotationDict[root])
        try:
            return (part1 + part2)/(len(g1.annotationDict[root]) + len(g2.annotationDict[root]))
        except:
            return 0
    def termSimilarity(self, t1, t2):
        numerator = 0.0
        #Intersect:
        if WangSimilarity.cache.has_key((t1,t2)):
            pass
        else:
            if len(t1.wang_s.keys()) > len(t2.wang_s.keys()):
                for ancestor in t1.wang_s.keys():
                    if t2.wang_s.has_key(ancestor):
                        numerator += t1.wang_s[ancestor] + t2.wang_s[ancestor]
            else:
                for ancestor in t2.wang_s.keys():
                    if t1.wang_s.has_key(ancestor):
                        numerator += t1.wang_s[ancestor] + t2.wang_s[ancestor]
            denominator = t1.wang_sv + t2.wang_sv
            WangSimilarity.cache[(t1,t2)] = numerator/denominator
        return WangSimilarity.cache[(t1,t2)]
        
            
    def termTermsSimilarity(self,t1,annotations):
        max = 0
        for a2 in annotations.keys():
            test = self.termSimilarity(t1,AnnotationSet.AnnotationSet.TERMS[annotations[a2].term_ID])
            if test > max:
                max = test
        return max