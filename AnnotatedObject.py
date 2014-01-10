'''
Created on Jun 13, 2013

MAy better be called Gene Product or just "Item" referring to genes and gene products.

@author: s-galvez
'''

#import GeneSimilarity
#import AnnotationSet
class AnnotatedObject(object):
    '''
    classdocs
    '''


    def __init__(self, object_ID, symbol, fileName, oneObjectAnnotationsCC, oneObjectAnnotationsMF = None, oneObjectAnnotationsBP = None):
        self.object_ID = object_ID
        if fileName == "gene_association.mgi":
            self.annotationDict = {"P" : oneObjectAnnotationsBP, "F": oneObjectAnnotationsMF, "C" : oneObjectAnnotationsCC} #geneSpecificAnnotations #Dict with 
        else:
            self.annotationDict = {"MP" : oneObjectAnnotationsCC}
        self.symbol = symbol
        self.termAncestors = {}
        self.methodStrings = ["jaccard", "jaccard_ext", "gic_ext", "wang"] #These strings are the indices of simMeasures
        self.resnikAvgMICAs = []
        self.resnikMaxMICA = ""
        self.disease = "None"