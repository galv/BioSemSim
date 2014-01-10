'''
Created on Jun 12, 2013

Contains data from individual row of data (therefore, for an individual gene or gene product) under "Gene Ontology Data" on this webpage: ftp://ftp.informatics.jax.org/pub/reports/index.html#go
Note that index 0 refers to column 1.

@author: s-galvez
'''

import AnnotatedObject
import AnnotationSet

class Annotation(object):
    '''
    classdocs
    '''

    #annotationStrings is list of strings containing all data for one annotation (one line of text file)
    def __init__(self, annotationStrings, annotationLoader):
    #See above to see what each row corresponds to.
        self.annotationStrings = annotationStrings
        self.object_ID = self.annotationStrings[1]
        self.term_ID = self.annotationStrings[4]
        self.ontology = self.annotationStrings[8]
        self.inferredFrom = self.annotationStrings[7]
        self.evidenceCode = self.annotationStrings[6]
        self.term = annotationLoader.terms[term_ID]
        #Create other shortcuts to values in annotation as needed.
    def __str__(self):
        return self.term_ID
class MPAnnotation(object):
    
    def __init__(self, annotationStrings):
        self.annotationStrings = annotationStrings
        self.object_ID = annotationStrings[0]
        self.object_symbol = annotationStrings[1]
        self.term_ID = annotationStrings[3]
        self.disease = None