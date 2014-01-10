'''
Created on Jun 12, 2013

AnnotationSet is a wrapper class for a list of annotation objects. The individual class was made to encapsulate the doAnnotations() function, which extracts the gene_association fiel into a list
whose objects are lists of strings. Each list of strings contains data on one annotation.

Update: AnnotationSet intializes the all Term, AnnotatedObject, and Annotation objects.

@author: s-galvez
'''

import Annotation
import AnnotatedObject
import math
import sets

class AnnotationSet(object):
    '''
    classdocs
    '''
    ALLANNOTATIONS = {} #Dict mapping term ID to Annotation object. Includes EVERY annotation in the ontology.
    TERMS = {} #Dict mapping id to Term object.
    #ontologyLengths = {"cellular_component" : 0, "biological_process": 0, "molecular_function": 0} #Dict mapping ontology name (see __init__ method) number of nodes in that ontology.

    def __init__(self, fileName, #fileName is the name of the file containing all annotations in the GO.
                 ontology, #A dictionary from DAG.py and Ontology.py  mapping Term objects to Tuples. See Ontology.py for details.
                 evidenceCodes =[], #List of evidence codes to be exlcluded
                 ):
        self.fileName = fileName
        self.annotations = {} #Dict mapping term ID to set of annotations.
        self.annotatedObjects = {} #
        self.evidenceCodes = evidenceCodes
        
        nodes = ontology.getNodes()
        
        for node in nodes:
            self.annotations[node.id] = sets.Set([])
        self.ontology = ontology
        self.terms = {}
        self.makeTermsPublic()
        self.doAnnotations()
        #self.defineICs()
        

    def makeTermsPublic(self): 
    #Change key from some value I don't know to the id (still mapping to Term object). Also makes Terms public variable.
        for term in self.ontology.getNodes(): 
            if not term.is_obsolete == ['true']: #For some reason, is_obselete is a one-index list containing a string stating "true" or "false" instead of a boolean
                AnnotationSet.TERMS[term.id] = term
                self.terms[term.id] = term
    def getAnnotationsForTerm(self,term_id):
    #Returns set of annotations for objects annotated to term.
        return self.annotations[term_id]
    def doAnnotations(self):
    #Creates all Annotation objects and puts them into annotations list.
        fp = open(self.fileName,'r')
        objects = fp.readlines()
        first = True
        lastObject_ID = ""
        oneObjectAnnotationsCC = {}
        oneObjectAnnotationsMF = {}
        oneObjectAnnotationsBP = {}
        #oneTermAnnotations = []
        for line in objects:
            if line.find("!", 0,1) == 0:
                pass
            else: 
                annotation_data = line.split("\t")
                column4List = annotation_data[3].split(" ")
                if(column4List[0] == "NOT" or not self.annotations.has_key(annotation_data[4]) or annotation_data[6] in self.evidenceCodes): #Remove all "NOT" annotations (see GO documentation) and annotations to obsolete nodes.
                    pass
                else:
                    annotation = Annotation.Annotation(annotation_data, self)
                    object_ID = annotation_data[1]
                    self.annotations[annotation.term_ID] |= sets.Set([annotation])
                    #except KeyError: #Must remove all obsolete terms
                        #pass
                        #self.annotations[annotation.term_ID] = sets.Set([annotation])                    
                    if(first):
                        first = False
                        lastObject_ID = annotation_data[1]
                        lastObject_symbol = annotation_data[2]
                    elif(object_ID != lastObject_ID):
                        #self.annotations[lastObject_ID] = annotation #oneObjectAnnotationsCC + oneObjectAnnotationsMF + oneObjectAnnotationsBP
                        self.annotatedObjects[lastObject_ID] = AnnotatedObject.AnnotatedObject(lastObject_ID, lastObject_symbol, self.fileName, oneObjectAnnotationsCC, 
                                                                                               oneObjectAnnotationsMF, oneObjectAnnotationsBP) #You should check this!
                        #AnnotationSet.ontologyLengths["cellular_component"] += len(oneObjectAnnotationsCC)
                        #AnnotationSet.ontologyLengths["biological_process"] += len(oneObjectAnnotationsBP)
                        #AnnotationSet.ontologyLengths["molecular_function"] += len(oneObjectAnnotationsMF)
                        lastObject_ID = object_ID
                        lastObject_symbol = annotation_data[2]
                        oneObjectAnnotationsCC = {}
                        oneObjectAnnotationsBP = {}
                        oneObjectAnnotationsMF = {}
                        if annotation.ontology == "C":
                            oneObjectAnnotationsCC[annotation.term_ID] = annotation
                        elif annotation.ontology == "F":
                            oneObjectAnnotationsMF[annotation.term_ID] = annotation
                        else:
                            oneObjectAnnotationsBP[annotation.term_ID] = annotation
                    else:
                        #Object_ID = annotation_data[4] #GO_ID in Gene Ontology
                        #self.annotations[object_ID] = annotation #Reference an annotation by its term ID, not its gene ID.
                        if annotation.ontology == "C":
                            oneObjectAnnotationsCC[annotation.term_ID] = annotation
                        elif annotation.ontology == "F":
                            oneObjectAnnotationsMF[annotation.term_ID] = annotation
                        else:
                            oneObjectAnnotationsBP[annotation.term_ID] = annotation
                        #lastObject_ID = object_ID
        object_symbol = annotation_data[2]
        self.annotatedObjects[object_ID] = AnnotatedObject.AnnotatedObject(object_ID, object_symbol, self.fileName, oneObjectAnnotationsCC, 
                                                                                               oneObjectAnnotationsMF, oneObjectAnnotationsBP)
        #AnnotationSet.ontologyLengths["cellular_component"] += len(oneObjectAnnotationsCC)
        #AnnotationSet.ontologyLengths["biological_process"] += len(oneObjectAnnotationsBP)
        #AnnotationSet.ontologyLengths["molecular_function"] += len(oneObjectAnnotationsMF)
    def makeAnnotationsPublic(self):
        AnnotationSet.ALLANNOTATIONS = self.annotations #key is object_ID
        
    def getAnnotationsOfGene(self,object_ID): #Returns list of annotation objects that belongs to genes.
        return AnnotationSet.ALLANNOTATIONS[object_ID]
    '''
    def getAnnotation(self, index):
        return self.annotations[index]
    '''
   
class MPAnnotationSet(AnnotationSet):
    ALLANNOTATIONS = {} #Dict mapping Object id to Annotation object. Includes EVERY annotation in the ontology.
    TERMS = {} #Dict mapping id to Term object.
    def __init__(self, fileName, #fileName is the name of the file containing all annotations in the GO.
                 ontology, #A dictionary from DAG.py and Ontology.py  mapping Term objects to Tuples. See Ontology.py for details.
                 diseaseFile,
                 evidenceCodes =[], #List of evidence codes to be exlcluded
                 ):
        super(MPAnnotationSet, self).__init__(fileName, ontology, evidenceCodes)
        #Upload diseases!
        fp = open(diseaseFile, 'r')
        diseases = fp.readlines()
        for disease_assocation in diseases:
            #disease_association = disease.split("\t")
            strain = disease_assocation[:12].replace("\t","")
            disease = disease_assocation[12:].replace("\n","")
            print "strain", strain
            print "disease", disease
            self.annotatedObjects[strain].disease = disease
    def doAnnotations(self):
        fp = open(self.fileName,'r')
        objects = fp.readlines()
        first = True
        lastObject_ID = ""
        oneObjectAnnotations = {}
        
        #self.convertEvidenceCodes(evidenceCodes)
        
        for line in objects:
            line_data = line.split("\t")
            if line_data[0] == "genotypeID" or line_data[5] == "normal" or line_data[2] == "no biological data":
                pass
            else:
                if True:# not self.annotations.has_key(line_data[3]) or line_data[2] in self.evidenceCodes:
                    annotation = Annotation.MPAnnotation(line_data)
                    self.annotations[annotation.term_ID] |= sets.Set([annotation])
                    object_ID = line_data[0]
                    if(first):
                        first = False
                        lastObject_ID = line_data[0]
                        lastObject_symbol = line_data[1]
                    elif(lastObject_ID != object_ID):
                        self.annotatedObjects[lastObject_ID] = AnnotatedObject.AnnotatedObject(lastObject_ID, lastObject_symbol.lstrip(" "), self.fileName, oneObjectAnnotations)
                        lastObject_ID = object_ID
                        lastObject_symbol = line_data[1]
                        oneObjectAnnotations = {}
                        oneObjectAnnotations[annotation.term_ID] = annotation
                    else:
                        oneObjectAnnotations[annotation.term_ID] = annotation
        object_ID = line_data[0]
        object_symbol = line_data[2]
        self.annotatedObjects[object_ID] = AnnotatedObject.AnnotatedObject(object_ID, object_symbol.lstrip(" "), self.fileName, oneObjectAnnotations)
class SubAnnotationSet(AnnotationSet):
    def __init__(self,fileName, #fileName is the name of the file containing all annotations in the GO.
                 ontology, #A dictionary from DAG.py and Ontology.py  mapping Term objects to Tuples. See Ontology.py for details.
                 keyStoneTerms,
                 evidenceCodes =[], #List of evidence codes to be exlcluded
                 ):
        self.keyStoneTerms = keyStoneTerms
        super(SubAnnotationSet, self).__init__(fileName, ontology, evidenceCodes)
        
    def doAnnotations(self):
        fp = open(self.fileName,'r')
        objects = fp.readlines()
        first = True
        lastObject_ID = ""
        oneObjectAnnotationsCC = {}
        oneObjectAnnotationsMF = {}
        oneObjectAnnotationsBP = {}
        testAnnotationsCC = {}
        testAnnotationsMF = {}
        testAnnotationsBP = {}
        tempEvidenceCodes =self.evidenceCodes
        #oneTermAnnotations = []
        for line in objects:
            if line.find("!", 0,1) == 0:
                pass
            else: 
                annotation_data = line.split("\t")
                column4List = annotation_data[3].split(" ")
                if(column4List[0] == "NOT" or not self.annotations.has_key(annotation_data[4]) or annotation_data[6] in self.evidenceCodes): #Remove all "NOT" annotations (see GO documentation) and annotations to obsolete nodes.
                    pass
                else:
                    annotation = Annotation.Annotation(annotation_data)
                    object_ID = annotation_data[1]
                    self.annotations[annotation.term_ID] |= sets.Set([annotation])
                    #except KeyError: #Must remove all obsolete terms
                        #pass
                        #self.annotations[annotation.term_ID] = sets.Set([annotation])                    
                    if(first):
                        first = False
                        lastObject_ID = annotation_data[1]
                        lastObject_symbol = annotation_data[2]
                    elif(object_ID != lastObject_ID):
                        relevantDomain = False
                        #print oneObjectAnnotationsCC.keys()
                        #print set(oneObjectAnnotationsCC.keys())
                        #print set(self.keyStoneTerms)
                        if len(set(oneObjectAnnotationsMF.keys()) & set(self.keyStoneTerms)) >= 1:
                            #oneObjectAnnotationsCC.has_key(self.keyStoneTersm)
                            self.annotatedObjects[lastObject_ID] = AnnotatedObject.AnnotatedObject(lastObject_ID, lastObject_symbol, self.fileName, oneObjectAnnotationsCC, 
                                                                                               oneObjectAnnotationsMF, oneObjectAnnotationsBP) #You should check this!
                            lastObject_ID = object_ID
                            lastObject_symbol = annotation_data[2]
                        oneObjectAnnotationsCC = {}
                        oneObjectAnnotationsBP = {}
                        oneObjectAnnotationsMF = {}
                        if annotation.ontology == "C":
                            #if testAnnotationsCC.has_key(annotation.term_ID):
                            oneObjectAnnotationsCC[annotation.term_ID] = annotation
                            #else:
                                #testAnnotationsCC[annotation.term_ID] = annotation
                        elif annotation.ontology == "F":
                            #if testAnnotationsMF.has_key(annotation.term_ID):
                            oneObjectAnnotationsMF[annotation.term_ID] = annotation
                            #else:
                                #testAnnotationsMF[annotation.term_ID] = annotation
                            #testAnnotationsMF[annotation.term_ID] = annotation
                        else:
                            #if testAnnotationsBP.has_key(annotation.term_ID):
                            oneObjectAnnotationsBP[annotation.term_ID] = annotation
                            #else:
                                #testAnnotationsBP[annotation.term_ID] = annotation
                            #testAnnotationsBP[annotation.term_ID] = annotation
                    else:
                        #Object_ID = annotation_data[4] #GO_ID in Gene Ontology
                        #self.annotations[object_ID] = annotation #Reference an annotation by its term ID, not its gene ID.
                        if annotation.ontology == "C":
                            oneObjectAnnotationsCC[annotation.term_ID] = annotation
                        elif annotation.ontology == "F":
                            oneObjectAnnotationsMF[annotation.term_ID] = annotation
                        else:
                            oneObjectAnnotationsBP[annotation.term_ID] = annotation
                        #lastObject_ID = object_ID
        object_symbol = annotation_data[2]
        relevantDomain = False
        if len(set(oneObjectAnnotationsMF.keys()) & set(self.keyStoneTerms)) >= 1:
                self.annotatedObjects[object_ID] = AnnotatedObject.AnnotatedObject(object_ID, object_symbol, self.fileName, oneObjectAnnotationsCC, 
                                                                               oneObjectAnnotationsMF, oneObjectAnnotationsBP)
        