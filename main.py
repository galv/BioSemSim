'''
Created on Jun 11, 2013

@author: s-galvez
'''

import AnnotationSet
import sys
import ObjectSimilarity

import Term
import DAG
import Ontology
import MyClosure
import os.path
import ReactionPathway
import DAGPrinter

import cProfile 
import subprocess
import AnnotatedObject


def newMain():
    evidenceCodesStrings = raw_input("What evidence codes would you like to remove? (Hit enter to use all. Separate evidence codes by commas. Spaces can be used. See http://www.geneontology.org/GO.evidence.shtml)")
    evidenceCodesStrings.replace(" ", "")
    evidenceCodes = evidenceCodesStrings.split(",")
    ontology = Ontology.load("gene_ontology.obo", nodeType = Term.Term, cullObsolete = True)    
    fileName = "gene_association.mgi" #hardcoded. Change later.
    annotationSet = AnnotationSet.AnnotationSet(fileName, ontology, evidenceCodes)
    object_ID = "MGI:87961"
    print annotationSet.annotatedObjects[object_ID].annotationDict["F"].keys()
    domainTerms = annotationSet.annotatedObjects[object_ID].annotationDict["F"].keys()
    c = MyClosure.ForwardClosure(annotationSet, ontology).go(ontology)#CFompute descendants.
    c_rev = MyClosure.ReverseClosure(annotationSet).go(ontology, reversed = True)#, allPaths = True) #Compute ancestors.
    roots = ["C","F", "P"]
    icType = "annotations"
    length = 25
    for root in roots:
            comparisonObject = annotationSet.annotatedObjects[object_ID]
            orderings = createSimilarityMatrix(annotationSet,comparisonObject, root, evidenceCodes, icType)
            comparisonData = compareMethodsSameObject(comparisonObject, orderings,25)
            comparisonData2 = compareMethodsSameObject(comparisonObject, orderings,5)
            printToFile(comparisonObject, orderings, length, root, evidenceCodes, icType,annotationSet,ontology, comparisonData, comparisonData2)
    
def main():
    pr = cProfile.Profile()
    pr.enable()
    ontoChoice = raw_input("Which ontology? 0 for GO, anything else for MP.")
    evidenceCodesStrings = raw_input("What evidence codes would you like to remove? (Hit enter to use all. Separate evidence codes by commas. Spaces can be used. See http://www.geneontology.org/GO.evidence.shtml)")
    evidenceCodesStrings.replace(" ", "")
    evidenceCodes = evidenceCodesStrings.split(",")
    if ontoChoice == "0":
        ontology = Ontology.load("gene_ontology.obo", nodeType = Term.Term, cullObsolete = True)
        fileName = "gene_association.mgi" #hardcoded. Change later.
        annotationSet = AnnotationSet.AnnotationSet(fileName, ontology, evidenceCodes)
        c = MyClosure.ForwardClosure(annotationSet, ontology).go(ontology)#Compute descendants.
        c_rev = MyClosure.ReverseClosure(annotationSet).go(ontology, reversed = True)#, allPaths = True) #Compute ancestors.
    else:
        ontology = Ontology.load("MPheno_OBO.obo", nodeType = Term.Term, cullObsolete = True)
        fileName = "MPannot-2013-07-16.txt"
        diseaseFile = "Geno_11_OMIM.txt"
        annotationSet = AnnotationSet.MPAnnotationSet(fileName, ontology, diseaseFile, evidenceCodes)
        c = MyClosure.ForwardClosure(annotationSet, ontology).go(ontology)#CFompute descendants.
        c_rev = MyClosure.ReverseClosure(annotationSet).go(ontology, reversed = True)#, allPaths = True) #Compute ancestors.
    countDict = {"C": 0, "F": 0, "P": 0}
    
    more = "y"
    while(more == "y"):
        if ontoChoice == "0":
            objectID = "MGI:96840"
            roots = ["C","F", "P"]
        else:
            objectID = raw_input("What genotype would you like to compare? Enter MGI ID: ")
            objectID = "MGI:3797573"
            roots = ["MP"]
        #icChoice = raw_input("Base IC on annotation frequency (Enter 0) or descendant frequency (Enter anything else)? ")
        icChoice = "0"
        if icChoice == "0":
            icType = "annotations"
        else:
            icType = "descendants"
        comparisonObject = annotationSet.annotatedObjects[objectID]
        print "entering createSimilarityMatrix"
        length = 25
        if ontoChoice == "0":
            object_IDs = ["MGI:94909", "MGI:87961","MGI:88057", "MGI:1859546","MGI:96840"]
            object_IDs = ["MGI:88057"]
            #object_IDs = ["MGI:88057", "MGI:1859546","MGI:96840"]
            #object_IDs = ["MGI:1099832"]
        else:
            object_IDs = []
            #disease = "Pulmonary Hypertension, Primary, 1; PPH1"
            #fp = open("Geno_11_OMIM.txt",'r')
            #diseases = fp.readlines()
            #lines = 0
            #for disease_data in diseases:
            #    if disease_data[12:].replace("\n","") == disease:
            #        if not disease_data[:12].replace("\t","") in object_IDs:
            #            object_IDs.append(disease_data[:12].replace("\t",""))
                #lines += 1
                #if len(object_IDs) == :
                #    break
                    
            #print len(object_IDs)
            object_IDs = ["MGI:2651594"]#, "MGI:3713721", "MGI:3714968", "MGI:3809472", "MGI:3814579"]
        #for root in roots:
        
        '''
        fakeID = "MGI:0000000000"
        #annotationStrings =["", fakeID, "", "", ]
        stompedOverUnpythonicAnnotations = annotationSet.annotations #Sets are annoying. Popping removes that element :( Implemented fake gene creation after 
        if ontoChoice == "0":
            annotationDict = {"C": {"GO:0031594": stompedOverUnpythonicAnnotations["GO:0031594"].pop(),
                                    #"GO:0007528": stompedOverUnpythonicAnnotations["GO:0007528"].pop(),
                                    #"GO:0007528": stompedOverUnpythonicAnnotations["GO:0007528"].pop(),
                                    #"GO:0007528": stompedOverUnpythonicAnnotations["GO:0007528"].pop(),
                                    #"GO:0007528": stompedOverUnpythonicAnnotations["GO:0007528"].pop()
                                    },
                              "F": {"GO:0035374": stompedOverUnpythonicAnnotations["GO:0035374"].pop(),
                                    "GO:0043167": stompedOverUnpythonicAnnotations["GO:0043167"].pop(),
                                    "GO:0030548": stompedOverUnpythonicAnnotations["GO:0030548"].pop(),
                                    #"GO:0007528": stompedOverUnpythonicAnnotations["GO:0007528"].pop(),
                                    #"GO:0007528": stompedOverUnpythonicAnnotations["GO:0007528"].pop()
                                    },
                              "P": {"GO:0043113": stompedOverUnpythonicAnnotations["GO:0043113"].pop(),
                                    "GO:0007528": stompedOverUnpythonicAnnotations["GO:0007528"].pop(),
                                    "GO:0048741": stompedOverUnpythonicAnnotations["GO:0048741"].pop(),
                                    "GO:0048747": stompedOverUnpythonicAnnotations["GO:0048747"].pop(),
                                    "GO:0050808": stompedOverUnpythonicAnnotations["GO:0050808"].pop()
                                    }}
            
            
            annotationDict = {"C":  {}, "F": {"GO:0000287": stompedOverUnpythonicAnnotations["GO:0003824"].pop()}, #{"GO:0003824" : stompedOverUnpythonicAnnotations["GO:0003824"].pop()},#, "GO:0005488": stompedOverUnpythonicAnnotations["GO:0005488"].pop(),
                                              #"GO:0016740": stompedOverUnpythonicAnnotations["GO:0016740"].pop(), "GO:0016787": stompedOverUnpythonicAnnotations["GO:0016787"].pop(),
                                              #"GO:0016787": stompedOverUnpythonicAnnotations["GO:0016787"].pop(), "GO:0003674": stompedOverUnpythonicAnnotations["GO:0003674"].pop()}, 
                              "P" : {"GO:0060326": stompedOverUnpythonicAnnotations["GO:0060326"].pop()}}#{"GO:0065007": stompedOverUnpythonicAnnotations["GO:0065007"].pop()}}#,"GO:0044763": stompedOverUnpythonicAnnotations["GO:0044763"].pop(),
                                     #"GO:0044699": stompedOverUnpythonicAnnotations["GO:0044699"].pop(),"GO:0050789": stompedOverUnpythonicAnnotations["GO:0050789"].pop(),
                                     #"GO:0008150": stompedOverUnpythonicAnnotations["GO:0008150"].pop(),"GO:0009987": stompedOverUnpythonicAnnotations["GO:0009987"].pop()}}
            
            
            comparisonObject = AnnotatedObject.AnnotatedObject(fakeID, "dummy", "gene_association.mgi",annotationDict["C"], annotationDict["F"], annotationDict["P"])
            
        else:
            #annotationDict = {#"MP:0000001": stompedOverUnpythonicAnnotations["MP:0000001"].pop()#, "MP:0003631": stompedOverUnpythonicAnnotations["MP:0003631"].pop(),"MP:0005385": stompedOverUnpythonicAnnotations["MP:0005385"].pop(),
            #                  #"MP:0001764": stompedOverUnpythonicAnnotations["MP:0001764"].pop(),"MP:0005376": stompedOverUnpythonicAnnotations["MP:0005376"].pop(),"MP:0003632": stompedOverUnpythonicAnnotations["MP:0003632"].pop(),
            #                  "MP:0005387": stompedOverUnpythonicAnnotations["MP:0005387"].pop()}
            annotationDict = {"MP:0006035": stompedOverUnpythonicAnnotations["MP:0006035"].pop()}
            comparisonObject = AnnotatedObject.AnnotatedObject(fakeID, "dummy", "MPanot-2013-07-16.txt",annotationDict)
        for root in roots:
            orderings = createSimilarityMatrix(annotationSet,comparisonObject, root, evidenceCodes, icType)
            comparisonData = compareMethodsSameObject(comparisonObject, orderings,25)
            printToFile(comparisonObject, orderings, length, root, evidenceCodes, icType,annotationSet,ontology, comparisonData, comparisonData)
        return
        #return
        '''
            
            
            
        '''
        diseaseOutput = ""
        object1 = annotationSet.annotatedObjects[annotationSet.annotatedObjects.keys()[0]]
        for methodName in object1.methodStrings:
            diseaseOutput += "\t" + methodName
        diseaseOutput += "\n"
        diseaseOutputList = []
        '''
        #fp = open("diseaseAnalysis.txt", "w")
        #fp.write(diseaseOutput)
        
        #strainCount  = 0
        for root in roots:
            for object_ID in object_IDs:
                comparisonObject = annotationSet.annotatedObjects[object_ID]
                orderings = createSimilarityMatrix(annotationSet,comparisonObject, root, evidenceCodes, icType)
                comparisonData = compareMethodsSameObject(comparisonObject, orderings,25)
                comparisonData2 = compareMethodsSameObject(comparisonObject, orderings,10)
                printToFile(comparisonObject, orderings, length, root, evidenceCodes, icType,annotationSet,ontology, comparisonData, comparisonData2)#, pathwayDict)
                #diseaseOutputList.append(printDiseaseResults(comparisonObject, orderings, annotationSet, length))
                #fp.write(printDiseaseResults(comparisonObject, orderings, annotationSet, length))
                #strainCount += 1
                #print 1.0*strainCount/len(object_IDs)
        #for stringData in diseaseOutputList:
        #    diseaseOutput += stringData
        
        
        #------
        #Using the DAGPrinter class.
        pr.print_stats()
        pr.disable()
        return
        if ontoChoice == "0":
            more = raw_input("Obtain metrics for another gene? y to continue. anything else to stop.")
        else:
            more = raw_input("Obtain metrics for another genotype? y to continue. anything else to stop.")
            
def geneEnrichment(object_IDs):#, threshold):
#Uses Vlad software to return statistically significant terms froma list of objects.
    
    fp = open('objects.txt', 'w')
    for object_ID in object_IDs:
        fp.write(object_ID + ",")
    subprocess.call(['vlad', '-qFile = objects.txt'])#, '-outDir = ' + os.path.dirname(os.path.realpath(__file__)), '-outRptFormat = text'])
    
def compareMethodsSameObject(comparisonObject, orderings, length):
    comparisons = {}
    for method1 in comparisonObject.methodStrings:
        comparisonOrdering = orderings[method1]
        for method2 in comparisonObject.methodStrings:
            secondOrdering = orderings[method2]
            if method1 == method2:
                comparisons[(method1,method2)] = 1.0*length/length
            else:
                count = 0
                for objectTuple1 in orderings[method1][0:length]:
                    for objectTuple2 in orderings[method2][0:length]:
                        #print objectTuple1[0], objectTuple2[0]
                        if objectTuple1[0] == objectTuple2[0]:
                            count += 1
                comparisons[(method1, method2)] = 1.0*count/length
                #comparisons[(method2, method1)] = 1.0*count/25 
    return comparisons
def printDiseaseResults(object1, orderings, annotationSet, length):
#Return number of genotypes all annnotated to the same disease
    returnString = object1.object_ID + " " + object1.disease
    for orderingName in orderings.keys():
        i = 0
        count = 0
        while i < length:
            if annotationSet.annotatedObjects[orderings[orderingName][i][0]].disease == object1.disease:
                count += 1
            i += 1
        returnString += "\t" + str(1.0*count/length)
    returnString += "\n"
    return returnString

#REturn distance in ranking between two objects according to two different rankings.
    
#Finds top length genes associated with input gene of interest.
def createSimilarityMatrix(annotationSet, object1, root, evidenceCodes, icType):
#Loops through all annotated objects and computes similarity scores with object1.
    simMeasures = {}
    '''
    jaccardSimilarity = ObjectSimilarity.JaccardSimilarity()
    simMeasures[object1.methodStrings[0]] = jaccardSimilarity
    diceSimilarity = ObjectSimilarity.DiceSimilarity()
    simMeasures[object1.methodStrings[1]] = diceSimilarity
    gicSimilarity = ObjectSimilarity.GICSimilarity(icType)
    simMeasures[object1.methodStrings[2]] = gicSimilarity
    jaccardExtSimilarity = ObjectSimilarity.JaccardSimilarityExtended()
    simMeasures[object1.methodStrings[3]] = jaccardExtSimilarity
    resnikSimilarityAvg = ObjectSimilarity.ResnikSimilarityAvg(icType)
    simMeasures[object1.methodStrings[4]] = resnikSimilarityAvg
    resnikSimilarityMax = ObjectSimilarity.ResnikSimilarityMax(icType)
    simMeasures[object1.methodStrings[5]] = resnikSimilarityMax
    wangSimilarity = ObjectSimilarity.WangSimilarity()
    simMeasures[object1.methodStrings[6]] = wangSimilarity
    gicExtSimilarity = ObjectSimilarity.GICSimilarityExtended(icType)
    simMeasures[object1.methodStrings[7]] = gicExtSimilarity
    '''
    
    jaccardSimilarity = ObjectSimilarity.JaccardSimilarity()
    simMeasures[object1.methodStrings[0]] = jaccardSimilarity
    jaccardExtSimilarity = ObjectSimilarity.JaccardSimilarityExtended()
    simMeasures[object1.methodStrings[1]] = jaccardExtSimilarity
    gicExtSimilarity = ObjectSimilarity.GICSimilarityExtended(icType)
    simMeasures[object1.methodStrings[2]] = gicExtSimilarity
    wangSimilarity = ObjectSimilarity.WangSimilarity()
    simMeasures[object1.methodStrings[3]] = wangSimilarity
    
    orderings = {}
    for measureName in object1.methodStrings:
        orderings[measureName] = []
    print "orderings declared"
    #count = 0
    #count = 0
    for key in annotationSet.annotatedObjects.keys():
        for methodName in object1.methodStrings:
            similarity = simMeasures[methodName].findSimilarity(object1, annotationSet.annotatedObjects[key], root)
            orderings[methodName].append((key, similarity, annotationSet.annotatedObjects[key].symbol)) #Appending a tuple here: current object's id, similarity value between object1 and current object, and current object's symbol.
        #count += 1
        #if count  == 20:
            #pr.print_stats
        #    return
        #print 1.0*count/len(annotationSet.annotatedObjects.keys())
    print "out of loop"
    for key in orderings:
        orderings[key].sort(key=lambda x: x[1], reverse = True)
    print "Done sorting"
    
    
    return orderings
    
def printToFile(object1, #the annotated object to which all other objects were compared.
                orderings, #A dict mapping methodStrings in the Annotated Object class to ranked lists of similar genes' tuples (their id, similarity, and symbol) according to the measures described by the strings
                length, #the first length genes are output from orderings (where smallest index corresponds to most similar object).
                root, #"F", "P", or "C" corresponding to one of the three ontologies in the GO. MP has root "MP".
                evidenceCodes, #Evidence codes removed during query.
                icType, 
                annotationSet, 
                ontology,
                comparisonData,
                comparisonData2):
    #Creates tab-delimited file for upload to excel or other program.
    modifier = root
    if root == "MP":
        fileName = object1.object_ID.replace("MGI:","") + "_" + modifier + ".txt"
    else:
        fileName = object1.symbol + "_" + modifier + ".txt"
    codeStr = ""
    for code in evidenceCodes: codeStr += code
    if codeStr == "":
        codeStr = "None"
    if root == "MP":
        dir = os.path.dirname(__file__) +"\Data\\" + object1.object_ID.replace("MGI:","") + "_" + root + "_" + icType + "IC_" + codeStr +"Removed"
    else:
        dir = os.path.dirname(__file__) +"\Data\\" + object1.symbol + "_" + icType + "IC_" + codeStr + "Removed" + "\\" + object1.symbol + "_" + root + "_" + icType + "IC_" + codeStr + "Removed"
    if not os.path.exists(dir):
        os.makedirs(dir)
    print dir
    completeName = os.path.join(dir,fileName)
    print fileName
    print completeName
    #completeName = os.path.join(, fileName)
    fp = open(completeName,'w')
    i = 0
    no_measurements = len(orderings)
    removedCodes = ""
    for code in evidenceCodes:
        removedCodes += code + "\t"
    fp.write(object1.object_ID + "\t" + object1.symbol + "\t" + root +"\t" + "Following codes removed: " + removedCodes + "\t\t" + "\n")
    measurementTypeLine = ""
    for key in object1.methodStrings:
        if root == "MP":
            measurementTypeLine += key + "\t\t\t\t\t\t\t\t\t"
        else:
            measurementTypeLine += key + "\t\t\t\t\t\t\t"
    measurementTypeLine += "\n"
    fp.write(measurementTypeLine)
    header = ""
    if root =="P" or root == "F" or root == "C":
        header += no_measurements*"Rank\tSymbol\tMGI ID\tscore\tNumber of Annotations\tOntology\tColon-Addressed Ontology\t" + "\n"
    else:
        header += no_measurements*"Rank\tSymbol\tMGI ID\tscore\tNumber of Annotations\tDisease\tSameDisease\tOntology\tColon-Addressed Ontology\t" + "\n"

        
    fp.write(header)
    
    
    #---------------------
    
    rankDict = {}
    for methodName in object1.methodStrings:
        rankDict[methodName] = 1
    firstTime = True
    '''
    allRanksExceded = False
    while not allRanksExceded: 
        line = ""
        start = True
        for methodName in object1.methodStrings:
            if not firstTime and not orderings[methodName][i][1] == orderings[methodName][i-1][1]:
                rankDict[methodName] += 1
            line += (str(rankDict[methodName]) + "\t" + orderings[methodName][i][2] + "\t" + orderings[methodName][i][0] + "\t" + str(orderings[methodName][i][1]) + "\t")
            
            if i == length and orderings[methodName][i][1] == orderings[methodName][i+1][1] and start:
                i = i - 1
                start = False
        line += "\n"
        fp.write(line)
        firstTime = False
        i = i + 1
        allRanksExceded = True
        for method in rankDict.keys():
            if rankDict[method] < length:
                allRanksExceded = False
                break
    '''
                
    while i < length:
        line = ""
        start = True
        for methodName in object1.methodStrings:
            #print object1.object_ID            
            if not firstTime and not orderings[methodName][i][1] == orderings[methodName][i-1][1]:
                rankDict[methodName] += 1
            if root == "MP":
                ontologyPath = "LEFT(CELL(\"filename\"),LEN(CELL(\"filename\"))-" + str(len(fileName) + 1) + ")" + "&\"" +"\\" + methodName + "\\" + object1.object_ID.replace("MGI:", "") + "_" + orderings[methodName][i][0].replace("MGI:","") + "_" + methodName + "_" + root + ".pdf\""
                secondOntologyPath = "LEFT(CELL(\"filename\"),LEN(CELL(\"filename\"))-" + str(len(fileName) + 1) + ")" + "&\"" +":" + methodName + ":" + object1.object_ID.replace("MGI:", "") + "_" + orderings[methodName][i][0].replace("MGI:","") + "_" + methodName + "_" + root + ".pdf\""
            else:
                ontologyPath =  "LEFT(CELL(\"filename\"),LEN(CELL(\"filename\"))-" + str(len(fileName) + 1) + ")" + "&\"" +"\\" + methodName + "\\" + object1.symbol + "_" + orderings[methodName][i][2] + "_" + methodName + "_" + root + ".pdf\""
                secondOntologyPath =  "LEFT(CELL(\"filename\"),LEN(CELL(\"filename\"))-" + str(len(fileName) + 1) + ")" + "&\"" +":" + methodName + ":" + object1.symbol + "_" + orderings[methodName][i][2] + "_" + methodName + "_" + root + ".pdf\""
            if root == "MP":
                line += (str(rankDict[methodName]) + "\t" + orderings[methodName][i][2] + "\t" + orderings[methodName][i][0] + "\t" + str(orderings[methodName][i][1]) + "\t" + str(len(annotationSet.annotatedObjects[orderings[methodName][i][0]].annotationDict[root])) + "\t" + annotationSet.annotatedObjects[orderings[methodName][i][0]].disease + "\t" + str(annotationSet.annotatedObjects[orderings[methodName][i][0]].disease == object1.disease) + "\t" + "=HYPERLINK(" + ontologyPath + "," + "\"" + object1.object_ID + " and " + orderings[methodName][i][0] + "\"" + ")" + "\t" + "=HYPERLINK(" + secondOntologyPath + "," + "\"" + object1.object_ID + " and " + orderings[methodName][i][0] + "\"" + ")"  + "\t")

            else:
                line += (str(rankDict[methodName]) + "\t" + orderings[methodName][i][2] + "\t" + orderings[methodName][i][0] + "\t" + str(orderings[methodName][i][1]) + "\t" + str(len(annotationSet.annotatedObjects[orderings[methodName][i][0]].annotationDict[root])) + "\t" + "=HYPERLINK(" + ontologyPath + "," + "\"" + object1.symbol + " and " + orderings[methodName][i][2] + "\"" + ")" + "\t" + "=HYPERLINK(" + secondOntologyPath + "," + "\"" + object1.symbol + " and " + orderings[methodName][i][2] + "\"" + ")"  + "\t")
            #print object1.object_ID
            if i == length and orderings[methodName][i][1] == orderings[methodName][i+1][1] and start:
                i = i - 1
                start = False
        line += "\n"
        fp.write(line)
        firstTime = False
        i = i + 1     
        
        #if i == length:
            #for method in orderings.keys():
                #print method
                #print orderings[method][i][1]
                #print orderings[method][i+1][1]
                #if orderings[method][i+1] != 0.0 and orderings[method][i+1][1] == orderings[method][i][1]:
                #    length += 1
    line = ""    
    line += "\t"
    
    for method in object1.methodStrings:
        line += method + "\t"
    line += "\n"                
    for method1 in object1.methodStrings:
        line += method1 +"\t"
        for method2 in object1.methodStrings:
            try:
                line += str(comparisonData[(method1, method2)]) + "\t"
            except:
                line += "0\t" 
        line += "\n"
    fp.write(line)

    line = ""    
    line += "\t"
    
    for method in object1.methodStrings:
        line += method + "\t"
    line += "\n"                
    for method1 in object1.methodStrings:
        line += method1 +"\t"
        for method2 in object1.methodStrings:
            try:
                line += str(comparisonData2[(method1, method2)]) + "\t"
            except:
                line += "0\t" 
        line += "\n"
    fp.write(line)
    
    fp.close()
    print "Created file called " + fileName
    
    dp = DAGPrinter.DAGPrinter()
    for method in object1.methodStrings:
        i = 0
        dir2 = dir
        path = os.path.join(dir2, method)
        print "path", path
        if not os.path.exists(path):
            os.makedirs(path)
        while i < length:
            textFile = dp.outputGraph(object1, annotationSet.annotatedObjects[orderings[method][i][0]], root, icType, ontology, annotationSet)
            fileName = "testing.dot"
            fp = open(fileName,'w')
            fp.write(textFile)
            fp.close()
            #print textFile
            if root == "MP":
                subprocess.call(['dot', '-Tpdf', fileName, '-o', os.path.join(path,object1.object_ID.replace("MGI:", "") + "_" + orderings[method][i][0].replace("MGI:","") + "_" + method + "_" + root + '.pdf')])
            else:
                subprocess.call(['dot', '-Tpdf', fileName, '-o', os.path.join(path,object1.symbol + "_" + orderings[method][i][2] + "_" + method + "_" + root + '.pdf')])
            i += 1
    
            
        
    #totalDir = "C://Users//s-galvez//Desktop"
    #totalDir = os.path.join(totalDir,fileName)
    #fp = open(totalDir,'w')
    #fp.write(file)
    #fp.close()
def readReactome(file): 
#Inputs the contents of NCBI_Reactome_MGI.txt
#Returns dict mapping MGI ID to pathways id's protein is involved in.
    reactions = {}
    fp = open(file,'r')
    lines = fp.readlines()
    i = 0
    for line in lines:
        line_data = line.split("\t")
        try:
            reactions[line_data[1]].append((line_data[2], line_data[3].replace("\n","")))
        except KeyError:
            reactions[line_data[1]] = []
            reactions[line_data[1]].append((line_data[2], line_data[3].replace("\n","")))
    #for reaction in reactions.keys():
    #    print reactions[reaction][0] + ", " + reactions[reaction][1]
    return reactions
    
def testMFBPResults():
    pass

def doInversions(list1, list2):
    pass
    #http://codereview.stackexchange.com/questions/12922/inversion-count-using-merge-sort

def compareDifferentMethodsLists(orderings, length):
#Find fraction of genes that are within the top length genes that two different genes have in common
    for method1 in orderings.keys():
        for method2 in orderings.keys():
            rankedList1 = orderings[method1]
            rankedList2 = orderings[method2]
            count  = 0
            i = 0
            while i < length:
                for geneTuple in rankedList2:
                    if rankedList1[i][2] == geneTuple[2]:
                        count += 1
                i += 1
            print method1, method2, 1.0*count/length
        
def findShallowObjects(ordering, annotationSet, length, root, method):
#Find the number of objects that are annotated only to the namespace in the top N similar objects for each method.
    #shallowObjects = {}
    #for method in orderings.keys():
    shallowObjects = {}
    index = 0
    while index < length:
        if len(annotationSet.annotatedObjects[ordering[index][0]].annotationDict[root]) == 1:
            for key in annotationSet.annotatedObjects[ordering[index][0]].annotationDict[root].keys():
                annotation = annotationSet.annotatedObjects[ordering[index][0]].annotationDict[root][key]
                if AnnotationSet.AnnotationSet.TERMS[annotation.term_ID].namespace == annotation.term_ID:
                    shallowObjects[annotation.term_ID] = annotationSet.annotatedObjects[ordering[index][0]]
        index += 1
    print method, 1.0*len(shallowObjects)/length 
    return shallowObjects
def findPoorlyAnnotatedObjects(ordering, threshold, annotationSet, length, root, method):
#Find the number of objects that are annotated only threshold times or less.
    poorlyAnnotatedObjects = {} #Dict where method name maps to dictionary mapping object ID to poorly annotated object.
    index = 0
    while index < length:
        if len(annotationSet.annotatedObjects[ordering[index][0]].annotationDict[root]) <= threshold:
            poorlyAnnotatedObjects[ordering[index][0]] = annotationSet.annotatedObjects[ordering[index][0]]
            print annotationSet.annotatedObjects[ordering[index][0]].annotationDict[root]
        index += 1
    print method, 1.0*len(poorlyAnnotatedObjects)/length 
    return poorlyAnnotatedObjects

def findAvgMaxResnikMICAs(orderings, length, root):
    avgResOrdering = orderings['resnikAvg']
    maxResOrdering = orderings['resnikMax']
    
    i = 0
    while i < length:
        termList = []
        objectTuple = avgResOrdering[i]
        
        i +=1


if __name__ == '__main__':
    main()
