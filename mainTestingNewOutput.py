'''
Created on Jun 11, 2013

@author: s-galvez
'''

import AnnotationSet
import sys
import ObjectSimilarity

#sys.path.insert(0, 'DagOntologyLibrary')
import Term
import DAG
import Ontology
import MyClosure

def main():
    exit = "y"
    while(exit == "y"):
        objectID = raw_input("What gene would you like to compare? Enter MGI ID. (See http://www.geneontology.org/) : ")
        #root = raw_input("Compare on only one GO ontology? (F for molecular function, P for biological process, C for cellular component.")
        ontology = Ontology.load("gene_ontology.obo", nodeType = Term.Term)
        c = MyClosure.MyClosure().go(ontology)#CFompute descendants.
        c_rev = MyClosure.MyClosure().go(ontology, reversed=True) #Compute ancestors.    
        fileName = "gene_association.mgi" #hardcoded. Change later.
        annotationSet = AnnotationSet.AnnotationSet(fileName, ontology.nodes)
        comparisonObject = annotationSet.annotatedObjects[objectID]
        createSimilarityMatrix(annotationSet,comparisonObject,11)
        print "Created file called " + comparisonObject.object_ID + ".sim"
        exit = raw_input("Obtain metrics for another gene? y to continue. anything else to stop.")
#Finds top length genes associated with input gene of interest.
def createSimilarityMatrix(annotationSet, object1, length):
    jaccardSimilarity = ObjectSimilarity.JaccardSimilarity()
    diceSimilarity = ObjectSimilarity.DiceSimilarity()
    gicSimilarity = ObjectSimilarity.GICSimilarity()
    jaccardExtendedSimilarity = ObjectSimilarity.JaccardSimilarityExtended()
    ccOrderings = {}
    bpOrderings = {}
    mfOrderings = {}
    ccOrderings["jaccard"] = []
    ccOrderings["dice"] = []
    ccOrderings["gic"] = []
    ccOrderings["resnikAvg"] = []
    ccOrderings["resnikMax"] = []
    
    mfOrderings["jaccard"] = []
    mfOrderings["dice"] = []
    mfOrderings["gic"] = []
    mfOrderings["resnikAvg"] = []
    mfOrderings["resnikMax"] = []
    
    bpOrderings["jaccard"] = []
    bpOrderings["dice"] = []
    bpOrderings["gic"] = []
    bpOrderings["resnikAvg"] = []
    bpOrderings["resnikMax"] = []
    
    for key in annotationSet.annotatedObjects.keys():
        ccRoot = "C"
        annotationSet.annotatedObjects[key].jaccardSim = jaccardSimilarity.findSimilarity(object1, annotationSet.annotatedObjects[key], ccRoot)
        ccOrderings["jaccard"].append((key, annotationSet.annotatedObjects[key].jaccardSim, annotationSet.annotatedObjects[key].symbol))
        
        annotationSet.annotatedObjects[key].diceSim = diceSimilarity.findSimilarity(object1, annotationSet.annotatedObjects[key], ccRoot)
        ccOrderings["dice"].append((key, annotationSet.annotatedObjects[key].diceSim, annotationSet.annotatedObjects[key].symbol))
        
        annotationSet.annotatedObjects[key].gicSim = gicSimilarity.findSimilarity(object1, annotationSet.annotatedObjects[key], ccRoot)
        ccOrderings["gic"].append((key, annotationSet.annotatedObjects[key].gicSim, annotationSet.annotatedObjects[key].symbol))
        
        annotationSet.annotatedObjects[key].resnikAvgSim = gicSimilarity.findSimilarity(object1, annotationSet.annotatedObjects[key], ccRoot)
        ccOrderings["resnikAvg"].append((key, annotationSet.annotatedObjects[key].gicSim, annotationSet.annotatedObjects[key].symbol))
        
        annotationSet.annotatedObjects[key].resnikMaxSim = gicSimilarity.findSimilarity(object1, annotationSet.annotatedObjects[key], ccRoot)
        ccOrderings["resnikMax"].append((key, annotationSet.annotatedObjects[key].gicSim, annotationSet.annotatedObjects[key].symbol))
        #-----------------------------------------------
        mfRoot = "F"
        annotationSet.annotatedObjects[key].jaccardSim = jaccardSimilarity.findSimilarity(object1, annotationSet.annotatedObjects[key], mfRoot)
        mfOrderings["jaccard"].append((key, annotationSet.annotatedObjects[key].jaccardSim, annotationSet.annotatedObjects[key].symbol))
        
        annotationSet.annotatedObjects[key].diceSim = diceSimilarity.findSimilarity(object1, annotationSet.annotatedObjects[key], mfRoot)
        mfOrderings["dice"].append((key, annotationSet.annotatedObjects[key].diceSim, annotationSet.annotatedObjects[key].symbol))
        
        annotationSet.annotatedObjects[key].gicSim = gicSimilarity.findSimilarity(object1, annotationSet.annotatedObjects[key], mfRoot)
        mfOrderings["gic"].append((key, annotationSet.annotatedObjects[key].gicSim, annotationSet.annotatedObjects[key].symbol))
        
        annotationSet.annotatedObjects[key].resnikAvgSim = gicSimilarity.findSimilarity(object1, annotationSet.annotatedObjects[key], mfRoot)
        mfOrderings["resnikAvg"].append((key, annotationSet.annotatedObjects[key].gicSim, annotationSet.annotatedObjects[key].symbol))
        
        annotationSet.annotatedObjects[key].resnikMaxSim = gicSimilarity.findSimilarity(object1, annotationSet.annotatedObjects[key], mfRoot)
        mfOrderings["resnikMax"].append((key, annotationSet.annotatedObjects[key].gicSim, annotationSet.annotatedObjects[key].symbol))
        
        #------------------------------------------------
        bpRoot = "P"
        annotationSet.annotatedObjects[key].jaccardSim = jaccardSimilarity.findSimilarity(object1, annotationSet.annotatedObjects[key], bpRoot)
        bpOrderings["jaccard"].append((key, annotationSet.annotatedObjects[key].jaccardSim, annotationSet.annotatedObjects[key].symbol))
        
        annotationSet.annotatedObjects[key].diceSim = diceSimilarity.findSimilarity(object1, annotationSet.annotatedObjects[key], bpRoot)
        bpOrderings["dice"].append((key, annotationSet.annotatedObjects[key].diceSim, annotationSet.annotatedObjects[key].symbol))
        
        annotationSet.annotatedObjects[key].gicSim = gicSimilarity.findSimilarity(object1, annotationSet.annotatedObjects[key], bpRoot)
        bpOrderings["gic"].append((key, annotationSet.annotatedObjects[key].gicSim, annotationSet.annotatedObjects[key].symbol))
        
        annotationSet.annotatedObjects[key].resnikAvgSim = gicSimilarity.findSimilarity(object1, annotationSet.annotatedObjects[key], bpRoot)
        bpOrderings["resnikAvg"].append((key, annotationSet.annotatedObjects[key].gicSim, annotationSet.annotatedObjects[key].symbol))
        
        annotationSet.annotatedObjects[key].resnikMaxSim = gicSimilarity.findSimilarity(object1, annotationSet.annotatedObjects[key], bpRoot)
        bpOrderings["resnikMax"].append((key, annotationSet.annotatedObjects[key].gicSim, annotationSet.annotatedObjects[key].symbol))
        #annotationSet.annotatedObjects[key].extendedJaccardSim = jaccardExtendedSimilarity.findSimilarity(object1, annotationSet.annotatedObjects[key])
        #jaccardExtendedOrdered.append((key, annotationSet.annotatedObjects[key].extendedJaccardSim, symbol))
        
    print "out of loop"
    for key in ccOrderings:
        ccOrderings[key].sort(key=lambda x: x[1], reverse = True)
    for key in mfOrderings:
        mfOrderings[key].sort(key=lambda x: x[1], reverse = True)
    for key in bpOrderings:
        bpOrderings[key].sort(key=lambda x: x[1], reverse = True)
    allOrderings = {}
    allOrderings["mf"] = mfOrderings
    allOrderings["cc"] = ccOrderings
    allOrderings["bp"] = bpOrderings
    #jaccardExtendedOrdered.sort(key=lambda x: x[1], reverse = True)
    print "Done sorting"
    #jaccardOrdered = jaccardSimilarity.selection(annotationSet, length)
    '''
    print "Molecular function, Jaccard:"
    count = 0
    for gene in allOrderings["mf"]["jaccard"]:
        count +=1
        print count, gene
    '''
    #Outputting data to tab-delimited file:
    printToFile(object1, allOrderings, length)
    
        
def printToFile(object1, allOrderings, length):
    fileName = object1.symbol + ".sim"
    fp = open(fileName,'w')
    i = 0
    no_measurements = len(allOrderings["cc"])
    fp.write(object1.object_ID + "\t" + object1.symbol + "Cellular Component Ontology" +"\t\t\t\t" + "\n")
    measurementTypeLine = "Jaccard\t\t\t\t\t\t\t\t\tDice\t\t\t\t\t\t\t\t\tgic\t\t\t\t\t\t\t\t\tResnik, Avg\t\t\t\t\t\t\t\t\tResnik, Max\t\t\t\t\t\t\t\t\tn"
    fp.write(measurementTypeLine)
    header = ""
    header += no_measurements*"Symbol\tMGI ID\tscore\t" +"\n"
    fp.write(header)
    while i < length:
        jaccard = (allOrderings["cc"]["jaccard"][i][2] + "\t" + allOrderings["cc"]["jaccard"][i][0] + "\t" + str(allOrderings["cc"]["jaccard"][i][1]) +"\t"
                   + allOrderings["bp"]["jaccard"][i][2] + "\t" + allOrderings["bp"]["jaccard"][i][0] + "\t" + str(allOrderings["bp"]["jaccard"][i][1]) +"\t"
                   + allOrderings["mf"]["jaccard"][i][2] + "\t" + allOrderings["mf"]["jaccard"][i][0] + "\t" + str(allOrderings["mf"]["jaccard"][i][1]))
        
        dice = (allOrderings["cc"]["dice"][i][2] + "\t" + allOrderings["cc"]["dice"][i][0] + "\t" + str(allOrderings["cc"]["dice"][i][1]) +"\t"
                   + allOrderings["bp"]["dice"][i][2] + "\t" + allOrderings["bp"]["dice"][i][0] + "\t" + str(allOrderings["bp"]["dice"][i][1]) +"\t"
                   + allOrderings["mf"]["dice"][i][2] + "\t" + allOrderings["mf"]["dice"][i][0] + "\t" + str(allOrderings["mf"]["dice"][i][1]))
        
        gic = (allOrderings["cc"]["gic"][i][2] + "\t" + allOrderings["cc"]["gic"][i][0] + "\t" + str(allOrderings["cc"]["gic"][i][1]) +"\t"
                   + allOrderings["bp"]["gic"][i][2] + "\t" + allOrderings["bp"]["gic"][i][0] + "\t" + str(allOrderings["bp"]["gic"][i][1]) +"\t"
                   + allOrderings["mf"]["gic"][i][2] + "\t" + allOrderings["mf"]["gic"][i][0] + "\t" + str(allOrderings["mf"]["gic"][i][1]))
        
        resnikAvg = (allOrderings["cc"]["resnikAvg"][i][2] + "\t" + allOrderings["cc"]["resnikAvg"][i][0] + "\t" + str(allOrderings["cc"]["resnikAvg"][i][1]) +"\t"
                   + allOrderings["bp"]["resnikAvg"][i][2] + "\t" + allOrderings["bp"]["resnikAvg"][i][0] + "\t" + str(allOrderings["bp"]["resnikAvg"][i][1]) +"\t"
                   + allOrderings["mf"]["resnikAvg"][i][2] + "\t" + allOrderings["mf"]["resnikAvg"][i][0] + "\t" + str(allOrderings["mf"]["resnikAvg"][i][1]))
        
        resnikMax = (allOrderings["cc"]["resnikMax"][i][2] + "\t" + allOrderings["cc"]["resnikMax"][i][0] + "\t" + str(allOrderings["cc"]["resnikMax"][i][1]) +"\t"
                   + allOrderings["bp"]["resnikMax"][i][2] + "\t" + allOrderings["bp"]["resnikMax"][i][0] + "\t" + str(allOrderings["bp"]["resnikMax"][i][1]) +"\t"
                   + allOrderings["mf"]["resnikMax"][i][2] + "\t" + allOrderings["mf"]["resnikMax"][i][0] + "\t" + str(allOrderings["mf"]["resnikMax"][i][1]))
        
        #jaccardExtended =
        line = jaccard + "\t" + dice + "\t" + gic + "\t" + resnikAvg + "\t" +resnikMax + "\n"
        fp.write(line)
        i = i + 1
    fp.close()
    #for key in jaccardOrdered:
    #    print annotationSet.annotatedObjects[key].jaccardSim
    #sorted(.items(), key=lambda x: x[1])
    #Sort
def doInversions(list1, list2):
    pass
        

if __name__ == '__main__':
    main()