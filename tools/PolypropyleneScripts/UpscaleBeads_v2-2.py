#Given a set of PMFs for a small bead, e.g. a polymer X-mer, produce a set of PMFs for a larger bead composed of N of these input beads at the specified density
#These can then be used with correction override 0 to make brushes from larger components to reduce the required computat

import numpy as np
import scipy.interpolate as scint
import os
import matplotlib.pyplot as plt


def loadPMF(target , verbose=True):
    PMFData = []
    try:
        pmfText = open(target , "r")
        for line in pmfText:
            if line[0] == "#":
                continue
            if "," in line:
                lineTerms = line.strip().split(",")
            else:
                lineTerms = line.split()
            PMFData.append([float(lineTerms[0]),float(lineTerms[1])])
        pmfText.close()
        foundPMF = 1
        if verbose==True:
            print("Loaded ", target)
        PMFData = np.array(PMFData)
        PMFData[:,1] = PMFData[:,1] - PMFData[-1,1] #set to 0 at end of PMF
        return PMFData
    except: 
        return [-1]
        
def buildRaspberry( innerRadius,outerRadius, beadRadius, makePlane=False):
    beadList = []
    gridPoints =3* (np.round(3.0*outerRadius/(2*np.sqrt(6)*beadRadius))+2).astype(int)
    for i in range(gridPoints):
        for j in range(gridPoints):
            for k in range(gridPoints):
                beadList.append([  beadRadius*(  2*i + ( (j+k)%2   )    ), 
                                 beadRadius*( np.sqrt(3) * (j + 1.0/3.0*(k%2))   ), 
                                 beadRadius*( 2*np.sqrt(6)/3.0 * k)   ])
    beadArray = np.array(beadList)
    beadArray = beadArray - outerRadius
    beadArray[:,0] = beadArray[:,0] - np.mean(beadArray[:,0] )
    beadArray[:,1] = beadArray[:,1] - np.mean(beadArray[:,1] )
    beadArray[:,2] = beadArray[:,2] - np.mean(beadArray[:,2] )
    if makePlane == True:
        beadArray[:,2] = beadArray[:,2] - np.amax( beadArray[:,2] )
        criteria =np.logical_and(   np.logical_and(  np.abs(beadArray[:,0]) < outerRadius+beadRadius  , np.abs(beadArray[:,1]) < outerRadius+beadRadius   )  , np.abs(beadArray[:,1]) < 2*( outerRadius+beadRadius ))
        return beadArray[criteria]
    
    return beadArray[
        np.logical_and(
            beadArray[:,0]**2 + beadArray[:,1]**2 + beadArray[:,2]**2 < outerRadius**2,
            beadArray[:,0]**2 + beadArray[:,1]**2 + beadArray[:,2]**2 >= innerRadius**2 
            )
        
        ]

def buildPlaneBrush( halfWidth, beadRadius, heightInBeads,firstOccupation,finalOccupation):
    firstOccupation = np.clip(firstOccupation, 0, 1.0)
    finalOccupation = np.clip(finalOccupation, 0, firstOccupation)
    beadList = []
    gridPoints =3* (np.round(3.0*outerRadius/(2*np.sqrt(6)*beadRadius))+2).astype(int)
    for i in range(gridPoints):
        for j in range(gridPoints):
            for k in [0]:
                beadList.append([  beadRadius*(  2*i + ( (j+k)%2   )    ), 
                                 beadRadius*( np.sqrt(3) * (j + 1.0/3.0*(k%2))   ), 
                                 beadRadius*( 2*np.sqrt(6)/3.0 * k)   ])
    brushBase = np.array(beadList)


def randomTranslate(beadArr, shiftMag):
    beadsWorking = np.copy(beadArr)
    beadsWorking[:,0] = beadsWorking[:,0] +  2*shiftMag*(np.random.random( len(beadsWorking[:,0] ) )- 0.5 ) #random number between -shiftMag and +shiftMag
    beadsWorking[:,1] = beadsWorking[:,1] +  2*shiftMag*(np.random.random( len(beadsWorking[:,1] ) )- 0.5 ) #random number between -shiftMag and +shiftMag
    beadsWorking[:,2] = beadsWorking[:,2] +  2*shiftMag*(np.random.random( len(beadsWorking[:,2] ) )- 0.5 ) #random number between -shiftMag and +shiftMag
    return beadsWorking

def getSpherePoints(phiDelta,thetaDelta):
    c1Range = np.arange( 0, 360.0 , phiDelta)*np.pi / 180.0
    c2Range = np.arange( thetaDelta, 180.0 - thetaDelta, thetaDelta)*np.pi / 180.0
    return c1Range,c2Range


def getArvoMatrix():
    piVal = np.pi;
    u1= np.random.random();
    u2= np.random.random();
    u3= np.random.random();
    #print(u1,u2,u3)
    
    rxx = -(np.cos(2*piVal*u1)*(1 - 2*u3*np.power(np.cos(2*piVal*u2),2))) - 2*u3*np.cos(2*piVal*u2)*np.sin(2*piVal*u1)*np.sin(2*piVal*u2)
    rxy = -((1 - 2*u3*np.power(np.cos(2*piVal*u2),2))*np.sin(2*piVal*u1)) +    2*u3*np.cos(2*piVal*u1)*np.cos(2*piVal*u2)*np.sin(2*piVal*u2)
    rxz = 2*np.sqrt(1 - u3)*np.sqrt(u3)*np.cos(2*piVal*u2)
    ryx = 2*u3*np.cos(2*piVal*u1)*np.cos(2*piVal*u2)*np.sin(2*piVal*u2) + np.sin(2*piVal*u1)*(1 - 2*u3*np.power(np.sin(2*piVal*u2),2))
    ryy = 2*u3*np.cos(2*piVal*u2)*np.sin(2*piVal*u1)*np.sin(2*piVal*u2) - np.cos(2*piVal*u1)*(1 - 2*u3*np.power(np.sin(2*piVal*u2),2))
    ryz = 2*np.sqrt(1 - u3)*np.sqrt(u3)*np.sin(2*piVal*u2)
    rzx = 2*np.sqrt(1 - u3)*np.sqrt(u3)*np.cos(2*piVal*u1)*np.cos(2*piVal*u2) -    2*np.sqrt(1 - u3)*np.sqrt(u3)*np.sin(2*piVal*u1)*np.sin(2*piVal*u2)
    rzy = 2*np.sqrt(1 - u3)*np.sqrt(u3)*np.cos(2*piVal*u2)*np.sin(2*piVal*u1) +    2*np.sqrt(1 - u3)*np.sqrt(u3)*np.cos(2*piVal*u1)*np.sin(2*piVal*u2)
    rzz = -1 + 2*(1 - u3)
    rotateMatrix = np.array( [[rxx,rxy,rxz],[ryx,ryy,ryz],[rzx,rzy,rzz]])
    #print(rotateMatrix)
    return rotateMatrix

def randomRotate(beads):
    beadsWorking = np.copy(beads)
    rotationMatrix = getArvoMatrix()
    beadsWorking[:,0] = beads[:,0]*rotationMatrix[0,0] + beads[:,1]*rotationMatrix[0,1] + beads[:,2]*rotationMatrix[0,2]
    beadsWorking[:,1] = beads[:,0]*rotationMatrix[1,0] + beads[:,1]*rotationMatrix[1,1] + beads[:,2]*rotationMatrix[1,2]
    beadsWorking[:,2] = beads[:,0]*rotationMatrix[2,0] + beads[:,1]*rotationMatrix[2,1] + beads[:,2]*rotationMatrix[2,2]
    beadsWorking[:,0] = beadsWorking[:,0] - np.mean( beadsWorking[:,0])
    beadsWorking[:,1] = beadsWorking[:,1] - np.mean( beadsWorking[:,1])
    beadsWorking[:,2] = beadsWorking[:,2] - np.mean( beadsWorking[:,2])
    return beadsWorking







inputPMFBeadRadius =  0.302 #radius of the sphere of volume occupied by the bead - ideally from a potential curve to get an actual realistic value, e.g. half-zero crossing distance with itself if this exists
#for a PEG trimer is this approx. 0.205 radius. The input PMFs are PEG-surface to AA-COM , so translating back to COM-COM requires adding this to all distances


#distance to add to the centre-centre distance to get whatever the distance in the PMF means, SCD = CCD + beadShift 
#if for example the PMF is centre-centre distance, this should be zero
#if the PMF is bead surface - centre distance, this should be -bead_radius
#inputBeadShift = -inputPMFBeadRadius
inputBeadShift = 0

makeSphere = False # if makeSphere is true, generate a spherical basis, if not planar
targetOutputRadius = 10 #radius for the output bead
beadDensityScale = 1.0 #probability that a potential bead is used in the final structure 
numSamples = 1 #number of random bead distributions to average over
targetPMFFolderBase = "../surface"
targetMaterial="PP"
pmfsFound = os.listdir(targetPMFFolderBase+"/"+targetMaterial)
inputPMFSet =  []
for pmf in pmfsFound:
    inputPMFSet.append( targetPMFFolderBase+"/"+targetMaterial+"/"+pmf )
inputPMFSet.sort()
#inputPMFSet = [targetPMFFolderBase+"/"+targetMaterial+"/ALA.dat"]
print(inputPMFSet)

#beadCenters = buildRasperry( 0, targetOutputRadius, inputPMFBeadRadius )


for i in range(numSamples):
    allBeadCenterList = []
    beadCenters = np.array( [[ 0, 0, 0  ]] )
    replacedArray = False
    if makeSphere == True:
        print("Making sphere")
        newBeads = randomRotate(  buildRaspberry( 0, targetOutputRadius, inputPMFBeadRadius) )
    else:
        print("Making plane")
        newBeads =  buildRaspberry( 0, targetOutputRadius, inputPMFBeadRadius, makePlane=True) 
    chosenIndex = np.random.random(len(newBeads)) < beadDensityScale
    #selectedBeads = np.reshape( newBeads[chosenIndex], (-1,3) )
    #print(chosenIndex)
    #print(newBeads[chosenIndex] )
    
    if len(newBeads[chosenIndex]) > 0:
        if replacedArray == True:
            beadCenters = np.concatenate( (beadCenters, newBeads[chosenIndex] ) ,axis=0)
        else:
            beadCenters = newBeads[chosenIndex]
            replacedArray = True
    if replacedArray == False:
        print("No beads were successfully generated, please check the density and try again")
        continue
    print("Beads generated: " , len(beadCenters)/numSamples, " per model")
    outputRadius = np.amax(  np.sqrt(beadCenters[:,0]**2 + beadCenters[:,1]**2 + beadCenters[:,2]**2) ) 
    print("Target outer radius: ", targetOutputRadius)
    print("Actual outer radius: ", outputRadius )
    phiRange,thetaRange = getSpherePoints(30,15)
    outputFolder = targetPMFFolderBase +"/"+targetMaterial+"-R"+str(round(outputRadius))+"-d"+ (str(beadDensityScale).replace(".", "p") )+"-n"+str(i)+"_planev2"
    if makeSphere == False:
        outputFolder = targetPMFFolderBase +"/"+targetMaterial+"-plane"+"-d"+ (str(beadDensityScale).replace(".", "p") )+"-n"+str(i)
    os.makedirs(outputFolder,exist_ok=True)



    atomNumRange = np.arange(len(beadCenters))
    xRange = np.arange(-0.5, 0.6, 0.1)
    yRange = np.arange(-0.5,0.6,0.1)
    print(np.amin(beadCenters[:,0]),  np.amax( beadCenters[:,0] ) )
    print(np.amin(beadCenters[:,1]),  np.amax( beadCenters[:,1] ) )
    print(np.amin(beadCenters[:,2]),  np.amax( beadCenters[:,2] ) )

    #print(  np.sqrt( (  beadCenters[1:,0] - beadCenters[0,0] )**2 + (  beadCenters[1:,1] - beadCenters[0,1] )**2 + (  beadCenters[1:,2] - beadCenters[0,2] )**2 ))
    #quit()
    c1grid,c2grid,atomIndexGrid = np.meshgrid(   xRange, yRange,atomNumRange) 
    pointWeights = np.ones_like(c2grid[:,:,0])
    conversionFactor = ( 8.314/1000.0) * 300

    #calculate the nominal surface via carbon probe
    if makeSphere == True:
        rRange = np.arange(outputRadius, outputRadius + 2.0 +0.01,0.01)
    else:
        rRange = np.arange( -0.5 , 1.5+0.01, 0.05 )
    
    print("Calculating surface offset from ALA.dat")
    cresList = []
    pmfDataRaw = loadPMF(targetPMFFolderBase+"/"+targetMaterial+"/ALA.dat")
    pmfData = np.copy( pmfDataRaw )
    print("Initial point: ", pmfData[0] )
    pmfData[:,0] = pmfData[:,0] - inputBeadShift #correct from surf-centre to centre-centre if needed
    print("Shifted point: ", pmfData[0] )
    pmfRes = pmfData[1,0]-pmfData[0,0]
    if pmfData[0,0] > 0.01:
        pmfBackFillRange = np.arange(0.01, pmfData[0,0], 0.01)
        pmfBackFillVal = pmfData[0,1] + (0.2/pmfBackFillRange)**10
        pmfBackFill = np.transpose( np.array( [ pmfBackFillRange, pmfBackFillVal] ) )
        #print(pmfBackFill)
        pmfData = np.concatenate( [pmfBackFill, pmfData] )
        #print(pmfData)
    pmfData2 = pmfData[  np.logical_and( pmfData[:,0] > 0 , pmfData[:,0] < 1.5 ) ]
    pmfData2[:,1] = pmfData2[:,1] - pmfData2[-1,1]
    pmfInterpolation = scint.interp1d( pmfData2[:,0], pmfData2[:,1] , bounds_error=False, fill_value=(50 ,0)  )
    for r in rRange:
        distArray = np.sqrt(   ( c1grid - beadCenters[atomIndexGrid,0])**2 + (c2grid - beadCenters[atomIndexGrid,1])**2 + (r - beadCenters[atomIndexGrid,2])**2  )
        #print(distArray)
        #print( np.amin(distArray), np.amax(distArray) )
        potVals = np.sum(  pmfInterpolation( distArray ), axis=-1 )
        distMinEnergy = np.amin(potVals)
        outputEnergy = distMinEnergy  -conversionFactor * np.log( np.sum(   pointWeights * np.exp( -(potVals-distMinEnergy) / conversionFactor) )   /np.sum(pointWeights)  )
        cresList.append([r, outputEnergy] )

        #distArray = np.sqrt(   ( c1grid - beadCenters[atomIndexGrid,0])**2 + (c2grid - beadCenters[atomIndexGrid,1])**2 + (r - beadCenters[atomIndexGrid,2])**2  )
        #potVals = np.sum(  pmfInterpolation( distArray ), axis=-1 )
        #distMinEnergy = np.amin(potVals)
        #outputEnergy = distMinEnergy  -conversionFactor * np.log( np.sum(   pointWeights * np.exp( -(potVals-distMinEnergy) / conversionFactor) )   /np.sum(pointWeights)  )
        #cresList.append([r, outputEnergy])
        #print(cresList[-1])
    alaSumPMF = np.array(cresList)
  
    #alaSumPMF[ alaSumPMF[:, 1] >  10 ] 

    #Find the point at which the alanine summation-PMF has an energy of 35 kJ/mol, which we then define to be the point at which r = 0.2
    targetEnergy = min(35, np.amax(alaSumPMF[:,1]))
    alignPointIndex = (np.where( np.diff(np.sign(alaSumPMF[:,1] - targetEnergy)))[-1])[0]
    #alignPointIndex =  (np.diff(np.sign(np.diff(alaSumPMF[:,1]))) < 0).nonzero()[0] + 1 
    targetEnergy = alaSumPMF[alignPointIndex, 1]
    print(" set to ", targetEnergy, " at ", alaSumPMF[alignPointIndex, 0] )
    alignPointLoc = alaSumPMF[ alignPointIndex, 0]
    print("Alanine alignment point at: ", alignPointLoc)
    #we want to set distance rprime such that U(rprime = 0.2) = 35 and know r for U(r) = 35
    offsetDist = 0.2 - alignPointLoc #sign convention: this will be negative if the original alignment was > 0.2 , positive if the alignment point was < 0.2 . thus we can add this value to future distances 
    #print(alaSumPMF)
    #print(offsetDist)
    offsetDist = 0.0 #override this 
    alaSumPMF[:,0] = alaSumPMF[:,0] + offsetDist
    #print(alaSumPMF)
    #quit()
    #offsetDist = 0.0
    for inputFile in inputPMFSet:
        pmfChemName = inputFile.split("/")[-1].split(".")[0]
        print("starting input file:", inputFile, pmfChemName)
        pmfDataRaw = loadPMF(inputFile)
        pmfData = np.copy( pmfDataRaw )





        pmfData[:,0] = pmfData[:,0] - inputBeadShift #correct from surf-centre to centre-centre if needed
        pmfRes = pmfData[1,0]-pmfData[0,0]

        if pmfData[0,0] > pmfRes:
            pmfBackFillRange = np.arange(pmfRes, pmfData[0,0], pmfRes)
            pmfBackFillVal = pmfData[0,1] + (0.2/pmfBackFillRange)**10
            pmfBackFill = np.transpose( np.array( [ pmfBackFillRange, pmfBackFillVal] ) )
            #print(pmfBackFill)
            pmfData = np.concatenate( [pmfBackFill, pmfData] )
        pmfData2 = pmfData[  np.logical_and( pmfData[:,0] > 0 , pmfData[:,0] < 1.5 ) ]
        pmfData2[:,1] = pmfData2[:,1] - pmfData2[-1,1]



        #pmfRes = pmfData[1,0]-pmfData[0,0]
        pmfInterpolation = scint.interp1d( pmfData2[:,0], pmfData2[:,1] , bounds_error=False, fill_value=(50 ,0)  )
        #construct the potential from the edge of the new bead out to 2 nm past the edge
        if makeSphere == True:
            rRange = np.arange(outputRadius, outputRadius + 2.0 +pmfRes, pmfRes)
        else:
            rRange = np.arange( -0.5 , 2.0+pmfRes, pmfRes )
        resList = []
        for r in rRange:
            numerator = 0
            denominator = 0
            if makeSphere == True:
                for phi in phiRange:
                    for theta in thetaRange:
                        targetPoint = [ r * np.cos(phi)*np.sin(theta) , r*np.sin(phi)*np.sin(theta), r*np.cos(theta) ]
                        potAtPoint = 0 
                        beadDists =  np.sqrt( np.sum(  (beadCenters - targetPoint)**2, axis=1)) #distance from bead centers to the target point
                        potVals = pmfInterpolation(beadDists )
                        numerator += np.sum(potVals)*np.sin(theta)
                        denominator += np.sin(theta)
                resList.append( [r - outputRadius, numerator/denominator] ) #map back to the surface - centre  distance for the output bead
            else:
                distArray = np.sqrt(   ( c1grid - beadCenters[atomIndexGrid,0])**2 + (c2grid - beadCenters[atomIndexGrid,1])**2 + (r - beadCenters[atomIndexGrid,2])**2  )
                potVals = np.sum(  pmfInterpolation( distArray ), axis=-1 )
                distMinEnergy = np.amin(potVals)
                outputEnergy = distMinEnergy  -conversionFactor * np.log( np.sum(   pointWeights * np.exp( -(potVals-distMinEnergy) / conversionFactor) )   /np.sum(pointWeights)  )
                resList.append([r, outputEnergy] )
                #print("original", resList[-1])
                #print("offset", [ r + offsetDist , outputEnergy])
        #print(np.array(resList))
        res = np.array(resList)    
        res[:,0] = res[:,0] + offsetDist
        res2 = res[   np.logical_and( res[:,0] > 0, res[:,0]< 1.51) ]
        res2[:,1] = res2[:,1] - res2[-1,1] 
        print("Writing to: "+outputFolder+"/"+pmfChemName+".dat" )
        pmfOut = open(outputFolder+"/"+pmfChemName+".dat","w")
        pmfOut.write("#"+targetMaterial+"_"+pmfChemName+"\n")
        pmfOut.write("#Upscaled from " + inputFile+" to radius: "+str(outputRadius)+" bead probs: "+str(beadDensityScale)+"\n")
        pmfOut.write("#surf-centre-dist[nm],U[kJ/mol]\n")
        for line in res2:     
            pmfOut.write(str(line[0])+","+str(line[1])+"\n") 
        pmfOut.close()

