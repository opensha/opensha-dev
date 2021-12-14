package scratch.ned.nshm23;

import static com.google.common.base.Preconditions.checkArgument;

import java.awt.geom.Area;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets.RupSetConfig;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.PolygonFaultGridAssociations;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.QuadSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.gui.infoTools.CalcProgressBar;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;
import scratch.UCERF3.enumTreeBranches.SpatialSeisPDF;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.ETAS.ETAS_SimAnalysisTools;
import scratch.UCERF3.erf.ETAS.SeisDepthDistribution;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.utils.MatrixIO;
import scratch.UCERF3.utils.RELM_RegionUtils;
import scratch.UCERF3.utils.U3FaultSystemIO;

/**
 * 
 * Questions: 
 * 
 * move getSectionPolygonRegion(*) method to faultSection; confirm that trace is not offset
 * 
 * 	 * TODO - confirm that trace is surface projection and not offset by DDW
 * 
 * move the following to a more general class: ETAS_SimAnalysisTools.writeMemoryUse()
 * 
 * Does rupSet.getMinMagForSection(s) really get the final minMag?
 * 
 * 
 * 
 * @author field
 *
 */
public class GridSourceProvider2023 {
	
	final static boolean D = true;
	
	final static double DEFAULT_MAX_FAULT_NUCL_DIST = 12d;		
	
	String defaultSectAtCubeCacheFilename = "/tmp/defaultSectAtCubeCache";
	String defaultSectDistForCubeCacheFilename = "/tmp/defaultSectDistForCubeCache";
	String defaultFracCubeUsedBySectCacheFilename = "/tmp/defaultFracCubeUsedBySectCache";
	
	double maxFaultNuclDist;
	
	CubedGriddedRegion cgr;
	
	double[] spatialPDF;
	FaultSystemSolution fss;
	FaultSystemRupSet rupSet;
	HistogramFunction depthNuclProbHist;
	
	IncrementalMagFreqDist totGriddedSeisMFD; // supplied as input

	SummedMagFreqDist totalSubSeisMFD, totalTrulyOffFaultMFD; // both computed

	List<int[]> sectAtCubeList;
	List<float[]> sectDistToCubeList;
	List<float[]> fracCubeUsedBySectList;

	// The following list contains, for each cube, a map of sections and their distance-fraction wts (where
	// the wts represent the fraction of seismicity assinged to the fault section below the min seismo mag).
	ArrayList<HashMap<Integer,Double>> sectDistFractWtMapList;
	// this is the total wt for each section summed from sectDistFractWtMapList (divide the wt directly above
	// by this value to get the nucleation fraction for the section in the associated cube) 
	double[] totSectDistFracWtArray;

	
	/**
	 * 
	 * @param fss
	 * @param griddedRegion
	 * @param spatialPDF
	 * @param totGriddedSeisMFD
	 */
	public GridSourceProvider2023(FaultSystemSolution fss, CubedGriddedRegion cgr, 
			double[] spatialPDF, IncrementalMagFreqDist totGriddedSeisMFD, HistogramFunction depthNuclProbHist) {

		this(fss, cgr, spatialPDF, totGriddedSeisMFD, depthNuclProbHist, DEFAULT_MAX_FAULT_NUCL_DIST);
	}
	
	
	/**
	 * 
	 * @param fss
	 * @param griddedRegion
	 * @param spatialPDF
	 * @param totGriddedSeisMFD
	 * @param maxDepth
	 * @param numCubeDepths
	 * @param numCubesPerGridEdge
	 */
	public GridSourceProvider2023(FaultSystemSolution fss, CubedGriddedRegion cgr, 
			double[] spatialPDF, IncrementalMagFreqDist totGriddedSeisMFD, HistogramFunction depthNuclProbHist, double maxFaultNuclDist) {
		
		this.fss = fss;
		this.cgr = cgr;
		this.rupSet = fss.getRupSet();
		this.spatialPDF = spatialPDF;
		this.totGriddedSeisMFD = totGriddedSeisMFD;
		this.depthNuclProbHist = depthNuclProbHist;
		this.maxFaultNuclDist = maxFaultNuclDist;
		
		// test that spatialPDF sums to 1.0
		double testSum=0;
		for(double val:spatialPDF) testSum += val;
		if(testSum>1.001 || testSum < 0.999)
			throw new RuntimeException("spatialPDF values must sum to 1.0; sum="+testSum);

		
		testSum = depthNuclProbHist.calcSumOfY_Vals();
		if(testSum>1.0001 || testSum < 0.9999)
			throw new RuntimeException("depthNuclProbHist y-axis values must sum to 1.0; sum=testSum");
		// could also check the exact x-axis discretization of depthNuclProbHist
			
		
		readOrGenerateCacheData();
		
		makeSectDistFractWtMapList();
		
		computeTotalOnAndOffFaultGriddedSeisMFDs();
		
		if(D) System.out.println("Done with constructor");
		
	}


	/**
	 * 
	 * @param faultSection
	 * @param distance
	 * @param accountForDip
	 * @return
	 */
	public static Region getSectionPolygonRegion(FaultSection faultSection, double distance, boolean accountForDip) {
		LocationList trace = faultSection.getFaultTrace();
		checkArgument(trace.size() > 1);
		double dipDir = faultSection.getDipDirection();
		double distPlusDip = distance;
		if(accountForDip)
			distPlusDip += faultSection.getOrigDownDipWidth() * Math.cos(faultSection.getAveDip() * Math.PI / 180d);
		LocationList locList = new LocationList();
		locList.add(trace.get(0));
		LocationVector v = new LocationVector(dipDir, distPlusDip, 0);

		for (int i = 0; i < trace.size(); i++) {
			locList.add(LocationUtils.location(trace.get(i), v));
		}
		locList.add(trace.get(trace.size()-1));
		double reverseDipDir = (dipDir + 180) % 360;
		v = new LocationVector(reverseDipDir, distance, 0);
		for (int i = trace.size()-1; i >= 0; i--) {
			locList.add(LocationUtils.location(trace.get(i), v));
		}
		return new Region(locList, BorderType.MERCATOR_LINEAR);
	}
	
	
	/**
	 * 
	 * @param sectIndex
	 */
	private void getCubeDistancesAndFractionsForFaultSection(int sectIndex, 
			HashMap<Integer,Double> cubeDistMap, 
			HashMap<Integer,Double> cubeFracUsedMap) {

		FaultSection fltSection = rupSet.getFaultSectionData(sectIndex);

		Region fltPolygon = getSectionPolygonRegion(fltSection, maxFaultNuclDist, true);
		
//		System.out.println(fltSection.getName()+"\nsectIndex = "+sectIndex+"\ndip = "+fltSection.getAveDip()
//		+"\ndipDir = "+fltSection.getDipDirection()+
//		"\nupSeisDep = "+fltSection.getOrigAveUpperDepth()+
//		"\nddw = "+fltSection.getOrigDownDipWidth()+"\nTrace:\n");
//		for(Location loc:fltSection.getFaultTrace())
//			System.out.println(loc);
//		System.out.println("\nPolygonRegion:\n");
//		for(Location loc:fltPolygon.getBorder())
//			System.out.println(loc);
		
//		System.out.println("\nGriddedSurface:\n");
//		RuptureSurface sectSurf = fltSection.getFaultSurface(1.0, false, false);
//		for(int i=0;i<sectSurf.getEvenlyDiscretizedNumLocs();i++) {
//			System.out.println(sectSurf.getEvenlyDiscretizedLocation(i));
//		}
		
		QuadSurface sectQuadSurf = new QuadSurface(fltSection,false);

		double subDiscFactor = 4; // how much to subdivide cubes for computing ave distance and fact cube used
		// Discretize polygon at 1/10 in lat, lon, and depth for computing average depth and faction of cube used
		GriddedRegion griddedPolygonReg = new GriddedRegion(fltPolygon, cgr.getCubeLatLonSpacing()/subDiscFactor, GriddedRegion.ANCHOR_0_0);
//		System.out.println("\nGriddedRegionPoints:\n");

		for(int i=0;i<griddedPolygonReg.getNumLocations();i++) {
			Location loc = griddedPolygonReg.getLocation(i);
			double depthDiscr = cgr.getCubeDepthDiscr() / subDiscFactor;
			for(double depth = depthDiscr/2;depth<cgr.getMaxDepth();depth+=depthDiscr) {
				Location loc2 = new Location(loc.getLatitude(),loc.getLongitude(),depth);
				double dist = LocationUtils.distanceToSurf(loc2, sectQuadSurf);
				int cubeIndex = cgr.getCubeIndexForLocation(loc2);
				if(dist>maxFaultNuclDist || cubeIndex==-1)
					continue;
				if(cubeDistMap.containsKey(cubeIndex)) {
					Double curVal = cubeDistMap.get(cubeIndex);
					cubeDistMap.replace(cubeIndex,curVal+dist);
					curVal = cubeFracUsedMap.get(cubeIndex);
					cubeFracUsedMap.replace(cubeIndex,curVal+1d);
				}
				else {
					cubeDistMap.put(cubeIndex,dist);
					cubeFracUsedMap.put(cubeIndex,1d);
				}
			}
		}
		for(int key:cubeDistMap.keySet()) {
			double num = cubeFracUsedMap.get(key);
			double aveDist = cubeDistMap.get(key)/num;
			cubeDistMap.replace(key, aveDist);
			cubeFracUsedMap.replace(key, num/(subDiscFactor*subDiscFactor*subDiscFactor));
		}
	}


	
	
	private void readOrGenerateCacheData() {
		
		File sectAtCubeCacheFile = new File(defaultSectAtCubeCacheFilename);
		File sectDistForCubeCacheFile = new File(defaultSectDistForCubeCacheFilename);
		File fracCubeUsedBySectCacheFile = new File(defaultFracCubeUsedBySectCacheFilename);
		
		// make cache files if they don't exist
		if (!sectAtCubeCacheFile.exists() || !sectDistForCubeCacheFile.exists() || !fracCubeUsedBySectCacheFile.exists()) { // read from file if it exists
			if(D) ETAS_SimAnalysisTools.writeMemoryUse("Memory before running generateAndWriteCacheDataToFiles()");
			generateAndWriteCacheDataToFiles();
			if(D) ETAS_SimAnalysisTools.writeMemoryUse("Memory after running generateAndWriteCacheDataToFiles()");
			System.gc(); // garbage collection
		}
		
		// now read them
		try {
			if(D) ETAS_SimAnalysisTools.writeMemoryUse("Memory before reading "+sectAtCubeCacheFile);
				sectAtCubeList = MatrixIO.intArraysListFromFile(sectAtCubeCacheFile);
			if(D) ETAS_SimAnalysisTools.writeMemoryUse("Memory before reading "+sectDistForCubeCacheFile);
				sectDistToCubeList = MatrixIO.floatArraysListFromFile(sectDistForCubeCacheFile);
			if(D) ETAS_SimAnalysisTools.writeMemoryUse("Memory before reading "+fracCubeUsedBySectCacheFile);
				fracCubeUsedBySectList = MatrixIO.floatArraysListFromFile(fracCubeUsedBySectCacheFile);
			if(D) ETAS_SimAnalysisTools.writeMemoryUse("Memory after reading isCubeInsideFaultPolygon");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	private void generateAndWriteCacheDataToFiles() {
		if(D) System.out.println("Starting "+this.getClass().getName()+".generateAndWriteListListDataToFile(); THIS WILL TAKE TIME AND MEMORY!");
		long st = System.currentTimeMillis();
		CalcProgressBar progressBar = null;
		if(D) {
			try {
				progressBar = new CalcProgressBar("Sections to process in generateAndWriteCacheDataToFiles()", "junk");
			} catch (Exception e1) {} // headless			
		}
		ArrayList<ArrayList<Integer>> sectAtCubeListTemp = new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Float>> sectDistToCubeListTemp = new ArrayList<ArrayList<Float>>();
		ArrayList<ArrayList<Float>> fracCubeUsedBySectListTemp = new ArrayList<ArrayList<Float>>();

		int numSect = rupSet.getNumSections();
		for(int i=0; i<cgr.getNumCubes();i++) {
			sectAtCubeListTemp.add(new ArrayList<Integer>());
			sectDistToCubeListTemp.add(new ArrayList<Float>());
			fracCubeUsedBySectListTemp.add(new ArrayList<Float>());
		}
		
		if (progressBar != null) progressBar.showProgress(true);
		for(int sectIndex=0;sectIndex<numSect;sectIndex++) {
			if (progressBar != null) progressBar.updateProgress(sectIndex, numSect);
			
			HashMap<Integer,Double> cubeDistMap = new HashMap<Integer,Double>();
			HashMap<Integer,Double> cubeFracUsedMap = new HashMap<Integer,Double>();
			
			getCubeDistancesAndFractionsForFaultSection(sectIndex, cubeDistMap, cubeFracUsedMap);

			if(cubeDistMap != null) {	// null if section is outside the region
				for(int cubeIndex:cubeDistMap.keySet()) {
					sectAtCubeListTemp.get(cubeIndex).add(sectIndex);
					sectDistToCubeListTemp.get(cubeIndex).add(new Float(cubeDistMap.get(cubeIndex)));
					fracCubeUsedBySectListTemp.get(cubeIndex).add(new Float(cubeFracUsedMap.get(cubeIndex)));
				}			
			}
		}
		
		ETAS_SimAnalysisTools.writeMemoryUse("Memory before writing files");
		File sectAtCubeCacheFile = new File(defaultSectAtCubeCacheFilename);
		File sectDistForCubeCacheFile = new File(defaultSectDistForCubeCacheFilename);
		File fracCubeUsedBySectCacheFile = new File(defaultFracCubeUsedBySectCacheFilename);
		try {
			MatrixIO.intListListToFile(sectAtCubeListTemp,sectAtCubeCacheFile);
			MatrixIO.floatListListToFile(sectDistToCubeListTemp, sectDistForCubeCacheFile);
			MatrixIO.floatListListToFile(fracCubeUsedBySectListTemp, fracCubeUsedBySectCacheFile);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
//System.exit(0);
		
		if (progressBar != null) progressBar.showProgress(false);
		
		if(D) System.out.println(this.getClass().getName()+".generateAndWriteListListDataToFile() took "+(System.currentTimeMillis()-st)/60000+ " min");

	}
	

	public void tempDoSection(int sectIndex) {

		HashMap<Integer,Double> cubeDistMap = new HashMap<Integer,Double>();
		HashMap<Integer,Double> cubeFracUsedMap = new HashMap<Integer,Double>();
		
		getCubeDistancesAndFractionsForFaultSection(sectIndex, cubeDistMap, cubeFracUsedMap);
		
		HashMap<Integer,Double> cubeNuclWtMap = new HashMap<Integer,Double>();

		double sum=0;
		for(int key:cubeDistMap.keySet()) {
			sum += ((maxFaultNuclDist-cubeDistMap.get(key))/maxFaultNuclDist) * cubeFracUsedMap.get(key);
		}
//		System.out.println("sum="+sum);
		for(int key:cubeDistMap.keySet()) {
			double wt = ((maxFaultNuclDist-cubeDistMap.get(key))/maxFaultNuclDist) * cubeFracUsedMap.get(key)/sum;
			cubeNuclWtMap.put(key, wt);
		}


//		for(int key:cubeDistMap.keySet()) {
//			Location loc = getCubeLocationForIndex(key);
//			System.out.println(loc.getLatitude()+"\t"+loc.getLongitude()+"\t"+loc.getDepth()+"\t"+
//			cubeDistMap.get(key).floatValue()+"\t"+cubeFracUsedMap.get(key).floatValue()+"\t"+cubeNuclWtMap.get(key).floatValue());
//		}
	}

	private double getDistWt(double dist) {
		return (maxFaultNuclDist-dist)/maxFaultNuclDist;
	}
	
	/**
	 * This list contains, for each cube, a map of the sections therein and their distance-fraction wts
	 */
	private void makeSectDistFractWtMapList() {
		
		sectDistFractWtMapList = new ArrayList<HashMap<Integer,Double>>();
		totSectDistFracWtArray = new double[rupSet.getNumSections()];
		
		for(int c=0;c<cgr.getNumCubes();c++) {
			HashMap<Integer,Double> sectWtMap = new HashMap<Integer,Double>();
			int numSect = sectAtCubeList.get(c).length;
			for(int i=0;i<numSect;i++) {
				int sectIndex = sectAtCubeList.get(c)[i];
				float dist = sectDistToCubeList.get(c)[i];
				float frac = fracCubeUsedBySectList.get(c)[i];
				double wt = getDistWt(dist)*frac/numSect;
				sectWtMap.put(sectIndex, wt);
				totSectDistFracWtArray[sectIndex] += wt;
			}
			sectDistFractWtMapList.add(sectWtMap);
		}
	}
	
	
	/**
	 * The computes how many different sections nucleate in each cube and then creates a
	 * histogram (how many have 0, 1, 2, etc sections in the cube)
	 */
	public void computeHistogramOfNumSectionsInCubes() {
		int[] numSectAtCubeList = new int[cgr.getNumCubes()];
		HistogramFunction numCubesWithNumSectHist = new HistogramFunction(0.0, 21,1.0);

		for(int c=0; c<cgr.getNumCubes(); c++) {
			numSectAtCubeList[c] = sectAtCubeList.get(c).length;
			numCubesWithNumSectHist.add(numSectAtCubeList[c], 1.0);
			if(numSectAtCubeList[c]==12) {
				System.out.println("\nCube "+c+ " has 12 sections; "+cgr.getCubeLocationForIndex(c));
				for(int i=0;i<sectAtCubeList.get(c).length;i++) {
					int s = sectAtCubeList.get(c)[i];
					float dist = sectDistToCubeList.get(c)[i];
					float frac = fracCubeUsedBySectList.get(c)[i];
					float wt = (float)getDistWt(dist)*frac;
					System.out.println(s+"\t"+dist+"\t"+frac+"\t"+wt+"\t"+rupSet.getFaultSectionData(s).getName());
				}
			}
		}
		System.out.println(numCubesWithNumSectHist);	
	}
	
	/**
	 * this creates a blank (zero y-axis values) MFD with the same discretization as the constructor supplied totGriddedSeisMFD.
	 * @return
	 */
	private SummedMagFreqDist getBlankMFD() {
		return new SummedMagFreqDist(totGriddedSeisMFD.getMinX(), totGriddedSeisMFD.size(),totGriddedSeisMFD.getDelta());
	}
	
	
	private void computeTotalOnAndOffFaultGriddedSeisMFDs() {
		
		totalSubSeisMFD = getBlankMFD();
		totalTrulyOffFaultMFD = getBlankMFD();
		
		double testWt = 0;
		for(int c=0;c<cgr.getNumCubes();c++) {
			SummedMagFreqDist mfd = getSubSeismoMFD_ForCube(c);
			if(mfd != null)
				totalSubSeisMFD.addIncrementalMagFreqDist(mfd);
		}
		
		for(int i=0;i<totGriddedSeisMFD.size();i++) {
			totalTrulyOffFaultMFD.add(i, totGriddedSeisMFD.getY(i) - totalSubSeisMFD.getY(i));
		}
		
		System.out.println("totGriddedSeisMFD:\n"+totGriddedSeisMFD);
		System.out.println("totSubSeisMFD:\n"+totalSubSeisMFD);
		System.out.println("totalTrulyOffFaultMFD:\n"+totalTrulyOffFaultMFD);
		System.out.println("testWt:\n"+testWt);
	}
	
	
	
	public SummedMagFreqDist getSubSeismoMFD_ForCube(int cubeIndex) {
		HashMap<Integer,Double> sectWtMap = sectDistFractWtMapList.get(cubeIndex);
		if(sectWtMap.size()==0) // no sections nucleate here
			return null;
		SummedMagFreqDist subSeisMFD = getBlankMFD();
		int gridIndex = cgr.getRegionIndexForCubeIndex(cubeIndex);
		int depIndex = cgr.getDepthIndexForCubeIndex(cubeIndex);
		for(int s:sectWtMap.keySet()) {
			double wt = sectWtMap.get(s)*spatialPDF[gridIndex]*depthNuclProbHist.getY(depIndex)/(cgr.getNumCubesPerGridEdge()*cgr.getNumCubesPerGridEdge());
			double minMag = rupSet.getMinMagForSection(s);
			double minMagIndex = totGriddedSeisMFD.getClosestXIndex(minMag);
			for(int i=0; i<minMagIndex;i++)
				subSeisMFD.add(i, wt*totGriddedSeisMFD.getY(i));
		}
		return subSeisMFD;
	}
	
	
	
	public SummedMagFreqDist getTrulyOffFaultMFD_ForCube(int cubeIndex) {
		
		double scaleFactor = totGriddedSeisMFD.getY(0)/totalTrulyOffFaultMFD.getY(0);
		
		HashMap<Integer,Double> sectWtMap = sectDistFractWtMapList.get(cubeIndex);
		double wtSum =0;
		for(int s:sectWtMap.keySet()) {
			wtSum+=sectWtMap.get(s);
		}
		SummedMagFreqDist trulyOffMFD = getBlankMFD();
		int gridIndex = cgr.getRegionIndexForCubeIndex(cubeIndex);
		int depIndex = cgr.getDepthIndexForCubeIndex(cubeIndex);
		double wt = (1d-wtSum)*scaleFactor*spatialPDF[gridIndex]*depthNuclProbHist.getY(depIndex)/(cgr.getNumCubesPerGridEdge()*cgr.getNumCubesPerGridEdge());
		
		for(int i=0; i<totalTrulyOffFaultMFD.size();i++)
			trulyOffMFD.add(i, wt*totalTrulyOffFaultMFD.getY(i));

		return trulyOffMFD;
	}
	
	
	public SummedMagFreqDist getSubSeismoMFD_ForGridCell(int gridIndex) {
		SummedMagFreqDist subSeisMFD = getBlankMFD();
		for(int c:cgr.getCubeIndicesForGridCell(gridIndex)) {
			SummedMagFreqDist mfd = getSubSeismoMFD_ForCube(c);
			if(mfd != null)
				subSeisMFD.addIncrementalMagFreqDist(mfd);
		}
		return subSeisMFD;
	}
	
	
	public SummedMagFreqDist getTrulyOffFaultMFD_ForCell(int gridIndex) {
		SummedMagFreqDist subSeisMFD = getBlankMFD();
		for(int c:cgr.getCubeIndicesForGridCell(gridIndex)) {
			subSeisMFD.addIncrementalMagFreqDist(getTrulyOffFaultMFD_ForCube(c));
		}
		return subSeisMFD;
	}


	public SummedMagFreqDist getGriddedSeisMFD_ForCell(int gridIndex) {
		SummedMagFreqDist gridSeisMFD = getBlankMFD();
		gridSeisMFD.addIncrementalMagFreqDist(getSubSeismoMFD_ForGridCell(gridIndex));
		gridSeisMFD.addIncrementalMagFreqDist(getTrulyOffFaultMFD_ForCell(gridIndex));
		return gridSeisMFD;
	}
	
	
	private void testTotalGriddedSeisMFD() {
		SummedMagFreqDist testMFD = getBlankMFD();

		for(int i=0;i<cgr.getGriddedRegion().getNumLocations();i++) {
			testMFD.addIncrementalMagFreqDist(getGriddedSeisMFD_ForCell(i));
		}
		System.out.println("testTotalGriddedSeisMFD():");
		for(int i=0;i<totGriddedSeisMFD.size();i++) {
			System.out.println(totGriddedSeisMFD.getX(i)+"\t"+totGriddedSeisMFD.getY(i)+"\t"+testMFD.getY(i)+"\t"+(float)(testMFD.getY(i)/totGriddedSeisMFD.getY(i)));
		}
		
	}

	private void testTotalTrulyOffFaultGriddedSeisMFD() {
		
		SummedMagFreqDist testMFD = getBlankMFD();
		
		for(int c=0;c<cgr.getNumCubes();c++) {
			SummedMagFreqDist mfd = getTrulyOffFaultMFD_ForCube(c);
			testMFD.addIncrementalMagFreqDist(mfd);
		}
		
		System.out.println("testTotalTrulyOffFaultGriddedSeisMFD():");
		for(int i=0;i<totalTrulyOffFaultMFD.size();i++) {
			System.out.println(totalTrulyOffFaultMFD.getX(i)+"\t"+totalTrulyOffFaultMFD.getY(i)+"\t"+testMFD.getY(i)+"\t"+(float)(testMFD.getY(i)/totalTrulyOffFaultMFD.getY(i)));
		}
		

		
	}

	
	
	public static void main(String[] args) {
		
		String fileName="/Users/field/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip";
		FaultSystemSolution fss;
		try {
			fss = U3FaultSystemIO.loadSol(new File(fileName));
		} catch (Exception e) {
			throw ExceptionUtils.asRuntimeException(e);
		}

		CaliforniaRegions.RELM_TESTING_GRIDDED griddedRegion = RELM_RegionUtils.getGriddedRegionInstance();

		SeisDepthDistribution seisDepthDistribution = new SeisDepthDistribution();
		double delta=2;
		HistogramFunction binnedDepthDistFunc = new HistogramFunction(1d, 12,delta);
		for(int i=0;i<binnedDepthDistFunc.size();i++) {
			double prob = seisDepthDistribution.getProbBetweenDepths(binnedDepthDistFunc.getX(i)-delta/2d,binnedDepthDistFunc.getX(i)+delta/2d);
			binnedDepthDistFunc.set(i,prob);
		}
		System.out.println("Total Depth Prob Sum: "+binnedDepthDistFunc.calcSumOfY_Vals());

		
		double[] spatialPDF = SpatialSeisPDF.UCERF3.getPDF();
		// this sums to 0.9994463999998295; correct it to 1.0
		double sum=0;
		for(double val:spatialPDF) sum+=val;
		for(int i=0;i<spatialPDF.length;i++)
			spatialPDF[i] = spatialPDF[i]/sum;
		
		// Get target total gridded seis MFD
		GridSourceProvider gridSrcProvider = fss.getGridSourceProvider();
		IncrementalMagFreqDist tempMFD = gridSrcProvider.getNodeMFD(0);
		SummedMagFreqDist totGriddedSeisMFD = new SummedMagFreqDist(tempMFD.getMinX(), tempMFD.size(),tempMFD.getDelta());
		for(int i=0;i<gridSrcProvider.size();i++)
			totGriddedSeisMFD.addIncrementalMagFreqDist(gridSrcProvider.getNodeMFD(i));		
//		System.out.println(totGriddedSeisMFD);
//		System.exit(0);
		
		CubedGriddedRegion cgr = new CubedGriddedRegion(griddedRegion);

		GridSourceProvider2023 gridProvider = new GridSourceProvider2023(fss, cgr, spatialPDF, totGriddedSeisMFD, binnedDepthDistFunc);
				
//		gridProvider.testGetCubeIndicesForGridCell();
		
		long startTime = System.currentTimeMillis();

		
		gridProvider.testTotalGriddedSeisMFD();
//		gridProvider.testTotalTrulyOffFaultGriddedSeisMFD();
		
//		gridProvider.computeHistogramOfNumSectionsInCubes();
		
		long runtime = System.currentTimeMillis()-startTime;
		double runtimeMin = runtime/60000d;
		System.out.println("Runtime = "+(float)runtimeMin+" min");
		
//		CalcProgressBar progressBar = null;
//		try {
//			progressBar = new CalcProgressBar("Sections to process", "junk");
//			progressBar.showProgress(true);
//		} catch (Exception e1) {} // headless
//
//		long startTime = System.currentTimeMillis();
//		int numSect = fss.getRupSet().getNumSections();
//		for(int s=0;s<numSect;s++) {
//			if (progressBar != null) progressBar.updateProgress(s, numSect);
//			gridProvider.tempDoSection(s);
//		}
//		long runtime = System.currentTimeMillis()-startTime;
//		double runtimeMin = runtime/60000d;
//		System.out.println("Runtime = "+(float)runtimeMin+" min");
//		
/*
		int sectIndex=0;
		for(int s=0;s<fss.getRupSet().getNumSections();s++) {
//			if(fss.getRupSet().getFaultSectionData(s).getOrigAveUpperDepth() > 3d) {
			if(fss.getRupSet().getFaultSectionData(s).getAveDip() > 89d) {
				sectIndex = s;
				break;
			}
		}

		gridProvider.doSection(sectIndex);

*/

	}

}
