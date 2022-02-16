package scratch.ned.nshm23;

import static com.google.common.base.Preconditions.checkArgument;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.GMT_MapGenerator;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.gui.GMT_MapGuiBean;
import org.opensha.commons.mapping.gmt.gui.ImageViewerWindow;
import org.opensha.commons.param.impl.CPTParameter;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FileUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.QuadSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.gui.infoTools.CalcProgressBar;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import scratch.UCERF3.analysis.GMT_CA_Maps;
import scratch.UCERF3.enumTreeBranches.SpatialSeisPDF;
import scratch.UCERF3.erf.ETAS.ETAS_SimAnalysisTools;
import scratch.UCERF3.erf.ETAS.SeisDepthDistribution;
import scratch.UCERF3.griddedSeismicity.AbstractGridSourceProvider;
import scratch.UCERF3.griddedSeismicity.GridReader;
import scratch.UCERF3.utils.MatrixIO;
import scratch.UCERF3.utils.RELM_RegionUtils;
import scratch.UCERF3.utils.U3FaultSystemIO;
import scratch.ned.FSS_Inversion2019.PlottingUtils;

/**
 * This class represents a grid source provider where, rather than in UCERF3 where a fault represented 
 * all supraseismogenic ruptures inside it's polygon, there is now a linear (ramp) transition between 
 * fault-based ruptures and gridded seismicity events. Thus, fault-based ruptures are most likely to
 * nucleate on the fault surface, but can also nucleate in the vicinity assuming a linear ramp with distance
 * out to the specified maxFaultNuclDist, and in 3D.  Thus, a fault is 50% less likely to nucleate at a distance
 * of maxFaultNuclDist/2 (from the nearest point on the fault surface), and has 0% likelihood of nucleating 
 * beyond maxFaultNuclDist.  Likewise, large gridded-seismicity events taper with the opposite trend near faults,
 * such that their rates are half at maxFaultNuclDist/2 and zero right on the fault surface.  Smaller events (those
 * less than the smallest supraseismogenic magnitudes on faults) have grid cell rates that match the given 
 * spatialPDF exactly (but internally they are apportioned between faults and off-fault ("unassociated") 
 * seismicity according to the same linear ramp. 
 * 
 * The class is backed by a CubedGriddedRegion for purposes of computing distances and associated nucleation rates.  
 * Sources could be given for each cube if higher resolution hazard calculations are desired.
 * 
 * The scaleAllMFDs(double[]) method only scales the gridded-region MFDs and sources (not the cube MFDs)
 * 
 * 
 * To Do:
 * 
 * 1) Finalize how cached files are handled
 *  * 
 * 2) Move getSectionPolygonRegion(*) method to faultSection? & confirm that trace is not offset
 *    when UpperSeisDepth != 0.
 *    
 * 3) Make sure computeLongTermSupraSeisMFD_OnSectArray() from fss has the final Mmin applied.
 * 
 * 4) QuadSurface - do we have analytical dist to a loc that is not on the earth surface?
 * 
 * 5) Move ETAS_SimAnalysisTools.writeMemoryUse() to a more general class/location
 * 
 * 6) improve input of fracStrikeSlip,fracNormal, and fracReverse (e.g., bundle and
 *    include fixed rake option?
 * 
 * 
 * 
 * @author field
 *
 */
public class GridSourceProvider2023 extends AbstractGridSourceProvider {
	
	final static boolean D = true;
	
	final static double DEFAULT_MAX_FAULT_NUCL_DIST = 12d;		
	
	String defaultSectAtCubeCacheFilename = "/Users/field/tmp/defaultSectAtCubeCache";
	String defaultSectDistForCubeCacheFilename = "/Users/field/tmp/defaultSectDistForCubeCache";
	String defaultTotDistWtsAtCubesForSectArrayFilename = "/Users/field/tmp/defaultTotDistWtsAtCubesForSectArrayCache";
	String defaultSectionsThatNucleateOutsideRegionListFilename = "/Users/field/tmp/defaultSectionsThatNucleateOutsideRegionListCache";
	
	double maxFaultNuclDist;
	
	CubedGriddedRegion cgr;
	
	double[] spatialPDF;
	FaultSystemSolution fss;
	FaultSystemRupSet rupSet;
	HistogramFunction depthNuclProbHist;
	GriddedRegion griddedRegion;
	
	IncrementalMagFreqDist totGriddedSeisMFD; // supplied as input

	SummedMagFreqDist totalSubSeisOnFaultMFD, totalTrulyOffFaultMFD, totalSupraSeisOnFaultMFD; // all computed
	
	SummedMagFreqDist[] subSeisOnFaultMFD_ForGridArray, unassociatedMFD_ForGridArray;

	List<int[]> sectAtCubeList;
	List<float[]> sectDistToCubeList;

	// The following list contains, for each cube, a map of sections and their distance-fraction wts (where
	// the wts represent the fraction of seismicity assinged to the fault section below the min seismo mag).
	ArrayList<HashMap<Integer,Double>> sectDistWtMapAtCubeList;
	
	// this is the total wt for each section summed from sectDistFractWtMapList (divide the wt directly above
	// by this value to get the nucleation fraction for the section in the associated cube) 
	double[] totDistWtsAtCubesForSectArray;
	
	IncrementalMagFreqDist[] longTermSupraSeisMFD_OnSectArray;

	ArrayList<Integer> sectionsThatNucleateOutsideRegionList;
	
	private static double[] fracStrikeSlip,fracNormal,fracReverse;

	
	/**
	 * 
	 * @param fss
	 * @param cgr
	 * @param spatialPDF
	 * @param totGriddedSeisMFD
	 * @param depthNuclProbHist
	 */
	public GridSourceProvider2023(FaultSystemSolution fss, CubedGriddedRegion cgr, 
			double[] spatialPDF, IncrementalMagFreqDist totGriddedSeisMFD, HistogramFunction depthNuclProbHist,
			double[] fracStrikeSlip, double[] fracNormal, double[] fracReverse) {

		this(fss, cgr, spatialPDF, totGriddedSeisMFD, depthNuclProbHist, fracStrikeSlip, fracNormal, fracReverse, DEFAULT_MAX_FAULT_NUCL_DIST);
	}
	
	/**
	 * 
	 * @param fss
	 * @param cgr
	 * @param spatialPDF
	 * @param totGriddedSeisMFD
	 * @param depthNuclProbHist
	 * @param maxFaultNuclDist
	 */
	public GridSourceProvider2023(FaultSystemSolution fss, CubedGriddedRegion cgr, 
			double[] spatialPDF, IncrementalMagFreqDist totGriddedSeisMFD, HistogramFunction depthNuclProbHist, 
			double[] fracStrikeSlip, double[] fracNormal, double[] fracReverse, double maxFaultNuclDist) {
		
		this.fss = fss;
		this.cgr = cgr;
		this.rupSet = fss.getRupSet();
		this.spatialPDF = spatialPDF;
		this.totGriddedSeisMFD = totGriddedSeisMFD;
		this.depthNuclProbHist = depthNuclProbHist;
		this.maxFaultNuclDist = maxFaultNuclDist;
		this.griddedRegion = cgr.getGriddedRegion();
		this.fracStrikeSlip = fracStrikeSlip;
		this.fracNormal = fracNormal;
		this.fracReverse = fracReverse;
		
		if(griddedRegion.getNodeCount() != spatialPDF.length)
			throw new RuntimeException("griddedRegion and spatialPDF have differe sizes: "+griddedRegion.getNodeCount()+" vs "+spatialPDF.length);
		
		sectionsThatNucleateOutsideRegionList = new ArrayList<Integer>();
		
		// test some things
		double testSum=0;
		for(double val:spatialPDF) testSum += val;
		if(testSum>1.001 || testSum < 0.999)
			throw new RuntimeException("spatialPDF values must sum to 1.0; sum="+testSum);
		
		if(spatialPDF.length != griddedRegion.getNumLocations())
			throw new RuntimeException("spatialPDF and griddedRegion must be the same size; they are "+spatialPDF.length+
					" and "+griddedRegion.getNumLocations()+",respectively");
		
		testSum = depthNuclProbHist.calcSumOfY_Vals();
		if(testSum>1.0001 || testSum < 0.9999)
			throw new RuntimeException("depthNuclProbHist y-axis values must sum to 1.0; sum=testSum");
		// could also check the exact x-axis discretization of depthNuclProbHist
			

		// generate or read the sections and their distances for cubes 
		long time = System.currentTimeMillis();
		readOrGenerateCacheData();
		long runtime = System.currentTimeMillis()-time;
		if(D) System.out.println("readOrGenerateCacheData() took "+(runtime/1000)+" seconds");
		
		// make the distWtMap for each cube
		time = System.currentTimeMillis();
		makeSectDistWtMapList();
		runtime = System.currentTimeMillis()-time;
		if(D) System.out.println("makeSectDistWtMapList() took "+(runtime/1000)+" seconds");
		
		// compute total MFDs
		time = System.currentTimeMillis();
		computeTotalOnAndOffFaultGriddedSeisMFDs();
		runtime = System.currentTimeMillis()-time;
		if(D) System.out.println("computeTotalOnAndOffFaultGriddedSeisMFDs() took "+(runtime/1000)+" seconds");

		
		// compute total MFDs
		time = System.currentTimeMillis();
		computeLongTermSupraSeisMFD_OnSectArray();
		runtime = System.currentTimeMillis()-time;
		if(D) System.out.println("computeLongTermSupraSeisMFD_OnSectArray() took "+(runtime/1000)+" seconds");
		
		// compute grid cell MFDs
		time = System.currentTimeMillis();
		subSeisOnFaultMFD_ForGridArray = new SummedMagFreqDist[spatialPDF.length];
		unassociatedMFD_ForGridArray = new SummedMagFreqDist[spatialPDF.length];
		for(int i=0;i<spatialPDF.length;i++) {
			subSeisOnFaultMFD_ForGridArray[i] = computeMFD_SubSeisOnFault(i);
			unassociatedMFD_ForGridArray[i] = computeMFD_Unassociated(i);
		}
		runtime = System.currentTimeMillis()-time;
		if(D) System.out.println("computing grid MFDs took "+(runtime/1000)+" seconds");
	
		if(D) System.out.println("Done with constructor");
		
	}


	/**
	 * This assumes the trace represents a surface projection of top edge of surface (and not projected
	 * up-dip when upper seismogenic depth is nonzero)
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
	 * @param cubeDistMap
	 * @return
	 */
	public double getCubeDistancesForFaultSection(int sectIndex, HashMap<Integer,Double> cubeDistMap) {

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
		
		RuptureSurface sectSurf = fltSection.getFaultSurface(0.25);
//		QuadSurface sectSurf = new QuadSurface(fltSection,false);

		GriddedRegion griddedPolygonReg = new GriddedRegion(fltPolygon, cgr.getCubeLatLonSpacing(), cgr.getCubeLocationForIndex(0));
		double totWt = 0;
		double testWt=0;
		for(int i=0;i<griddedPolygonReg.getNumLocations();i++) {
			Location loc = griddedPolygonReg.getLocation(i);
			double depthDiscr = cgr.getCubeDepthDiscr();
			for(double depth = depthDiscr/2;depth<cgr.getMaxDepth();depth+=depthDiscr) {
				Location loc2 = new Location(loc.getLatitude(),loc.getLongitude(),depth);
				double dist = LocationUtils.distanceToSurf(loc2, sectSurf);
				if(dist<=maxFaultNuclDist) { 
					totWt += getDistWt(dist); // this will includes cubes outside the region where the section could nucleate
					int cubeIndex = cgr.getCubeIndexForLocation(loc2);
					if(cubeIndex>=0) {// make sure it's in the region
						cubeDistMap.put(cubeIndex,dist);
						testWt += getDistWt(dist);
					}
				}
			}
		}

		float ratio = (float)testWt/(float)totWt;
		if(ratio != 1.0) {
			sectionsThatNucleateOutsideRegionList.add(sectIndex);
			if(D) System.out.println((1f-ratio)+" of "+rupSet.getFaultSectionData(sectIndex).getName()+ " nucleates outside the region");
		}
		return totWt;
	}


	
	
	private void readOrGenerateCacheData() {
		
		File sectAtCubeCacheFile = new File(defaultSectAtCubeCacheFilename);
		File sectDistForCubeCacheFile = new File(defaultSectDistForCubeCacheFilename);
		File totDistWtsAtCubesForSectArrayCacheFile = new File(defaultTotDistWtsAtCubesForSectArrayFilename);
		File sectionsThatNucleateOutsideRegionListFilename = new File(defaultSectionsThatNucleateOutsideRegionListFilename);

		
		// make cache files if they don't exist
		if (!sectAtCubeCacheFile.exists() || !sectDistForCubeCacheFile.exists() || !totDistWtsAtCubesForSectArrayCacheFile.exists()) { // read from file if it exists
			if(D) ETAS_SimAnalysisTools.writeMemoryUse("Memory before running generateAndWriteCacheDataToFiles()");
			generateAndWriteCacheDataToFiles();
			if(D) ETAS_SimAnalysisTools.writeMemoryUse("Memory after running generateAndWriteCacheDataToFiles()");
			System.gc(); // garbage collection
		}
		
		// now read them
		int[] tempIntArray = null;
		try {
			if(D) ETAS_SimAnalysisTools.writeMemoryUse("Memory before reading files "+sectAtCubeCacheFile);
			sectAtCubeList = MatrixIO.intArraysListFromFile(sectAtCubeCacheFile);
			sectDistToCubeList = MatrixIO.floatArraysListFromFile(sectDistForCubeCacheFile);
			totDistWtsAtCubesForSectArray = MatrixIO.doubleArrayFromFile(totDistWtsAtCubesForSectArrayCacheFile);
			tempIntArray = MatrixIO.intArrayFromFile(sectionsThatNucleateOutsideRegionListFilename);
			if(D) ETAS_SimAnalysisTools.writeMemoryUse("Memory after reading files");
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		sectionsThatNucleateOutsideRegionList = new ArrayList<Integer>();
		for(int i:tempIntArray)
			sectionsThatNucleateOutsideRegionList.add(i);
		
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
		double[] totDistWtsForSectArray = new double[rupSet.getNumSections()];

		for(int i=0; i<cgr.getNumCubes();i++) {
			sectAtCubeListTemp.add(new ArrayList<Integer>());
			sectDistToCubeListTemp.add(new ArrayList<Float>());
		}
		
		if (progressBar != null) progressBar.showProgress(true);
		int numSect = rupSet.getNumSections();
		for(int sectIndex=0;sectIndex<numSect;sectIndex++) {
			if (progressBar != null) progressBar.updateProgress(sectIndex, numSect);
			
			HashMap<Integer,Double> cubeDistMap = new HashMap<Integer,Double>();
			
			// this will fill in cubeDistMap:
			totDistWtsForSectArray[sectIndex] = getCubeDistancesForFaultSection(sectIndex, cubeDistMap);
// System.out.println(+totDistWtsForSectArray[sectIndex]+" for "+rupSet.getFaultSectionData(sectIndex).getName());

			if(cubeDistMap != null) {	// null if section is outside the region
				for(int cubeIndex:cubeDistMap.keySet()) {
					sectAtCubeListTemp.get(cubeIndex).add(sectIndex);
					sectDistToCubeListTemp.get(cubeIndex).add(new Float(cubeDistMap.get(cubeIndex)));
				}			
			}
		}
		
		ETAS_SimAnalysisTools.writeMemoryUse("Memory before writing files");
		
		File sectAtCubeCacheFile = new File(defaultSectAtCubeCacheFilename);
		File sectDistForCubeCacheFile = new File(defaultSectDistForCubeCacheFilename);
		File totDistWtsAtCubesForSectArrayCacheFile = new File(defaultTotDistWtsAtCubesForSectArrayFilename);
		File sectionsThatNucleateOutsideRegionListFilename = new File(defaultSectionsThatNucleateOutsideRegionListFilename);
		
		int[] tempIntArray = new int[sectionsThatNucleateOutsideRegionList.size()];
		for(int i=0;i<tempIntArray.length;i++)
			tempIntArray[i] = sectionsThatNucleateOutsideRegionList.get(i);
		try {
			MatrixIO.intListListToFile(sectAtCubeListTemp,sectAtCubeCacheFile);
			MatrixIO.floatListListToFile(sectDistToCubeListTemp, sectDistForCubeCacheFile);
			MatrixIO.doubleArrayToFile(totDistWtsForSectArray, totDistWtsAtCubesForSectArrayCacheFile);
			MatrixIO.intArrayToFile(tempIntArray, sectionsThatNucleateOutsideRegionListFilename);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
//System.exit(0);
		
		if (progressBar != null) progressBar.showProgress(false);
		
		if(D) System.out.println(this.getClass().getName()+".generateAndWriteListListDataToFile() took "+(System.currentTimeMillis()-st)/60000+ " min");

	}
	


	public double getDistWt(double dist) {
		return (maxFaultNuclDist-dist)/maxFaultNuclDist;
	}
	

	
	
	/**
	 * This list contains, for each cube, a map of the sections therein and their distance wts.  
	 * If multiple sections exist in a cube and their wts exceed 1.0, they are all normalized by
	 * the sum so they add up to 1.0.
	 */
	private void makeSectDistWtMapList() {
		
		sectDistWtMapAtCubeList = new ArrayList<HashMap<Integer,Double>>();
		
		for(int c=0;c<cgr.getNumCubes();c++) {
			HashMap<Integer,Double> sectWtMap = new HashMap<Integer,Double>();
			int numSect = sectAtCubeList.get(c).length;
			double[] prelimWts = new double[numSect];
			double wtSum = 0;
			for(int s=0;s<numSect;s++) {
				float dist = sectDistToCubeList.get(c)[s];
				prelimWts[s] = getDistWt(dist);
				wtSum += prelimWts[s];
			}
			double[] finalWts = prelimWts;
			if(wtSum>1.0) {
				finalWts = new double[numSect];
				for(int s=0;s<numSect;s++) {
					finalWts[s] = prelimWts[s]/wtSum;
					int sectIndex = sectAtCubeList.get(c)[s];
					totDistWtsAtCubesForSectArray[sectIndex] += (finalWts[s]-prelimWts[s]); // reduce this for the wt reduction here
				}
			}

			for(int s=0;s<numSect;s++) {
				int sectIndex = sectAtCubeList.get(c)[s];
				double wt = finalWts[s];
				sectWtMap.put(sectIndex, wt);
			}
			sectDistWtMapAtCubeList.add(sectWtMap);
		}
	}

		
	
	/**
	 * this creates a blank (zero y-axis values) MFD with the same discretization as the constructor supplied totGriddedSeisMFD.
	 * @return
	 */
	public SummedMagFreqDist initSummedMFD() {
		return new SummedMagFreqDist(totGriddedSeisMFD.getMinX(), totGriddedSeisMFD.size(),totGriddedSeisMFD.getDelta());
	}
	
	
	private void computeTotalOnAndOffFaultGriddedSeisMFDs() {
		
		totalSubSeisOnFaultMFD = initSummedMFD();
		totalSubSeisOnFaultMFD.setName("totalSubSeisMFD");
		totalTrulyOffFaultMFD = initSummedMFD();
		totalTrulyOffFaultMFD.setName("totalTrulyOffFaultMFD");
		
		for(int c=0;c<cgr.getNumCubes();c++) {
			SummedMagFreqDist mfd = getSubSeismoMFD_ForCube(c);
			if(mfd != null)
				totalSubSeisOnFaultMFD.addIncrementalMagFreqDist(mfd);
		}
		
		for(int i=0;i<totGriddedSeisMFD.size();i++) {
			totalTrulyOffFaultMFD.add(i, totGriddedSeisMFD.getY(i) - totalSubSeisOnFaultMFD.getY(i));
		}
//		if(D) {
//			System.out.println("totGriddedSeisMFD:\n"+totGriddedSeisMFD);
//			System.out.println("totGriddedSeisMFD Cumulative::\n"+totGriddedSeisMFD.getCumRateDistWithOffset());
//			System.out.println("totSubSeisMFD:\n"+totalSubSeisOnFaultMFD);
//			System.out.println("totalTrulyOffFaultMFD:\n"+totalTrulyOffFaultMFD);
//		}
	}
	
	/**
	 * Need to figure out how to compute this when fss has the module: ModSectMinMags.
	 * 
	 * Note that this includes any fault sections that are outside the region
	 */
	private void computeLongTermSupraSeisMFD_OnSectArray() {
		
		SummedMagFreqDist mfd = initSummedMFD();

		// this didn't work so use ERF to get section mdfs
//		ModSectMinMags mod = fss.getModule(ModSectMinMags.class);
		
		longTermSupraSeisMFD_OnSectArray = new IncrementalMagFreqDist[rupSet.getNumSections()];
		for(int s=0;s<rupSet.getNumSections();s++) {
			IncrementalMagFreqDist nuclMFD = fss.calcNucleationMFD_forSect(s, mfd.getMinX(), mfd.getMaxX(), mfd.size());
			longTermSupraSeisMFD_OnSectArray[s] = nuclMFD;
			mfd.addIncrementalMagFreqDist(nuclMFD);
		}
		mfd.setName("totalSupraSeisOnFaultMFD");
		totalSupraSeisOnFaultMFD=mfd;
		
	}

	
	
	
	public SummedMagFreqDist getSubSeismoMFD_ForCube(int cubeIndex) {
		HashMap<Integer,Double> sectWtMap = sectDistWtMapAtCubeList.get(cubeIndex);
		if(sectWtMap.size()==0) // no sections nucleate here
			return null;
		SummedMagFreqDist subSeisMFD = initSummedMFD();
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
	
	
	
	public SummedMagFreqDist getUnassociatedMFD_ForCube(int cubeIndex) {
		
		double scaleFactor = totGriddedSeisMFD.getY(0)/totalTrulyOffFaultMFD.getY(0);
		
		HashMap<Integer,Double> sectWtMap = sectDistWtMapAtCubeList.get(cubeIndex);
		double wtSum =0;
		for(int s:sectWtMap.keySet()) {
			wtSum+=sectWtMap.get(s);
		}
		SummedMagFreqDist trulyOffMFD = initSummedMFD();
		int gridIndex = cgr.getRegionIndexForCubeIndex(cubeIndex);
		int depIndex = cgr.getDepthIndexForCubeIndex(cubeIndex);
		double wt = (1d-wtSum)*scaleFactor*spatialPDF[gridIndex]*depthNuclProbHist.getY(depIndex)/(cgr.getNumCubesPerGridEdge()*cgr.getNumCubesPerGridEdge());
		
		for(int i=0; i<totalTrulyOffFaultMFD.size();i++)
			trulyOffMFD.add(i, wt*totalTrulyOffFaultMFD.getY(i));

		return trulyOffMFD;
	}
	
	public SummedMagFreqDist getGriddedSeisMFD_ForCube(int cubeIndex) {
		SummedMagFreqDist cubeMFD = initSummedMFD();
		SummedMagFreqDist mfd = getSubSeismoMFD_ForCube(cubeIndex);
		if(mfd != null)
			cubeMFD.addIncrementalMagFreqDist(mfd);
		mfd = getUnassociatedMFD_ForCube(cubeIndex);
		if(mfd != null)
			cubeMFD.addIncrementalMagFreqDist(mfd);
		return cubeMFD;
	}
	
	
	
	public SummedMagFreqDist getSupraSeisMFD_ForCube(int cubeIndex) {

		HashMap<Integer,Double> sectWtMap = sectDistWtMapAtCubeList.get(cubeIndex);
		if(sectWtMap.size()==0) // no sections nucleate here
			return null;
		SummedMagFreqDist supraMFD = initSummedMFD();
		for(int s:sectWtMap.keySet()) {
			double wt = sectWtMap.get(s)/totDistWtsAtCubesForSectArray[s];
			IncrementalMagFreqDist mfd = longTermSupraSeisMFD_OnSectArray[s].deepClone();
			mfd.scale(wt);
			supraMFD.addIncrementalMagFreqDist(mfd);
		}
		return supraMFD;
	}
	
	public SummedMagFreqDist getTotalMFD_ForCube(int cubeIndex) {
		SummedMagFreqDist cubeMFD = initSummedMFD();
		SummedMagFreqDist mfd = getGriddedSeisMFD_ForCube(cubeIndex);
		if(mfd != null)
			cubeMFD.addIncrementalMagFreqDist(mfd);
		mfd = getSupraSeisMFD_ForCube(cubeIndex);
		if(mfd != null)
			cubeMFD.addIncrementalMagFreqDist(mfd);
		return cubeMFD;
	}



	
	@Override
	public SummedMagFreqDist getMFD_SubSeisOnFault(int gridIndex) {
		return subSeisOnFaultMFD_ForGridArray[gridIndex];
	}
	
	@Override
	public SummedMagFreqDist getMFD_Unassociated(int gridIndex) {
		return unassociatedMFD_ForGridArray[gridIndex];
	}
	
	
	private SummedMagFreqDist computeMFD_SubSeisOnFault(int gridIndex) {
		SummedMagFreqDist summedMFD = initSummedMFD();
		for(int c:cgr.getCubeIndicesForGridCell(gridIndex)) {
			SummedMagFreqDist mfd = getSubSeismoMFD_ForCube(c);
			if(mfd != null)
				summedMFD.addIncrementalMagFreqDist(mfd);
		}
		if(summedMFD.getTotalIncrRate() < 1e-10)
			summedMFD=null;
		return summedMFD;
	}
	
	
	private SummedMagFreqDist computeMFD_Unassociated(int gridIndex) {
		SummedMagFreqDist summedMFD = initSummedMFD();
		for(int c:cgr.getCubeIndicesForGridCell(gridIndex)) {
			summedMFD.addIncrementalMagFreqDist(getUnassociatedMFD_ForCube(c));
		}
		if(summedMFD.getTotalIncrRate() < 1e-10)
			summedMFD=null;
		return summedMFD;
	}

	
	


//	public SummedMagFreqDist getMFD(int gridIndex) {
//		SummedMagFreqDist gridSeisMFD = initSummedMFD();
//		gridSeisMFD.addIncrementalMagFreqDist(getMFD_SubSeisOnFault(gridIndex));
//		gridSeisMFD.addIncrementalMagFreqDist(getMFD_Unassociated(gridIndex));
//		return gridSeisMFD;
//	}
	
	/**
	 * this returns the total sub-seismogenic on-fault MFD
	 * @return
	 */
	public SummedMagFreqDist getTotalSubSeisOnFaultMFD() {
		return totalSubSeisOnFaultMFD;
	}
	
	/**
	 * This returns the total unassociated (truly off-fault) MFD
	 * @return
	 */
	public SummedMagFreqDist getTotalUnassociatedMFD() {
		return totalTrulyOffFaultMFD;
	}

	/**
	 * This returns the total supraseismogenic on-fault MFD
	 * @return
	 */
	public SummedMagFreqDist getTotalSupraSeisOnFaultMFD() {
		return totalSupraSeisOnFaultMFD;
	}
	
	/**
	 * This returns the total gridded seismicity MFD (the total sub-seismogenic 
	 * on-fault MFD plus the total unassociated MFD)
	 * @return
	 */
	public IncrementalMagFreqDist getTotalGriddedSeisMFD() {
		return totGriddedSeisMFD;
	}
	
	
	@Override
	public double getFracStrikeSlip(int gridIndex) {
		return fracStrikeSlip[gridIndex];
	}


	@Override
	public double getFracReverse(int gridIndex) {
		return fracReverse[gridIndex];
	}


	@Override
	public double getFracNormal(int gridIndex) {
		return fracNormal[gridIndex];
	}
	

	@Override
	public GriddedRegion getGriddedRegion() {
		return griddedRegion;
	}
	
	
	@Override
	public String getName() {
		return "GridSourceProvider2023";
	}	
	
	/**
	 * This test that the total gridded seismicity summed over all cubes equals the target
	 */
	private void testTotalGriddedSeisMFD() {
		SummedMagFreqDist testMFD = initSummedMFD();

		for(int i=0;i<cgr.getGriddedRegion().getNumLocations();i++) {
			testMFD.addIncrementalMagFreqDist(getMFD(i));
		}
		for(int i=0;i<totGriddedSeisMFD.size();i++) {
			if(totGriddedSeisMFD.getY(i)>0.0 ) {
				double ratio = testMFD.getY(i)/totGriddedSeisMFD.getY(i);
				if(ratio>1.00001 || ratio < 0.99999)
					throw new RuntimeException("testTotalGriddedSeisMFD() failed");
			}
			else if(testMFD.getY(i)>0.0)
				throw new RuntimeException("testTotalGriddedSeisMFD() failed at index "+i+"; "+testMFD.getY(i)+" vs "+totGriddedSeisMFD.getY(i));
//			System.out.println(totGriddedSeisMFD.getX(i)+"\t"+totGriddedSeisMFD.getY(i)+"\t"+testMFD.getY(i)+"\t"+(float)(testMFD.getY(i)/totGriddedSeisMFD.getY(i)));
		}	
		if(D) System.out.println("testTotalGriddedSeisMFD() succeeded.");
	}

	/**
	 * This test that the total truly off fault gridded seismicity summed over all cubes equals the target
	 */
	private void testTotalTrulyOffFaultGriddedSeisMFD() {
		
		SummedMagFreqDist testMFD = initSummedMFD();
		
		for(int c=0;c<cgr.getNumCubes();c++) {
			SummedMagFreqDist mfd = getUnassociatedMFD_ForCube(c);
			testMFD.addIncrementalMagFreqDist(mfd);
		}
		
		for(int i=0;i<totalTrulyOffFaultMFD.size();i++) {
			if(totalTrulyOffFaultMFD.getY(i)>0.0 ) {
				double ratio = testMFD.getY(i)/totalTrulyOffFaultMFD.getY(i);
				if(ratio>1.00001 || ratio < 0.99999)
					throw new RuntimeException("testTotalTrulyOffFaultGriddedSeisMFD() failed");
			}
			else if(testMFD.getY(i)>0.0)
				throw new RuntimeException("testTotalTrulyOffFaultGriddedSeisMFD() failed at index "+i+"; "+testMFD.getY(i)+" vs "+totalTrulyOffFaultMFD.getY(i));
//			System.out.println(totalTrulyOffFaultMFD.getX(i)+"\t"+totalTrulyOffFaultMFD.getY(i)+"\t"+testMFD.getY(i)+"\t"+(float)(testMFD.getY(i)/totalTrulyOffFaultMFD.getY(i)));
		}
		
		if(D) System.out.println("testTotalTrulyOffFaultGriddedSeisMFD() succeeded.");

	}

	
	/**
	 * This test fails at M=5.05 because of a low min supra mag for Quien Sabe (5.07)
	 */
	private void testTotalMgt4_RatesInCells() {
		
		double[] ratioArray = new double[spatialPDF.length];
		double totRate = totGriddedSeisMFD.getY(4.05);
		double ave=0, min=Double.MAX_VALUE, max=-Double.MAX_VALUE;
		
		for(int i=0;i<spatialPDF.length;i++) {
			double ratio = getMFD(i).getY(4.05)/(totRate*spatialPDF[i]);
			ratioArray[i]=ratio;
			ave+=ratio;
			if(min>ratio) min=ratio;
			if(max<ratio) max=ratio;
		}
		ave /= spatialPDF.length;
		
//		System.out.println("testTotalMgt4_RatesInCells(): \nave="+(float)ave+"\nmin="+
//		(float)min+"\nmax="+(float)max+"\n");
		
		if(min<0.9999 || max > 1.0001)
			throw new RuntimeException("testTotalMgt4_RatesInCells() failure: \nave="+(float)ave+"\nmin="+
					(float)min+"\nmax="+(float)max+"\n");
		
		if(D) System.out.println("testTotalMgt4_RatesInCells() succeeded.");


//		for(int i=0;i<spatialPDF.length;i++) {
//			Location loc = this.griddedRegion.getLocation(i);
//			System.out.println(loc.getLongitude()+"\t"+loc.getLatitude()+"\t"+ratioArray[i]);
//		}
	}
	
	
	/**
	 * This tests that sectDistFractWtMapAtCubeList and  totSectDistFracWtAtCubeArray are correct
	 * (they are for all fault sections inside the RELM region)
	 */
	private void testSectDistFractWtMapList() {
		
		double[] testArray = new double[rupSet.getNumSections()];
		
		for(int c=0;c<cgr.getNumCubes();c++) {
			HashMap<Integer,Double> sectWtMap = sectDistWtMapAtCubeList.get(c);
			for(int sectIndex: sectWtMap.keySet()) {
				testArray[sectIndex] += sectWtMap.get(sectIndex)/totDistWtsAtCubesForSectArray[sectIndex];
			}
		}
		
		for(int s=0;s<testArray.length;s++) {
			if(testArray[s]<0.9999 || testArray[s]>1.0001) {
				if(!sectionsThatNucleateOutsideRegionList.contains(s))
					throw new RuntimeException("testSectDistFractWtMapList() failure for section index "+s+"; "+rupSet.getFaultSectionData(s).getName());
			}
		}
		
		if(D) System.out.println("testSectDistFractWtMapList() succeeded.");
	}
	

	
	

	
	public static void main(String[] args) {
		
		// example instantiation using UCERF3 ingredients
		
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

		double[] spatialPDF = SpatialSeisPDF.UCERF3.getPDF();  // this sums to 0.9994463999998295; correct it to 1.0
		double sum=0;
		for(double val:spatialPDF) sum+=val;
		for(int i=0;i<spatialPDF.length;i++)
			spatialPDF[i] = spatialPDF[i]/sum;
		
		GridSourceProvider gridSrcProviderU3 = fss.getGridSourceProvider();
		IncrementalMagFreqDist tempMFD = gridSrcProviderU3.getMFD(0);
		SummedMagFreqDist totGriddedSeisMFD = new SummedMagFreqDist(tempMFD.getMinX(), tempMFD.size(),tempMFD.getDelta());
		for(int i=0;i<gridSrcProviderU3.size();i++) {
			totGriddedSeisMFD.addIncrementalMagFreqDist(gridSrcProviderU3.getMFD(i));	
		}
		
		CubedGriddedRegion cgr = new CubedGriddedRegion(griddedRegion);
		
		double[] fracStrikeSlip, fracReverse, fracNormal;
		fracStrikeSlip = new GridReader("StrikeSlipWts.txt").getValues();
		fracReverse = new GridReader("ReverseWts.txt").getValues();
		fracNormal = new GridReader("NormalWts.txt").getValues();
		
		long startTime = System.currentTimeMillis();

		GridSourceProvider2023 gridProvider = new GridSourceProvider2023(fss, cgr, spatialPDF, totGriddedSeisMFD, 
				binnedDepthDistFunc, fracStrikeSlip, fracReverse, fracNormal);
				
		long runtime = System.currentTimeMillis()-startTime;
		System.out.println("Runtime = "+(float)(runtime/60000d)+" min");
	}

}
