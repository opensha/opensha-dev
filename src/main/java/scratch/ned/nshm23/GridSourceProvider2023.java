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
import org.opensha.sha.gui.infoTools.CalcProgressBar;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import scratch.UCERF3.analysis.GMT_CA_Maps;
import scratch.UCERF3.enumTreeBranches.SpatialSeisPDF;
import scratch.UCERF3.erf.ETAS.ETAS_SimAnalysisTools;
import scratch.UCERF3.erf.ETAS.SeisDepthDistribution;
import scratch.UCERF3.utils.MatrixIO;
import scratch.UCERF3.utils.RELM_RegionUtils;
import scratch.UCERF3.utils.U3FaultSystemIO;
import scratch.ned.FSS_Inversion2019.PlottingUtils;

/**
 * 
 * Questions: 
 * 
 * move getSectionPolygonRegion(*) method to faultSection; confirm that trace is not offset
 * 
 * 	 * TODO - confirm that trace is surface projection and not offset by DDW
 * 
 * 	have computeLongTermSupraSeisMFD_OnSectArray() from fss apply the  if it's not already done
 * for the supplied fss.  Does rupSet.getMinMagForSection(s) really get the final minMag?
 * 
 * move the following to a more general class: ETAS_SimAnalysisTools.writeMemoryUse()
 * 
 * 
 * @author field
 *
 */
public class GridSourceProvider2023 {
	
	final static boolean D = true;
	
	final static double DEFAULT_MAX_FAULT_NUCL_DIST = 12d;		
	
	String defaultSectAtCubeCacheFilename = "/Users/field/tmp/defaultSectAtCubeCache";
	String defaultSectDistForCubeCacheFilename = "/Users/field/tmp/defaultSectDistForCubeCache";
	
	double maxFaultNuclDist;
	
	CubedGriddedRegion cgr;
	
	double[] spatialPDF;
	FaultSystemSolution fss;
	FaultSystemRupSet rupSet;
	HistogramFunction depthNuclProbHist;
	GriddedRegion griddedRegion;
	
	IncrementalMagFreqDist totGriddedSeisMFD; // supplied as input

	SummedMagFreqDist totalSubSeisOnFaultMFD, totalTrulyOffFaultMFD, totalSupraSeisOnFaultMFD; // all computed

	List<int[]> sectAtCubeList;
	List<float[]> sectDistToCubeList;

	// The following list contains, for each cube, a map of sections and their distance-fraction wts (where
	// the wts represent the fraction of seismicity assinged to the fault section below the min seismo mag).
	ArrayList<HashMap<Integer,Double>> sectDistWtMapAtCubeList;
	
	// this is the total wt for each section summed from sectDistFractWtMapList (divide the wt directly above
	// by this value to get the nucleation fraction for the section in the associated cube) 
	double[] totSectDistWtAtCubeArray;
	
	IncrementalMagFreqDist[] longTermSupraSeisMFD_OnSectArray;


	
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
		this.griddedRegion = cgr.getGriddedRegion();
		
		if(griddedRegion.getNodeCount() != spatialPDF.length)
			throw new RuntimeException("griddedRegion and spatialPDF have differe sizes: "+griddedRegion.getNodeCount()+" vs "+spatialPDF.length);
		
		// test that spatialPDF sums to 1.0
		double testSum=0;
		for(double val:spatialPDF) testSum += val;
		if(testSum>1.001 || testSum < 0.999)
			throw new RuntimeException("spatialPDF values must sum to 1.0; sum="+testSum);

		
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
			HashMap<Integer,Double> cubeDistMap) {

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

		GriddedRegion griddedPolygonReg = new GriddedRegion(fltPolygon, cgr.getCubeLatLonSpacing(), cgr.getCubeLocationForIndex(0));

		for(int i=0;i<griddedPolygonReg.getNumLocations();i++) {
			Location loc = griddedPolygonReg.getLocation(i);
			double depthDiscr = cgr.getCubeDepthDiscr();
			for(double depth = depthDiscr/2;depth<cgr.getMaxDepth();depth+=depthDiscr) {
				Location loc2 = new Location(loc.getLatitude(),loc.getLongitude(),depth);
				double dist = LocationUtils.distanceToSurf(loc2, sectQuadSurf);
				int cubeIndex = cgr.getCubeIndexForLocation(loc2);
				if(dist>maxFaultNuclDist || cubeIndex==-1)
					continue;
				cubeDistMap.put(cubeIndex,dist);
			}
		}
	}


	
	
	private void readOrGenerateCacheData() {
		
		File sectAtCubeCacheFile = new File(defaultSectAtCubeCacheFilename);
		File sectDistForCubeCacheFile = new File(defaultSectDistForCubeCacheFilename);
		
		// make cache files if they don't exist
		if (!sectAtCubeCacheFile.exists() || !sectDistForCubeCacheFile.exists()) { // read from file if it exists
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
			if(D) ETAS_SimAnalysisTools.writeMemoryUse("Memory after reading files");
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

		for(int i=0; i<cgr.getNumCubes();i++) {
			sectAtCubeListTemp.add(new ArrayList<Integer>());
			sectDistToCubeListTemp.add(new ArrayList<Float>());
		}
		
		if (progressBar != null) progressBar.showProgress(true);
		int numSect = rupSet.getNumSections();
		for(int sectIndex=0;sectIndex<numSect;sectIndex++) {
			if (progressBar != null) progressBar.updateProgress(sectIndex, numSect);
			
			HashMap<Integer,Double> cubeDistMap = new HashMap<Integer,Double>();
			
			getCubeDistancesAndFractionsForFaultSection(sectIndex, cubeDistMap);

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
		try {
			MatrixIO.intListListToFile(sectAtCubeListTemp,sectAtCubeCacheFile);
			MatrixIO.floatListListToFile(sectDistToCubeListTemp, sectDistForCubeCacheFile);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
//System.exit(0);
		
		if (progressBar != null) progressBar.showProgress(false);
		
		if(D) System.out.println(this.getClass().getName()+".generateAndWriteListListDataToFile() took "+(System.currentTimeMillis()-st)/60000+ " min");

	}
	


	private double getDistWt(double dist) {
		return (maxFaultNuclDist-dist)/maxFaultNuclDist;
	}
	

	
	
	/**
	 * This list contains, for each cube, a map of the sections therein and their distance wts.  
	 * If multiple sections exist in a cube and their wts exceed 1.0, they are all normalized by
	 * the sum so they add up to 1.0.
	 */
	private void makeSectDistWtMapList() {
		
		sectDistWtMapAtCubeList = new ArrayList<HashMap<Integer,Double>>();
		totSectDistWtAtCubeArray = new double[rupSet.getNumSections()];
		
		for(int c=0;c<cgr.getNumCubes();c++) {
			HashMap<Integer,Double> sectWtMap = new HashMap<Integer,Double>();
			int numSect = sectAtCubeList.get(c).length;
			double[] prelimWts = new double[numSect];
			double wtSum = 0;
			for(int i=0;i<numSect;i++) {
				float dist = sectDistToCubeList.get(c)[i];
				prelimWts[i] = getDistWt(dist);
				wtSum += prelimWts[i];
			}
			if(wtSum>1.0) {
				for(int i=0;i<numSect;i++) {
					prelimWts[i] /= wtSum;
				}
			}

			for(int i=0;i<numSect;i++) {
				int sectIndex = sectAtCubeList.get(c)[i];
				double wt = prelimWts[i];
				sectWtMap.put(sectIndex, wt);
				totSectDistWtAtCubeArray[sectIndex] += wt;
			}
			sectDistWtMapAtCubeList.add(sectWtMap);
		}
	}

		
	
	/**
	 * this creates a blank (zero y-axis values) MFD with the same discretization as the constructor supplied totGriddedSeisMFD.
	 * @return
	 */
	private SummedMagFreqDist getBlankMFD() {
		return new SummedMagFreqDist(totGriddedSeisMFD.getMinX(), totGriddedSeisMFD.size(),totGriddedSeisMFD.getDelta());
	}
	
	
	private void computeTotalOnAndOffFaultGriddedSeisMFDs() {
		
		totalSubSeisOnFaultMFD = getBlankMFD();
		totalSubSeisOnFaultMFD.setName("totalSubSeisMFD");
		totalTrulyOffFaultMFD = getBlankMFD();
		totalTrulyOffFaultMFD.setName("totalTrulyOffFaultMFD");
		
		for(int c=0;c<cgr.getNumCubes();c++) {
			SummedMagFreqDist mfd = getSubSeismoMFD_ForCube(c);
			if(mfd != null)
				totalSubSeisOnFaultMFD.addIncrementalMagFreqDist(mfd);
		}
		
		for(int i=0;i<totGriddedSeisMFD.size();i++) {
			totalTrulyOffFaultMFD.add(i, totGriddedSeisMFD.getY(i) - totalSubSeisOnFaultMFD.getY(i));
		}
		if(D) {
			System.out.println("totGriddedSeisMFD:\n"+totGriddedSeisMFD);
			System.out.println("totGriddedSeisMFD Cumulative::\n"+totGriddedSeisMFD.getCumRateDistWithOffset());
			System.out.println("totSubSeisMFD:\n"+totalSubSeisOnFaultMFD);
			System.out.println("totalTrulyOffFaultMFD:\n"+totalTrulyOffFaultMFD);
		}
	}
	
	/**
	 * Need to figure out how to compute this when fss has the module: ModSectMinMags
	 */
	private void computeLongTermSupraSeisMFD_OnSectArray() {
		
		SummedMagFreqDist mfd = getBlankMFD();

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
		
		HashMap<Integer,Double> sectWtMap = sectDistWtMapAtCubeList.get(cubeIndex);
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
	
	public SummedMagFreqDist getGriddedSeisMFD_ForCube(int cubeIndex) {
		SummedMagFreqDist cubeMFD = getBlankMFD();
		SummedMagFreqDist mfd = getSubSeismoMFD_ForCube(cubeIndex);
		if(mfd != null)
			cubeMFD.addIncrementalMagFreqDist(mfd);
		mfd = getTrulyOffFaultMFD_ForCube(cubeIndex);
		if(mfd != null)
			cubeMFD.addIncrementalMagFreqDist(mfd);
		return cubeMFD;
	}
	
	public SummedMagFreqDist getSupraSeisMFD_ForCube(int cubeIndex) {
		if(longTermSupraSeisMFD_OnSectArray == null)
			computeLongTermSupraSeisMFD_OnSectArray();
		HashMap<Integer,Double> sectWtMap = sectDistWtMapAtCubeList.get(cubeIndex);
		if(sectWtMap.size()==0) // no sections nucleate here
			return null;
		SummedMagFreqDist supraMFD = getBlankMFD();
		for(int s:sectWtMap.keySet()) {
			double wt = sectWtMap.get(s)/totSectDistWtAtCubeArray[s];
			IncrementalMagFreqDist mfd = longTermSupraSeisMFD_OnSectArray[s].deepClone();
			mfd.scale(wt);
			supraMFD.addIncrementalMagFreqDist(mfd);
		}
		return supraMFD;
	}
	
	public SummedMagFreqDist getTotalMFD_ForCube(int cubeIndex) {
		SummedMagFreqDist cubeMFD = getBlankMFD();
		SummedMagFreqDist mfd = getGriddedSeisMFD_ForCube(cubeIndex);
		if(mfd != null)
			cubeMFD.addIncrementalMagFreqDist(mfd);
		mfd = getSupraSeisMFD_ForCube(cubeIndex);
		if(mfd != null)
			cubeMFD.addIncrementalMagFreqDist(mfd);
		return cubeMFD;
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
	
	/**
	 * The computes how many different sections nucleate in each cube and then creates a
	 * histogram (how many have 0, 1, 2, etc sections in the cube)
	 */
	public HistogramFunction computeHistogramOfNumSectionsInCubes() {
		HistogramFunction numCubesWithNumSectHist = new HistogramFunction(0.0, 21,1.0);

		int maxNum = 0;
		for(int c=0; c<cgr.getNumCubes(); c++) {
			int num = sectAtCubeList.get(c).length;
			numCubesWithNumSectHist.add(num, 1.0);
			if(maxNum<num)
				maxNum = num;
		}
		
		// write out those with the max muber:
		if(D) {
			for(int c=0; c<cgr.getNumCubes(); c++) {
				int num = sectAtCubeList.get(c).length;
				if(num==maxNum) {
					System.out.println("\nCube "+c+ " has "+ maxNum+" sections; "+cgr.getCubeLocationForIndex(c));
					for(int i=0;i<sectAtCubeList.get(c).length;i++) {
						int s = sectAtCubeList.get(c)[i];
						float dist = sectDistToCubeList.get(c)[i];
						float wt = (float)getDistWt(dist);
						System.out.println("\t"+s+"\t"+dist+"\t"+wt+"\t"+rupSet.getFaultSectionData(s).getName());
					}
				}
			}
			System.out.println(numCubesWithNumSectHist);	
		}
		return numCubesWithNumSectHist;
	}

	
	/**
	 * This test that the total gridded seismicity summed over all cubes equals the target
	 */
	private void testTotalGriddedSeisMFD() {
		SummedMagFreqDist testMFD = getBlankMFD();

		for(int i=0;i<cgr.getGriddedRegion().getNumLocations();i++) {
			testMFD.addIncrementalMagFreqDist(getGriddedSeisMFD_ForCell(i));
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
	}

	/**
	 * This test that the total truly off fault gridded seismicity summed over all cubes equals the target
	 */
	private void testTotalTrulyOffFaultGriddedSeisMFD() {
		
		SummedMagFreqDist testMFD = getBlankMFD();
		
		for(int c=0;c<cgr.getNumCubes();c++) {
			SummedMagFreqDist mfd = getTrulyOffFaultMFD_ForCube(c);
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
	}

	
	/**
	 * This test fails at M=5.05 because of a low min supra mag for Quien Sabe (5.07)
	 */
	private void testTotalMgt4_RatesInCells() {
		
		double[] ratioArray = new double[spatialPDF.length];
		double totRate = totGriddedSeisMFD.getY(4.05);
		double ave=0, min=Double.MAX_VALUE, max=-Double.MAX_VALUE;
		
		for(int i=0;i<spatialPDF.length;i++) {
			double ratio = getGriddedSeisMFD_ForCell(i).getY(4.05)/(totRate*spatialPDF[i]);
			ratioArray[i]=ratio;
			ave+=ratio;
			if(min>ratio) min=ratio;
			if(max<ratio) max=ratio;
		}
		ave /= spatialPDF.length;
		
		System.out.println("testTotalMgt5_RatesInCells(): \nave="+(float)ave+"\nmin="+
		(float)min+"\nmax="+(float)max+"\n");

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
				testArray[sectIndex] += sectWtMap.get(sectIndex)/totSectDistWtAtCubeArray[sectIndex];
			}
		}
		
		System.out.println("testSectDistFractWtMapList():");
		
		double ave=0, min=Double.MAX_VALUE, max=-Double.MAX_VALUE;
		for(int i=0;i<testArray.length;i++) {
			ave+=testArray[i];
			if(min>testArray[i]) min=testArray[i];
			if(max<testArray[i]) max=testArray[i];
			if(testArray[i]<0.99)
				System.out.println(testArray[i]+" for "+rupSet.getFaultSectionData(i).getName());
		}
		
		ave /= testArray.length;
		
		System.out.println("nave="+(float)ave+"\nmin="+(float)min+"\nmax="+(float)max+"\n");
		
	}
	
	/**
	 * This plots the event rates above the specified magnitude for cubes at the given depth
	 * @param depth
	 * @param dirName
	 * @return
	 */
	public String plotRateAtDepthMap(double depth, double mag, String dirName) {
		
		GriddedRegion gridRegForCubes = cgr.getGridRegForCubes();
		GriddedGeoDataSet xyzDataSet = new GriddedGeoDataSet(gridRegForCubes, true);
		int depthIndex = cgr.getCubeDepthIndex(depth);
		int numCubesAtDepth = xyzDataSet.size();
		CalcProgressBar progressBar = new CalcProgressBar("Looping over all points", "junk");
		progressBar.showProgress(true);
		
		int magIndex = totGriddedSeisMFD.getClosestXIndex(mag);

		for(int i=0; i<numCubesAtDepth;i++) {
			progressBar.updateProgress(i, numCubesAtDepth);
			int cubeIndex = cgr.getCubeIndexForRegAndDepIndices(i, depthIndex);
			SummedMagFreqDist mfd = getGriddedSeisMFD_ForCube(cubeIndex);
			double rate = 0.0;
			if(mfd != null)
				rate = mfd.getCumRate(magIndex);
			if(rate == 0.0)
				rate = 1e-16;
			xyzDataSet.set(i, rate);
//			Location loc = xyzDataSet.getLocation(i);
//			System.out.println(loc.getLongitude()+"\t"+loc.getLatitude()+"\t"+xyzDataSet.get(i));
		}
		progressBar.showProgress(false);

		
		GMT_MapGenerator mapGen = GMT_CA_Maps.getDefaultGMT_MapGenerator();
		
		CPTParameter cptParam = (CPTParameter )mapGen.getAdjustableParamsList().getParameter(GMT_MapGenerator.CPT_PARAM_NAME);
		cptParam.setValue(GMT_CPT_Files.MAX_SPECTRUM.getFileName());
		cptParam.getValue().setBelowMinColor(Color.WHITE);
		
		mapGen.setParameter(GMT_MapGenerator.MIN_LAT_PARAM_NAME,gridRegForCubes.getMinGridLat());
		mapGen.setParameter(GMT_MapGenerator.MAX_LAT_PARAM_NAME,gridRegForCubes.getMaxGridLat());
		mapGen.setParameter(GMT_MapGenerator.MIN_LON_PARAM_NAME,gridRegForCubes.getMinGridLon());
		mapGen.setParameter(GMT_MapGenerator.MAX_LON_PARAM_NAME,gridRegForCubes.getMaxGridLon());
		mapGen.setParameter(GMT_MapGenerator.GRID_SPACING_PARAM_NAME, gridRegForCubes.getLatSpacing());	// assume lat and lon spacing are same
		mapGen.setParameter(GMT_MapGenerator.GMT_SMOOTHING_PARAM_NAME, false);
		mapGen.setParameter(GMT_MapGenerator.GRD_VIEW_PARAM_NAME,true);

		
		mapGen.setParameter(GMT_MapGenerator.LOG_PLOT_NAME,true);
//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_FROMDATA);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_MANUALLY);
//		double maxZ = Math.ceil(Math.log10(xyzDataSet.getMaxZ()))+0.5;
//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,maxZ-5);
//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,maxZ);
		
		if(mag<5) {
			mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,-5d);
			mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,-1d);			
		}
		else {
			mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,-11d);
			mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,-6d);
		}

		String metadata = "Map from calling plotRateAtDepthMap(*) method";
		
		try {
				String url = mapGen.makeMapUsingServlet(xyzDataSet, "M≥"+mag+" Rates at "+depth+" km depth", metadata, dirName);
				metadata += GMT_MapGuiBean.getClickHereHTML(mapGen.getGMTFilesWebAddress());
				ImageViewerWindow imgView = new ImageViewerWindow(url,metadata, true);		
				
				File downloadDir = new File(dirName);
				if (!downloadDir.exists())
					downloadDir.mkdir();
				File zipFile = new File(downloadDir, "allFiles.zip");
				// construct zip URL
				String zipURL = url.substring(0, url.lastIndexOf('/')+1)+"allFiles.zip";
				FileUtils.downloadURL(zipURL, zipFile);
				FileUtils.unzipFile(zipFile, downloadDir);

//			System.out.println("GMT Plot Filename: "+name);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return "For rates at depth above mag map: "+mapGen.getGMTFilesWebAddress()+" (deleted at midnight)";
	}

	
	/**
	 * This plots the event rates above the specified magnitude for cubes at the given depth
	 * @param depth
	 * @param dirName
	 * @return
	 */
	public void plotCubeFractionOnAndOffFaultAtDepth(double depth) {
		
		GriddedRegion gridRegForCubes = cgr.getGridRegForCubes();
		String dirNameOn = "CubeFractionOnFaultAtDepth"+depth+"km";
		String dirNameOff = "CubeFractionOffFaultAtDepth"+depth+"km";
		GriddedGeoDataSet xyzDataSetFractOn = new GriddedGeoDataSet(gridRegForCubes, true);
		GriddedGeoDataSet xyzDataSetFractOff = new GriddedGeoDataSet(gridRegForCubes, true);

		int depthIndex = cgr.getCubeDepthIndex(depth);
		int numCubesAtDepth = xyzDataSetFractOn.size();
		CalcProgressBar progressBar = new CalcProgressBar("Looping over all points", "junk");
		progressBar.showProgress(true);
		
		double max = -1, min = 2;;
		for(int i=0; i<numCubesAtDepth;i++) {
			progressBar.updateProgress(i, numCubesAtDepth);
			int cubeIndex = cgr.getCubeIndexForRegAndDepIndices(i, depthIndex);
			
			HashMap<Integer,Double> sectWtMap = sectDistWtMapAtCubeList.get(cubeIndex);
			double wtSum =0;
			for(int s:sectWtMap.keySet()) {
				wtSum+=sectWtMap.get(s);
			}
			xyzDataSetFractOn.set(i, wtSum);
			xyzDataSetFractOff.set(i, 1.0-wtSum);
			
			if(max<wtSum) max=wtSum;
			if(min>wtSum) min=wtSum;
		}
		
		if(D) System.out.println("plotCubeFractionOnAndOffFaultAtDepth:"+"\n\tmin="+min+"\n\tmax="+max);
		
		progressBar.showProgress(false);

		
		GMT_MapGenerator mapGen = GMT_CA_Maps.getDefaultGMT_MapGenerator();
		
		CPTParameter cptParam = (CPTParameter )mapGen.getAdjustableParamsList().getParameter(GMT_MapGenerator.CPT_PARAM_NAME);
		cptParam.setValue(GMT_CPT_Files.GMT_POLAR.getFileName());
		cptParam.getValue().setBelowMinColor(Color.WHITE);
		
		mapGen.setParameter(GMT_MapGenerator.MIN_LAT_PARAM_NAME,gridRegForCubes.getMinGridLat());
		mapGen.setParameter(GMT_MapGenerator.MAX_LAT_PARAM_NAME,gridRegForCubes.getMaxGridLat());
		mapGen.setParameter(GMT_MapGenerator.MIN_LON_PARAM_NAME,gridRegForCubes.getMinGridLon());
		mapGen.setParameter(GMT_MapGenerator.MAX_LON_PARAM_NAME,gridRegForCubes.getMaxGridLon());
		mapGen.setParameter(GMT_MapGenerator.GRID_SPACING_PARAM_NAME, gridRegForCubes.getLatSpacing());	// assume lat and lon spacing are same
		mapGen.setParameter(GMT_MapGenerator.GMT_SMOOTHING_PARAM_NAME, false);
		mapGen.setParameter(GMT_MapGenerator.GRD_VIEW_PARAM_NAME,true);

		
		mapGen.setParameter(GMT_MapGenerator.LOG_PLOT_NAME,false);
//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_FROMDATA);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_MANUALLY);
//		double maxZ = Math.ceil(Math.log10(xyzDataSet.getMaxZ()))+0.5;
//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,maxZ-5);
//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,maxZ);
		
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,-1.0);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,1.0);			

		String metadata = "Map from calling plotCubeFractionOnAndOffFaultAtDepth(*) method";
		
		try {
				String url = mapGen.makeMapUsingServlet(xyzDataSetFractOn, "Fraction On Fault at "+depth+" km depth", metadata, dirNameOn);
				metadata += GMT_MapGuiBean.getClickHereHTML(mapGen.getGMTFilesWebAddress());
				ImageViewerWindow imgView = new ImageViewerWindow(url,metadata, true);		
				
				File downloadDir = new File(dirNameOn);
				if (!downloadDir.exists())
					downloadDir.mkdir();
				File zipFile = new File(downloadDir, "allFiles.zip");
				// construct zip URL
				String zipURL = url.substring(0, url.lastIndexOf('/')+1)+"allFiles.zip";
				FileUtils.downloadURL(zipURL, zipFile);
				FileUtils.unzipFile(zipFile, downloadDir);

//			System.out.println("GMT Plot Filename: "+name);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		try {
			String url = mapGen.makeMapUsingServlet(xyzDataSetFractOff, "Fraction Off Fault at "+depth+" km depth", metadata, dirNameOff);
			metadata += GMT_MapGuiBean.getClickHereHTML(mapGen.getGMTFilesWebAddress());
			ImageViewerWindow imgView = new ImageViewerWindow(url,metadata, true);		
			
			File downloadDir = new File(dirNameOff);
			if (!downloadDir.exists())
				downloadDir.mkdir();
			File zipFile = new File(downloadDir, "allFiles.zip");
			// construct zip URL
			String zipURL = url.substring(0, url.lastIndexOf('/')+1)+"allFiles.zip";
			FileUtils.downloadURL(zipURL, zipFile);
			FileUtils.unzipFile(zipFile, downloadDir);

//		System.out.println("GMT Plot Filename: "+name);
	} catch (Exception e) {
		e.printStackTrace();
	}

//		return "For rates at depth above mag map: "+mapGen.getGMTFilesWebAddress()+" (deleted at midnight)";
	}
	
	
	
	
	

	
	public void plotCharFactorAtDepth(double depth, double bVal) {
		
		if(longTermSupraSeisMFD_OnSectArray == null)
			computeLongTermSupraSeisMFD_OnSectArray();
		
		GriddedRegion gridRegForCubes = cgr.getGridRegForCubes();
		GriddedGeoDataSet xyzDataSet = new GriddedGeoDataSet(gridRegForCubes, true);

		int depthIndex = cgr.getCubeDepthIndex(depth);
		int numCubesAtDepth = xyzDataSet.size();
		CalcProgressBar progressBar = new CalcProgressBar("Looping over all points", "junk");
		progressBar.showProgress(true);
		
		for(int i=0; i<numCubesAtDepth;i++) {
			progressBar.updateProgress(i, numCubesAtDepth);
			int cubeIndex = cgr.getCubeIndexForRegAndDepIndices(i, depthIndex);
			
			HashMap<Integer,Double> sectWtMap = sectDistWtMapAtCubeList.get(cubeIndex);
			
			double aveMinSupraMag=0;
			if(sectWtMap.size()==0) { // no sections nucleate here
				aveMinSupraMag=6.35;
			}
			else {
				double totWt=0;
				for(int s:sectWtMap.keySet()) {
					IncrementalMagFreqDist mfd = longTermSupraSeisMFD_OnSectArray[s];
					double minMag = mfd.getMinMagWithNonZeroRate();
					double wt = mfd.getTotalIncrRate()*sectWtMap.get(s)/totSectDistWtAtCubeArray[s];
					aveMinSupraMag += wt*minMag;
					totWt+=wt;
				}
				aveMinSupraMag /= totWt;			
			}
			
			SummedMagFreqDist totalMFD = getTotalMFD_ForCube(cubeIndex);
			
			int index = totalMFD.getClosestXIndex(aveMinSupraMag);
			aveMinSupraMag = totalMFD.getX(index);
			
			double minMag = totalMFD.getMinMagWithNonZeroRate();
			double maxMag = totalMFD.getMaxMagWithNonZeroRate();
			
			GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(totalMFD.getMinX(), totalMFD.size(), totalMFD.getDelta(),
					minMag, maxMag, 1.0, bVal);
			gr.scaleToIncrRate(minMag, totalMFD.getY(minMag));
			
			double charFact = totalMFD.getCumRate(aveMinSupraMag)/gr.getCumRate(aveMinSupraMag);
			xyzDataSet.set(i, charFact);
			
		}
		
		progressBar.showProgress(false);
		
		String dirName = "CharFactorAtDepth"+depth+"km";
		
		GMT_MapGenerator mapGen = GMT_CA_Maps.getDefaultGMT_MapGenerator();
		
		CPTParameter cptParam = (CPTParameter )mapGen.getAdjustableParamsList().getParameter(GMT_MapGenerator.CPT_PARAM_NAME);
		cptParam.setValue(GMT_CPT_Files.UCERF3_RATIOS.getFileName());
//		cptParam.getValue().setBelowMinColor(Color.WHITE);
		
		mapGen.setParameter(GMT_MapGenerator.MIN_LAT_PARAM_NAME,gridRegForCubes.getMinGridLat());
		mapGen.setParameter(GMT_MapGenerator.MAX_LAT_PARAM_NAME,gridRegForCubes.getMaxGridLat());
		mapGen.setParameter(GMT_MapGenerator.MIN_LON_PARAM_NAME,gridRegForCubes.getMinGridLon());
		mapGen.setParameter(GMT_MapGenerator.MAX_LON_PARAM_NAME,gridRegForCubes.getMaxGridLon());
		mapGen.setParameter(GMT_MapGenerator.GRID_SPACING_PARAM_NAME, gridRegForCubes.getLatSpacing());	// assume lat and lon spacing are same
		mapGen.setParameter(GMT_MapGenerator.GMT_SMOOTHING_PARAM_NAME, false);
		mapGen.setParameter(GMT_MapGenerator.GRD_VIEW_PARAM_NAME,true);
		mapGen.setParameter(GMT_MapGenerator.LOG_PLOT_NAME,true);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_MANUALLY);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,-3d);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,3d);
		
		
		String metadata = "Map from calling plotCharFactorAtDepth(*) method";
		
		try {
				String url = mapGen.makeMapUsingServlet(xyzDataSet, "Bulge at "+depth+" km depth", metadata, dirName);
				metadata += GMT_MapGuiBean.getClickHereHTML(mapGen.getGMTFilesWebAddress());
				ImageViewerWindow imgView = new ImageViewerWindow(url,metadata, true);		
				
				File downloadDir = new File(dirName);
				if (!downloadDir.exists())
					downloadDir.mkdir();
				File zipFile = new File(downloadDir, "allFiles.zip");
				// construct zip URL
				String zipURL = url.substring(0, url.lastIndexOf('/')+1)+"allFiles.zip";
				FileUtils.downloadURL(zipURL, zipFile);
				FileUtils.unzipFile(zipFile, downloadDir);

//			System.out.println("GMT Plot Filename: "+name);
		} catch (Exception e) {
			e.printStackTrace();
		}


		
	}
	
	/**
	 * This plots the classic MFDs for the paper
	 */
	public void plotMFDs() {
		
		totGriddedSeisMFD.setName("totGriddedSeisMFD");
		
		SummedMagFreqDist testMFD = getBlankMFD();
		testMFD.setName("Test totGriddedSeisMFD");
		testMFD.addIncrementalMagFreqDist(totalSubSeisOnFaultMFD);
		testMFD.addIncrementalMagFreqDist(totalTrulyOffFaultMFD);
		
		SummedMagFreqDist totalMFD = getBlankMFD();
		totalMFD.setName("totalMFD");
		totalMFD.addIncrementalMagFreqDist(totGriddedSeisMFD);
		totalMFD.addIncrementalMagFreqDist(totalSupraSeisOnFaultMFD);
		
		ArrayList<XY_DataSet> funcs = new  ArrayList<XY_DataSet>();
		funcs.add(totGriddedSeisMFD);
		funcs.add(totalSubSeisOnFaultMFD);
		funcs.add(totalTrulyOffFaultMFD);
		funcs.add(testMFD);
		funcs.add(totalSupraSeisOnFaultMFD);
		funcs.add(totalMFD);


		
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 4f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GREEN));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));

		
		Range xRange = new Range(5.0,9.0);
		Range yRange = new Range(1e-7, 1);
		
		PlottingUtils.writeAndOrPlotFuncs(funcs, plotChars, null, "Magnitude", "Rate (per yr)", xRange, yRange, 
				false, true, 3.5, 3.0, null, true);	

	}
	
	
	/**
	 * 
	 * To plot this result in 3D (rotatable) in Igor64 paste the data in and then paste 
	 * this into the command window (all together, but without the leading "	 * "):
	 * 
	 * Concatenate {lon,lat,depth}, tripletWave
	 * NewGizmo
	 * AppendToGizmo DefaultScatter= root:tripletWave
	 * ModifyGizmo makeTripletColorWave = {wt, grays, 1 }
	 * ModifyGizmo stopUpdates
	 * ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ scatterColorType,1}
	 * ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ colorWave,root:wt_C}
	 * ModifyGizmo resumeUpdates
	 * ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ Shape,2}
	 * ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ size,0.5}
	 * 
	 * @param parSectID
	 * @param fileName
	 */
	public void writeCubeLocsAndWtsForForParentSect(int parSectID, String fileName) {
		
		
//		for(int s=0;s<rupSet.getNumSections();s++) {
//			String name = rupSet.getFaultSectionData(s).getName();
//			if(name.contains("Mojave")) {
//				System.out.println(rupSet.getFaultSectionData(s).getParentSectionId()+"\t"+rupSet.getFaultSectionData(s).getParentSectionName());
//			}
//		}
//		System.exit(-1);
		
				
		System.out.println("writeLocAndFactionOnForCubes:\n");
		boolean didIt = false;
		
		FileWriter fw;
		try {
			fw = new FileWriter(fileName);
			fw.write("lon\tlat\tdepth\twt\n");
			for(int s=0;s<rupSet.getNumSections();s++) {
				if(rupSet.getFaultSectionData(s).getParentSectionId() == parSectID) {
					if(!didIt) {
						System.out.println(rupSet.getFaultSectionData(s).getParentSectionName());
						didIt=true;
					}
					
					HashMap<Integer,Double> cubeDistMap = new HashMap<Integer,Double>();
					getCubeDistancesAndFractionsForFaultSection(s, cubeDistMap);
					if(cubeDistMap != null) {	// null if section is outside the region
						for(int cubeIndex:cubeDistMap.keySet()) {
							HashMap<Integer,Double> sectWtMap = sectDistWtMapAtCubeList.get(cubeIndex);
							double wt = sectWtMap.get(s);
							Location loc = cgr.getCubeLocationForIndex(cubeIndex);
							fw.write((float)loc.getLongitude()+"\t"+(float)loc.getLatitude()+
									"\t"+(float)-loc.getDepth()+"\t"+(float)wt+"\n");
						}			
					}
				}
			}
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
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
		
//		// make spatialPDF constant for testing
//		for(int i=0;i<spatialPDF.length;i++)
//			spatialPDF[i] = 1.0/spatialPDF.length;

		
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
		
		// this should not exactly equal Figure 9c in U3ETAS BSSA paper (due to how char factors for mult sections in cube are handled)
//		gridProvider.plotCharFactorAtDepth(7d,1d);

		// the following matches the U3ETAS map (Figure 5a)
//		gridProvider.plotRateAtDepthMap(7d,2.55,"RatesAboveM2pt5_AtDepth7km");

//		gridProvider.plotCubeFractionOnAndOffFaultAtDepth(7d);
		
//		gridProvider.writeCubeLocsAndWtsForForParentSect(301, "CubeLocWtsForSAF_MojaveSouth.txt");
//		gridProvider.writeCubeLocsAndWtsForForParentSect(141, "CubeLocWtsForGreatValley14_KH.txt");
		
//		gridProvider.testTotalGriddedSeisMFD();
//		gridProvider.testTotalTrulyOffFaultGriddedSeisMFD();

		gridProvider.plotMFDs();
		
//		gridProvider.testTotalMgt4_RatesInCells();
//		gridProvider.testSectDistFractWtMapList();
		
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
