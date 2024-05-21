package scratch.ned.GK_Declustering;

import java.awt.Color;
import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.GregorianCalendar;
import java.util.List;
import java.util.Random;
import java.util.TimeZone;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.AbstractDiscretizedFunc;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc_3D;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.ParameterList;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.params.PtSrcDistanceCorrectionParam;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.observedEarthquake.Declustering.GardnerKnopoffDeclustering;
import org.opensha.sha.earthquake.observedEarthquake.parsers.USGS_NSHMP_CatalogParser;
import org.opensha.sha.earthquake.param.AleatoryMagAreaStdDevParam;
import org.opensha.sha.earthquake.param.ApplyGardnerKnopoffAftershockFilterParam;
import org.opensha.sha.earthquake.param.BPTAveragingTypeOptions;
import org.opensha.sha.earthquake.param.BPTAveragingTypeParam;
import org.opensha.sha.earthquake.param.BackgroundRupParam;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.faultSurface.utils.PtSrcDistCorr;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncLevelParam;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncTypeParam;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.math.Quantiles;

import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO.ETAS_Catalog;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_SimAnalysisTools;
import scratch.UCERF3.erf.ETAS.FaultSystemSolutionERF_ETAS;
import scratch.UCERF3.erf.ETAS.launcher.ETAS_Launcher;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.UCERF3.griddedSeismicity.AbstractGridSourceProvider;
import scratch.UCERF3.utils.U3FaultSystemIO;
import scratch.UCERF3.utils.GardnerKnopoffAftershockFilter;
import scratch.ned.FSS_Inversion2019.FaultSystemRuptureRateInversion;
import scratch.ned.FSS_Inversion2019.PlottingUtils;
import scratch.ned.GK_Declustering.MakeFigures;

public class U3ETAS_LossSimulationAnalysis {
	
	static final boolean D = true; // debug flag
	
//	static String fssFileName = "/Users/field/workspace/git/opensha-ucerf3/src/scratch/UCERF3/data/scratch/InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip";
	static String fssFileName = "/Users/field/FilesFromOldComputerTransfer/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip";
	
	//	static String catalogsFileName = "/Users/field/Field_Other/CEA_WGCEP/UCERF3/DeclusteringAnalysis/Data/U3ETAS_SimulationData/results_m5_preserve_chain_FullTD_totRateScaleFactor1pt14.bin"; 
	static String catalogsFileName = "/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/CEA_WGCEP/UCERF3/DeclusteringOnHazardAnalysisPaper/OtherStuff/Data/U3ETAS_SimulationData/results_m5_preserve_chain_FullTD_totRateScaleFactor1pt14.bin";
	
	static double catalogDuration=500; // hard coded
	static double catalogStartYear=2012;
	static long millisPerYr = (long)(1000*60*60*24*365.25);

	
	ArrayList<ObsEqkRupList> catalogList;
	ArrayList<ObsEqkRupList> catalogDeclusteredList;
	ArrayList<ObsEqkRupList> catalogRandomizedList;
	static double minMag = 5.05;
	static double deltaMag = 0.1;
	static int numMag = 40;
	
	static double lossCurveLnMin = Math.log(0.001);
	static double lossCurveLnMax = Math.log(10);
	static int lossCurveNum = 100; //40;
	static double lossCurveDelta = (lossCurveLnMax-lossCurveLnMin)/(double)(lossCurveNum-1);
	
	
	public U3ETAS_LossSimulationAnalysis() {
		
		try {
			this.catalogList = loadCatalogs(new File(fssFileName), new File(catalogsFileName), 5.0);
		} catch (IOException | DocumentException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		System.out.println("num Catalogs = "+catalogList.size());
		
		catalogDeclusteredList = getGK_DeclusteredCatalog(catalogList);
		if(D) System.out.println("Done Declustering");
		

	}
	
	
	public static ArrayList<ObsEqkRupList> getGK_DeclusteredCatalog(ArrayList<ObsEqkRupList> catalogList) {
		
		ArrayList<ObsEqkRupList> declusteredList = new ArrayList<ObsEqkRupList>();
		for(ObsEqkRupList rupList: catalogList) {
			ObsEqkRupList declusteredCatalog = GardnerKnopoffDeclustering.getDeclusteredCatalog(rupList);
			declusteredList.add(declusteredCatalog);
		}
		return declusteredList;
	}
	
	
	/**
	 * This declusters according to whether it was spontaneous in U3ETAS
	 * @param catalogList
	 * @return
	 */
	public static ArrayList<ObsEqkRupList> getSpontaneousEventsCatalog(ArrayList<ObsEqkRupList> catalogList) {
		ArrayList<ObsEqkRupList> declusteredList = new ArrayList<ObsEqkRupList>();
		for(ObsEqkRupList rupList: catalogList) {
			ObsEqkRupList filteredCatalog = new ObsEqkRupList();
			for(ObsEqkRupture rup:rupList) {
				if(((ETAS_EqkRupture)rup).getGeneration() == 0)
							filteredCatalog.add(rup);
			}
			declusteredList.add(filteredCatalog);
		}
		
		return declusteredList;
	}

	
	
	/**
	 * This declusters using the UCERF GK filter
	 * @param catalogList
	 * @return
	 */
	public static ArrayList<ObsEqkRupList> getU3_GK_FilteredCatalog(ArrayList<ObsEqkRupList> catalogList) {
		int numGriddedRups=0;
		int numFSS_Rups=0;
		GardnerKnopoffAftershockFilter gkFilter = new GardnerKnopoffAftershockFilter(0.05, 9.95, 100);
		ArrayList<ObsEqkRupList> declusteredList = new ArrayList<ObsEqkRupList>();
		for(ObsEqkRupList rupList: catalogList) {
			ObsEqkRupList filteredCatalog = new ObsEqkRupList();
			for(ObsEqkRupture rup:rupList) {
				if(rup.getRuptureSurface().isPointSurface()) {
					numGriddedRups+=1;
					double probKeep = gkFilter.getInterpolatedY(rup.getMag());
						if(Math.random()<probKeep)
							filteredCatalog.add(rup);
				}
				else {
					numFSS_Rups+=1;
					if(Math.random()<0.97) // 3% chance of removal
						filteredCatalog.add(rup);
				}
			}
			declusteredList.add(filteredCatalog);
		}
		if(D)
			System.out.println("numGriddedRups="+numGriddedRups+"\nnumFSS_Rups"+numFSS_Rups);
		
		return declusteredList;
	}

	
	
	public ArrayList<ObsEqkRupList> getRandomizedCatalogs(Random random) {
		if(catalogRandomizedList==null) {
			if(D) System.out.println("Making randomized catalogs");
			if(random==null) random = new Random();
			catalogRandomizedList = getRandomizedCatalogs(this.getCatalogs(), random);
		}
		return catalogRandomizedList;
	}

	

	/**
	 * This not only randomizes event times within each catalog, but also among catalogs
	 * so there isn't a high variance of total rates among catalogs.  This clones the
	 * ruptures and modifies the origin times.
	 */
	public static ArrayList<ObsEqkRupList> getRandomizedCatalogs(ArrayList<ObsEqkRupList> catalogList, Random random) {
		if(random==null)
			random = new Random();
		ArrayList<ObsEqkRupList> catalogRandomizedList = new ArrayList<ObsEqkRupList>();
		for(int i=0;i<catalogList.size();i++)
			catalogRandomizedList.add(new ObsEqkRupList());
		for(ObsEqkRupList catalog:catalogList) {
			for(ObsEqkRupture rup:catalog) {
				int ithCatalog = (int)Math.floor(random.nextDouble()*catalogList.size());
				double randomTimeYr = catalogStartYear+random.nextDouble()*catalogDuration;
				if(randomTimeYr<catalogStartYear || randomTimeYr>catalogStartYear+catalogDuration)
					throw new RuntimeException("ERROR: "+randomTimeYr);
//				rup.setOriginTime((long)((double)millisPerYr*(randomTimeYr-1970)));
				ObsEqkRupture newRup = (ObsEqkRupture) rup.clone();
				newRup.setOriginTime((long)((double)millisPerYr*(randomTimeYr-1970)));
//				catalogRandomizedList.get(ithCatalog).add(rup);
				catalogRandomizedList.get(ithCatalog).add(newRup);
			}
		}
		for(int i=0;i<catalogRandomizedList.size();i++) {
			catalogRandomizedList.get(i).sortByOriginTime();
//			System.out.println(catalogRandomizedList.get(i).size()+" events in randomized catalog "+i);
		}
		return catalogRandomizedList;
	}
	
	

	/**
	 * First three
	 * @return - first three functions (index 0, 1, 2) are orig, GK declustered, and spontaneous only MDFs,
	 * then come the cumulative MFDs of these in the same order (index 3, 4, 5), index 6 is the ratio of 
	 * GK declustered to original incremental MFD, index 7 is spontaneous over orig incr MFD, and index 
	 * 8 is for the UCERF3 GK filter.
	 */
	public ArrayList<XY_DataSet> makeMFDs(boolean mkPopUpPlots) {
		
		IncrementalMagFreqDist mfd = makeMFD(catalogList);
		double bVal = Math.log10(mfd.getY(5.05)/mfd.getY(6.95))/(6.95-5.05);
		mfd.setName("Full TD Catalog");
		mfd.setInfo("bVal = "+(float)bVal);
		
		IncrementalMagFreqDist mfdDeclustered = makeMFD(catalogDeclusteredList);
		bVal = Math.log10(mfdDeclustered.getY(5.05)/mfdDeclustered.getY(6.95))/(6.95-5.05);
		mfdDeclustered.setInfo("GK Declustered");
		mfdDeclustered.setInfo("bVal = "+(float)bVal);
		
//		IncrementalMagFreqDist mfdDeclustered_U3_GKfilter = makeMFD(getU3_GK_FilteredCatalog(catalogList));
//		bVal = Math.log10(mfdDeclustered_U3_GKfilter.getY(5.05)/mfdDeclustered_U3_GKfilter.getY(6.95))/(6.95-5.05);
//		mfdDeclustered_U3_GKfilter.setName("Declustered with U3 GK Filter");
//		mfdDeclustered_U3_GKfilter.setInfo("bVal = "+(float)bVal);

		IncrementalMagFreqDist mfdSpontaneousEventsOnly = makeMFD(getSpontaneousEventsCatalog(catalogList));
		bVal = Math.log10(mfdSpontaneousEventsOnly.getY(5.05)/mfdSpontaneousEventsOnly.getY(6.95))/(6.95-5.05);
		mfdSpontaneousEventsOnly.setName("Spontaneous Events Only");
		mfdSpontaneousEventsOnly.setInfo("bVal = "+(float)bVal);


				// cumulative MFDs
		EvenlyDiscretizedFunc cumMFD = mfd.getCumRateDistWithOffset();
		bVal = Math.log10(cumMFD.getClosestYtoX(5.0)/cumMFD.getClosestYtoX(6.9))/(6.9-5.0);
		cumMFD.setInfo("bVal = "+(float)bVal);
		
		EvenlyDiscretizedFunc cumMFD_Declustered = mfdDeclustered.getCumRateDistWithOffset();
		bVal = Math.log10(cumMFD_Declustered.getClosestYtoX(5.0)/cumMFD_Declustered.getClosestYtoX(6.9))/(6.9-5.0);
		cumMFD_Declustered.setInfo("bVal = "+(float)bVal);

//		EvenlyDiscretizedFunc cumMFD_Declustered_U3_GKfilter = mfdDeclustered_U3_GKfilter.getCumRateDistWithOffset();
//		bVal = Math.log10(cumMFD_Declustered_U3_GKfilter.getClosestYtoX(5.0)/cumMFD_Declustered_U3_GKfilter.getClosestYtoX(6.9))/(6.9-5.0);
//		cumMFD_Declustered_U3_GKfilter.setInfo("bVal = "+(float)bVal);


		EvenlyDiscretizedFunc cumMFD_SpontaneousEventsOnly = mfdSpontaneousEventsOnly.getCumRateDistWithOffset();
		bVal = Math.log10(cumMFD_SpontaneousEventsOnly.getClosestYtoX(5.0)/cumMFD_SpontaneousEventsOnly.getClosestYtoX(6.9))/(6.9-5.0);
		cumMFD_SpontaneousEventsOnly.setInfo("bVal = "+(float)bVal);
		
		// ratio plot
		IncrementalMagFreqDist ratio = makeMFD();
		IncrementalMagFreqDist ratio2 = makeMFD();
		for(int i=0;i<ratio.size();i++) {
			ratio.add(i, mfdDeclustered.getY(i)/mfd.getY(i));
			ratio2.add(i, mfdSpontaneousEventsOnly.getY(i)/mfd.getY(i));
		}
		ratio.setName("MFD ratio of GK declustered over Orig");
		ratio2.setName("MFD ratio of Spontaneous only over Orig");

		GardnerKnopoffAftershockFilter gkFilter = new GardnerKnopoffAftershockFilter(0.05, 9.95, 100);

		
		if(mkPopUpPlots) {
			ArrayList<XY_DataSet> funcs = new  ArrayList<XY_DataSet>();
			funcs.add(mfd);
			funcs.add(cumMFD);
			funcs.add(mfdDeclustered);
			funcs.add(cumMFD_Declustered);
//			funcs.add(mfdDeclustered_U3_GKfilter);
//			funcs.add(cumMFD_Declustered_U3_GKfilter);
			funcs.add(mfdSpontaneousEventsOnly);
			funcs.add(cumMFD_SpontaneousEventsOnly);
			
			ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLUE));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.RED));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLACK));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GREEN));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GREEN));
			Range xRange = new Range(5, 9);
			Range yRange = new Range(1e-5, 20.);
			PlottingUtils.writeAndOrPlotFuncs(funcs, plotChars, null, "Magnitude", "Num", xRange, yRange, 
					false, true, 3.5, 3.0, null, true);
			
			

			ArrayList<XY_DataSet> funcs2 = new  ArrayList<XY_DataSet>();
			funcs2.add(ratio);
			funcs2.add(gkFilter);
			ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();
			plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			plotChars2.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLACK));
			Range yRange2 = new Range(0, 1.1);
			PlottingUtils.writeAndOrPlotFuncs(funcs2, plotChars2, null, "Magnitude", "Num", xRange, yRange2, 
					false, false, 3.5, 3.0, null, true);
		}
		
		ArrayList<XY_DataSet> returnList = new ArrayList<XY_DataSet>();
		returnList.add(mfd);
		returnList.add(mfdDeclustered);
		returnList.add(mfdSpontaneousEventsOnly);
		returnList.add(cumMFD);
		returnList.add(cumMFD_Declustered);
		returnList.add(cumMFD_SpontaneousEventsOnly);
		returnList.add(ratio);
		returnList.add(ratio2);
		returnList.add(gkFilter);
		
		return returnList;
	}
	
	
	public static ArrayList<ObsEqkRupList> loadCatalogs(File fssFile, File catalogsFile, double minMag) throws IOException, DocumentException {
		

		FaultSystemSolution sol = U3FaultSystemIO.loadSol(fssFile);;
		List<ETAS_Catalog> catalogs = ETAS_CatalogIO.loadCatalogsBinary(catalogsFile, minMag); 
		
		// temporary hack
		AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF = 2.55;
		
		FaultSystemSolutionERF_ETAS erf = ETAS_Launcher.buildERF(sol,false, 1d, 2012);
		erf.updateForecast();
		
		System.out.println("Sol has RupMFDs? "+sol.hasModule(RupMFDsModule.class));
		System.out.println("ERF totNumRups: "+erf.getTotNumRups());
		System.out.println("ERF getTotNumRupsFromFaultSystem: "+erf.getTotNumRupsFromFaultSystem());
		System.out.println("ERF getNumFaultSystemSources: "+erf.getNumFaultSystemSources());
		System.out.println("ERF getNumSources: "+erf.getNumSources());

		// set the finite rupture surfaces
		String errorMessage = null;
		for (ETAS_Catalog catalog : catalogs) {
			for (ETAS_EqkRupture rup : catalog) {
				int nth = rup.getNthERF_Index();
				EqkRupture nthRup = erf.getNthRupture(nth);
				
				// verify that magnitudes are the same
				try {
					Preconditions.checkState((float)rup.getMag() == (float)nthRup.getMag(), "Nth rupture mag=%s doesn't match ETAS mag=%s",
							nthRup.getMag(), rup.getMag());
				} catch (Exception e) {
					if(errorMessage == null)
						errorMessage = "ERROR: one or more magnitudes did not match";
//					System.out.println("Error: "+(float)rup.getMag()+"\t"+(float)nthRup.getMag()+"\t"+
//							rup.getHypocenterLocation()+"\t"+
//							nthRup.getHypocenterLocation());
//					e.printStackTrace();
//					System.exit(0);
				}
				
				// replace surface if it's a FSS rupture
				if(rup.getFSSIndex() != -1) {
					rup.setRuptureSurface(nthRup.getRuptureSurface());
					rup.setAveRake(nthRup.getAveRake());
				}
				else {
					rup.setAveRake(0.);
					rup.setPointSurface(rup.getHypocenterLocation(), 0., 90.);
				}
			}
		}
		
		if(D && errorMessage != null)
			System.out.println(errorMessage);

				
		// convert to a list of ObsEqkRupList objects (and filter M≤5 events out)
		ArrayList<ObsEqkRupList> obsEqkRupListList = new ArrayList<ObsEqkRupList>();
		
		for (ETAS_Catalog catalog : catalogs) {
			ObsEqkRupList obsEqkRupList = new ObsEqkRupList();
			for (ETAS_EqkRupture rup : catalog) {
				if(rup.getMag()>=minMag)
					obsEqkRupList.add(rup);
			}
			obsEqkRupListList.add(obsEqkRupList);
		}
		
		return obsEqkRupListList;
	}
	
	
	
	private static IncrementalMagFreqDist makeMFD() {
		IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(minMag, numMag, deltaMag);
		mfd.setTolerance(deltaMag);
		return mfd;
	}
	
	
	private static IncrementalMagFreqDist makeMFD(ArrayList<ObsEqkRupList> catalogList) {
		IncrementalMagFreqDist mfd = makeMFD();
		for(ObsEqkRupList rupList: catalogList)
			addToMFD(rupList, mfd);
		mfd.scale(1.0/(catalogList.size()*catalogDuration));
		return mfd;
	}

	
	

	
	private static IncrementalMagFreqDist makeMagNumDist(ObsEqkRupList rupList) {
		IncrementalMagFreqDist mfd = makeMFD();
		for(ObsEqkRupture rup:rupList)
			mfd.add(rup.getMag(), 1.0);
		return mfd;
	}

	
	


	private static void addToMFD(ObsEqkRupList rupList, IncrementalMagFreqDist mfd) {
		for (ObsEqkRupture rup : rupList) {
			if(rup.getMag()>1)
				mfd.add(rup.getMag(), 1.0);
		}
	}
	
	
	
	
	

	


	
	public void computeNumEventDistribution(ArrayList<ObsEqkRupList> catalogList, double duration, double[] magArray) {
		
		// get subcatalogs
		ArrayList<ObsEqkRupList> eqkRupListList = getSubcatalogList(catalogList, duration);
		
		EvenlyDiscretizedFunc tempFunc = makeMagNumDist(eqkRupListList.get(0)).getCumRateDistWithOffset();
		ArbDiscrEmpiricalDistFunc_3D cumMFD_FromAllCatalogsFunc_3D = new ArbDiscrEmpiricalDistFunc_3D(tempFunc.getMinX(), tempFunc.size(), tempFunc.getDelta());
		
		int showProgressAt = eqkRupListList.size()/10;
		for(int i=0;i<eqkRupListList.size();i++) {
			if(i>showProgressAt && D) {
				if(D) System.out.println(i);
				showProgressAt+=eqkRupListList.size()/10;
			}
			IncrementalMagFreqDist func = makeMagNumDist(eqkRupListList.get(i));
			cumMFD_FromAllCatalogsFunc_3D.set(func.getCumRateDistWithOffset(), 1.0);
		}		
		
		EvenlyDiscretizedFunc meanExpNumDist = cumMFD_FromAllCatalogsFunc_3D.getMeanCurve();
		
		System.out.println(meanExpNumDist);

		
		for(double mag:magArray) {
//			System.out.println("working on M="+mag);
			int magIndex = meanExpNumDist.getClosestXIndex(mag);
			double expNum = meanExpNumDist.getY(magIndex);
//			System.out.println("Exp Num at M="+mag+" is: "+(float)meanExpNumDist.getY(magIndex));
			ArbDiscrEmpiricalDistFunc numDist = cumMFD_FromAllCatalogsFunc_3D.getArbDiscrEmpDistFuncArray()[magIndex].getNormalizedDist();
//			System.out.println(numDist);

			int maxNum = (int)Math.round(numDist.getMaxX());
			EvenlyDiscretizedFunc poissNumDist = new EvenlyDiscretizedFunc(0.,maxNum+1,1.0);
			PoissonDistribution poisDist = new PoissonDistribution(expNum);
			for(int i=0;i<poissNumDist.size();i++) {
				poissNumDist.set(i,poisDist.probability(i));
			}
			numDist.setName("Catalog Num Distribution for M="+mag);
			numDist.setInfo("ExpNum = : "+expNum);
			poissNumDist.setName("Poisson Num Distribution for M="+mag);
			poissNumDist.setInfo("ExpNum = : "+expNum);
			ArrayList<XY_DataSet> plottingFuncsArray = new ArrayList<XY_DataSet>();
			plottingFuncsArray.add(numDist);
			plottingFuncsArray.add(poissNumDist);
			this.quickPlot(plottingFuncsArray,"Num","Probability","M="+mag+"; Duration="+duration);

		}

	}

	
	/**
	 * This divides the catalogs up into subcatalogs of the given duration
	 * @param catalogList
	 * @param duration - yrs
	 * @return
	 */
	public static ArrayList<ObsEqkRupList> getSubcatalogList(ArrayList<ObsEqkRupList> catalogList, double duration) {
		
		int numSubCatPerCat = (int)Math.floor((catalogDuration+0.003)/duration);	// add a day (0.03) to make sure we get the last window
		if(D) System.out.println("numSubCatPerCat = "+numSubCatPerCat);
		
		ArrayList<ObsEqkRupList> eqkRupListList = new ArrayList<ObsEqkRupList>();
		for(ObsEqkRupList catalog : catalogList) {
			long endEpoch = (long)((catalogStartYear-1970)*millisPerYr) + (long)(duration*millisPerYr);
			ObsEqkRupList currEqkList = new ObsEqkRupList();
			eqkRupListList.add(currEqkList);
			for(ObsEqkRupture rup: catalog) {
				if(rup.getOriginTime()<endEpoch)
					currEqkList.add(rup);
				else {
					currEqkList = new ObsEqkRupList();
					eqkRupListList.add(currEqkList);
					endEpoch += (long)(duration*millisPerYr);
				}
			}
		}
		if(D) System.out.println("eqkRupListList.size() = "+eqkRupListList.size()+"; it should be "+numSubCatPerCat*catalogList.size());
		
		return eqkRupListList;
	}
	
	
	
	
	/**
	 * This plots a 500-year cumulative number of events versus times for the least, most, and median catalogs in terms
	 * total number of events (at 500 years).
	 * @param catalogList
	 * @param dirName
	 * @param popupWindow
	 * @param plotTitle
	 */
	public void writeAndOrPlotRateVersusTime(ArrayList<ObsEqkRupList> catalogList, String dirName, boolean popupWindow, String plotTitle) {
		
		double[] totNumArray = new double[catalogList.size()];
		
		ArbitrarilyDiscretizedFunc[] allFuncsArray = new ArbitrarilyDiscretizedFunc[catalogList.size()];

		for(int i=1;i<catalogList.size();i++ ) {
			ArbitrarilyDiscretizedFunc cumNumMgt5 = new ArbitrarilyDiscretizedFunc();
			long firstTime = catalogList.get(i).get(0).getOriginTime();
			int num = 1;
			for(ObsEqkRupture rup : catalogList.get(i)) {
				double yrs = ((double)(rup.getOriginTime()-firstTime))/(double)millisPerYr;
				cumNumMgt5.set(yrs,num);
				num+=1;
			}
			allFuncsArray[i] = cumNumMgt5;
			totNumArray[i] = cumNumMgt5.getY(cumNumMgt5.size()-1);
		}
		
		double median = Quantiles.median().compute(totNumArray);
		double min = Double.MAX_VALUE;
		double max = Double.NEGATIVE_INFINITY;
		int medianIndex=-1, minIndex=-1, maxIndex=-1;
		double minMedianDiff = Double.MAX_VALUE;
		
		for(int i=1;i<catalogList.size();i++ ) {
			double medianDiff = Math.abs(totNumArray[i]-median);
			if(medianDiff < minMedianDiff) {
				minMedianDiff = medianDiff;
				medianIndex=i;
			}
			if(totNumArray[i]<min) {
				min = totNumArray[i];
				minIndex=i;
			}
			if(totNumArray[i]>max) {
				max = totNumArray[i];
				maxIndex=i;
			}
		}		
		
		if(D) {
			double maxMag = 0;
			long ot = 0;
			for(ObsEqkRupture rup : catalogList.get(maxIndex)) {
				if(rup.getMag()>maxMag) {
					maxMag = rup.getMag();
					ot = rup.getOriginTime();
				}
			}
			double yrs = ((double)ot - (catalogStartYear-1970.)*millisPerYr)/millisPerYr;
			System.out.println("maxMag = "+maxMag+" at "+yrs+" yrs");
		}


		
		ArrayList<XY_DataSet> plottingFuncsArray = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
		
		plottingFuncsArray.add(allFuncsArray[medianIndex]);
		plottingFuncsArray.add(allFuncsArray[minIndex]);
		plottingFuncsArray.add(allFuncsArray[maxIndex]);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		
		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/cumNumEventsVsTime"+plotTitle;
		String xAxisLabel = "Time (years)";
		String yAxisLabel = "Cumulative Number of Events";
		Range xAxisRange = null;
		Range yAxisRange = null;
		boolean logX = false;
		boolean logY = false;

		PlottingUtils.writeAndOrPlotFuncs(plottingFuncsArray, plotChars, plotTitle, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);
		
	}

	
	public void quickPlot(ArrayList<XY_DataSet> plottingFuncsArray, String xAxisLabel, String yAxisLabel, String title) {
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
		
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		
		Range xAxisRange = null;
		Range yAxisRange = null;
		boolean logX = false;
		boolean logY = false;

		PlottingUtils.writeAndOrPlotFuncs(plottingFuncsArray, plotChars, title, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, null, true);

	}
	
	public static void writeAndOrPlotLossCurves(ArrayList<UncertainArbDiscFunc[]> dataSetsArray, ArrayList<XY_DataSet> funcsArray, double duration, 
			double saPeriod, String dirName, boolean popupWindow, String plotTitle) {
		
		Color[] colorArray = {Color.BLUE, Color.RED, Color.BLACK, Color.MAGENTA, Color.CYAN};
		int colorIndex = 0;
		
		ArrayList<XY_DataSet> plottingFuncsArray = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	

		for(UncertainArbDiscFunc[] dataSets:dataSetsArray) {
			plottingFuncsArray.add(dataSets[1]); // for solid line
			plottingFuncsArray.add(dataSets[1]); // for shaded region
			plottingFuncsArray.add(dataSets[0].getLower());
			plottingFuncsArray.add(dataSets[0].getUpper());
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1f, colorArray[colorIndex]));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, colorArray[colorIndex]));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, colorArray[colorIndex]));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, colorArray[colorIndex]));
			colorIndex += 1;
		}
		
		
		// add other curves
		for(XY_DataSet func:funcsArray) {
			plottingFuncsArray.add(func);
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, colorArray[colorIndex]));
			colorIndex += 1;
		}
		
		
		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/"+plotTitle+"_lossCurve";
		String xAxisLabel = "FixLabel";
		String yAxisLabel = "Probability (in "+duration+" yr)";
		Range xAxisRange = new Range(1e-2,10);
		Range yAxisRange = new Range(1e-4,1.0);
		boolean logX = true;
		boolean logY = true;

		PlottingUtils.writeAndOrPlotFuncs(plottingFuncsArray, plotChars, null, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);

	}
	
	
	

	

	
	
	

	
	public static  FaultSystemSolutionERF getTimeIndERF_Instance(double duration, boolean gkFilterOn) {
		FaultSystemSolution sol=null;
		try {
			sol = U3FaultSystemIO.loadSol(new File(fssFileName));
		} catch (Exception e) {
			e.printStackTrace();
		}
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		// set parameters
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.INCLUDE);
		erf.setParameter(BackgroundRupParam.NAME, BackgroundRupType.POINT);
		erf.setParameter(ApplyGardnerKnopoffAftershockFilterParam.NAME, gkFilterOn);
		erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.POISSON);
		erf.setParameter(AleatoryMagAreaStdDevParam.NAME, 0.0);
		erf.getTimeSpan().setDuration(duration);
		erf.updateForecast();
		
		// write out parameter values
		if(D) {
			ParameterList paramList = erf.getAdjustableParameterList();
			for(int i=0;i<paramList.size(); i++) {
				Parameter<?> param = paramList.getByIndex(i);
				System.out.println(param.getName()+"\t"+param.getValue());
			}			
		}
		if(D) System.out.println("NumFaultSystemSources = "+erf.getNumFaultSystemSources());
		
		double minMag = erf.getSource(erf.getNumSources()-1).getRupture(0).getMag();
		if(D) System.out.println("minMag = "+minMag);

		return erf;

	}
	
	
	public ArrayList<ObsEqkRupList> getCatalogs() {
		return catalogList;
	}
	
	public ArrayList<ObsEqkRupList> getCatalogsDeclustered() {
		return catalogDeclusteredList;
	}
	
	


	private ArrayList<ObsEqkRupList> makeRandomCatalogListFromERF(FaultSystemSolutionERF erf) {
		ArrayList<ObsEqkRupList> catalogList = new ArrayList<ObsEqkRupList>();
		erf.getTimeSpan().setDuration(catalogDuration);
		erf.updateForecast();
		int numCatalogs = 1000;
		for(int i=0;i<numCatalogs; i++) {
			List<EqkRupture> rupList = erf.drawRandomEventSet();
			ObsEqkRupList obsRupList = new ObsEqkRupList();
			for(EqkRupture rup:rupList) {
				if(rup.getMag()<5)
					continue;
				ObsEqkRupture obsRup = new ObsEqkRupture();
				obsRup.setAveRake(rup.getAveRake());
				obsRup.setMag(rup.getMag());
				obsRup.setRuptureSurface(rup.getRuptureSurface());
				double randYrs = (catalogStartYear-1970.) + Math.random()*catalogDuration;
				long ot = (long)(randYrs*this.millisPerYr);
				obsRup.setOriginTime(ot);
				obsRupList.add(obsRup);
			}
			obsRupList.sortByOriginTime();
			catalogList.add(obsRupList);
		}
				
				return catalogList;
	}
	
	
	
	
	
	public static void main(String[] args) throws IOException, DocumentException {

		

		U3ETAS_LossSimulationAnalysis analysis = new U3ETAS_LossSimulationAnalysis();
		Random random = new Random(102553864); // for reproducibility; change argument to get different results
		analysis.getRandomizedCatalogs(random);  // make and save the randomized catalogs

		
		// STUFF FROM U3ETAS_SimulationAnalysis
		
		// FIGURE 2 (MFDs)
//		ArrayList<XY_DataSet> funcList = analysis.makeMFDs(false);
//		MakeFigures.makeFigure2_Parts(funcList);

		
		
//		// FIGURE 3 (Hazard Curves)
//		filePrefix = "Figure3";
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogs(), loc, duration, saPeriod, false, imr));
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogsDeclustered(), loc, duration, saPeriod, false, imr));
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, false, imr));
//		funcsArray.add(computeHazardCurvesFromCatalogsPoisson(analysis.getCatalogs(), loc, duration, saPeriod, imr));
//		// add random IML for Full TD
//		funcsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogs(), loc, duration, saPeriod, true, imr)[1]);
//		String dirName = "/Users/field/Field_Other/CEA_WGCEP/UCERF3/DeclusteringAnalysis/FiguresFromEclipse/Figure3";
//		MakeFigures.makeFigure3_Parts(dataSetsArray, funcsArray, duration, saPeriod, dirName, true, filePrefix);

		
		
		// FIGURE 9 (1-yr Hazard Curves)
//		duration = 1;
//		filePrefix = "Figure9";
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogs(), loc, duration, saPeriod, false, imr));
//		funcsArray.add(computeHazardCurvesFromCatalogsPoisson(analysis.getCatalogs(), loc, duration, saPeriod, imr));
//		String dirName = "/Users/field/Field_Other/CEA_WGCEP/UCERF3/DeclusteringAnalysis/FiguresFromEclipse/Figure9";
//		MakeFigures.makeFigure9_Parts(dataSetsArray, funcsArray, duration, saPeriod, dirName, true, filePrefix);

		
		// FIGURES 4 and 5
//		// THIS LOOKS AT NUM EXCEEDANCE DISTRIBUTIONS AT THE SPECIFIED IMLS USING RANDOM IML SAMPLES
//		double[] imlArray = new double[]{0.423, 1.67}; // the first value is where the probability difference is maximum
//		ArrayList<ArrayList<XY_DataSet>> listList1 = analysis.computeNumHazardExceedDistRandomIML(analysis.getCatalogs(), loc, duration, saPeriod, imlArray, imr, "Full TD Cats & Rand IML", random);
//		ArrayList<ArrayList<XY_DataSet>> listList2 = analysis.computeNumHazardExceedDistRandomIML(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, imlArray, imr, "Full Randomized Cats & Rand IML", random);
//		MakeFigures.makeNumHazExceedanceRandIML_Figures(imlArray, listList1, listList2);
//		// THIS LOOKS AT NUM EXCEEDANCE DISTRIBUTIONS AT THE SPECIFIED IMLS
//		ArrayList<ArrayList<XY_DataSet>> listList3 = analysis.computeNumHazardExceedDist(analysis.getCatalogs(), loc, duration, saPeriod, imlArray, imr, "Full TD Cats");
//		ArrayList<ArrayList<XY_DataSet>> listList4 = analysis.computeNumHazardExceedDist(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, imlArray, imr, "Full Randomized Cats");
//		MakeFigures.makeNumHazExceedanceFigures(imlArray, listList3, listList4, null);

		// THIS IS TO LOOK AT HAZARD CURVES WHERE THE FULL TD AND RANDOMIZED CASES ARE MOST DIFFERENT (SE CORNER OF REGION)
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogs(), loc, duration, saPeriod, false, imr));
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, false, imr));
//		loc = new Location(33.9,-113.5);
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogs(), loc, duration, saPeriod, false, imr));
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, false, imr));
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogsDeclustered(), loc, duration, saPeriod, false, imr));
//		String dirName = "/Users/field/Field_Other/CEA_WGCEP/UCERF3/DeclusteringAnalysis/FiguresFromEclipse/???????";
//		U3ETAS_SimulationAnalysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, "Lat=33.9, Lon=-113.5");
//		double[] imlArray = new double[]{0.10,0.17};
//		ArrayList<ArrayList<XY_DataSet>> listList1 = analysis.computeNumHazardExceedDistRandomIML(analysis.getCatalogs(), loc, duration, saPeriod, imlArray, imr, "Full TD Cats & Rand IML", random);
//		ArrayList<ArrayList<XY_DataSet>> listList2 = analysis.computeNumHazardExceedDistRandomIML(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, imlArray, imr, "Full Randomized Cats & Rand IML", random);
//		MakeFigures.makeNumHazExceedanceRandIML_Figures(imlArray, listList1, listList2, "maxDiffNumHazExceedanceRandIML_Test");
//		ArrayList<ArrayList<XY_DataSet>> listList3 = analysis.computeNumHazardExceedDist(analysis.getCatalogs(), loc, duration, saPeriod, imlArray, imr, "Full TD Cats");
//		ArrayList<ArrayList<XY_DataSet>> listList4 = analysis.computeNumHazardExceedDist(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, imlArray, imr, "Full Randomized Cats");
//		MakeFigures.makeNumHazExceedanceFigures(imlArray, listList3, listList4, "maxDiffNumHazExceedanceTest");

		
		// THIS IS TO LOOK AT HAZARD CURVES WHERE THE FULL TD AND RANDOMIZED CASES ARE 2nd MOST DIFFERENT (NEAR BIG LAGOON - BALD MTN FAULT)
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogs(), loc, duration, saPeriod, false, imr));
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, false, imr));
//		loc = new Location(42.0,-124.6);
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogs(), loc, duration, saPeriod, false, imr));
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, false, imr));
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogsDeclustered(), loc, duration, saPeriod, false, imr));
//		String dirName = "/Users/field/Field_Other/CEA_WGCEP/UCERF3/DeclusteringAnalysis/FiguresFromEclipse/???????";
//		U3ETAS_SimulationAnalysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, "Lat=33.9, Lon=-113.5");
//		double[] imlArray = new double[]{0.65,0.80};
//		ArrayList<ArrayList<XY_DataSet>> listList1 = analysis.computeNumHazardExceedDistRandomIML(analysis.getCatalogs(), loc, duration, saPeriod, imlArray, imr, "Full TD Cats & Rand IML", random);
//		ArrayList<ArrayList<XY_DataSet>> listList2 = analysis.computeNumHazardExceedDistRandomIML(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, imlArray, imr, "Full Randomized Cats & Rand IML", random);
//		MakeFigures.makeNumHazExceedanceRandIML_Figures(imlArray, listList1, listList2, "maxDiffNumHazExceedanceRandIML_Test_BigLagoon");
//		ArrayList<ArrayList<XY_DataSet>> listList3 = analysis.computeNumHazardExceedDist(analysis.getCatalogs(), loc, duration, saPeriod, imlArray, imr, "Full TD Cats");
//		ArrayList<ArrayList<XY_DataSet>> listList4 = analysis.computeNumHazardExceedDist(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, imlArray, imr, "Full Randomized Cats");
//		MakeFigures.makeNumHazExceedanceFigures(imlArray, listList3, listList4, "maxDiffNumHazExceedanceTest_BigLagoon");

		
		// REMOVE THIS METHOD ????????????????????? (RUNS TAKE ~90 MINUTES)
//		makeERF_GKfilter_HazMapData(0.2, false);
//		System.exit(0);
		
		// FOR KEVIN **********************
		
		// loop over sites in RELM gridded region
//		Location loc = new Location(34.05,-118.25);
//		Location loc = new Location(34.1,-118.2); // test loc to confirm gridded region curve
		// loop over saPeriods of 0, 0.2, 1, 5
//		double saPeriod = 0; // 0 = PGA
//		
//		ScalarIMR imr = AttenRelRef.CB_2014.instance(null);
//		
//		double duration = 50;
		
//		Random random = new Random();
//		
//		ArrayList<ObsEqkRupList> catalogsList = U3ETAS_SimulationAnalysis.loadCatalogs(new File(fssFileName), new File(catalogsFileName));
//		UncertainArbDiscDataset[] datasetsArray1 = U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(catalogsList, loc, duration, saPeriod, false, imr);
//		ArbitrarilyDiscretizedFunc poissonCurve = U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogsPoisson(catalogsList, loc, duration, saPeriod, imr);
//		ArrayList<ObsEqkRupList> declusteredCatalogsList = U3ETAS_SimulationAnalysis.getGK_DeclusteredCatalog(catalogsList);
//		UncertainArbDiscDataset[] datasetsArray2 = U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(declusteredCatalogsList, loc, duration, saPeriod, false, imr);
//		ArrayList<ObsEqkRupList> catalogsRandmizedList = U3ETAS_SimulationAnalysis.getRandomizedCatalogs(catalogsList, random);
//		UncertainArbDiscDataset[] datasetsArray3 = U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(catalogsRandmizedList, loc, duration, saPeriod, false, imr);
		
		// write out contents of poissonCurve, datasetsArray1, datasetsArray2, and datasetsArray3,
		
		// END FOR KEVIN **********************
		
		
		
		
//		// TEST PLOT THE ABOVE:
//		ArrayList<UncertainArbDiscDataset[]> dataSetsArray = new ArrayList<UncertainArbDiscDataset[]>();
//		ArrayList<XY_DataSet> funcsArray = new ArrayList<XY_DataSet>();
//		
//		String locName = "Los Angeles";
		
//		dataSetsArray.add(datasetsArray1);
//		dataSetsArray.add(datasetsArray2);
//		dataSetsArray.add(datasetsArray3);
//		funcsArray.add(poissonCurve);
//		U3ETAS_SimulationAnalysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, locName);

	
		
//		U3ETAS_SimulationAnalysis analysis = new U3ETAS_SimulationAnalysis();
//		
//		// TEST OF OUR PRIMARY COMPARISON
//		funcsArray.add(computeHazardCurvesFromCatalogsPoisson(analysis.getCatalogs(), loc, duration, saPeriod, imr)); // Magenta
//		funcsArray.add(computeHazardCurvesFromCatalogsPoisson(getU3_GK_FilteredCatalog(analysis.getCatalogs()), loc, duration, saPeriod, imr)); // Magenta
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogs(), loc, duration, saPeriod, false, imr)); // Blue
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogsDeclustered(), loc, duration, saPeriod, false, imr));	// Red

//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getRandomizedCatalogs(), loc, duration, saPeriod, false, imr)); // Black
//		U3ETAS_SimulationAnalysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, locName);

		
//		analysis.doMFDs();
		
//		// COMPARE ERF-TI CURVES WITH AND WITHOUT GK FILTER
//		FaultSystemSolutionERF erf = getTimeIndERF_Instance(duration, true);
//		ArbitrarilyDiscretizedFunc curveWith = computeHazardCurveFromERF(erf, loc, saPeriod, duration, imr);
//		curveWith.setName("ERF-TI with With GK Filter");
//		funcsArray.add(curveWith); 
//		erf.setParameter(ApplyGardnerKnopoffAftershockFilterParam.NAME, false);
//		erf.updateForecast();
//		ArbitrarilyDiscretizedFunc curveWithOut = computeHazardCurveFromERF(erf, loc, saPeriod, duration, imr);
//		curveWithOut.setName("ERF-TI with Without GK Filter");
//		funcsArray.add(curveWithOut); 
//		writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, locName);
		
//		// THIS COMPARES HAZARD CURVES USING THE TI ERF WITH RANDOM CATALOGS FROM THE ERF
//		// THE MATCH IS PERFECT
//		FaultSystemSolutionERF erf = analysis.getTimeIndERF_Instance(duration);
//		funcsArray.add(analysis.computeHazardCurveFromERF(erf, loc, saPeriod, duration, imr)); 
//		System.out.println("starting Rand Cats");
//		ArrayList<ObsEqkRupList> randCats = analysis.makeRandomCatalogListFromERF(erf);
//		System.out.println("done with Rand Cats");
//		dataSetsArray.add(analysis.computeHazardCurvesFromCatalogs(randCats, loc, duration, saPeriod, false, imr));
//		funcsArray.add(analysis.computeHazardCurvesFromCatalogsPoisson(randCats, loc, duration, saPeriod, imr));
//		analysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, locName);
		
//		// THIS COMPARES THE USGS LA HAZARD CURVE TO THAT COMPUTED WITH THE ERF FOR 4 DIFF IMRs
//		// MATCH IS GOOD AT 2 IN 50 ; 10% OFF AT LOWER AND HIGHES END
//		// NO ADDED GMM EPISTEMIC UNCERTAINTY AND ONE SMOOTHED SIES BRANCH HERE SO SOME DIFFS EXPECTED
//		FaultSystemSolutionERF erf = analysis.getTimeIndERF_Instance(duration);
//		erf.setParameter(ApplyGardnerKnopoffAftershockFilterParam.NAME, true);
//		erf.updateForecast();
//		imr = AttenRelRef.CB_2014.instance(null);
//		funcsArray.add(analysis.computeHazardCurveFromERF(erf, loc, saPeriod, duration, imr)); // HAD TO SKIP M≤5 EVENTS IN THE HAZARD CURVE CALCULATOR (HARD CODED) FOR THIS TO WORK
//		imr = AttenRelRef.ASK_2014.instance(null);
//		funcsArray.add(analysis.computeHazardCurveFromERF(erf, loc, saPeriod, duration, imr)); // HAD TO SKIP M≤5 EVENTS IN THE HAZARD CURVE CALCULATOR (HARD CODED) FOR THIS TO WORK
//		imr = AttenRelRef.BSSA_2014.instance(null);
//		funcsArray.add(analysis.computeHazardCurveFromERF(erf, loc, saPeriod, duration, imr)); // HAD TO SKIP M≤5 EVENTS IN THE HAZARD CURVE CALCULATOR (HARD CODED) FOR THIS TO WORK
//		imr = AttenRelRef.CY_2014.instance(null);
//		funcsArray.add(analysis.computeHazardCurveFromERF(erf, loc, saPeriod, duration, imr)); // HAD TO SKIP M≤5 EVENTS IN THE HAZARD CURVE CALCULATOR (HARD CODED) FOR THIS TO WORK
//		funcsArray.add(analysis.getUSGS_LosAngelesPGA_2in50_HazardCurve());
//		analysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, locName);

		
//		// THIS COMPARES ERF HAZARD CURVE TO U3ETAS RANDOMIZED CURVE; DIFFERENCE IS 16% AT 2 IN 50, WHICH SEEMS HIGH
//		// BUT IS WITHIN UNCERTAINTIES AT LA FOR U3-TI SUPPLEMENTAL DATA
//		// MFD COMPARISON SHOWS THE PROBLEM THAT U3ETAS IS UNDERPREDICTING (by 26%) EVENTS AROUND M≥7, WHICH IS ABOUT THE SAME
//		// AMOUNT THE BSSA PAPER SHOWS.  NEED TO BUMP UP THE totRateScaleFactor BY 31% (1.14*1.31=1.5) TO 
//		// BALANCE MFD DISCREPENCIES. NoERT CASE SHOWS LESS, BUT STILL A DISCREPANCY (14% ON MFD; 12% ON 2IN50).
//		FaultSystemSolutionERF erf = analysis.getTimeIndERF_Instance(duration);
//		funcsArray.add(analysis.computeHazardCurveFromERF(erf, loc, saPeriod, duration, imr));
//		funcsArray.add(analysis.computeHazardCurvesFromCatalogsPoisson(analysis.getCatalogs(), loc, duration, saPeriod, imr));
//		funcsArray.add(analysis.getUSGS_LosAngelesPGA_2in50_HazardCurve());
//		analysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, locName);
//		// COMPARE MFDS
//		ArrayList<XY_DataSet> funcsArray2 = new ArrayList<XY_DataSet>();
//		SummedMagFreqDist erfMFD = ERF_Calculator.getTotalMFD_ForERF(erf, minMag, 8.95, numMag, true);
//		funcsArray2.add(erfMFD);
//		funcsArray2.add(erfMFD.getCumRateDistWithOffset());
//		IncrementalMagFreqDist catMFD = analysis.makeMFD(analysis.getCatalogs());
//		funcsArray2.add(catMFD);
//		funcsArray2.add(catMFD.getCumRateDistWithOffset());
//		double aveLnRatio = 0;
//		int num=0;
//		for(int i=0;i<erfMFD.size();i++)
//			if(catMFD.getY(i) > 0.0) {
//				aveLnRatio += Math.log((erfMFD.getY(i)/catMFD.getY(i)));
//				num+=1;
//			}
//		aveLnRatio /=num;
//		System.out.println("aveRatio="+Math.exp(aveLnRatio));
//		analysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray2, duration, saPeriod, null, true, locName);


//		// THIS LOOKS AT NUM EXCEEDANCE DISTRIBUTIONS AT THE SPECIFIED IMLS USING RANDOM IML SAMPLES
//		analysis.computeNumHazardExceedDistRandomIML(analysis.getCatalogs(), loc, duration, saPeriod, new double[]{0.02,0.2,0.7,2.0}, imr);
//		analysis.computeNumHazardExceedDistRandomIML(analysis.getRandomizedCatalogs(), loc, duration, saPeriod, new double[]{0.02,0.2,0.7,2.0}, imr);

//		// THIS LOOKS AT NUM EXCEEDANCE DISTRIBUTIONS AT THE SPECIFIED IMLS
//		analysis.computeNumHazardExceedDist(analysis.getCatalogs(), loc, duration, saPeriod, new double[]{0.02,0.2,0.7,2.0}, imr);
//		analysis.computeNumHazardExceedDist(analysis.getRandomizedCatalogs(), loc, duration, saPeriod, new double[]{0.02,0.2,0.7,2.0}, imr);

		
		
		
//		// THIS LOOKS AT THE NUM DISTRIBUTION OF EVENTS ABOVE VARIOUS MAG THRESHOLDS
//		// THINGS LOOK POISSONIAN FOR M≥7 ISH
//		duration = 50;
//		analysis.computeNumEventDistribution(analysis.getCatalogs(), duration, new double[]{5.0,6.0, 7.0, 8.0});
//		analysis.computeNumEventDistribution(analysis.getCatalogsDeclustered(), duration, new double[]{5.0,6.0, 7.0, 8.0});
//		analysis.computeNumEventDistribution(analysis.getRandomizedCatalogs(), duration, new double[]{5.0,6.0, 7.0, 8.0});
		
//		// THIS PLOTS CUMULATIVE RATES VERSUS TIME
//		analysis.writeAndOrPlotRateVersusTime(analysis.getCatalogs(), null, true, "Test U3ETAS");
//		analysis.writeAndOrPlotRateVersusTime(analysis.getRandomizedCatalogs(), null, true, "Test Poisson");
		
//		// THIS TESTS THAT CUTTING SUMULATIONS IN HALF PRODUCES CONSISENT RESULTS
//		// cut catalogs in half:
//		ArrayList<ObsEqkRupList> allCatalogs = analysis.getCatalogs();
//		ArrayList<ObsEqkRupList> firstHalfCats = new ArrayList<ObsEqkRupList> ();
//		ArrayList<ObsEqkRupList> secondHalfCats = new ArrayList<ObsEqkRupList> ();
//		for(int i=0;i<500;i++) {
//			firstHalfCats.add(allCatalogs.get(i));
//			secondHalfCats.add(allCatalogs.get(i+500));
//		}
//		dataSetsArray.add(analysis.computeHazardCurvesFromCatalogs(firstHalfCats, loc, duration, saPeriod,false, imr));
//		dataSetsArray.add(analysis.computeHazardCurvesFromCatalogs(secondHalfCats, loc, duration, saPeriod,false, imr));
//		U3ETAS_SimulationAnalysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, locName);
		
//		// THIS TESTS THAT A RANDOMIZED CATALOG PRODUCES THE SAME CURVE AS THAT FROM A POISSON MODEL USING THE AVERAGE NUMBER OF EXCEEDANCES
//		funcsArray.add(analysis.computeHazardCurvesFromCatalogsPoisson(analysis.getCatalogs(), loc, duration, saPeriod, imr));
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getRandomizedCatalogs(), loc, duration, saPeriod, false, imr));	
//		U3ETAS_SimulationAnalysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, locName);
		
//		// THIS TESTS THAT RANDOM IML SAMPLES PRODUCE THE SAME RESULT
//		boolean randIML = false;
//		dataSetsArray.add(analysis.computeHazardCurvesFromCatalogs(analysis.getCatalogs(), loc, duration, saPeriod, randIML, imr));
//		randIML = true;
//		dataSetsArray.add(analysis.computeHazardCurvesFromCatalogs(analysis.getCatalogs(), loc, duration, saPeriod, randIML, imr));
//		analysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, locName);
		
		// THIS TESTS THAT RANDOM IMLs FROM THE IMR REPRODUCE THE DISTRIBUTION, INCLUDING TRUNCATION
//		analysis.testRandomSamplesFromIMR(analysis.getCatalogs().get(0).get(0), loc, saPeriod, imr);


	}

}
