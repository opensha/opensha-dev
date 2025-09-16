package scratch.ned.nshm23;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.geo.json.FeatureProperties;
import org.opensha.commons.geo.json.Geometry;
import org.opensha.commons.geo.json.Geometry.Polygon;
import org.opensha.sha.earthquake.FocalMechanism;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import gov.usgs.earthquake.nshmp.model.NshmErf;
import gov.usgs.earthquake.nshmp.model.NshmSource;


/**
 * The creates all the NSHM 2023 source zones used in CEUS.
 * Mean rates here match those in Peter's 2023 csv files within 2% (exception thrown in debug more if not).  
 * Differences come from other things, such as rounding errors or whether mean is from 5-point distribution 
 * or from the true pdf/cdf; we use the latter here because it's more accurate.
 * 
 * Debug mode (D=true) is quite slow because Peter's ERF wrapper must be instantiated.  Keep this off once finalized.
 * 
 * This creates sources for the following:
 * 
 *	AR/Crowleys Ridge (south)/	slip rate
 *	AR/Crowleys Ridge (west)/	slip rate
 *	AR/Joiner Ridge/			slip rate
 *	AR/Marianna/				N in T
 *	AR/Saline River/			N in T
 *	IL/Wabash Valley/			N in T
 *	ME/Charlevoix/				N in T
 *	SC/Charleston/				N in T		3 sources
 *	VA/Central Virginia/		N in T		2 sources
 * 
 * @author field
 *
 */
public class CEUS_FaultZones_creator {
	
	final static boolean D=true;
	
	final static double mMin = 6.45;
	final static double mMax = 7.9;
	final static double deltaMag = 0.05;
	final static int numMag = (int) Math.round((mMax-mMin)/deltaMag+1.0);
	
	final static double upperSesDepth = 5;
	final static double lowerSeisDepth = 22;
	
	public static ArrayList<AreaSourceNSHM23> getFaultSourceZoneList(String nshmModelDirPath, double duration, boolean lineSource) {
		
		ArrayList<AreaSourceNSHM23> sourceList = new ArrayList<AreaSourceNSHM23>();
		
		// hard coded list of source zone IDs
		ArrayList<Integer> srcZoneID_List = getSourceZoneID_List();
		
		NshmErf erf = null;
		
		if(D) {
			// get the ERf for tests
			erf = getNshmERF(nshmModelDirPath);
			
			ArrayList<String> skippedSrcNames = new ArrayList<String>();
		    System.out.println("NSHM ERF NumSources: " + erf.getNumSources());
		    ArrayList<Integer> testIdList = new ArrayList<Integer>();
		    for(int s=0;s<erf.getNumSources();s++) {
		    	NshmSource src = (NshmSource)erf.getSource(s);
		    	int id = src.getNSHM_ID();
		    	if(!srcZoneID_List.contains(id)) { // skip other sources
		    		if(!skippedSrcNames.contains(src.getName()))
		    			skippedSrcNames.add(src.getName());
		    		continue;
		    	}
		    	if(!testIdList.contains(id))
		    		testIdList.add(id);
		    }
		    
		    // these should be the fault-based sources
		    System.out.println("Skipped source names:");
		    for(String nm:skippedSrcNames)
			    System.out.println("\t"+nm);
		    System.out.println("\n");
		    
		    if(testIdList.size() != srcZoneID_List.size())
		    	throw new RuntimeException("testIdList.size() != srcZoneID_List.size()");

		    System.out.println("Unique IDs:");
		    for(int id:testIdList)
			    System.out.println("\t"+id);
		    System.out.println("\n");
		}

		
//		System.exit(0);
		
		// make list of grid source locations 	THIS WON'T WORK WHEN PETER IMPLEMENTS THE FINITE SURFACES
//		System.out.println("All point source locatations:");
//		if(erf != null && D) {
//		    for(int s=0;s<erf.getNumSources();s++) {
//		    	NshmSource src = (NshmSource)erf.getSource(s);
//		    	if(CEUS_FSS_creator.isPointSource(src)) {
//		    	Location loc = src.getRupture(0).getRuptureSurface().getEvenlyDiscretizedLocation(0);
//				System.out.println(loc.getLatitude()+"\t"+loc.getLongitude());
//		    	}
//		    }		
//		}
//		System.out.println("--------------------------------------");

//		System.exit(0);

				
		// WABASH VALLEY 
		// discrepancy with ERF rate is due to two of Peter's points being excluded due to diffs in region.contains(loc)
		String srcName = "Wabash Valley Zone Source";
		ArrayList<NinT_Data> nInTDataList = new ArrayList<NinT_Data>();
		nInTDataList.add(new NinT_Data(1.0, 1, 11000d, 13000d));
		String geojsonPath = nshmModelDirPath+"stable-crust/zone/IL/Wabash Valley/wabash.geojson";
		String csvPath = nshmModelDirPath+"stable-crust/zone/IL/Wabash Valley/wabash.csv";
		double activityWeight = 1.0;
		double netBranchWeight = 1.0;	// what's in any source-tree.json file(s)
		sourceList.add(processNinT_SourceZone(srcName,  geojsonPath,  csvPath, nInTDataList, activityWeight, netBranchWeight, erf, duration, lineSource));
		
		
		//CHARLEVOIX
		// Since strike=NaN, a random strike is being applied until Peter tells me how to handle these
		// 5-pt rates come from CEUS SSCn (2012) Table H-5.1-2; weights (0.2, 0.6, 0.2) come from Petersen et al (2014) figure.
		srcName = "Charlevoix Zone Source";
		nInTDataList = new ArrayList<NinT_Data>();
		// First one needs a T adjustment to match the mean of CEUS SSCn (2012) Table H-5.1-2 5-point rates;
		// note also that N here is the one interval not the number of quakes
		double[] targetCharlFivePoint = {7.7E-04,2.2E-03,4.2E-03,6.7E-03,9.3E-03};  // column 1 of CEUS SSCn (2012)Table 6.1.2-4
		double origT = 348;
		double meanRate = PoissonRateFromNinT_Calc.getMeanFromFivePointRates(targetCharlFivePoint);
		int n=1; // intervals, not quakes
		double t=(n+1)/meanRate; // the latter is the target
		if(D) System.out.println("targetMean="+(float)meanRate+"; new T = "+(float)t+"; Orig T = "+origT);
		nInTDataList.add(new NinT_Data(0.2, n, t, Double.NaN));  // Charlevoix 1 interval in 348
		nInTDataList.add(new NinT_Data(0.6, 3, 6e3, 7e3));  // Charlevoix 3 in 6k-7k
		nInTDataList.add(new NinT_Data(0.2, 4, 9.5e3, 10.2e3));  // Charlevoix 4 in 9.5k-10.2k
		geojsonPath = nshmModelDirPath+"stable-crust/zone/ME/Charlevoix/charlevoix.geojson";
		csvPath = nshmModelDirPath+"stable-crust/zone/ME/Charlevoix/charlevoix.csv";
		activityWeight = 1.0;
		netBranchWeight = 1.0;	// what's in any source-tree.json file(s)
		sourceList.add(processNinT_SourceZone(srcName,  geojsonPath,  csvPath, nInTDataList, activityWeight, netBranchWeight, erf, duration, lineSource));
		
		// MARIANNA
		srcName = "Marianna Zone Source";
		nInTDataList = new ArrayList<NinT_Data>();
		nInTDataList.add(new NinT_Data(0.5, 2, 9.6e3, 10.2e3));  
		nInTDataList.add(new NinT_Data(0.5, 3, 9.6e3, 10.2e3));  // 
		geojsonPath = nshmModelDirPath+"stable-crust/zone/AR/Marianna/active/marianna.geojson";
		csvPath = nshmModelDirPath+"stable-crust/zone/AR/Marianna/active/marianna.csv";
		activityWeight = 1.0;
		netBranchWeight = 0.5;	// what's in any source-tree.json file(s)
		sourceList.add(processNinT_SourceZone( srcName,  geojsonPath,  csvPath, nInTDataList, activityWeight, netBranchWeight, erf, duration, lineSource));

		
		// SALINE RIVER
		// for this one the 5-point mean is ~2% different from the true mean
		srcName = "Saline River Zone Source";
		nInTDataList = new ArrayList<NinT_Data>();
		nInTDataList.add(new NinT_Data(1.0, 1, 322, 5100));  
		geojsonPath = nshmModelDirPath+"stable-crust/zone/AR/Saline River/active/saline_river.geojson";
		csvPath = nshmModelDirPath+"stable-crust/zone/AR/Saline River/active/saline-river.csv";
		activityWeight = 1.0;
		netBranchWeight = 0.5;	// what's in any source-tree.json file(s)
		sourceList.add(processNinT_SourceZone( srcName,  geojsonPath,  csvPath, nInTDataList, activityWeight, netBranchWeight, erf, duration, lineSource));

		// CENTRAL VIRGINIA
		// these three apply to both branches
		nInTDataList = new ArrayList<NinT_Data>();
		nInTDataList.add(new NinT_Data(1.0, 1, 1800, 2800));  
		activityWeight = 1.0; 
		srcName = "Central Virginia Local Zone Source";
		geojsonPath = nshmModelDirPath+"stable-crust/zone/VA/Central Virginia/active/local/local.geojson";
		csvPath = nshmModelDirPath+"stable-crust/zone/VA/Central Virginia/active/local/local.csv";
		netBranchWeight = 0.5*0.5;	// what's in any source-tree.json file(s)
		sourceList.add(processNinT_SourceZone( srcName,  geojsonPath,  csvPath, nInTDataList, activityWeight, netBranchWeight, erf, duration, lineSource));
		srcName = "Central Virginia Regional Zone Source";
		geojsonPath = nshmModelDirPath+"stable-crust/zone/VA/Central Virginia/active/regional/regional.geojson";
		csvPath = nshmModelDirPath+"stable-crust/zone/VA/Central Virginia/active/regional/regional.csv";
		netBranchWeight = 0.5*0.5;	// what's in any source-tree.json file(s)
		sourceList.add(processNinT_SourceZone( srcName,  geojsonPath,  csvPath, nInTDataList, activityWeight, netBranchWeight, erf, duration, lineSource));
		

		
		
		// CHARLESTON LOCAL
		// note that N below are num intervals not num events, and T has to be solved 
		// for because CEUS SSCn (2012) monte-carlo sampled them over event-date uncertainties
		// for the oldest event
		// five-point rate targets come from CEUS SSCn (2012)Table 6.1.2-4
		// the following and (activityWeight) can be used for all three Charleston source zones
		nInTDataList = new ArrayList<NinT_Data>();
		// Post-2000 years Earthquakes 1886, A, B, and C;
		double[] targetBranchRates1 = {6.8E-04,1.3E-03,2.1E-03,3.1E-03,4.7E-03};  // column 1 of CEUS SSCn (2012)Table 6.1.2-4
		meanRate = PoissonRateFromNinT_Calc.getMeanFromFivePointRates(targetBranchRates1);
		n=3;
		t = (n+1)/meanRate;
		nInTDataList.add(new NinT_Data(0.80, n, t, Double.NaN));
		// Post-5,500 years Earthquakes 1886, A, B, and C; 
		double[] targetBranchRates2 = {6.8E-04,1.3E-03,2.1E-03,3.1E-03,4.7E-03};
		meanRate = PoissonRateFromNinT_Calc.getMeanFromFivePointRates(targetBranchRates2);
		n=3; // num intervals, not num events
		t = (n+1)/meanRate;
		nInTDataList.add(new NinT_Data(0.04, n, t, Double.NaN));
		// Post-5,500 years Earthquakes 1886, A, B, C and D; 
		double[] targetBranchRates3 = {5.0E-04,8.8E-04,1.3E-03,1.9E-03,2.7E-03};
		meanRate = PoissonRateFromNinT_Calc.getMeanFromFivePointRates(targetBranchRates3);
		n=4; // num intervals, not num events
		t = (n+1)/meanRate;
		nInTDataList.add(new NinT_Data(0.06, n, t, Double.NaN));
		// Post-5,500 years Earthquakes 1886, A, B, C and E; 
		double[] targetBranchRates4 = {3.4E-04,6.4E-04,9.2E-04,1.3E-03,1.9E-03};
		meanRate = PoissonRateFromNinT_Calc.getMeanFromFivePointRates(targetBranchRates4);
		n=4; // num intervals, not num events
		t = (n+1)/meanRate;
		nInTDataList.add(new NinT_Data(0.04, n, t, Double.NaN));
		// Post-5,500 years Earthquakes 1886, A, B, C, D and E; 
		double[] targetBranchRates5 = {4.6E-04,7.8E-04,1.1E-03,1.5E-03,2.2E-03};
		meanRate = PoissonRateFromNinT_Calc.getMeanFromFivePointRates(targetBranchRates5);
		n=5; // num intervals, not num events
		t = (n+1)/meanRate;
		nInTDataList.add(new NinT_Data(0.06, n, t, Double.NaN));

		// the 0.9 activity weight for all three branches:
		activityWeight = 0.9;

		srcName = "Charleston Local Zone Source";
		geojsonPath = nshmModelDirPath+"stable-crust/zone/SC/Charleston/local/local.geojson";
		csvPath = nshmModelDirPath+"stable-crust/zone/SC/Charleston/local/local.csv";
		netBranchWeight = 0.5;	// what's in any source-tree.json file(s)
		sourceList.add(processNinT_SourceZone( srcName,  geojsonPath,  csvPath, nInTDataList, activityWeight, netBranchWeight, erf, duration, lineSource));
		
		srcName = "Charleston Regional Zone Source";
		geojsonPath = nshmModelDirPath+"stable-crust/zone/SC/Charleston/regional/regional.geojson";
		csvPath = nshmModelDirPath+"stable-crust/zone/SC/Charleston/regional/regional.csv";
		netBranchWeight = 0.2;	// what's in any source-tree.json file(s)
		sourceList.add(processNinT_SourceZone( srcName,  geojsonPath,  csvPath, nInTDataList, activityWeight, netBranchWeight, erf, duration, lineSource));

		srcName = "Charleston Narrow Zone Source";
		geojsonPath = nshmModelDirPath+"stable-crust/zone/SC/Charleston/narrow/narrow.geojson";
		csvPath = nshmModelDirPath+"stable-crust/zone/SC/Charleston/narrow/narrow.csv";
		netBranchWeight = 0.3;	// what's in any source-tree.json file(s)
		sourceList.add(processNinT_SourceZone( srcName,  geojsonPath,  csvPath, nInTDataList, activityWeight, netBranchWeight, erf, duration, lineSource));



		// SLIP-RATE BASED BELOW
		
		// CROWLEY'S RIDGE (SOUTH)
		srcName = "Crowleys Ridge (South) Zone Source";
		activityWeight = 1.0;
		double length = 84.0; // km
		double ddw = 15.0; //km
		double slipRate = 0.005; // cm/yr
		geojsonPath = nshmModelDirPath+"stable-crust/zone/AR/Crowleys Ridge (south)/active/crowleys_ridge_south.geojson";
		csvPath = nshmModelDirPath+"stable-crust/zone/AR/Crowleys Ridge (south)/active/crowleys-ridge-south.csv";
		activityWeight = 1.0;
		netBranchWeight = 0.5;	// what's in any source-tree.json file(s)
		sourceList.add(processSlipRateSourceZone( srcName,  geojsonPath,  csvPath, length, ddw, slipRate, activityWeight, netBranchWeight, erf, duration, lineSource));

		// CROWLEY'S RIDGE (WEST)
		// discrepancy with ERF rate is due to two of Peter's points being excluded due to diffs in region.contains(loc)
		srcName = "Crowleys Ridge (West) Zone Source";
		activityWeight = 1.0;
		length = 53.0; // km
		ddw = 15.0; //km
		slipRate = 0.009; // cm/yr
		geojsonPath = nshmModelDirPath+"stable-crust/zone/AR/Crowleys Ridge (west)/active/crowleys_ridge_west.geojson";
		csvPath = nshmModelDirPath+"stable-crust/zone/AR/Crowleys Ridge (west)/active/crowleys-ridge-west.csv";
		activityWeight = 1.0;
		netBranchWeight = 0.5;	// what's in any source-tree.json file(s)
		sourceList.add(processSlipRateSourceZone( srcName,  geojsonPath,  csvPath, length, ddw, slipRate, activityWeight, netBranchWeight, erf, duration, lineSource));

		// JOINER RIDGE
		// 
		srcName = "Joiner Ridge Zone Source";
		activityWeight = 1.0;
		length = 35.0; // km
		ddw = 15.0; //km
		slipRate = 0.03; 
		geojsonPath = nshmModelDirPath+"stable-crust/zone/AR/Joiner Ridge/active/joiner_ridge.geojson";
		csvPath = nshmModelDirPath+"stable-crust/zone/AR/Joiner Ridge/active/joiner-ridge.csv";
		activityWeight = 1.0;
		netBranchWeight = 0.5;	// what's in any source-tree.json file(s)
		sourceList.add(processSlipRateSourceZone( srcName,  geojsonPath,  csvPath, length, ddw, slipRate, activityWeight, netBranchWeight, erf, duration, lineSource));
	
		
		return sourceList;

	}
	
	
	/**
	 * This is for slip-rate based sources
	 * @param srcName
	 * @param geojsonPath
	 * @param csvPath
	 * @param length - km
	 * @param ddw - km
	 * @param slipRate - cm/yr
	 * @param activityWeight - this is a weight applied before csv grid rate file created (only applicable to Charleston)
	 * @param netBranchWt - this is the weight applied to values in csv grid rate file (specified in Peter's source-tree.json files)
	 * @param erf - Peter's erf wrapper to test results 
	 */
	private static AreaSourceNSHM23 processSlipRateSourceZone(String srcName, String geojsonPath, String csvPath, 
			double length, double ddw, double slipRate, double activityWeight, double netBranchWt, 
			NshmErf erf, double duration, boolean lineSource) {

		// in the following there should really be a separate rate for each mag rather than using the average rate for all mags
		double moRate = FaultMomentCalc.getMoment(length*ddw*1e6, slipRate*1e-2); // SI units
		double meanRate=0;

		NSHM_SrcZoneData data = new NSHM_SrcZoneData(geojsonPath, csvPath);
		IncrementalMagFreqDist mpd = data.magProbDist;

//		System.out.println("\n\tmagnitude\tmoRate\trate");
		for(int i=0;i<mpd.size();i++) {
			if(mpd.getY(i) > 0) {
				double rate = activityWeight*moRate/MagUtils.magToMoment(mpd.getX(i));
				meanRate += mpd.getY(i)*rate;
//				System.out.println("\t"+(float)mfd.getX(i)+"\t"+(float)moRate+"\t"+(float)rate);
			}
		}
		
		return processSourceZone(srcName, data, meanRate, activityWeight, netBranchWt, erf, duration, lineSource);
		
	}

	
	
	
	/**
	 * 
	 * @param srcName
	 * @param geojsonPath - 2023 NSHM geojson file
	 * @param csvPath - 2023 NSHM csv file for grid point rates
	 * @param nInTDataList
	 * @param activityWeight - this is a weight applied before csv grid rate file created (only applicable to Charleston)
	 * @param netBranchWt - this is the weight applied to values in csv grid rate file (specified in Peter's source-tree.json files)
	 * @param erf - Peter's erf wrapper to test results 
	 */
	private static AreaSourceNSHM23 processNinT_SourceZone(String srcName, String geojsonPath, String csvPath, 
			List<NinT_Data> nInTDataList, double activityWeight, double netBranchWt, 
			NshmErf erf, double duration, boolean lineSource) {

		EvenlyDiscretizedFunc cdf = PoissonRateFromNinT_Calc.getCumRateDistribution(nInTDataList);
		double meanRate = PoissonRateFromNinT_Calc.getMeanFromCDF(cdf);
		// the following not good for CDFs that are sums of others (five point mean not a good approx of actual mean); 
		// actually, means can also differ for single NinT_Data lists (~2% for saline river)
//		double[] fivePointRates = PoissonRateFromNinT_Calc.getFivePointRatesFromDistribution(cdf);
//		double meanRate5pt = PoissonRateFromNinT_Calc.getMeanFromFivePointRates(fivePointRates);
//		System.out.println("meanRate5pt="+(float)meanRate5pt);
//		PoissonRateFromNinT_Calc.plotRateDistributions(cdf);
		
		meanRate *= activityWeight; // this should match csv file total
		
		// load info from files
		NSHM_SrcZoneData data = new NSHM_SrcZoneData(geojsonPath, csvPath);
		
		return processSourceZone(srcName, data, meanRate, activityWeight, netBranchWt, erf, duration, lineSource);

	}
	
	
	/**
	 * 
	 * @param srcName
	 * @param data
	 * @param meanRate
	 * @param activityWeight
	 * @param netBranchWt
	 * @param erf
	 * @param duration
	 * @param lineSource
	 * @return
	 */
	private static AreaSourceNSHM23 processSourceZone(String srcName, NSHM_SrcZoneData data, double meanRate,
			double activityWeight, double netBranchWt, NshmErf erf, double duration, boolean lineSource) {

		String nshmName = data.name;
		double csvTotRate = data.totalRate;
//		double strike = data.strike;
		Region zonePolygon = data.region;
		IncrementalMagFreqDist mpd = data.magProbDist;
		int id = data.id;
		FocalMechanism focalMech = new FocalMechanism(data.strike,90d,0d); // vertical strike slip
		
		if(mpd == null) {
			// these come from a different file named "stable-crust/zone/SC/Charleston/mfd-map.json"
			double[] charlestonMags = {6.7, 6.9, 7.1, 7.3, 7.5};
			double[] charlestonProbs = {0.1, 0.25, 0.3, 0.25, 0.1};
			mpd = new IncrementalMagFreqDist(mMin, mMax,numMag);
			for(int i=0;i<charlestonMags.length;i++)
				mpd.set(charlestonMags[i], charlestonProbs[i]);
		}

		double meanBranchRate = meanRate*netBranchWt; // use this one for the source

		if(D) {
			System.out.println(srcName+"\n\tID: "+id+"\n\tName: "+nshmName+"\n\tRate:"+(float)meanRate+"\n\tTestRate:"+(float)csvTotRate+
					"\n\tRate test:"+(float)(meanRate/csvTotRate)+"\n\tStrike = "+data.strike);

			System.out.println("\tMag Prob Distribution:");
			for(int i=0; i<mpd.size();i++) {
				if(mpd.getY(i)>0.0) {
					System.out.println("\t"+(float)mpd.getX(i)+"\t"+(float)mpd.getY(i));
				}
			}
			System.out.println("\tZone Polygon:");
			for(Location loc: zonePolygon.getBorder()) {
				System.out.println("\t"+loc.getLatitude()+"\t"+loc.getLongitude());
			}
			System.out.println("\n");
		}
		
		if(Math.abs(meanRate/csvTotRate - 1.0) > 0.01) {
			String message = "rate discrepancy more then 1% for "+srcName+"; rate = "+meanRate+"; rate from Peter's file = "+csvTotRate+"; ratio = "+(float)(meanRate/csvTotRate);
//			throw new RuntimeException(message);
			System.out.flush(); // this forces immediate writing
			System.err.println("WARNING: "+message);
			System.err.println();
			System.err.flush();

		}
		// test against ERF
		LocationList erfLocList = new LocationList();
		if(erf != null && D) {
			double erfRate=0;
		    for(int s=0;s<erf.getNumSources();s++) {
		    	NshmSource src = (NshmSource)erf.getSource(s);
		    	if(src.getNSHM_ID() == id) {
		    		// this only works while Peter treats these as point sources
			    	Location loc = src.getRupture(0).getRuptureSurface().getEvenlyDiscretizedLocation(0);
		    		erfLocList.add(loc);
			    	for(int r=0;r<src.getNumRuptures();r++)
			    		erfRate += src.getRupture(r).getMeanAnnualRate(1.0);
		    	}
		    }		
		    double erfRateTest = Math.abs(1.0-erfRate/(csvTotRate*netBranchWt));
		    if(erfRateTest>1e-4)
		    	throw new RuntimeException("csvTotRate inconsistent with ERF rate for "+nshmName+"; erfRate="+erfRate+" & csvTotRate*netBranchWt="+(csvTotRate*netBranchWt));
		    
			System.out.println("\terfRate:"+(float)erfRate+"\tTestRatio:"+(float)(meanBranchRate/erfRate)+"\n");
		}
		
		IncrementalMagFreqDist mfd = mpd.deepClone();
		mfd.scale(meanBranchRate);
		AreaSourceNSHM23 src = new AreaSourceNSHM23(id, nshmName, zonePolygon, 0.1, 
				mfd, upperSesDepth, lowerSeisDepth, focalMech, duration, lineSource);
		src.setTectonicRegionType(TectonicRegionType.STABLE_SHALLOW);

		GriddedRegion gridRegion = src.getGriddedRegion();
		if(D && (gridRegion.getNumLocations() != data.locList.size())) {
			double ratio = (double)gridRegion.getNumLocations()/(double)data.locList.size();
			System.out.println("\tDiff number of locs: here="+gridRegion.getNumLocations()+
					"; nshm="+data.locList.size()+"\tRatio="+(float)ratio+"\n");
		}
		
		return src;
		
	}


	
	 /**
	  * This returns the ERF with timespan duration set to 1.0
	  * @param nshmModelDirPath
	  * @return
	  */
	private static NshmErf getNshmERF(String nshmModelDirPath) {
		long startTime = System.currentTimeMillis();
	    Set<TectonicRegionType> trts = EnumSet.of(TectonicRegionType.STABLE_SHALLOW);
	    NshmErf erf = new NshmErf(Path.of(nshmModelDirPath), trts, IncludeBackgroundOption.EXCLUDE);
	    erf.getTimeSpan().setDuration(1.0);
	    erf.updateForecast();
	    long runtime = System.currentTimeMillis()-startTime;
	    float runTimeMin = (float)runtime/(float)(1000*60);
	    if(D) System.out.println("ERF instantiation took "+runTimeMin +" minutes");
	    return erf;
	}
	
	

		
		
//		if (geometry instanceof MultiLineString && ((MultiLineString)geometry).lines.size() == 2) {
//			// subduction interface, convert to simple fault data (TODO: actually store the lower trace and do it right)
//			MultiLineString multiGeom = (MultiLineString)geometry;
//			LocationList upperTrace = multiGeom.lines.get(0);
//			FaultTrace upperAsFaultTrace = new FaultTrace(null);
//			upperAsFaultTrace.addAll(upperTrace);
//			LocationList lowerTrace = multiGeom.lines.get(1);
//			FaultTrace lowerAsFaultTrace = new FaultTrace(null);
//			lowerAsFaultTrace.addAll(lowerTrace);
//			ApproxEvenlyGriddedSurface approxSurf = new ApproxEvenlyGriddedSurface(
//					upperAsFaultTrace, lowerAsFaultTrace, 1d);
//			
//			mappedProps.set(DIP, approxSurf.getAveDip());
//			if (origProps.containsKey("rake"))
//				mappedProps.set(RAKE, origProps.require("rake", Double.class));
//			mappedProps.set(UPPER_DEPTH, approxSurf.getAveRupTopDepth());
//			mappedProps.set(LOW_DEPTH, approxSurf.getAveRupBottomDepth());
//			mappedProps.set(DIP_DIR, approxSurf.getAveDipDirection());
//			
//			// convert geometry to simple line string with only upper trace
//			geometry = new LineString(upperTrace);
//			
//			// set lower trace as separate property
//			mappedProps.put(LOWER_TRACE, new LineString(lowerTrace));
//		} else {
//			mappedProps.set(DIP, origProps.require("dip", Double.class));
//			mappedProps.set(RAKE, origProps.require("rake", Double.class));
//			mappedProps.set(UPPER_DEPTH, origProps.require("upper-depth", Double.class));
//			mappedProps.set(LOW_DEPTH, origProps.require("lower-depth", Double.class));
//		}
//		
//		// optional ones to carry forward
//		if (origProps.containsKey("state"))
//			mappedProps.set("PrimState", origProps.require("state", String.class));
//		if (origProps.containsKey("references"))
//			mappedProps.set("references", origProps.get("references"));
	
	
	/**
	 * This class holds the info needed from NSHM geojson and csv files
	 * @author field
	 *
	 */
	static class NSHM_SrcZoneData {

		String name;
		protected int id;
		protected double strike;
		protected Region region;
		protected IncrementalMagFreqDist magProbDist;
		protected double totalRate;
		protected LocationList locList;
		protected ArrayList<Double> rateList;
		protected int numLocs;

		public NSHM_SrcZoneData(String geojsonPath, String csvPath) {
			Feature feature=null;
			try {
				feature = Feature.read(new File(geojsonPath));
			} catch (IOException e) {
				System.out.println("Problem with input file: "+geojsonPath);
				e.printStackTrace();
			}

			id = Integer.valueOf(feature.id.toString());

			FeatureProperties props = feature.properties;
			name = null;
			if (props.containsKey("name"))
				name = props.require("name", String.class);
			
			strike = Double.NaN;
			if (props.containsKey("strike"))
				strike = props.require("strike", Double.class);
			
			Geometry geometry = feature.geometry;
			region = null;
			if (geometry instanceof Polygon) {
				Polygon polygon = (Polygon)geometry;
				region = polygon.asRegion();
			}
			
			magProbDist = null;
	        // read in the 'mfd-tree' as a list of FeatureProperties instances, which are basically fancy maps
	        // this will return null if an error occurs
	        List<FeatureProperties> treePropsList = feature.properties.getPropertiesList("mfd-tree");
	        int numGR_branches=0; // only one branch currently supported
	        if(treePropsList != null) {
	        	magProbDist = new IncrementalMagFreqDist(mMin, mMax,numMag);
		        for (FeatureProperties treeProps : treePropsList) {
		            String tree_id = treeProps.getString("id");
		            // need to supply a default value for primitives as null can't be returned if it's not found, thus the Double.NaN
		            double weight = treeProps.getDouble("weight", Double.NaN);
		 //           System.out.println("\tid="+id+" with weight="+weight);
		            // 'value' is a map, so we'll load it in as its own FeatureProperties
		            FeatureProperties value = treeProps.getProperties("value");
		            String type = value.getString("type");
//System.out.println("type="+type);
		            if(type.equals("SINGLE")) {
			            double m = value.getDouble("m", Double.NaN);
			   		 //           System.out.println("\t\ttype: "+type+", m="+m);
			   		            magProbDist.set(m,weight);		            	
		            }
		            else {
//		                "a": 0.0,
//		                "b": 0.8,
//		                "Δm": 0.1,
//		                "mMin": 6.5,
//		                "mMax": 7.9
		            	if(numGR_branches>0)
		            		throw new RuntimeException("Code only handels one GR branch; need to fix this");
		            	double b = value.getDouble("b", Double.NaN);
		            	double deltaM = value.getDouble("Δm", Double.NaN);
		            	double mMin = value.getDouble("mMin", Double.NaN)+deltaM/2d; // make bin centered
		            	double mMax = value.getDouble("mMax", Double.NaN)-deltaM/2d;
						int numMag = (int)Math.round((mMax-mMin)/deltaM+1);
//		            	System.out.println(b+"; "+deltaM+"; "+mMin+"; "+mMax);

						magProbDist = new GutenbergRichterMagFreqDist(mMin, numMag, deltaM, mMin, mMax, 1.0, b);
		            	double totRate=magProbDist.getTotalIncrRate();
		            	magProbDist.scale(weight/totRate);	// normalize to PDF
		            	numGR_branches+=1;
		            }

		        }
		        double testTotWt = Math.abs(magProbDist.getTotalIncrRate()-1.0);
		        if(testTotWt > 0.001)
		        	throw new RuntimeException("Mag weights don't sum to 1.0 for "+name+"; sum = "+(float)magProbDist.getTotalIncrRate());
	        }
	        
	        // read csv rate file
			File file = new File(csvPath);
			locList = new LocationList();
			rateList = new ArrayList<Double>();
			try {
				CSVFile<String> csvFile = CSVFile.readFile(file, true);
				for(int i=1; i<csvFile.getNumRows(); i++) {
					locList.add(new Location(csvFile.getDouble(i, 1),csvFile.getDouble(i, 0)));
					rateList.add(csvFile.getDouble(i, 2));
					totalRate += csvFile.getDouble(i, 2);
				}
				numLocs = rateList.size();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	public static ArrayList<Integer> getSourceZoneID_List() {
		// hard coded list of source zone IDs
		ArrayList<Integer> srcZoneID_List = new ArrayList<Integer>();

		srcZoneID_List.add(3520);
		srcZoneID_List.add(3530);
		srcZoneID_List.add(3540);
		srcZoneID_List.add(3510);
		srcZoneID_List.add(3550);
		srcZoneID_List.add(3710);
		srcZoneID_List.add(3711);
		srcZoneID_List.add(3712);
		srcZoneID_List.add(3560);
		srcZoneID_List.add(3500);
		srcZoneID_List.add(3571);
		srcZoneID_List.add(3572);

// OLD models:
//		srcZoneID_List.add(10001);
//		srcZoneID_List.add(10002);
//		srcZoneID_List.add(3710);
//		srcZoneID_List.add(3711);
//		srcZoneID_List.add(3712);
//		srcZoneID_List.add(3500);
//		srcZoneID_List.add(30000);
//		srcZoneID_List.add(3300);
//		srcZoneID_List.add(40000);
//		srcZoneID_List.add(50000);
//		srcZoneID_List.add(20000);
//		srcZoneID_List.add(3400);
		
		return srcZoneID_List;

	}




	public static void main(String[] args) {
		
		String nshmModelDirPath = "/Users/field/nshm-haz_data/nshm-conus-6.1.2/";

		// the following will no longer work:
//		String nshmModelDirPath = "/Users/field/nshm-haz_data/oldVersions/nshm-conus-6.0.0/";
//		String nshmModelDirPath = "/Users/field/nshm-haz_data/nshm-conus-6.b.4/";
//		String nshmModelDirPath = "/Users/field/nshm-haz_data/nshm-conus-6.0.0/";
		
		double duration = 50;
		boolean lineSource = false;
		ArrayList<AreaSourceNSHM23> srcList = getFaultSourceZoneList(nshmModelDirPath, duration, lineSource);
		
		System.out.println("Names from srcList:");
		for(AreaSourceNSHM23 src:srcList)
			System.out.println(src.id+"\t"+src.getName());
//		System.exit(0);
		
		
		

		
		
		
		
		
		

		
//		readGeojson(nshmModelDirPath+"stable-crust/zone/AR/Crowleys Ridge (south)/active/crowleys_ridge_south.geojson");
		
//		String rateTestString = nshmModelDirPath+"stable-crust/zone/AR/Crowleys Ridge (south)/active/crowleys-ridge-south.csv";
//		double rate = getTotRateFromNshmHazCSV_file(rateTestString);
//		
//		System.out.println("rate="+rate);
//		
//
//		NshmErf erf = getNshmERF(nshmModelDirPath);
//	    System.out.println("NSHM ERF NumSources: " + erf.getNumSources());
//	    ArrayList<String> nameList = new ArrayList<String>();
//	    ArrayList<String> idList = new ArrayList<String>();
//	    double mMin=10d;
//	    Location targetLoc = new Location(38.0,-78); // Cnetal Virginia fault zone loc
//	    for(int s=0;s<erf.getNumSources();s++) {
//		    double mMinSrc=10;
//		    double mMaxSrc=0;
//	    	NshmSource src = (NshmSource)erf.getSource(s);
//    	String id = Integer.toString(src.getNSHM_ID());
//    	if(!srcZoneID_List.contains(id)) // skip other sources
//    		continue;
//	    	for(int r=0;r<src.getNumRuptures();r++) {
//	    		double mag = src.getRupture(r).getMag();
//	    		if(mMin>mag)  mMin=mag;
//	    		if(mMinSrc>mag)  mMinSrc=mag;
//	    		if(mMaxSrc<mag)  mMaxSrc=mag;
//	    	}
//	    	RuptureSurface rupSurf = src.getRupture(0).getRuptureSurface();
//	    	Location loc = rupSurf.getEvenlyDiscretizedLocation(0);
//    		System.out.println(s+"\t"+id+"\t"+src.getNumRuptures()+"\t"+mMinSrc+"\t"+mMaxSrc+"\t"+
//    		loc.toString()+"\t"+src.getName());
//
//////	    	System.out.println(id+"\t"+src.getName());
////	    	// -78	38
////	    	if(Math.abs(loc.getLatitude()-targetLoc.getLatitude())<0.0001 && Math.abs(loc.getLongitude()-targetLoc.getLongitude())<0.0001) {
//////	    		gov.usgs.earthquake.nshmp.model.NshmSource.Point
//////	    		double strike = src.getDelegate().get(0).surface().strike();
////	    		int numSurfPts = src.getRupture(0).getRuptureSurface().getEvenlyDiscretizedNumLocs();
////	    		System.out.println(s+"\t"+id+"\t"+src.getNumRuptures()+"\t"+mMinSrc+"\t"+mMaxSrc+"\t"+
////	    		loc.toString()+"\t"+src.getName()+"\t"+numSurfPts);
////	    	}
//	    	
//	    	
//// each pt source has different name and all IDs = -1
////	    	if(!nameList.contains(src.getName()))
////	    		nameList.add(src.getName());
////	    	if(!idList.contains(id))
////	    		idList.add(id);
//	    }
//		System.out.println("mMin = "+mMin);
//	    
////	    System.out.println("Unique Source Names:");
////	    for(String str:nameList)
////		    System.out.println("\t"+str);
////	    System.out.println("Unique IDs:");
////	    for(String str:idList)
////		    System.out.println("\t"+str);
	}

}
