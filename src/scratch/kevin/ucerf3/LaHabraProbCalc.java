package scratch.kevin.ucerf3;

import gov.usgs.earthquake.event.EventQuery;
import gov.usgs.earthquake.event.EventWebService;
import gov.usgs.earthquake.event.Format;
import gov.usgs.earthquake.event.JsonEvent;

import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.net.MalformedURLException;
import java.net.URL;
import java.sql.Date;
import java.text.DecimalFormat;
import java.util.GregorianCalendar;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.observedEarthquake.parsers.UCERF3_CatalogParser;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.analysis.CompoundFSSPlots.RupInRegionsCache;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.aftershockStatistics.ComcatAccessor;

public class LaHabraProbCalc {
	
	private static ObsEqkRupList fetchComcat(Region reg) {
		EventWebService service;
		try {
			service = new EventWebService(new URL("http://earthquake.usgs.gov/fdsnws/event/1/"));
		} catch (MalformedURLException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		EventQuery query = new EventQuery();
		
		query.setMinDepth(new BigDecimal(0d));
		query.setMaxDepth(new BigDecimal(100d));
		
		query.setStartTime(new GregorianCalendar(1900, 0, 1).getTime());
		query.setEndTime(new GregorianCalendar(2015, 10, 1).getTime());
		
		query.setMinLatitude(new BigDecimal(reg.getMinLat()));
		query.setMaxLatitude(new BigDecimal(reg.getMaxLat()));
		query.setMinLongitude(new BigDecimal(reg.getMinLon()));
		query.setMaxLongitude(new BigDecimal(reg.getMaxLon()));
		
		query.setMinMagnitude(new BigDecimal(5d));
		query.setMaxMagnitude(new BigDecimal(9d));
		
		try {
			System.out.println(service.getUrl(query, Format.GEOJSON));
		} catch (MalformedURLException e) {
			e.printStackTrace();
		}
		List<JsonEvent> events;
		try {
			events = service.getEvents(query);
		} catch (Exception e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		ObsEqkRupList rups = new ObsEqkRupList();
		for (JsonEvent event : events) {
			ObsEqkRupture rup = ComcatAccessor.eventToObsRup(event);
			rups.add(rup);
		}
		
		if (!reg.isRectangular()) {
			System.out.println("Fetched "+rups.size()+" events before region filtering");
			for (int i=rups.size(); --i>=0;)
				if (!reg.contains(rups.get(i).getHypocenterLocation()))
					rups.remove(i);
		}
		
		System.out.println("Returning "+rups.size()+" eqs");
		
		return rups;
	}
	
	private static void checkCatalog(Region reg, boolean comcat) throws IOException {
		ObsEqkRupList histQkList;
		if (comcat) {
			System.out.println("Fetching from ComCat");
			histQkList = fetchComcat(reg);
		} else {
			System.out.println("Loading UCERF3 file");
			File file = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/EarthquakeCatalog/ofr2013-1165_EarthquakeCat.txt");
			histQkList = UCERF3_CatalogParser.loadCatalog(file);
		}
		histQkList = histQkList.getRupsAboveMag(5d).getRupsInside(reg);
		histQkList.sortByOriginTime();
		
		long threshold_millis = (long)(3d*ProbabilityModelsCalc.MILLISEC_PER_YEAR);
		
		for (ObsEqkRupture rup : histQkList)
			System.out.println("M"+(float)rup.getMag()+" rup on "+new Date(rup.getOriginTime()));
		
		DecimalFormat df = new DecimalFormat("0.00");
		
		for (int i=histQkList.size(); --i>=1;) {
//		for (int i=1; i<histQkList.size(); i++) {
			long prevTime = histQkList.get(i-1).getOriginTime();
			long curTime = histQkList.get(i).getOriginTime();
			
			long delta = curTime - prevTime;
			
			if (delta > threshold_millis) {
				double timeYears = (double)delta/ProbabilityModelsCalc.MILLISEC_PER_YEAR;
				System.out.println(df.format(timeYears)+" yr period ending "+new Date(curTime));
			}
		}
	}

	public static void main(String[] args) throws IOException, DocumentException {
		Region reg = new Region(new Location(33.9225, -117.9352), 100d);
		checkCatalog(reg, true);
		checkCatalog(reg, false);
//		System.exit(0);
		
		FaultSystemSolution sol = FaultSystemIO.loadSol(new File(
				"/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
//				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_2_MEAN_BRANCH_AVG_SOL.zip"));
		
		// 1 month, 1 year, 3 year
		String[] durStrs = { "1 Mo", "1 Yr", "3 Yr" };
		double[] durations = { 1/12d, 1d, 3d };
		double[] mags = { 5d, 6d, 7d, 8d };
		
		boolean nucleation = true;
		
		double[][] table = new double[durations.length][mags.length];
		
		RupInRegionsCache cache = null;
		if (!nucleation)
			cache = new RupInRegionsCache();
//		RupInRegionsCache cache = null;
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
//		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.ONLY);
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
		erf.getTimeSpan().setStartTime(2014);
//		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		
		double minX = 5.05;
		int num = 50;
		double delta = 0.1;
		
		for (int i=0; i<durations.length; i++) {
			System.out.println("Calculating "+durStrs[i]);
			double duration = durations[i];
			erf.getTimeSpan().setDuration(duration);
			erf.updateForecast();
			
			SummedMagFreqDist incrMFD;
			if (nucleation)
				incrMFD = ERF_Calculator.getMagFreqDistInRegionFaster(erf, reg, minX, num, delta, true);
			else
				incrMFD = ERF_Calculator.getParticipationMagFreqDistInRegion(erf, reg, minX, num, delta, true, cache);
//			System.out.println(incrPartMFD);
			
			EvenlyDiscretizedFunc cmlMFD = incrMFD.getCumRateDistWithOffset();
//			System.out.println(cmlPartMFD);
			EvenlyDiscretizedFunc cmlMPD =
					FaultSysSolutionERF_Calc.calcProbsFromSummedMFD(cmlMFD, duration);
			
			for (int m=0; m<mags.length; m++) {
//				int closest = cmlPartMFD.getClosestXIndex(mags[m]);
//				System.out.println(mags[m]+" Index: "+cmlPartMFD.getXIndex(mags[m]));
//				System.out.println("Closest: "+closest+": "+cmlPartMFD.getX(closest));
				table[i][m] = cmlMPD.getY(cmlMPD.getClosestXIndex(mags[m]));
			}
		}
		
		DecimalFormat df = new DecimalFormat("0.00%");
		
		System.out.print("Mag");
		for (String durStr : durStrs)
			System.out.print("\t"+durStr);
		System.out.println();
		for (int m=0; m<mags.length; m++) {
//			System.out.print("M≥"+(float)mags[m]);
			System.out.print("M>"+(float)mags[m]);
			for (int d=0; d<durations.length; d++)
				System.out.print("\t"+df.format(table[d][m]));
			System.out.println();
		}
	}

}
