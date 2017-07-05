package scratch.peter.ucerf3.calc;

import static com.google.common.base.Charsets.US_ASCII;
import static org.opensha.nshmp.NEHRP_TestCity.SANTA_ROSA;
import static org.opensha.nshmp2.util.Period.GM0P00;
import static org.opensha.nshmp2.util.Period.GM0P20;
import static org.opensha.nshmp2.util.Period.GM1P00;
import static org.opensha.sha.earthquake.param.IncludeBackgroundOption.EXCLUDE;
import static scratch.peter.curves.ProbOfExceed.PE10IN50;
import static scratch.peter.curves.ProbOfExceed.PE2IN50;

import java.io.File;
import java.io.IOException;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.io.IOUtils;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.util.DataUtils;
import org.opensha.nshmp.NEHRP_TestCity;
import org.opensha.nshmp2.calc.ERF_ID;
import org.opensha.nshmp2.calc.HazardCalc;
import org.opensha.nshmp2.calc.HazardResult;
import org.opensha.nshmp2.erf.NSHMP2008;
import org.opensha.nshmp2.util.Period;
import org.opensha.nshmp2.util.SourceIMR;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.EpistemicListERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.imr.ScalarIMR;

import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.peter.curves.CurveUtilsUC33;

import com.google.common.base.CharMatcher;
import com.google.common.base.Functions;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.Collections2;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Files;

/**
 * Add comments here
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class UC33_NicoUtils {
	private static Splitter SPLIT = Splitter.on(CharMatcher.WHITESPACE).omitEmptyStrings();
	private static Joiner JOIN = Joiner.on(',');
	private static final String LF = IOUtils.LINE_SEPARATOR;
	private static Set<NEHRP_TestCity> cities = NEHRP_TestCity.getCA();
//	private static Set<NEHRP_TestCity> cities = EnumSet.of(LOS_ANGELES); //, LOS_ANGELES, CENTURY_CITY, NORTHRIDGE);
	
	private static Set<Period> periods = EnumSet.of(GM0P00, GM0P20, GM1P00);
	private static String solDir = "tmp/UC33/src/bravg/FM-DM/";

	// use 8 branch averaged solutions to compute curves for NEHRP cities
	// in CA
	private static void buildCityCurves() throws IOException {
		Map<String, Double> fmWts = Maps.newHashMap();
		fmWts.put("FM31", 0.5);
		fmWts.put("FM32", 0.5);
		Map<String, Double> dmWts = Maps.newHashMap();
		dmWts.put("ABM", 0.1);
		dmWts.put("GEOL", 0.3);
		dmWts.put("NEOK", 0.3);
		dmWts.put("ZENGBB", 0.3);

		Map<Period, Map<NEHRP_TestCity, DiscretizedFunc>> map = Maps.newHashMap();
		
		for (String fm : fmWts.keySet()) {
			for (String dm : dmWts.keySet()) {
				double wt = fmWts.get(fm) * dmWts.get(dm);
				String solPath = solDir + "UC33brAvg_" + fm + "_" + dm + ".zip";
				
				// init erf
				FaultSystemSolutionERF erf = UC3_CalcUtils.getUC3_ERF(solPath,
					IncludeBackgroundOption.INCLUDE, false, true, 1.0);
				erf.updateForecast();
				EpistemicListERF wrappedERF = ERF_ID.wrapInList(erf);
				SourceIMR imr = SourceIMR.WUS_FAULT_14;
				System.out.println("Starting erf: " + solPath);
				
				for (Period p : periods) {
					System.out.println("   Starting period: " + p);
					Map<NEHRP_TestCity, DiscretizedFunc> cityMap = map.get(p);
					if (cityMap == null) {
						cityMap = Maps.newHashMap();
						map.put(p, cityMap);
					}
					for (NEHRP_TestCity city : cities) {
						System.out.println("      Starting city: " + city);
						Site site = new Site(city.location());
						HazardCalc hc = HazardCalc.create(wrappedERF, imr, site, p, true, false);
						HazardResult result = hc.call();
						
						DiscretizedFunc addFunc = result.curve().deepClone();
						addFunc.scale(wt); // scale by erf weight
						
						DiscretizedFunc mapFunc = cityMap.get(city);
						if (mapFunc == null) {
							cityMap.put(city, addFunc);
						} else {
							for (int i=0; i<mapFunc.size(); i++) {
								mapFunc.set(i, mapFunc.getY(i) + addFunc.getY(i));
							}
						}
					}
				}
			}
		}
		
		// write curve files
		for (Period p : map.keySet()) {
			File outFile = new File("tmp/forNico", p + "curves.csv");
			String header = createHeader(p);
			Files.write(header, outFile, US_ASCII);

			Map<NEHRP_TestCity, DiscretizedFunc> cityMap = map.get(p);
			for (NEHRP_TestCity city : cityMap.keySet()) {
				DiscretizedFunc cityCurve = cityMap.get(city);
				
				// add cascadia curve
				DiscretizedFunc cascFunc = readCurve(city, p, "casc");
				for (int i=0; i<cascFunc.size(); i++) {
					cityCurve.set(i, cityCurve.getY(i) + cascFunc.getY(i));
				}
					
				// add deep curve
				DiscretizedFunc deepFunc = readCurve(city, p, "deep");
				for (int i=0; i<deepFunc.size(); i++) {
					cityCurve.set(i, cityCurve.getY(i) + deepFunc.getY(i));
				}
				
				double pe2in50 = CurveUtilsUC33.getPE(cityCurve, PE2IN50);
				double pe10in50 = CurveUtilsUC33.getPE(cityCurve, PE10IN50);
				double rtgm = CurveUtilsUC33.getRTGM(cityCurve, p);
				
				String curveLine = createResult(city, pe2in50, pe10in50, rtgm, cityCurve);
				Files.append(curveLine, outFile, US_ASCII);				
			}
		}
	}
	
	private static void buildCityCurveUCERF2w2104() throws IOException {

		Map<Period, Map<NEHRP_TestCity, DiscretizedFunc>> map = Maps.newHashMap();
		
		// init erf
		EpistemicListERF erf = NSHMP2008.createCalifornia();
		erf.updateForecast();
		SourceIMR imr = SourceIMR.WUS_FAULT_14;
		System.out.println("Starting erf: " + erf.getName());
		
		ExecutorService ex = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		ExecutorCompletionService<HazardResult> ecs = new ExecutorCompletionService<HazardResult>(ex);

		for (Period p : periods) {
			System.out.println("    Starting period: " + p);
			Map<NEHRP_TestCity, DiscretizedFunc> cityMap = map.get(p);
			if (cityMap == null) {
				cityMap = Maps.newHashMap();
				map.put(p, cityMap);
			}
			
			for (NEHRP_TestCity city : cities) {
				System.out.println("        Submitted city: " + city);
				Site site = new Site(city.location());
				HazardCalc hc = HazardCalc.create(erf, imr, site, p, true, false);
				ecs.submit(hc);
			}
			
			for (int i = 0; i < cities.size(); i++) {
				try {
					HazardResult result = ecs.take().get();
					NEHRP_TestCity city = NEHRP_TestCity.forLocation(result.location());
					System.out.println("        Completed city: " + city);
					cityMap.put(city, result.curve());
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
		
		ex.shutdown();

		
		// write curve files
		for (Period p : map.keySet()) {
			File outFile = new File("tmp/forNico", p + "curves.csv");
			String header = createHeader(p);
			Files.write(header, outFile, US_ASCII);

			Map<NEHRP_TestCity, DiscretizedFunc> cityMap = map.get(p);
			for (NEHRP_TestCity city : cities) { //cityMap.keySet()) {
				DiscretizedFunc cityCurve = cityMap.get(city);
				
				// add cascadia curve
				DiscretizedFunc cascFunc = readCurve(city, p, "casc");
				for (int i=0; i<cascFunc.size(); i++) {
					cityCurve.set(i, cityCurve.getY(i) + cascFunc.getY(i));
				}
					
				// add deep curve
				DiscretizedFunc deepFunc = readCurve(city, p, "deep");
				for (int i=0; i<deepFunc.size(); i++) {
					cityCurve.set(i, cityCurve.getY(i) + deepFunc.getY(i));
				}
				
				double pe2in50 = CurveUtilsUC33.getPE(cityCurve, PE2IN50);
				double pe10in50 = CurveUtilsUC33.getPE(cityCurve, PE10IN50);
				double rtgm = CurveUtilsUC33.getRTGM(cityCurve, p);
				
				String curveLine = createResult(city, pe2in50, pe10in50, rtgm, cityCurve);
				Files.append(curveLine, outFile, US_ASCII);				
			}
		}
	}

	
	private static String createResult(NEHRP_TestCity city, double pe2in50,
			double pe10in50, double rtgm, DiscretizedFunc curve) {
		List<String> data = Lists.newArrayList(city.name());
		data.add(Double.toString(pe2in50));
		data.add(Double.toString(pe10in50));
		data.add(Double.toString(rtgm));
		Iterable<String> values = Collections2.transform(curve.yValues(),
			Functions.toStringFunction());
		String result = JOIN.join(Iterables.concat(data, values)) + LF;
		return result;
	}

//	private static void addFunc()
	private static String createHeader(Period period) {
		List<String> headerVals = Lists.newArrayList("loc", "2in50", "10in50", "rtgm");
		for (Double d : period.getIMLs()) {
			headerVals.add(d.toString());
		}
		return JOIN.join(headerVals) + LF;
	}

	
	private static String sol = "tmp/UC33/src/tree/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip";
	private static String FM31 = "FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3";
	private static String FM32 = "FM3_2_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3";
	
	private static void buildDeterministicLists() throws IOException {
		// both fault models
		String[] branches = new String[] {FM31, FM32};
		for (Period p : periods) {

			System.out.println("   Starting period: " + p);
			ScalarIMR imr = SourceIMR.WUS_FAULT_14.instance(p);

			File out = new File("tmp/forNico/det/" + p.name() + "top10.text");
			Files.write(LF, out, US_ASCII);

			Map<NEHRP_TestCity, Double> maxMean = Maps.newHashMap();
			Map<NEHRP_TestCity, Double> max_rRup = Maps.newHashMap(); // rRup for maxMean
			Map<NEHRP_TestCity, Double> max_mag = Maps.newHashMap(); // mag for maxMean
			
			for (String branch : branches) {
				ERF erf = UC3_CalcUtils.getUC3_ERF(sol, branch, EXCLUDE, false, true, 1.0);
				erf.updateForecast();
				
				for (NEHRP_TestCity city : cities) {
					System.out.println("      Starting city: " + city);
					Files.append(city.name() + "  " + branch + LF, out, US_ASCII);
					List<Double> means = Lists.newArrayList();
					List<Double> rRups = Lists.newArrayList();
					List<Double> mags = Lists.newArrayList();
					List<String> names = Lists.newArrayList();
					
					Site site = new Site(city.location());
					HazardCalc.initSite(site);
					imr.setSite(site);
					
					for (ProbEqkSource src : erf) {
						for (ProbEqkRupture rup : src) {
							imr.setEqkRupture(rup);
							means.add(imr.getMean());
							rRups.add(rup.getRuptureSurface().getDistanceRup(city.location()));
							mags.add(rup.getMag());
							names.add(src.getName());
						}
					}
					
					System.out.println("        ...processing results");
					List<Integer> indices = DataUtils.sortedIndices(means, false);
					
					// write top 10
					for (int i=0; i<10; i++) {
						int idx = indices.get(i);
						String line = Math.exp(means.get(idx)) + "\t" + 
									  mags.get(idx) + "\t" + 
									  rRups.get(idx) + "\t" +
									  names.get(idx) + LF;
						Files.append(line, out, US_ASCII);
						
					}
					Files.append(LF, out, US_ASCII);

					// save max mean and rRup
					int top = indices.get(0);
					Double max1 = maxMean.get(city); // already exp'd
					Double max2 = Math.exp(means.get(top));
					Double rRup = rRups.get(top);
					Double mag = mags.get(top);
					if (max1 == null || max2 > max1) {
						maxMean.put(city, max2);
						max_rRup.put(city, rRup);
						max_mag.put(city, mag);
					}
				}
			}
			
			// write max vals
			File maxOut = new File("tmp/forNico/det/" + p.name() + "maxVals.csv");
			Files.write("city, maxMedian,mag,rRup" + LF, maxOut, US_ASCII);
			for (NEHRP_TestCity city : cities) {
				Files.append(city + ", " + maxMean.get(city)
					              + ", " + max_mag.get(city)
					              + ", " + max_rRup.get(city) + LF, maxOut, US_ASCII);
			}
		}		
		
	}
	
	
	
	
	public static void main(String[] args) throws IOException {
//		DiscretizedFunc f = readCurve(NEHRP_TestCity.LOS_ANGELES, GM0P20, "casc");
//		System.out.println(f);
		
//		buildCityCurves();
//		buildCityCurveUCERF2w2104();
//		buildDeterministicLists();
		
		System.out.println(readCurve(SANTA_ROSA, GM0P00, "deep"));
	}
	
	private static String curveDir = "/Users/pmpowers/projects/NSHMP/tmp/UC3";
	
	// load Cascadia and CAdeep grids
	private static DiscretizedFunc readCurve(NEHRP_TestCity city, Period p,
			String id) throws IOException {
		String pString = (p == GM0P00) ? "pga" : (p == GM0P20) ? "5hz" : "1hz";
		String fName = city.name() + "." + id + "." + pString;
		File cFile = new File(curveDir, fName);
		List<String> lines = Files.readLines(cFile, US_ASCII);
		DiscretizedFunc f = new ArbitrarilyDiscretizedFunc();
		for (String line : lines) {
			if (line.trim().isEmpty()) continue;
			if (line.startsWith("#")) continue;
			Iterator<String> vals = SPLIT.split(line).iterator();
			double x = Double.valueOf(vals.next());
			double y = Double.valueOf(vals.next());
			f.set(x, y);
		}
		return f;
	}
}
