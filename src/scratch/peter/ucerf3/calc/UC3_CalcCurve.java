package scratch.peter.ucerf3.calc;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.util.ClassUtils;
import org.opensha.nshmp2.calc.ERF_ID;
import org.opensha.nshmp2.calc.HazardResultWriterSites;
import org.opensha.nshmp2.calc.ThreadedHazardCalc;
import org.opensha.nshmp2.util.Period;
import org.opensha.sha.earthquake.EpistemicListERF;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;

import scratch.UCERF3.erf.FaultSystemSolutionERF;

import com.google.common.base.Enums;
import com.google.common.base.Splitter;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

/**
 * Curves at sites for a single solution.
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class UC3_CalcCurve {

	private static final String S = File.separator;

	public UC3_CalcCurve(String solSetPath, String sitePath, String outDir,
		List<Period> periods, boolean epiUncert)
			throws IOException, InterruptedException, ExecutionException {

		FaultSystemSolutionERF erf = UC3_CalcUtils.getUC3_ERF(solSetPath,
				IncludeBackgroundOption.INCLUDE, false, true, 1.0);
		erf.updateForecast();
		EpistemicListERF wrappedERF = ERF_ID.wrapInList(erf);
		Map<String, Location> siteMap = UC3_CalcUtils.readSiteFile(sitePath);
		LocationList locs = new LocationList();
		for (Location loc : siteMap.values()) {
			locs.add(loc);
		}
		
		for (Period period : periods) {
			String outPath = outDir + S + erf.getName();
			System.out.println(outPath);
			HazardResultWriterSites writer = new HazardResultWriterSites(outPath,
				siteMap);
			writer.writeHeader(period);
			ThreadedHazardCalc thc = new ThreadedHazardCalc(wrappedERF, null, locs,
				period, epiUncert, writer, false);
			thc.calculate(null);
		}
	}


	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length != 4) {
			System.out.println("USAGE: " +
					ClassUtils.getClassNameWithoutPackage(UC3_CalcCurve.class) +
					" <filepath> <sitefile> <periods> <outDir>");
			System.exit(1);
		}
		
		String solSetPath = args[0];
		String siteFile = args[1];
		List<Period> periods = readArgAsList(args[2], Period.class);
		String outDir = args[3];
		boolean epiUnc = false;

		try {
			new UC3_CalcCurve(solSetPath, siteFile, outDir, periods, epiUnc);
		} catch (Exception ioe) {
			ioe.printStackTrace();
		}
	}
	
	private static final Splitter SPLIT = Splitter.on(',');

	private static <T extends Enum<T>> List<T> readArgAsList(String arg,
			Class<T> clazz) {
		Iterable<T> it = Iterables.transform(SPLIT.split(arg),
			Enums.stringConverter(clazz));
		return Lists.newArrayList(it);
	}


	
}
