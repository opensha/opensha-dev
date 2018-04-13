package scratch.kevin.simCompare;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.compress.utils.Lists;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

import com.google.common.base.Preconditions;

import scratch.kevin.util.MarkdownUtils;

public abstract class SiteHazardCurveComarePageGen<E> {
	private SimulationRotDProvider<E> simProv;
	SimulationHazardCurveCalc<E> simCalc;
	private String simName;
	
	private static double[] gmpe_truncs = { 3d, 2d, 1d };
	private static double[] gmpe_fixed_sigmas = { 0.3, 0d };
	
	private static ExecutorService exec;
	
	private List<SimulationRotDProvider<?>> compSimProvs;
	private List<SimulationHazardCurveCalc<?>> compCurveCals;
	
	private static Map<AttenRelRef, LinkedList<ScalarIMR>> gmpesInstancesCache;

	public SiteHazardCurveComarePageGen(SimulationRotDProvider<E> simProv, String simName) {
		this(simProv, simName, new ArrayList<>());
	}

	public SiteHazardCurveComarePageGen(SimulationRotDProvider<E> simProv, String simName,
			SimulationRotDProvider<?>... compSimProvs) {
		this(simProv, simName, toList(compSimProvs));
	}

	public SiteHazardCurveComarePageGen(SimulationRotDProvider<E> simProv, String simName,
			List<SimulationRotDProvider<?>> compSimProvs) {
		super();
		this.simProv = simProv;
		simCalc = new SimulationHazardCurveCalc<>(simProv);
		this.simName = simName;
		if (compSimProvs == null)
			compSimProvs = new ArrayList<>();
		this.compSimProvs = compSimProvs;
		compCurveCals = new ArrayList<>();
		for (SimulationRotDProvider<?> compSimProv : compSimProvs)
			compCurveCals.add(new SimulationHazardCurveCalc<>(compSimProv));
		
		gmpesInstancesCache = new HashMap<>();
	}
	
	private static List<SimulationRotDProvider<?>> toList(SimulationRotDProvider<?>... compSimProvs) {
		if (compSimProvs == null || compSimProvs.length == 0)
			return new ArrayList<>();
		List<SimulationRotDProvider<?>> list = new ArrayList<>();
		for (SimulationRotDProvider<?> simProv : compSimProvs)
			list.add(simProv);
		return list;
	}
	
	public void generateSitePage(Site site, List<? extends RuptureComparison<E>> comps, File outputDir, List<String> headerLines,
			double[] periods, AttenRelRef gmpeRef) throws IOException {
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		List<String> lines = new LinkedList<>();
		if (headerLines != null && !headerLines.isEmpty()) {
			lines.addAll(headerLines);
			if (!lines.get(lines.size()-1).isEmpty())
				lines.add("");
		}
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		lines.add("## Hazard Curves");
		lines.add(topLink); lines.add("");
		
		double curveDuration = 1d;
		
		List<String> curveNames = new ArrayList<>();
		curveNames.add(simName);
		for (SimulationRotDProvider<?> compProv : compSimProvs)
			curveNames.add(compProv.getName());
		
		lines.addAll(MultiRupGMPE_ComparePageGen.getCurveLegend(curveNames, gmpeRef.getShortName(), gmpe_truncs, gmpe_fixed_sigmas));
		lines.add("");
		
		for (double period : periods) {
			System.out.println("Calculating primary hazard curve");
			List<DiscretizedFunc> simCurves = new ArrayList<>();
			
			DiscretizedFunc simCurve = simCalc.calc(site, period, curveDuration);
			simCurves.add(simCurve);
			
			for (int i=0; i<compCurveCals.size(); i++) {
				SimulationHazardCurveCalc<?> calc = compCurveCals.get(i);
				System.out.println("Calculating for "+compSimProvs.get(i).getName());
				simCurves.add(calc.calc(site, period, curveDuration));
			}
			
			String prefix = site.getName().replaceAll(" ", "_")+"_curves_"+(float)period+"s_"+gmpeRef.getShortName();
			
			File curvePlot = MultiRupGMPE_ComparePageGen.plotHazardCurve(simCurves, comps, simCalc.getXVals(), site, period,
					curveDuration, gmpeRef, gmpe_truncs, gmpe_fixed_sigmas, resourcesDir, prefix);
			
			lines.add("### "+optionalDigitDF.format(period)+"s Hazard Curves");
			lines.add(topLink); lines.add("");
			lines.add("![Hazard Curve]("+resourcesDir.getName()+"/"+curvePlot.getName()+")");
			lines.add("");
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2, 4));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	protected synchronized static ExecutorService getExec() {
		if (exec == null)
			exec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		return exec;
	}
	
	protected static  ScalarIMR checkOutGMPE(AttenRelRef gmpeRef) {
		synchronized (SiteHazardCurveComarePageGen.class) {
			if (gmpesInstancesCache == null)
				gmpesInstancesCache = new HashMap<>();
		}
		synchronized (gmpesInstancesCache) {
			LinkedList<ScalarIMR> gmpes = gmpesInstancesCache.get(gmpeRef);
			if (gmpes == null) {
				gmpes = new LinkedList<>();
				gmpesInstancesCache.put(gmpeRef, gmpes);
			}
			if (!gmpes.isEmpty())
				return gmpes.pop();
		}
		ScalarIMR gmpe = gmpeRef.instance(null);
		gmpe.setParamDefaults();
		gmpe.setIntensityMeasure(SA_Param.NAME);
		return gmpe;
	}
	
	protected static void checkInGMPE(AttenRelRef gmpeRef, ScalarIMR gmpe) {
		synchronized (gmpesInstancesCache) {
			gmpesInstancesCache.get(gmpeRef).push(gmpe);
		}
	}
	
	static final DecimalFormat optionalDigitDF = new DecimalFormat("0.##");
}
