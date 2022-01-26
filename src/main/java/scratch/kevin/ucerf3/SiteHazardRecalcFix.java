package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipException;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.nshmp2.util.Period;
import org.opensha.sha.calc.hazardMap.BinaryHazardCurveReader;
import org.opensha.sha.calc.hazardMap.components.BinaryCurveArchiver;
import org.opensha.sha.calc.hazardMap.components.CurveMetadata;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.U3CompoundFaultSystemSolution;
import scratch.UCERF3.U3FaultSystemSolutionFetcher;
import scratch.UCERF3.analysis.CompoundFSSPlots;
import scratch.UCERF3.analysis.CompoundFSSPlots.ERFBasedSiteHazardHistPlot;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.logicTree.U3APrioriBranchWeightProvider;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;

public class SiteHazardRecalcFix {
	
	public static void main(String[] args) throws Exception {
		final U3FaultSystemSolutionFetcher origFetch = U3CompoundFaultSystemSolution.fromZipFile(
				new File("/tmp/comp_plots/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip"));
		List<U3LogicTreeBranch> origBranches = Lists.newArrayList(origFetch.getBranches());
		Collections.sort(origBranches);
		int replaceIndex = 95;
		final U3LogicTreeBranch fixBranch = origBranches.get(replaceIndex);
		
		boolean calc = false;
		
		
		U3FaultSystemSolutionFetcher fetch = new U3FaultSystemSolutionFetcher() {
			
			@Override
			public Collection<U3LogicTreeBranch> getBranches() {
				List<U3LogicTreeBranch> branches = Lists.newArrayList();
				branches.add(fixBranch);
				return branches;
			}
			
			@Override
			protected InversionFaultSystemSolution fetchSolution(U3LogicTreeBranch branch) {
				return origFetch.getSolution(fixBranch);
			}
		};
		
		File tempCurveDir = new File("/tmp/comp_plots/curvefix");
		if (!tempCurveDir.exists())
			tempCurveDir.mkdir();
		
		if (calc) {
			ERFBasedSiteHazardHistPlot plot = new ERFBasedSiteHazardHistPlot(
					new U3APrioriBranchWeightProvider(), tempCurveDir, 1);
			
			List<CompoundFSSPlots> plots = Lists.newArrayList();
			plots.add(plot);
			CompoundFSSPlots.batchPlot(plots, fetch);
		}
		
		File origCurveDir = new File("/tmp/comp_plots/site_hazard_curve_cache");
		
		List<Period> periods = ERFBasedSiteHazardHistPlot.getPeriods();
		
		for (Site site : ERFBasedSiteHazardHistPlot.getSites()) {
			System.out.println("Processing "+site.getName());
			File origDir = new File(origCurveDir, site.getName());
			File newDir = new File(origCurveDir, site.getName());
			
			Map<String, DiscretizedFunc> xValsMap = Maps.newHashMap();
			fileLoop:
			for (File file : newDir.listFiles()) {
				String name = file.getName();
				if (!name.endsWith(".bin"))
					continue;
				
				for (Period p : periods) {
					if (name.contains(p.getLabel())) {
						String prefix = name.replaceAll(".bin", "");
						xValsMap.put(prefix, p.getFunction());
						continue fileLoop;
					}
				}
			}
			
			BinaryCurveArchiver archiver = new BinaryCurveArchiver(origDir, origBranches.size(), xValsMap);
			
			for (File file : newDir.listFiles()) {
				String name = file.getName();
				if (!name.endsWith(".bin"))
					continue;
				BinaryHazardCurveReader read = new BinaryHazardCurveReader(file.getAbsolutePath());
				
				String prefix = name.replaceAll(".bin", "");
				
				ArbitrarilyDiscretizedFunc curve = read.nextCurve();
				
				archiver.archiveCurve(curve, new CurveMetadata(site, replaceIndex, null, prefix));
			}
			
			archiver.close();
		}
	}

}
