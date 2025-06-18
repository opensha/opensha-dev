package scratch.kevin.nshm23.devinSlipRateTests;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;

import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;

public class MultiFaultRateDebug {

	public static void main(String[] args) throws IOException {
		FaultSystemSolution sol1 = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_11_15-nshm23_branches-devin_tapered_slip_tests-WUS_FM_v3-GEO_AVG_FROM_DEVIN-UNIFORM/"
				+ "results_WUS_FM_v3_UNIFORM_branch_averaged.zip"));
		String name1 = "Average Slip Rates, Uniform Dsr";
		FaultSystemSolution sol2 = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_11_15-nshm23_branches-devin_tapered_slip_tests-WUS_FM_v3-GEO_FROM_DEVIN-TAPER_OVERRIDE_INDIVIDUAL/"
				+ "results_WUS_FM_v3_TAPER_OVERRIDE_INDIVIDUAL_branch_averaged.zip"));
		String name2 = "High Res Slip Rates, Multi-Rainbow Dsr";
		File outputDir = new File("/tmp/devin_mfds");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir(),
                "Output directory doesn't exist and could not be created: %s", outputDir.getAbsolutePath());
		
		FaultSystemRupSet rupSet = sol1.getRupSet();
		
		int id1 = FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Cucamonga");
		int id2 = FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Andreas", "Mojave", "south");
		String fault1 = null;
		String fault2 = null;
		for (FaultSection sect : rupSet.getFaultSectionDataList()) {
			if (sect.getParentSectionId() == id1)
				fault1 = sect.getParentSectionName();
			if (sect.getParentSectionId() == id2)
				fault2 = sect.getParentSectionName();
		}
		
		HashSet<Integer> rups1 = new HashSet<>(rupSet.getRupturesForParentSection(id1));
		HashSet<Integer> rups2 = new HashSet<>(rupSet.getRupturesForParentSection(id2));
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(5.01d, 8.45);
		
		HashSet<Integer> rupsBoth = new HashSet<>(rups1);
		rupsBoth.retainAll(rups2);
		
		IncrementalMagFreqDist[] fault1targets = new IncrementalMagFreqDist[2];
		IncrementalMagFreqDist[] fault2targets = new IncrementalMagFreqDist[2];
		IncrementalMagFreqDist[] fault1solPartics = new IncrementalMagFreqDist[2];
		IncrementalMagFreqDist[] fault2solPartics = new IncrementalMagFreqDist[2];
		IncrementalMagFreqDist[] fault1solNucls = new IncrementalMagFreqDist[2];
		IncrementalMagFreqDist[] fault2solNucls = new IncrementalMagFreqDist[2];
		IncrementalMagFreqDist[] combSolPartics = new IncrementalMagFreqDist[2];
		IncrementalMagFreqDist[] combSolNucls1 = new IncrementalMagFreqDist[2];
		IncrementalMagFreqDist[] combSolNucls2 = new IncrementalMagFreqDist[2];
		
		for (int n=0; n<2; n++) {
			FaultSystemSolution sol = n == 0 ? sol1 : sol2;
			double rate1 = calcRate(rups1, sol, 0d);
			double rate1_m7 = calcRate(rups1, sol, 7d);
			double moRate1 = calcMomentRate(rups1, sol, 0d);
			double rate2 = calcRate(rups2, sol, 0d);
			double rate2_m7 = calcRate(rups2, sol, 7d);
			double moRate2 = calcMomentRate(rups2, sol, 0d);
			double rateBoth = calcRate(rupsBoth, sol, 0d);
			double rateBothM7 = calcRate(rupsBoth, sol, 7d);
			double moRateBoth = calcMomentRate(rupsBoth, sol, 0d);
			
			DecimalFormat pDF = new DecimalFormat("0.00%");
			
			if (n == 0)
				System.out.println(name1);
			else
				System.out.println(name2);
			
			System.out.println("Rates "+fault1);
			System.out.println("\tSupra:\t"+(float)rate1);
			System.out.println("\tM>7:\t"+(float)rate1_m7);
			System.out.println("\tMo:\t"+(float)moRate1);
			System.out.println("Rates "+fault2);
			System.out.println("\tSupra:\t"+(float)rate2);
			System.out.println("\tM>7:\t"+(float)rate2_m7);
			System.out.println("\tMo:\t"+(float)moRate2);
			System.out.println("Rates Both");
			System.out.println("\tSupra:\t"+(float)rateBoth);
			System.out.println("\t\t"+pDF.format(rateBoth/rate1)+" of 1");
			System.out.println("\t\t"+pDF.format(rateBoth/rate2)+" of 2");
			System.out.println("\tM>7:\t"+(float)rateBothM7);
			System.out.println("\t\t"+pDF.format(rateBothM7/rate1_m7)+" of 1");
			System.out.println("\t\t"+pDF.format(rateBothM7/rate2_m7)+" of 2");
			System.out.println("\tMo:\t"+(float)moRateBoth);
			System.out.println("\t\t"+pDF.format(moRateBoth/moRate1)+" of 1");
			System.out.println("\t\t"+pDF.format(moRateBoth/moRate2)+" of 2");
			System.out.println();
			
			InversionTargetMFDs targetMFDs = sol.getRupSet().requireModule(InversionTargetMFDs.class);
			List<? extends IncrementalMagFreqDist> sectMFDs = targetMFDs.getOnFaultSupraSeisNucleationMFDs();
			
			SummedMagFreqDist target1 = new SummedMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
			SummedMagFreqDist target2 = new SummedMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
			for (FaultSection sect : rupSet.getFaultSectionDataList()) {
				if (sect.getParentSectionId() == id1)
					target1.addIncrementalMagFreqDist(sectMFDs.get(sect.getSectionId()));
				else if (sect.getParentSectionId() == id2)
					target2.addIncrementalMagFreqDist(sectMFDs.get(sect.getSectionId()));
			}
			
			fault1targets[n] = target1;
			fault2targets[n] = target2;
			
			fault1solPartics[n] = sol.calcParticipationMFD_forParentSect(id1, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
			fault2solPartics[n] = sol.calcParticipationMFD_forParentSect(id2, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
			fault1solNucls[n] = sol.calcNucleationMFD_forParentSect(id1, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
			fault2solNucls[n] = sol.calcNucleationMFD_forParentSect(id2, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
			
			combSolPartics[n] = sol.calcParticipationMFD_forRups(rupsBoth, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
			
			SummedMagFreqDist combSolNucl1 = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
			SummedMagFreqDist combSolNucl2 = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
			
			for (int sectIndex=0; sectIndex<rupSet.getNumSections(); sectIndex++) {
				int parent = rupSet.getFaultSectionData(sectIndex).getParentSectionId();
				if (parent == id1)
					combSolNucl1.addIncrementalMagFreqDist(sol.calcNucleationMFD_forSect(
							sectIndex, rupsBoth, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size()));
				else if (parent == id2)
					combSolNucl2.addIncrementalMagFreqDist(sol.calcNucleationMFD_forSect(
							sectIndex, rupsBoth, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size()));
			}
			combSolNucls1[n] = combSolNucl1;
			combSolNucls2[n] = combSolNucl2;
		}
		
		for (int f=0; f<2; f++) {
			IncrementalMagFreqDist[] targets = f == 0 ? fault1targets : fault2targets;
			IncrementalMagFreqDist[] partics = f == 0 ? fault1solPartics : fault2solPartics;
			IncrementalMagFreqDist[] nucls = f == 0 ? fault1solNucls : fault2solNucls;
			IncrementalMagFreqDist[] combNucls = f == 0 ? combSolNucls1 : combSolNucls2;
			String faultName = f == 0 ? fault1 : fault2;
			for (int i=0; i<3; i++) {
				IncrementalMagFreqDist[] mfds;
				IncrementalMagFreqDist[] coruptureMFDs;
				String title = faultName;
				String prefix;
				if (i == 0) {
					title += " Target Nucleation";
					prefix = "target_nucl";
					mfds = targets;
					coruptureMFDs = null;
				} else if (i == 1) {
					title += " Solution Participation";
					prefix = "sol_partic";
					mfds = partics;
					coruptureMFDs = combSolPartics;
				} else {
					title += " Solution Nucleation";
					prefix = "sol_nucl";
					mfds = nucls;
					coruptureMFDs = combNucls;
				}
				prefix += "_"+(f+1);
				
				List<IncrementalMagFreqDist> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();

				mfds[0].setName(name1);
				funcs.add(mfds[0]);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_orange));
				
				if (coruptureMFDs != null) {
					coruptureMFDs[0].setName("Corupture w/ "+(f == 0 ? fault2 : fault1));
					funcs.add(coruptureMFDs[0]);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Colors.tab_lightorange));
				}

				mfds[1].setName(name2);
				funcs.add(mfds[1]);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_blue));
				
				if (coruptureMFDs != null) {
					coruptureMFDs[1].setName(coruptureMFDs[0].getName());
					funcs.add(coruptureMFDs[1]);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Colors.tab_lightblue));
				}
				
				PlotSpec plot = new PlotSpec(funcs, chars, title, "Magnitude", "Incremental Rate (1/yr)");
				plot.setLegendInset(true);
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				
				gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
				
				gp.drawGraphPanel(plot, false, true, new Range(5d, 8.5d), new Range(1e-7, 1e-2));
				
				PlotUtils.writePlots(outputDir, prefix, gp, 800, 750, true, false, false);
				
				List<EvenlyDiscretizedFunc> cmlFuncs = new ArrayList<>();
				for (IncrementalMagFreqDist mfd : funcs)
					cmlFuncs.add(mfd.getCumRateDistWithOffset());
				
				plot = new PlotSpec(cmlFuncs, chars, title, "Magnitude", "Cumulative Rate (1/yr)");
				plot.setLegendInset(true);
				
				gp.drawGraphPanel(plot, false, true, new Range(5d, 8.5d), new Range(1e-7, 1e-2));
				
				PlotUtils.writePlots(outputDir, prefix+"_cml", gp, 800, 750, true, false, false);
			}
		}
	}
	
	private static double calcRate(Collection<Integer> rups, FaultSystemSolution sol, double minMag) {
		double sum = 0d;
		for (int rupIndex : rups)
			if (sol.getRupSet().getMagForRup(rupIndex) >= minMag)
				sum += sol.getRateForRup(rupIndex);
		return sum;
	}
	
	private static double calcMomentRate(Collection<Integer> rups, FaultSystemSolution sol, double minMag) {
		double sum = 0d;
		for (int rupIndex : rups) {
			double mag = sol.getRupSet().getMagForRup(rupIndex);
			if (mag >= minMag)
				sum += MagUtils.magToMoment(mag)*sol.getRateForRup(rupIndex);
		}
		return sum;
	}

}
