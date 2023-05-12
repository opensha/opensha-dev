package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.UncertainDataConstraint.SectMappedUncertainDataConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.PaleoseismicConstraintData;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionSlipRates;

import com.google.common.base.Preconditions;

public class SlipRateScattersPaleoSects {

	public static void main(String[] args) throws IOException {
		File mainDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		
		File solFile = new File(mainDir, "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip");
		
		File evenFile = new File(mainDir, "node_branch_averaged/PaleoUncert_EvenFitPaleo.zip");
		File underFile = new File(mainDir, "node_branch_averaged/PaleoUncert_UnderFitPaleo.zip");
		File overFile = new File(mainDir, "node_branch_averaged/PaleoUncert_OverFitPaleo.zip");
		
		File outputDir = new File(mainDir, "misc_plots/slip_rate_scatters_paleo_sects");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		List<String> prefixes = new ArrayList<>();
		List<String> names = new ArrayList<>();
		List<FaultSystemSolution> sols = new ArrayList<>();
		
		prefixes.add("full_ba");
		names.add("Full Branch Average");
		sols.add(FaultSystemSolution.load(solFile));
		
		prefixes.add("even_fit");
		names.add("Even-Fit Paleo Branch");
		sols.add(FaultSystemSolution.load(evenFile));
		
		prefixes.add("under_fit");
		names.add("Under-Fit Paleo Branch");
		sols.add(FaultSystemSolution.load(underFile));
		
		prefixes.add("over_fit");
		names.add("Over-Fit Paleo Branch");
		sols.add(FaultSystemSolution.load(overFile));
		
		for (int i=0; i<prefixes.size(); i++) {
			String prefix = prefixes.get(i);
			FaultSystemSolution sol = sols.get(i);
			FaultSystemRupSet rupSet = sol.getRupSet();
			String name = names.get(i);
			
			CSVFile<String> csv = new CSVFile<>(false);
			csv.addLine("Section Index", "Section Name", "Target Slip Rate (mm/yr)", "Solution Slip Rate (mm/yr)", "Parent Section Mapped Paleo Site(s)");
			
			SectSlipRates targets = rupSet.getSectSlipRates();
			SolutionSlipRates solRates = sol.requireModule(SolutionSlipRates.class);
			
			DefaultXY_DataSet scatter = new DefaultXY_DataSet();
			DefaultXY_DataSet scatterNoPaleo = new DefaultXY_DataSet();
			DefaultXY_DataSet scatterWithPaleo = new DefaultXY_DataSet();
			
			PaleoseismicConstraintData paleoData = rupSet.requireModule(PaleoseismicConstraintData.class);
			Map<Integer, List<SectMappedUncertainDataConstraint>> paleoParents = new HashMap<>();
			for (SectMappedUncertainDataConstraint constr : paleoData.getPaleoRateConstraints()) {
				int parentID = rupSet.getFaultSectionData(constr.sectionIndex).getParentSectionId();
				List<SectMappedUncertainDataConstraint> sectConstraints = paleoParents.get(parentID);
				if (sectConstraints == null) {
					sectConstraints = new ArrayList<>();
					paleoParents.put(parentID, sectConstraints);
				}
				sectConstraints.add(constr);
			}
			
			for (int s=0; s<rupSet.getNumSections(); s++) {
				double target = targets.getSlipRate(s)*1e3;
				double solSlip = solRates.get(s)*1e3;
				
				scatter.set(target, solSlip);
				
				List<SectMappedUncertainDataConstraint> sectConstraints = paleoParents.get(rupSet.getFaultSectionData(s).getParentSectionId());
				if (sectConstraints == null)
					scatterNoPaleo.set(target, solSlip);
				else
					scatterWithPaleo.set(target, solSlip);
				
				List<String> line = new ArrayList<>();
				line.add(s+"");
				line.add(rupSet.getFaultSectionData(s).getSectionName());
				line.add((float)target+"");
				line.add((float)solSlip+"");
				if (sectConstraints != null)
					for (SectMappedUncertainDataConstraint constr : sectConstraints)
						line.add(constr.name);
				
				csv.addLine(line);
			}
			
			csv.writeToFile(new File(outputDir, prefix+".csv"));
			
			for (boolean log : new boolean[] {false,true}) {
				Range range = log ? new Range(1e-3, 1e2) : new Range(0, 35);
				for (int j=0; j<3; j++) {
					DefaultXY_DataSet myScatter;
					String myPrefix, title;
					if (j == 0) {
						myScatter = scatter;
						myPrefix = prefix+"_all_sects";
						title = name+", All Subsections";
					} else if (j == 1) {
						myScatter = scatterNoPaleo;
						myPrefix = prefix+"_no_paleo";
						title = name+", Subsections Without Paleo Data";
					} else {
						myScatter = scatterWithPaleo;
						myPrefix = prefix+"_with_paleo";
						title = name+", Subsections With Paleo Data";
					}
					
					if (log)
						myPrefix += "_log";
					
					List<XY_DataSet> funcs = new ArrayList<>();
					List<PlotCurveCharacterstics> chars = new ArrayList<>();
					
					DefaultXY_DataSet oneToOne = new DefaultXY_DataSet();
					oneToOne.set(range.getLowerBound(), range.getLowerBound());
					oneToOne.set(range.getUpperBound(), range.getUpperBound());
					
					funcs.add(oneToOne);
					chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.GRAY));
					
					funcs.add(myScatter);
					chars.add(new PlotCurveCharacterstics(PlotSymbol.X, 3f, Color.BLACK));
					
					PlotSpec spec = new PlotSpec(funcs, chars, title, "Target Slip Rate (mm/yr)", "Solution Slip Rate (mm/yr)");
					
					HeadlessGraphPanel gp = PlotUtils.initHeadless();
					
					gp.drawGraphPanel(spec, log, log, range, range);
					
					PlotUtils.writePlots(outputDir, myPrefix, gp, 850, false, true, true, false);
				}
			}
		}
	}

}
