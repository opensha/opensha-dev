package scratch.kevin.ucerf3;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;

import scratch.UCERF3.griddedSeismicity.FaultPolyMgr;
import scratch.UCERF3.griddedSeismicity.GridSourceProvider;
import scratch.UCERF3.inversion.InversionTargetMFDs;
import scratch.UCERF3.inversion.U3InversionTargetMFDs;
import scratch.UCERF3.utils.U3FaultSystemIO;

public class FaultMFDCalc {

	public static void main(String[] args) throws IOException, DocumentException {
		File outputDir = new File("/tmp/parsons");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		FaultSystemSolution sol = U3FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		FaultSystemRupSet rupSet = sol.getRupSet();
		GridSourceProvider gridProv = sol.getGridSourceProvider();
		FaultGridAssociations polyMgr = FaultPolyMgr.create(rupSet.getFaultSectionDataList(), U3InversionTargetMFDs.FAULT_BUFFER);
		
		Map<String, int[]> parentIDsMap = new HashMap<>();
		
		boolean participation = true;
		
		parentIDsMap.put("San Andreas", new int[] {287, 300, 285, 295, 658, 301, 654, 653, 32, 655, 282, 283, 284, 657});
		parentIDsMap.put("Hayward Calaveras Rogers Creek", new int[] {639, 638, 637, 602, 601, 621, 603, 651});
		parentIDsMap.put("Garlock", new int[] {641, 48, 49});
		
		for (String faultName : parentIDsMap.keySet()) {
			System.out.println(faultName);
			int[] parentIDs = parentIDsMap.get(faultName);
			HashSet<Integer> parentIDsSet = new HashSet<>();
			for (int parentID : parentIDs)
				parentIDsSet.add(parentID);
			
			HashSet<Integer> subSects = new HashSet<>();
			for (FaultSection sect : rupSet.getFaultSectionDataList())
				if (parentIDsSet.contains(sect.getParentSectionId()))
					subSects.add(sect.getSectionId());
			
			HashSet<Integer> rups = new HashSet<>();
			for (int s : subSects)
				rups.addAll(rupSet.getRupturesForSection(s));
			
			IncrementalMagFreqDist faultMFD = new IncrementalMagFreqDist(5.05, 41, 0.1);
			for (int r : rups) {
				double rate = sol.getRateForRup(r);
				if (!participation) {
					// convert to fault nucleation rate
					double areaOn = 0;
					double areaOff = 0;
					for (Integer s : rupSet.getSectionsIndicesForRup(r)) {
						double area = rupSet.getAreaForSection(s);
						if (subSects.contains(s))
							areaOn += area;
						else
							areaOff += area;
					}
					rate *= areaOn / (areaOn + areaOff);
				}
				faultMFD.add(faultMFD.getClosestXIndex(rupSet.getMagForRup(r)), rate);
			}
			faultMFD.setName("Supra-Seimogenic");
			
			IncrementalMagFreqDist offMFD = new IncrementalMagFreqDist(5.05, 41, 0.1);
			offMFD.setName("Sub-Seimogenic");
			int numMatches = 0;
			double totFractMatch = 0;
			for (int i=0; i<gridProv.size(); i++) {
				IncrementalMagFreqDist subSeisMFD = gridProv.getNodeSubSeisMFD(i);
				if (subSeisMFD == null)
					continue;
				double fractMatch = 0d;
				Map<Integer, Double> sectFracsOnNode = polyMgr.getSectionFracsOnNode(i);
				for (Integer sectIndex : sectFracsOnNode.keySet())
					if (subSects.contains(sectIndex))
						fractMatch += sectFracsOnNode.get(sectIndex);
				if (fractMatch > 0) {
					numMatches++;
					totFractMatch += fractMatch;
					
					for (Point2D pt : subSeisMFD)
						if (pt.getX() > offMFD.getMinX() - 0.5*offMFD.getDelta())
							offMFD.add(pt.getX(), pt.getY()*fractMatch);
				}
			}
			System.out.println(numMatches+" node matches, "+(float)totFractMatch+" fractional");
			
			SummedMagFreqDist totalMFD = new SummedMagFreqDist(faultMFD.getMinX(), faultMFD.getMaxX(), faultMFD.size());
			totalMFD.setName("Total MFD");
			totalMFD.addIncrementalMagFreqDist(faultMFD);
			totalMFD.addIncrementalMagFreqDist(offMFD);
			
			for (boolean cumulative : new boolean [] {true, false}) {
				List<DiscretizedFunc> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				String yAxisLabel;
				if (participation)
					yAxisLabel = "Participation Rate (1/yr)";
				else
					yAxisLabel = "Nucleation Rate (1/yr)";
				
				if (cumulative) {
					funcs.add(offMFD.getCumRateDistWithOffset());
					funcs.add(faultMFD.getCumRateDistWithOffset());
					funcs.add(totalMFD.getCumRateDistWithOffset());
					yAxisLabel = "Cumulative "+yAxisLabel;
				} else {
					funcs.add(offMFD);
					funcs.add(faultMFD);
					funcs.add(totalMFD);
					yAxisLabel = "Incremental "+yAxisLabel;
				}
				
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
				
				PlotSpec spec = new PlotSpec(funcs, chars, faultName, "Magnitude", yAxisLabel);
				spec.setLegendVisible(true);
				
				HeadlessGraphPanel gp = new HeadlessGraphPanel();
				gp.setBackgroundColor(Color.WHITE);
				gp.setTickLabelFontSize(18);
				gp.setAxisLabelFontSize(20);
				gp.setPlotLabelFontSize(21);
				
				String prefix = new File(outputDir, faultName.replaceAll(" ", "_")).getAbsolutePath();
				if (participation)
					prefix += "_participation";
				else
					prefix += "_nucleation";
				if (cumulative)
					prefix += "_cumulative";
				else
					prefix += "_incremental";
				
				gp.drawGraphPanel(spec, false, true, new Range(5d, 9d), new Range(1e-8, 1d));
				gp.getChartPanel().setSize(1000, 800);
				gp.saveAsPNG(prefix+".png");
				gp.saveAsPDF(prefix+".pdf");
				gp.saveAsTXT(prefix+".txt");
			}
		}
	}

}
