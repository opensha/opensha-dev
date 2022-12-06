package scratch.kevin.nshm23;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.uncertainty.BoundedUncertainty;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.geo.LocationUtils;
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

import com.google.common.base.Preconditions;

public class PaleoSlipUncertScatterPlot {
	
	public static void main(String[] args) throws IOException {
		File inputDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/2022_11_22-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		
		File solFile = new File(inputDir, "results_NSHM23_v2_CoulombRupSet_branch_averaged.zip");
		
		File outputDir = new File(inputDir, "misc_plots");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		FaultSystemSolution sol = FaultSystemSolution.load(solFile);
		
		boolean parents = true;
		
		PaleoseismicConstraintData paleoData = sol.getRupSet().requireModule(PaleoseismicConstraintData.class);
		
		Map<Integer, List<SectMappedUncertainDataConstraint>> riDataMap = sectMapped(sol.getRupSet(), parents,
				paleoData.getPaleoRateConstraints());
		List<? extends SectMappedUncertainDataConstraint> origSlip = paleoData.getPaleoSlipConstraints();
		List<SectMappedUncertainDataConstraint> inferredSlip = PaleoseismicConstraintData.inferRatesFromSlipConstraints(
				sol.getRupSet(), origSlip, true);
		Map<SectMappedUncertainDataConstraint, SectMappedUncertainDataConstraint> inferredToOrig = new HashMap<>();
		for (int i=0; i<inferredSlip.size(); i++)
			inferredToOrig.put(inferredSlip.get(i), origSlip.get(i));
		Map<Integer, List<SectMappedUncertainDataConstraint>> slipDataMap = sectMapped(sol.getRupSet(), parents,
				inferredSlip);
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		DefaultXY_DataSet scatter = new DefaultXY_DataSet();
		
		for (int parentID : slipDataMap.keySet()) {
			List<SectMappedUncertainDataConstraint> slipDatas = slipDataMap.get(parentID);
			List<SectMappedUncertainDataConstraint> riDatas = riDataMap.get(parentID);
			if (riDatas == null)
				continue;
			for (SectMappedUncertainDataConstraint slipData : slipDatas) {
				double minDist = Double.POSITIVE_INFINITY;
				SectMappedUncertainDataConstraint closestRI = null;
				for (SectMappedUncertainDataConstraint riData : riDatas) {
					double dist = LocationUtils.horzDistanceFast(slipData.dataLocation, riData.dataLocation);
					if (dist < minDist) {
						minDist = dist;
						closestRI = riData;
					}
				}
				double x = closestRI.bestEstimate;
				double y = slipData.bestEstimate;
				plotWhiskers(x, y, true, closestRI, funcs, chars);
				plotWhiskers(x, y, false, slipData, funcs, chars);
				scatter.set(x, y);
				System.out.println("Paleo Slip Site: "+slipData.getName());
				System.out.println("\tSlip Sect Mapping: "+slipData.sectionName);
				System.out.println("\tSlip Sect Rate (RI): "+(float)slipData.bestEstimate+" ("+(int)(1d/slipData.bestEstimate + 0.5)+" yrs)");
				System.out.println("\tSlip Fractional SD: "+(float)(slipData.getPreferredStdDev()/slipData.bestEstimate));
				System.out.println("\tRate Sect Mapping: "+closestRI.sectionName);
				System.out.println("\tRate Sect Rate (RI): "+(float)closestRI.bestEstimate+" ("+(int)(1d/closestRI.bestEstimate + 0.5)+" yrs)");
				System.out.println("\tRate Fractional SD: "+(float)(closestRI.getPreferredStdDev()/closestRI.bestEstimate));
//				System.out.println("\tProxy Slip Inferrences:");
//				SectMappedUncertainDataConstraint orig
//				System.out.println("\t\t")
			}
		}
		
//		double max = 0d;
//		double min = Double.POSITIVE_INFINITY;
//		for (XY_DataSet xy : funcs) {
//			for (Point2D pt : xy) {
//				if (pt.getX() > 0) {
//					min = Math.min(pt.getX(), min);
//					max = Math.max(pt.getX(), max);
//				}
//				if (pt.getY() > 0) {
//					min = Math.min(pt.getY(), min);
//					max = Math.max(pt.getY(), max);
//				}
//			}
//		}
//		min = Math.max(min, 1e-4);
//		Range range = new Range(Math.pow(10, Math.floor(Math.log10(min))), Math.pow(10, Math.ceil(Math.log10(max))));
		Range range = new Range(1e-4, 2e-2);
		
		DefaultXY_DataSet oneToOne = new DefaultXY_DataSet();
		oneToOne.set(range.getLowerBound(), range.getLowerBound());
		oneToOne.set(range.getUpperBound(), range.getUpperBound());
		funcs.add(oneToOne);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		
		funcs.add(scatter);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 4f, Color.BLACK));
		
		String title;
		if (parents)
			title = "Parent Section Paleo Data Comparison";
		else
			title = "Subsection Paleo Data Comparison";
		PlotSpec spec = new PlotSpec(funcs, chars, title, "Paleo RI Data Rate (/yr)", "Paleo Slip Data Proxy Rate (/yr)");
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec, true, true, range, range);
		
		if (parents)
			PlotUtils.writePlots(outputDir, "parent_sect_paleo_data_scatter", gp, 1000, false, true, true, false);
		else
			PlotUtils.writePlots(outputDir, "sect_paleo_data_scatter", gp, 1000, false, true, true, false);
	}
	
	private static final double logWhiskerDelta95 = 0.03;
	private static final double logWhiskerDelta68 = 0.04;
	
	private static void plotWhiskers(double x, double y, boolean horizontal, SectMappedUncertainDataConstraint constr,
			List<XY_DataSet> funcs, List<PlotCurveCharacterstics> chars) {
		BoundedUncertainty conf68 = constr.estimateUncertaintyBounds(UncertaintyBoundType.CONF_68);
		BoundedUncertainty conf95 = constr.estimateUncertaintyBounds(UncertaintyBoundType.CONF_95);
		double lower68 = conf68.lowerBound;
		double upper68 = conf68.upperBound;
		double lower95 = conf95.lowerBound;
		double upper95 = conf95.upperBound;
		
		double whiskerAbove95, whiskerBelow95, whiskerAbove68, whiskerBelow68;
		DefaultXY_DataSet confRange95 = new DefaultXY_DataSet();
		DefaultXY_DataSet confRange68 = new DefaultXY_DataSet();
		
		if (horizontal) {
			whiskerAbove95 = Math.pow(10, Math.log10(y)+logWhiskerDelta95);
			whiskerBelow95 = Math.pow(10, Math.log10(y)-logWhiskerDelta95);
			confRange95.set(lower95, whiskerAbove95);
			confRange95.set(lower95, whiskerBelow95);
			confRange95.set(lower95, y);
			confRange95.set(upper95, y);
			confRange95.set(upper95, whiskerAbove95);
			confRange95.set(upper95, whiskerBelow95);
			
			whiskerAbove68 = Math.pow(10, Math.log10(y)+logWhiskerDelta68);
			whiskerBelow68 = Math.pow(10, Math.log10(y)-logWhiskerDelta68);
			confRange68.set(lower68, whiskerAbove68);
			confRange68.set(lower68, whiskerBelow68);
			confRange68.set(lower68, y);
			confRange68.set(upper68, y);
			confRange68.set(upper68, whiskerAbove68);
			confRange68.set(upper68, whiskerBelow68);
		} else {
			whiskerAbove95 = Math.pow(10, Math.log10(x)+logWhiskerDelta95);
			whiskerBelow95 = Math.pow(10, Math.log10(x)-logWhiskerDelta95);
			confRange95.set(whiskerAbove95, lower95);
			confRange95.set(whiskerBelow95, lower95);
			confRange95.set(x, lower95);
			confRange95.set(x, upper95);
			confRange95.set(whiskerAbove95, upper95);
			confRange95.set(whiskerBelow95, upper95);
			
			whiskerAbove68 = Math.pow(10, Math.log10(x)+logWhiskerDelta68);
			whiskerBelow68 = Math.pow(10, Math.log10(x)-logWhiskerDelta68);
			confRange68.set(whiskerAbove68, lower68);
			confRange68.set(whiskerBelow68, lower68);
			confRange68.set(x, lower68);
			confRange68.set(x, upper68);
			confRange68.set(whiskerAbove68, upper68);
			confRange68.set(whiskerBelow68, upper68);
		}
		
		funcs.add(confRange95);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.LIGHT_GRAY));
		
//		funcs.add(confRange68);
//		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
	}
	
	private static Map<Integer, List<SectMappedUncertainDataConstraint>> sectMapped(FaultSystemRupSet rupSet,
			boolean parents, List<? extends SectMappedUncertainDataConstraint> datas) {
		Map<Integer, List<SectMappedUncertainDataConstraint>> ret = new HashMap<>();
		
		for (SectMappedUncertainDataConstraint data : datas) {
			int id;
			if (parents)
				id = rupSet.getFaultSectionData(data.sectionIndex).getParentSectionId();
			else
				id = data.sectionIndex;
			List<SectMappedUncertainDataConstraint> parentDatas = ret.get(id);
			if (parentDatas == null) {
				parentDatas = new ArrayList<>();
				ret.put(id, parentDatas);
			}
			parentDatas.add(data);
		}
		
		return ret;
	}

}
