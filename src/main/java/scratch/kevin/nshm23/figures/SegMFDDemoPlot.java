package scratch.kevin.nshm23.figures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.RupSetScalingRelationship;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets.FullySegmentedRupSetConfig;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SegmentationMFD_Adjustment;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs.Builder;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.SectNucleationMFD_Estimator;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import net.mahdilamb.colormap.Colors;

public class SegMFDDemoPlot {

	public static void main(String[] args) throws IOException {
		FaultTrace trace1 = new FaultTrace();
		trace1.add(new Location(0d, 0d));
		trace1.add(LocationUtils.location(trace1.first(), 0d, 60d));
		GeoJSONFaultSection sect1 = new GeoJSONFaultSection.Builder(0, "Fault 1", trace1)
				.lowerDepth(15d)
				.upperDepth(0d)
				.rake(0d)
				.dip(90)
				.slipRate(10d)
				.slipRateStdDev(1)
				.aseismicity(0d)
				.build();
		
		FaultTrace trace2 = new FaultTrace();
		trace2.add(LocationUtils.location(trace1.last(), 0d, 5d));
		trace2.add(LocationUtils.location(trace2.first(), 0d, 60d));
		GeoJSONFaultSection sect2 = new GeoJSONFaultSection.Builder(1, "Fault 2", trace2)
				.lowerDepth(15d)
				.upperDepth(0d)
				.rake(0d)
				.dip(90)
				.slipRate(10d)
				.slipRateStdDev(1)
				.aseismicity(0d)
				.build();
		
		double subSectLen = sect1.getOrigDownDipWidth()*0.25;
		int minSubSectsPerRup = 4;
		
		List<FaultSection> subSects = new ArrayList<>();
		subSects.addAll(sect1.getSubSectionsList(subSectLen, 0, 2));
		subSects.addAll(sect2.getSubSectionsList(subSectLen, subSects.size(), 2));
		
		RupSetScalingRelationship scale = NSHM23_ScalingRelationships.LOGA_C4p2;
		
		RuptureSets.SimpleAzimuthalRupSetConfig config = new RuptureSets.SimpleAzimuthalRupSetConfig(subSects, scale);
		config.setMaxJumpDist(10d);
		config.setMinSectsPerParent(minSubSectsPerRup);
		FaultSystemRupSet rupSet = config.build(1);
		
		double bVal = 1d;
		
		Builder mfdBuilder = new SupraSeisBValInversionTargetMFDs.Builder(rupSet, bVal);
		mfdBuilder.magDepDefaultRelStdDev(M->0.1*Math.max(1, Math.pow(10, bVal*0.5*(M-6))));
		SupraSeisBValInversionTargetMFDs targetMFDs = mfdBuilder.build();
		
		IncrementalMagFreqDist unsegTarget = targetMFDs.getTotalOnFaultSupraSeisMFD();
		
		JumpProbabilityCalc segModel = NSHM23_SegmentationModels.MID.getModel(rupSet, rupSet.getModule(LogicTreeBranch.class));
		
		SectNucleationMFD_Estimator adjustment = SegmentationMFD_Adjustment.REL_GR_THRESHOLD_AVG.getAdjustment(segModel);
		mfdBuilder.adjustTargetsForData(adjustment);
		
		targetMFDs = mfdBuilder.build();
		
		IncrementalMagFreqDist segTarget = targetMFDs.getTotalOnFaultSupraSeisMFD();
		
		segModel = NSHM23_SegmentationModels.CLASSIC.getModel(rupSet, rupSet.getModule(LogicTreeBranch.class));
		
		adjustment = SegmentationMFD_Adjustment.REL_GR_THRESHOLD_AVG.getAdjustment(segModel);
		mfdBuilder.clearTargetAdjustments().adjustTargetsForData(adjustment);
		
		targetMFDs = mfdBuilder.build();
		
		IncrementalMagFreqDist classicTarget = targetMFDs.getTotalOnFaultSupraSeisMFD();
		
		for (boolean incr : new boolean[] {true,false}) {
			List<DiscretizedFunc> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			DiscretizedFunc segMFD, classicMFD, unsegMFD;
			if (incr) {
				segMFD = segTarget;
				classicMFD = classicTarget;
				unsegMFD = unsegTarget;
			} else {
				segMFD = segTarget.getCumRateDistWithOffset();
				unsegMFD = unsegTarget.getCumRateDistWithOffset();
				classicMFD = classicTarget.getCumRateDistWithOffset();
			}
			unsegMFD.setName("Unsegmented");
			funcs.add(unsegMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_green));
			
			segMFD.setName("Partially-segmented");
			funcs.add(segMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_blue));
			
			classicMFD.setName("Fully-segmented");
			funcs.add(classicMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_orange));
			
			PlotSpec plot = new PlotSpec(funcs, chars, " ", "Magnitude", incr ? "Incremental rate (1/yr)" : "Cumulative rate (1/yr)");
			plot.setLegendInset(true);
			
			HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
			
			gp.drawGraphPanel(plot, false, true, new Range(6.5d, 7.5d), new Range(1e-4, 1e-1));
			
			String prefix = incr ? "mfds_incr" : "mfds_cml";
			PlotUtils.writePlots(new File("/tmp"), prefix, gp, 800, 800, true, true, false);
		}
		
		
	}

}
