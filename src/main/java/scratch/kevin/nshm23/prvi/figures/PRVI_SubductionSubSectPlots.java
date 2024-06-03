package scratch.kevin.nshm23.prvi.figures;

import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.PRVI25_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionScalingRelationships;
import org.opensha.sha.faultSurface.FaultSection;

public class PRVI_SubductionSubSectPlots {

	public static void main(String[] args) throws IOException {
		LogicTreeBranch<LogicTreeNode> branch = PRVI25_LogicTreeBranch.DEFAULT_SUBDUCTION_INTERFACE;
		PRVI25_SubductionScalingRelationships scale = PRVI25_SubductionScalingRelationships.LOGA_C4p0;
		branch = branch.copy();
		branch.setValue(scale);
		PRVI25_InvConfigFactory factory = new PRVI25_InvConfigFactory();
		
		File outputDir = new File("/tmp");
		
		Font font = new Font(Font.SANS_SERIF, Font.BOLD, 24);
		DecimalFormat magDF = new DecimalFormat("0.00");
		
		CPT magCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(7.5d, 9d);
		CPT slipCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, 5d);
		CPT slipUncertCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, 2d);
		CPT rakeCPT = GMT_CPT_Files.SEQUENTIAL_NAVIA_UNIFORM.instance().rescale(0d, 90d);
		for (PRVI25_SubductionFaultModels fm : PRVI25_SubductionFaultModels.values()) {
			if (fm.getNodeWeight(branch) == 0d)
				continue;
			branch = branch.copy();
			branch.setValue(fm);
			branch.requireValue(PRVI25_SubductionDeformationModels.class).build(fm);
			FaultSystemRupSet rupSet = factory.buildRuptureSet(branch, FaultSysTools.defaultNumThreads());
			
			List<? extends FaultSection> sects = rupSet.getFaultSectionDataList();
			GeographicMapMaker mapMaker = new GeographicMapMaker(sects);
			mapMaker.setFillSurfaces(true);
			
			List<Double> minMags = new ArrayList<>();
			List<Double> maxMags = new ArrayList<>();
			for (int s=0; s<sects.size(); s++) {
				minMags.add(rupSet.getMinMagForSection(s));
				maxMags.add(rupSet.getMaxMagForSection(s));
			}
			
			double minMin = minMags.stream().mapToDouble(D->D).min().getAsDouble();
			double maxMin = minMags.stream().mapToDouble(D->D).max().getAsDouble();
			double minMax = maxMags.stream().mapToDouble(D->D).min().getAsDouble();
			double maxMax = maxMags.stream().mapToDouble(D->D).max().getAsDouble();
			
			Range xRange = mapMaker.getXRange();
			Range yRange = mapMaker.getYRange();
			
			double annX = xRange.getLowerBound() + 0.97*xRange.getLength();
			double annY = yRange.getLowerBound() + 0.95*yRange.getLength();

			mapMaker.plotSectScalars(minMags, magCPT, "Minimum Magnitude ("+scale.getShortName()+")");
			XYTextAnnotation rangeAnn = new XYTextAnnotation("["+magDF.format(minMin)+", "+magDF.format(maxMin)+"]", annX, annY);
			rangeAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
			rangeAnn.setFont(font);
			mapMaker.addAnnotation(rangeAnn);
			mapMaker.plot(outputDir, "subduction_min_mag_"+fm.getFilePrefix(), " ");
			
			mapMaker.plotSectScalars(maxMags, magCPT, "Maximum Magnitude ("+scale.getShortName()+")");
			mapMaker.clearAnnotations();
			rangeAnn = new XYTextAnnotation("["+magDF.format(minMax)+", "+magDF.format(maxMax)+"]", annX, annY);
			rangeAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
			rangeAnn.setFont(font);
			mapMaker.addAnnotation(rangeAnn);
			mapMaker.plot(outputDir, "subduction_max_mag_"+fm.getFilePrefix(), " ");
			
			mapMaker.clearAnnotations();
			
			for (PRVI25_SubductionDeformationModels dm : PRVI25_SubductionDeformationModels.values()) {
				if (dm.getNodeWeight(branch) == 0d)
					continue;
				sects = dm.build(fm);
				
				mapMaker.setFaultSections(sects);
				List<Double> slips = new ArrayList<>(sects.size());
				List<Double> slipUncerts = new ArrayList<>(sects.size());
				List<Double> rakes = new ArrayList<>(sects.size());
				for (FaultSection sect : sects) {
					slips.add(sect.getOrigAveSlipRate());
					slipUncerts.add(sect.getOrigSlipRateStdDev());
					rakes.add(sect.getAveRake());
				}
				
				mapMaker.plotSectScalars(slips, slipCPT, dm.getShortName()+" Slip Rate (mm/yr)");
				mapMaker.plot(outputDir, "subduction_slip_"+fm.getFilePrefix()+"_"+dm.getFilePrefix(), " ");
				
				mapMaker.plotSectScalars(slipUncerts, slipUncertCPT, dm.getShortName()+" Slip Rate Uncertainty (mm/yr)");
				mapMaker.plot(outputDir, "subduction_slip_uncert_"+fm.getFilePrefix()+"_"+dm.getFilePrefix(), " ");
				
				mapMaker.plotSectScalars(rakes, rakeCPT, dm.getShortName()+" Rake");
				mapMaker.plot(outputDir, "subduction_rake_"+fm.getFilePrefix()+"_"+dm.getFilePrefix(), " ");
				
			}
		}
	}

}
