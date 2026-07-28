package scratch.kevin.nshm27.figures;

import static scratch.kevin.nshm27.figures.NSHM27_PaperPaths.*;

import java.io.File;
import java.io.IOException;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.modules.ModuleContainer;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.nshmp.erf.nshm27.NSHM27_InvConfigFactory;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceObsSeisDMAdjustment;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_LogicTree;
import gov.usgs.earthquake.nshmp.erf.nshm27.util.NSHM27_RegionLoader.NSHM27_SeismicityRegions;

public class InterfaceObsSeisAdjFigures {

	public static void main(String[] args) throws IOException {
		ModuleContainer.VERBOSE_DEFAULT = false;
		File outputDir = new File(FIGURES_DIR, "interface_obs_seis_adj");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		NSHM27_InvConfigFactory factory = new NSHM27_InvConfigFactory();
		
		CPT redCPT = GMT_CPT_Files.SEQUENTIAL_LAJOLLA_UNIFORM.instance().reverse().rescale(0d, 1d);
		
		for (NSHM27_SeismicityRegions seisReg : NSHM27_SeismicityRegions.values()) {
			LogicTreeBranch<LogicTreeNode> branch = NSHM27_LogicTree.buildDefault(seisReg, TectonicRegionType.SUBDUCTION_INTERFACE, false);
			System.out.println("Default branch: "+branch);
			
			CSVFile<String> csv = new CSVFile<>(true);
			csv.addLine("Adjustment", "Overall Moment Reduction Factor", "Section-Averaged Reduction Factor",
					"Maximum Reduction Factor", "Median Reduction Factor", "Average Slip Rate (mm/yr)", "Maximum Slip Rate (mm/yr)");
			
			branch.setValue(NSHM27_InterfaceObsSeisDMAdjustment.NONE);
			
			FaultSystemRupSet rs = factory.buildRuptureSet(branch, FaultSysTools.defaultNumThreads());
			SectSlipRates fullSlips = rs.requireModule(SectSlipRates.class);
			
			GeographicMapMaker mapMaker = new GeographicMapMaker(rs.getFaultSectionDataList());
			
			double fullMoment = fullSlips.calcTotalMomentRate();
			
			for (NSHM27_InterfaceObsSeisDMAdjustment adj : NSHM27_InterfaceObsSeisDMAdjustment.values()) {
				SectSlipRates slips;
				if (adj == NSHM27_InterfaceObsSeisDMAdjustment.NONE) {
					slips = fullSlips;
				} else {
					branch.setValue(adj);
					rs = factory.buildRuptureSet(branch, FaultSysTools.defaultNumThreads());
					slips = rs.requireModule(SectSlipRates.class);
				}
				MinMaxAveTracker slipTrack = new MinMaxAveTracker();
				
				double branchMoment = slips.calcTotalMomentRate();
				
				double[] factors = new double[slips.size()];
				for (int s=0; s<slips.size(); s++) {
					double fullSlip = fullSlips.getSlipRate(s);
					double redSlip = slips.getSlipRate(s);
					factors[s] = (fullSlip - redSlip)/fullSlip;
					slipTrack.addValue(slips.getSlipRate(s)*1e3);
				}
				
				mapMaker.plotSectScalars(factors, redCPT, "Slip-deficit rate reduction factor");
				
				mapMaker.plot(outputDir, seisReg.name()+"_"+adj.name(), adj.getShortName());
				
				double overall = (fullMoment - branchMoment)/fullMoment;
				csv.addLine(adj.getShortName(),
						twoDF.format(overall),
						twoDF.format(StatUtils.mean(factors)),
						twoDF.format(StatUtils.max(factors)),
						twoDF.format(DataUtils.median(factors)),
						twoDF.format(slipTrack.getAverage()),
						twoDF.format(slipTrack.getMax()));
			}
			
			csv.writeToFile(new File(outputDir, seisReg.name()+"_factors.csv"));
		}
	}

}
