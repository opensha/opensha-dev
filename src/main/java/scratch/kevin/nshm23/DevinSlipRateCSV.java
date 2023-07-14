package scratch.kevin.nshm23;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionSlipRates;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;

public class DevinSlipRateCSV {

	public static void main(String[] args) throws IOException {
		File solDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
//				+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/");
				+ "2023_06_29-nshm23_branches-NSHM23_v2-CoulombRupSet-DsrTap-TotNuclRate-NoRed-ThreshAvgIterRelGR/");
//		File solFile = new File(solDir, "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip");
		File solFile = new File(solDir, "results_NSHM23_v2_CoulombRupSet_branch_averaged.zip");
		FaultSystemSolution sol = FaultSystemSolution.load(solFile);
		
		SolutionSlipRates solSlips = sol.requireModule(SolutionSlipRates.class);
		FaultSystemRupSet rupSet = sol.getRupSet();
		SectSlipRates slipRates = rupSet.requireModule(SectSlipRates.class);
		
		List<List<? extends FaultSection>> dmSubSectsList = new ArrayList<>();
		List<String> header = new ArrayList<>();
		header.add("Subection Index");
		header.add("Subection Name");
		header.add("Start Latitude");
		header.add("Start Longitude");
		header.add("End Latitude");
		header.add("End Longitude");
		NSHM23_FaultModels fm = NSHM23_FaultModels.NSHM23_v2;
		for (NSHM23_DeformationModels dm : NSHM23_DeformationModels.values()) {
			if (dm.getNodeWeight(null) > 0) {
				dmSubSectsList.add(dm.build(fm));
				header.add(dm.getShortName()+" Slip Rate (mm/yr)");
			}
		}
		header.add("Branch-Averagea Target Slip Rate (mm/yr)");
		header.add("Branch-Averagea Solution Slip Rate (mm/yr)");
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine(header);
		
		for (int s=0; s<rupSet.getNumSections(); s++) {
			FaultSection sect = rupSet.getFaultSectionData(s);
			List<String> line = new ArrayList<>(header.size());
			line.add(s+"");
			line.add(sect.getSectionName());
			FaultTrace trace = sect.getFaultTrace();
			Location first = trace.first();
			Location last = trace.last();
			line.add((float)first.lat+"");
			line.add((float)first.lon+"");
			line.add((float)last.lat+"");
			line.add((float)last.lon+"");
			for (List<? extends FaultSection> dmSubSects : dmSubSectsList) {
				FaultSection dmSect = dmSubSects.get(s);
				line.add((float)dmSect.getOrigAveSlipRate()+"");
			}
			double target = slipRates.getSlipRate(s)*1e3;
			line.add((float)target+"");
			double solution = solSlips.get(s)*1e3;
			line.add((float)solution+"");
			csv.addLine(line);
		}
		
		csv.writeToFile(new File("/tmp/slips_for_devin.csv"));
		
		// now write out segmentation-model-dependent slips
		File nodesDir = new File(solDir, "node_branch_averaged");
		
		List<SolutionSlipRates> segModelSlips = new ArrayList<>();
		List<NSHM23_SegmentationModels> segModels = new ArrayList<>();
		
		for (NSHM23_SegmentationModels segModel : NSHM23_SegmentationModels.values()) {
			if (segModel.getNodeWeight(null) == 0d)
				continue;
			File segModelSolFile = new File(nodesDir, "SegModel_"+segModel.getFilePrefix()+".zip");
			FaultSystemSolution segModelSol = FaultSystemSolution.load(segModelSolFile);
			
			segModels.add(segModel);
			segModelSlips.add(segModelSol.requireModule(SolutionSlipRates.class));
		}
		
		header = new ArrayList<>();
		header.add("Subection Index");
		header.add("Subection Name");
		header.add("Start Latitude");
		header.add("Start Longitude");
		header.add("End Latitude");
		header.add("End Longitude");
		header.add("Branch-Averagea Target Slip Rate (mm/yr)");
		header.add("Branch-Averagea Solution Slip Rate (mm/yr)");
		for (NSHM23_SegmentationModels segModel : segModels)
			header.add(segModel.getShortName()+" Solution Slip Rate (mm/yr)");
		csv = new CSVFile<>(true);
		csv.addLine(header);
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		CPT segCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, segModels.size()-1d);
		
		for (int s=0; s<rupSet.getNumSections(); s++) {
			FaultSection sect = rupSet.getFaultSectionData(s);
			List<String> line = new ArrayList<>(header.size());
			line.add(s+"");
			line.add(sect.getSectionName());
			FaultTrace trace = sect.getFaultTrace();
			Location first = trace.first();
			Location last = trace.last();
			line.add((float)first.lat+"");
			line.add((float)first.lon+"");
			line.add((float)last.lat+"");
			line.add((float)last.lon+"");
			double target = slipRates.getSlipRate(s)*1e3;
			line.add((float)target+"");
			double solution = solSlips.get(s)*1e3;
			line.add((float)solution+"");
			for (SolutionSlipRates segSlips : segModelSlips)
				line.add((float)(segSlips.get(s)*1e3)+"");
			csv.addLine(line);
			
			if (sect.getName().contains("Cucamonga") || sect.getName().contains("Madre")) {
				boolean label = funcs.isEmpty();
				
				DefaultXY_DataSet targetXY = new DefaultXY_DataSet();
				targetXY.set(first.lon, target);
				targetXY.set(last.lon, target);
				DefaultXY_DataSet solXY = new DefaultXY_DataSet();
				solXY.set(first.lon, solution);
				solXY.set(last.lon, solution);
				
				if (label) {
					targetXY.setName("Target");
					solXY.setName("Solution");
				}
				
				funcs.add(targetXY);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.CYAN));
				
				funcs.add(solXY);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
				
				for (int i=0; i<segModels.size(); i++) {
					DefaultXY_DataSet segXY = new DefaultXY_DataSet();
					double val = segModelSlips.get(i).get(s)*1e3;
					segXY.set(first.lon, val);
					segXY.set(last.lon, val);
					if (label)
						segXY.setName(segModels.get(i).getShortName()+" Solution");
					
					funcs.add(segXY);
					chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, segCPT.getColor((float)i)));
				}
			}
		}
		
		csv.writeToFile(new File("/tmp/seg_slips_for_devin.csv"));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Sierra Madre/Cucamonga Slips", "Longitude", "Slip Rate (mm/yr)");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec);
		
		PlotUtils.writePlots(new File("/tmp"), "seg_slips_for_devin", gp, 1500, 800, true, false, false);
	}

}
