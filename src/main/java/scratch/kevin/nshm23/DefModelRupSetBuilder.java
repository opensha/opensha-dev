package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen.PlotLevel;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_ConstraintBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;

import com.google.common.base.Preconditions;

public class DefModelRupSetBuilder {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/home/kevin/OpenSHA/nshm23/def_models");
		
		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
		
		LogicTreeBranch<LogicTreeNode> defaultBranch = NSHM23_LogicTreeBranch.DEFAULT;
		defaultBranch = defaultBranch.copy();
		defaultBranch.setValue(NSHM23_ScalingRelationships.AVERAGE);
		
		NSHM23_DeformationModels defaultDM = defaultBranch.getValue(NSHM23_DeformationModels.class);
		
		PlotLevel plotLevelDefaultChoice = PlotLevel.FULL;
		PlotLevel plotLevelOtherChoices = PlotLevel.DEFAULT;
		
		NSHM23_FaultModels fm = defaultBranch.requireValue(NSHM23_FaultModels.class);
		
		outputDir = new File(outputDir, fm.getFilePrefix());
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		FaultSystemRupSet rawRupSet = null;
		
		int threads = FaultSysTools.defaultNumThreads();
		
		for (NSHM23_DeformationModels dm : NSHM23_DeformationModels.values()) {
			if (dm.isApplicableTo(fm) && (dm.getNodeWeight(null) > 0d || dm == defaultDM)) {
				LogicTreeBranch<LogicTreeNode> branch = defaultBranch.copy();
				branch.setValue(dm);
				
				File rupSetFile = new File(outputDir, dm.getFilePrefix()+".zip");
				if (rawRupSet == null && rupSetFile.exists())
					rawRupSet = FaultSystemRupSet.load(rupSetFile);
				
				FaultSystemRupSet rupSet;
				if (rawRupSet == null) {
					rupSet = factory.buildRuptureSet(branch, threads);
					rawRupSet = rupSet;
				} else {
					rupSet = factory.updateRuptureSetForBranch(rawRupSet, branch);
				}
				
				rupSet.write(rupSetFile);
				
				File reportDir = new File(outputDir, dm.getFilePrefix());
				
				PlotLevel level = dm == defaultDM ? plotLevelDefaultChoice : plotLevelOtherChoices;
				
				ReportPageGen report = new ReportPageGen(rupSet, null, dm.getName(), reportDir,
						ReportPageGen.getDefaultRupSetPlots(level));
				report.setReplot(true);
				
				report.generatePage();
			}
		}
		
		// now parkfield
		NSHM23_ConstraintBuilder constrBuilder = new NSHM23_ConstraintBuilder(rawRupSet, 0.5d);
		List<Integer> parkfieldRups = constrBuilder.findParkfieldRups();
		Preconditions.checkState(!parkfieldRups.isEmpty());
		CSVFile<String> parkfieldCSV = new CSVFile<>(true);
		List<String> header = new ArrayList<>();
		header.add("Rupture");
		List<MinMaxAveTracker> scaleTracks = new ArrayList<>();
		List<FaultSystemRupSet> scaleRupSets = new ArrayList<>();
		for (NSHM23_ScalingRelationships scale : NSHM23_ScalingRelationships.values()) {
			header.add(scale.getShortName());
			scaleTracks.add(new MinMaxAveTracker());
			scaleRupSets.add(FaultSystemRupSet.buildFromExisting(rawRupSet).forScalingRelationship(scale).build());
		}
		parkfieldCSV.addLine(header);
		for (int rupIndex : parkfieldRups) {
			List<String> line = new ArrayList<>();
			line.add(rupIndex+"");
			for (int s=0; s<scaleRupSets.size(); s++) {
				double mag = scaleRupSets.get(s).getMagForRup(rupIndex);
				line.add((float)mag+"");
				scaleTracks.get(s).addValue(mag);
			}
			parkfieldCSV.addLine(line);
		}
		List<String> line = new ArrayList<>();
		line.add("AVERAGE");
		for (MinMaxAveTracker track : scaleTracks)
			line.add((float)track.getAverage()+"");
		parkfieldCSV.addLine(line);
		line = new ArrayList<>();
		line.add("MINIMUM");
		for (MinMaxAveTracker track : scaleTracks)
			line.add((float)track.getMin()+"");
		parkfieldCSV.addLine(line);
		line = new ArrayList<>();
		line.add("MAXIMUM");
		for (MinMaxAveTracker track : scaleTracks)
			line.add((float)track.getMax()+"");
		parkfieldCSV.addLine(line);
		
		parkfieldCSV.writeToFile(new File(outputDir, "parkfield_rups.csv"));
	}

}
