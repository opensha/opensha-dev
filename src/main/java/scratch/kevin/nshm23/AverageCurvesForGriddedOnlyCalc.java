package scratch.kevin.nshm23;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;

import com.google.common.base.Preconditions;

public class AverageCurvesForGriddedOnlyCalc {

	public static void main(String[] args) throws IOException {
		File origResults = new File(args[0]);
		File gridResults = new File(args[1]);
		String dirName = args[2];
		
		File origResultsZip = new File(origResults.getParentFile(), origResults.getName()+".zip");
//		File gridResultsZip = new File(gridResults.getParentFile(), gridResults.getName()+".zip");
		
		SolutionLogicTree origSLT = SolutionLogicTree.load(origResultsZip);
//		SolutionLogicTree gridSLT = SolutionLogicTree.load(gridResultsZip);
		
		// TODO: this does not handle cases with multiple branch averages
		
		LogicTreeBranch<?> commonBranch = origSLT.getLogicTree().getBranch(0).copy();
		for (LogicTreeBranch<?> branch : origSLT.getLogicTree()) {
			for (int i=0; i<commonBranch.size(); i++) {
				LogicTreeNode commonVal = commonBranch.getValue(i);
				if (commonVal != null) {
					LogicTreeNode branchVal = branch.getValue(i);
					if (!commonVal.equals(branchVal))
						commonBranch.clearValue(i);
				}
			}
		}
		
		List<LogicTreeLevel<? extends LogicTreeNode>> faultLevels = new ArrayList<>();
		List<LogicTreeNode> faultNodes = new ArrayList<>();
		for (int i=0; i<commonBranch.size(); i++) {
			LogicTreeLevel<?> level = commonBranch.getLevel(i);
			LogicTreeNode value = commonBranch.getValue(i);
			if (value != null && level.affects(FaultSystemSolution.RATES_FILE_NAME, true)) {
				faultLevels.add(level);
				faultNodes.add(value);
			}
		}
		LogicTreeBranch<LogicTreeNode> subBranch = new LogicTreeBranch<>(faultLevels, faultNodes);
		
		System.out.println("Common sub-branch: "+subBranch);
		
		File subBranchDir = new File(gridResults, subBranch.buildFileName());
		System.out.println("Destination sub-branch dir: "+subBranchDir.getAbsolutePath());
		Preconditions.checkState(subBranchDir.exists() || subBranchDir.mkdir());
		
		File subBranchHazardDir = new File(subBranchDir, dirName);
		System.out.println("Destination sub-branch hazard dir: "+subBranchHazardDir.getAbsolutePath());
		Preconditions.checkState(subBranchHazardDir.exists() || subBranchHazardDir.mkdir());
		
		Map<String, DiscretizedFunc[]> averages = null;
		double sumWeight = 0d;
		
		LocationList locs = null;
		
		for (LogicTreeBranch<?> branch : origSLT.getLogicTree()) {
			System.out.println("Processing original branch: "+branch);
			File solDir = new File(origResults, branch.buildFileName());
			Preconditions.checkState(solDir.exists(), "dir doesn't exist: %s", solDir.getAbsolutePath());
			
			File hazardDir = new File(solDir, dirName);
			Preconditions.checkState(hazardDir.exists(), "dir doesn't exist: %s", hazardDir.getAbsolutePath());
			
			if (averages == null) {
				// first time
				averages = new HashMap<>();
				for (File file : hazardDir.listFiles()) {
					if (file.getName().endsWith(".csv") || file.getName().endsWith(".csv.gz")) {
						averages.put(file.getName(), new DiscretizedFunc[0]);
					}
				}
			}
			
			double weight = origSLT.getLogicTree().getBranchWeight(branch);
			sumWeight += weight;
			
			for (String fileName : averages.keySet()) {
				File hazardFile = new File(hazardDir, fileName);
				Preconditions.checkState(hazardFile.exists(), "file doesn't exist: %s", hazardFile.getAbsolutePath());
				
				CSVFile<String> csv = CSVFile.readFile(hazardFile, true);
				DiscretizedFunc[] curves = SolHazardMapCalc.loadCurvesCSV(csv, null);
				
				if (locs == null) {
					locs = new LocationList();
					for (int row=1; row<csv.getNumRows(); row++) {
						double lat = csv.getDouble(row, 1);
						double lon = csv.getDouble(row, 2);
						locs.add(new Location(lat, lon));
					}
				}
				
				DiscretizedFunc[] avgCurves = averages.get(fileName);
				if (avgCurves.length == 0) {
					// first time
					avgCurves = new DiscretizedFunc[curves.length];
					averages.put(fileName, avgCurves);
				} else {
					Preconditions.checkState(avgCurves.length == curves.length);
				}
				
				for (int i=0; i<curves.length; i++) {
					if (avgCurves[i] == null) {
						avgCurves[i] = new ArbitrarilyDiscretizedFunc();
						for (Point2D pt : curves[i])
							avgCurves[i].set(pt.getX(), 0d);
					}
					Preconditions.checkState(curves[i].size() == avgCurves[i].size());
					for (int j=0; j<curves[i].size(); j++) {
						Preconditions.checkState((float)curves[i].getX(j) == (float)avgCurves[i].getX(j));
						avgCurves[i].set(j, avgCurves[i].getY(j) + weight*curves[i].getY(j));
					}
				}
			}
		}
		
		System.out.println("Writing averages with sumWeight="+sumWeight);
		for (String fileName : averages.keySet()) {
			File hazardFile = new File(subBranchHazardDir, fileName);
			System.out.println("Writing "+hazardFile.getAbsolutePath());
			
			DiscretizedFunc[] avgCurves = averages.get(fileName);
			for (int i=0; i<avgCurves.length; i++)
				avgCurves[i].scale(1d/sumWeight);
			
			SolHazardMapCalc.writeCurvesCSV(hazardFile, avgCurves, locs);
		}
		
		System.out.println("DONE");
	}

}
