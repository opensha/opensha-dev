package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;

import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.util.ClassUtils;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.griddedSeismicity.GridSourceProvider;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.logicTree.LogicTreeBranchNode;
import scratch.UCERF3.utils.FaultSystemIO;

public class GridSourcesNodeWriter {

	public static void main(String[] args) throws IOException, DocumentException {
		FaultSystemSolution sol = FaultSystemIO.loadSol(new File(
				"/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		
		GridSourceProvider grid = sol.getGridSourceProvider();
		
		GriddedRegion reg = grid.getGriddedRegion();
		
		File outputFile = new File("/tmp/relm_gridded_region.csv");
		
		CSVFile<String> csv = new CSVFile<String>(true);
		
		csv.addLine("Index", "Latitude", "Longitude");
		
		LocationList nodeList = reg.getNodeList();
		
		for (int i=0; i<reg.getNodeCount(); i++) {
			Location loc = nodeList.get(i);
			csv.addLine(i+"", loc.getLatitude()+"",loc.getLongitude()+"");
		}
		
		csv.writeToFile(outputFile);
		
		for (Class<? extends LogicTreeBranchNode<?>> clazz : LogicTreeBranch.getLogicTreeNodeClasses()) {
			System.out.println("== "+ClassUtils.getClassNameWithoutPackage(clazz)+" ==");
			System.out.println("||= Prefix =||= Description =||");
			for (LogicTreeBranchNode<?> node : clazz.getEnumConstants()) {
				if (node.getRelativeWeight(LogicTreeBranch.DEFAULT.getValue(InversionModels.class)) <= 0)
					continue;
				String str = "||";
				if (LogicTreeBranch.DEFAULT.getValueUnchecked(clazz) == node)
					str += "'''"+node.encodeChoiceString()+"'''";
				else
					str += node.encodeChoiceString();
				str += "||"+node.getName()+"||";
				System.out.println(str);
			}
			System.out.println();
		}
	}

}
