package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.zip.ZipException;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.griddedSeismicity.FaultPolyMgr;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.logicTree.LogicTreeBranch;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class PolygonFileWriter {

	public static void main(String[] args) throws ZipException, IOException {
		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip"));
		
		File outputDir = new File("/tmp/polygons");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		LogicTreeBranch branch = LogicTreeBranch.DEFAULT;
		
		FaultModels[] fms = { FaultModels.FM3_1, FaultModels.FM3_2 };
		
		for (FaultModels fm : fms) {
			LogicTreeBranch fmBranch = (LogicTreeBranch)branch.clone();
			fmBranch.setValue(fm);
			
			InversionFaultSystemSolution sol = cfss.getSolution(fmBranch);
			
			// parent sections
			writePolygons(new File(outputDir, fm.name()+"_parents_geologic.csv"), fm.fetchFaultSections(), null);
			writePolygons(new File(outputDir, fm.name()+"_sub_sect_hazard.csv"), sol.getRupSet().getFaultSectionDataList(),
					sol.getRupSet().getInversionTargetMFDs().getGridSeisUtils().getPolyMgr());
		}
	}
	
	private static void writePolygons(File outputFile, List<FaultSectionPrefData> sects, FaultPolyMgr polys)
			throws IOException {
		CSVFile<String> csv = new CSVFile<String>(true);
		
		csv.addLine("Section ID", "Name", "Parent ID", "Parent Name", "Latitude", "Longitude");
		
		for (FaultSectionPrefData sect : sects) {
			List<String> line = Lists.newArrayList(sect.getSectionId()+"", sect.getName());
			if (sect.getParentSectionId() >= 0) {
				line.add(sect.getParentSectionId()+"");
				line.add(sect.getParentSectionName());
			} else {
				line.add("");
				line.add("");
			}
			Region poly;
			if (polys != null)
				poly = polys.getPoly(sect.getSectionId());
			else
				poly = sect.getZonePolygon();
			if (poly == null) {
				System.out.println("NULL polygon for "+sect.getSectionId()+": "+sect.getSectionName());
				continue;
			}
			for (Location loc : poly.getBorder()) {
				List<String> subLine = Lists.newArrayList(line);
				subLine.add(loc.getLatitude()+"");
				subLine.add(loc.getLongitude()+"");
				csv.addLine(subLine);
			}
		}
		
		csv.writeToFile(outputFile);
	}

}
