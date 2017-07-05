package scratch.kevin;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.griddedSeismicity.GridSourceProvider;
import scratch.UCERF3.utils.FaultSystemIO;

public class IrisFocalParse {

	public static void main(String[] args) throws IOException, DocumentException {
		CSVFile<String> csv = CSVFile.readFile(new File(
				"/home/kevin/Documents/Geol 551/eq presentation/focals.iris.csv"), true);
		FileWriter fw = new FileWriter(new File("/home/kevin/Documents/Geol 551/eq presentation/focals.iris.sf"));
		for (int row=1; row<csv.getNumRows(); row++) {
			List<String> line = csv.getLine(row);
			String lon = line.get(5);
			String lat = line.get(4);
			String dep = line.get(6);
			String mag = line.get(11);
			if (line.get(19).trim().isEmpty() || mag.trim().isEmpty() || mag.toLowerCase().contains("nan"))
				continue;
			int strike = (int)Math.round(Double.parseDouble(line.get(19)));
			int dip = (int)Math.round(Double.parseDouble(line.get(20)));
			int rake = (int)Math.round(Double.parseDouble(line.get(21)));
			
			
			
			fw.write(lon+" "+lat+" "+dep+" "+strike+" "+dip+" "+rake+" "+mag+"\n");
		}
		fw.close();
		
		// now find UCERF3 M5 rate
		FaultSystemSolution fss = FaultSystemIO.loadSol(new File("/home/kevin/workspace/OpenSHA/dev/"
				+ "scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		GridSourceProvider prov = fss.getGridSourceProvider();
		GriddedRegion reg = prov.getGriddedRegion();
		Location eqLoc = new Location(40.969, -124.747);
		Location closestNode = null;
		double minDist = Double.POSITIVE_INFINITY;
		for (Location loc : reg.getNodeList()) {
			double dist = LocationUtils.horzDistanceFast(eqLoc, loc);
			if (dist < minDist) {
				minDist = dist;
				closestNode = loc;
			}
		}
		System.out.println("Closest node is "+minDist+" km away");
		IncrementalMagFreqDist mfd = prov.getNodeUnassociatedMFD(reg.indexForLocation(closestNode));
//		IncrementalMagFreqDist mfd = prov.getNodeUnassociatedMFD(reg.indexForLocation(eqLoc));
		System.out.println("MFD at eq loc:");
		for (Point2D pt : mfd)
			System.out.println("Mag: "+(float)pt.getX()+"\tRate: "+(float)pt.getY()+"\tRI: "+(float)(1d/pt.getY()));
		System.out.println("Cumulative MFD at eq loc:");
		for (Point2D pt : mfd.getCumRateDistWithOffset())
			System.out.println("Mag: "+(float)pt.getX()+"\tRate: "+(float)pt.getY()+"\tRI: "+(float)(1d/pt.getY()));
		
		FaultSectionPrefData salmonData = FaultModels.FM3_1.fetchFaultSectionsMap().get(16);
		System.out.println("Salmon Distance: "+salmonData.getStirlingGriddedSurface(1d).getDistanceRup(eqLoc)+" km");
	}

}
