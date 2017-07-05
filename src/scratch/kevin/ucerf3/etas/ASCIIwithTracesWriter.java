package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.geo.Location;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.RuptureSurface;

import com.google.common.base.Preconditions;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.utils.FaultSystemIO;

public class ASCIIwithTracesWriter {

	public static void main(String[] args) throws IOException, DocumentException {
		File dir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2016_02_17-spontaneous-1000yr-scaleMFD1p14-full_td-subSeisSupraNucl-gridSeisCorr");
		File inputFile = new File(dir, "results_m4.bin");
		
		File outputDir = new File(dir, "results_ascii_m4_with_traces");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		FaultSystemRupSet rupSet = FaultSystemIO.loadRupSet(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		
		int index = 0;
		for (List<ETAS_EqkRupture> catalog : ETAS_CatalogIO.getBinaryCatalogsIterable(inputFile, 4d)) {
			File file = new File(outputDir, "catalog_"+(index++)+".txt");
			System.out.println("Writing "+file.getAbsolutePath());
			FileWriter fw = new FileWriter(file);
			
			String header = ETAS_CatalogIO.EVENT_FILE_HEADER;
			header += "\t[NumTraceLocs\tLat1\tLon1\t..\tLatN\tLonN]";
			
			fw.write("% "+header+"\n");
			
			for (ETAS_EqkRupture rup : catalog) {
				String line = ETAS_CatalogIO.getEventFileLine(rup);
				if (rup.getFSSIndex() >= 0) {
					// add fault trace
					RuptureSurface surf = rupSet.getSurfaceForRupupture(rup.getFSSIndex(), 1d, false);
					
					FaultTrace upper = surf.getUpperEdge();
					
					line += "\t"+upper.size();
					for (Location loc : upper)
						line += "\t"+(float)loc.getLatitude()+"\t"+(float)loc.getLongitude();
				}
				
				fw.write(line+"\n");
			}
			
			fw.close();
		}
	}

}
