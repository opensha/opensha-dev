package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;

public class SantaClaraAftershockCatalogSearch {

	public static void main(String[] args) throws IOException {
		// centroids, approx
		Location cuppertinoLoc = new Location(37.320086, -122.048582);
		Location paloAltoLoc = new Location(37.429103, -122.128197);
		
		double minMag = 6;
		
		double bufferKM = 5;
		Region cuppertinoReg = new Region(cuppertinoLoc, bufferKM);
		Region paloAltoReg = new Region(paloAltoLoc, bufferKM);
		
//		Region santaClara = loadSantaClara(new File("/tmp/tl_2010_06085_nowater_nad83.csv"));
		Region santaClara = coarseSantaClara();
		
		File catFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2016_06_15-haywired_m7-10yr-full_td-no_ert-combined/results_descendents_m5.bin");
		if (args.length > 0)
			catFile = new File(args[0]);
		Iterable<List<ETAS_EqkRupture>> catalogs = ETAS_CatalogIO.getBinaryCatalogsIterable(catFile, 0d);
		
		List<List<ETAS_EqkRupture>> cuppertinoCatalogs = Lists.newArrayList();
		List<List<ETAS_EqkRupture>> paloAltoCatalogs = Lists.newArrayList();
		List<List<ETAS_EqkRupture>> bothCatalogs = Lists.newArrayList();
		List<List<ETAS_EqkRupture>> scCatalogs = Lists.newArrayList();
		
		int count = 0;
		
		for (List<ETAS_EqkRupture> catalog : catalogs) {
			if (count % 10000 == 0)
				System.out.println("Catalog "+count);
			count++;
			boolean hasCuppertino = false;
			boolean hasPaloAlto = false;
			boolean hasSC = false;
			
			for (ETAS_EqkRupture rup : catalog) {
				if (rup.getMag() < minMag)
					continue;
				Location hypo = rup.getHypocenterLocation();
				if (cuppertinoReg.contains(hypo))
					hasCuppertino = true;
				if (paloAltoReg.contains(hypo))
					hasPaloAlto = true;
				if (santaClara.contains(hypo))
					hasSC = true;
			}
			
			if (hasCuppertino && hasPaloAlto) {
				bothCatalogs.add(catalog);
			} else {
				if (hasCuppertino)
					cuppertinoCatalogs.add(catalog);
				if (hasPaloAlto)
					paloAltoCatalogs.add(catalog);
			}
			if (hasSC)
				scCatalogs.add(catalog);
		}
		
		File baseDir = new File("/tmp/santa_clara_catalogs");
		Preconditions.checkState(baseDir.exists() || baseDir.mkdir());
		writeCatalogs(cuppertinoCatalogs, "cuppertino_only", new File(baseDir, "cuppertino_only"));
		writeCatalogs(paloAltoCatalogs, "palo_alto_only", new File(baseDir, "palo_alto_only"));
		writeCatalogs(bothCatalogs, "both", new File(baseDir, "both"));
		writeCatalogs(scCatalogs, "santa_clara", new File(baseDir, "santa_clara"));
	}
	
	private static void writeCatalogs(List<List<ETAS_EqkRupture>> catalogs, String prefix, File outputDir)
			throws IOException {
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		System.out.println("Writing "+catalogs.size()+" catalogs to "+outputDir.getAbsolutePath());
		int digits = ((catalogs.size()-1)+"").length();
		for (int i=0; i<catalogs.size(); i++) {
			String numStr = i+"";
			while (numStr.length() < digits)
				numStr = "0"+numStr;
			File outFile = new File(outputDir, prefix+"_"+numStr+".txt");
			ETAS_CatalogIO.writeEventDataToFile(outFile, catalogs.get(i));
		}
	}
	
	private static Region loadSantaClara(File file) throws IOException {
		CSVFile<String> csv = CSVFile.readFile(file, true);
		LocationList border = new LocationList();
		for (int row=1; row<csv.getNumRows(); row++) {
			double lat = Double.parseDouble(csv.get(row, 2));
			double lon = Double.parseDouble(csv.get(row, 3));
			border.add(new Location(lat, lon));
		}
//		border.add(border.get(0));
		return new Region(border, BorderType.GREAT_CIRCLE);
	}
	
	private static Region coarseSantaClara() {
		LocationList border = new LocationList();
		border.add(new Location(37.399569, -122.191787));
		border.add(new Location(37.431191, -122.189359));
		border.add(new Location(37.457404, -122.152459));
		border.add(new Location(37.453550, -122.122356));
		border.add(new Location(37.465497, -122.115558));
		border.add(new Location(37.445070, -122.060208));
		border.add(new Location(37.460873, -122.032047));
		border.add(new Location(37.461644, -121.911635));
		border.add(new Location(37.484379, -121.865024));
		border.add(new Location(37.481682, -121.482910));
		border.add(new Location(37.310416, -121.404740));
		border.add(new Location(37.281447, -121.458148));
		border.add(new Location(37.150370, -121.397457));
		border.add(new Location(37.183644, -121.361042));
		border.add(new Location(37.157335, -121.236260));
		border.add(new Location(36.962057, -121.215382));
		border.add(new Location(36.960505, -121.418335));
		border.add(new Location(36.987270, -121.448923));
		border.add(new Location(36.982228, -121.489222));
		border.add(new Location(36.894522, -121.575647));
		border.add(new Location(37.016739, -121.736358));
		border.add(new Location(37.022941, -121.731017));
		border.add(new Location(37.091524, -121.835892));
		border.add(new Location(37.144565, -121.991748));
		border.add(new Location(37.274879, -122.150031));
		border.add(new Location(37.358672, -122.203439));
		border.add(border.get(0));
		
		return new Region(border, BorderType.GREAT_CIRCLE);
	}

}
