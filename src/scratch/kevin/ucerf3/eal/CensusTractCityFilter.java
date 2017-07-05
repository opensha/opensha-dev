package scratch.kevin.ucerf3.eal;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.HashSet;

import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.sra.gui.portfolioeal.Asset;
import org.opensha.sra.gui.portfolioeal.Portfolio;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

public class CensusTractCityFilter {

	public static void main(String[] args) throws IOException {
		File ealDir = new File("/home/kevin/OpenSHA/UCERF3/eal");
		File cityDir = new File(ealDir, "city_specific");
		// go here: https://www.census.gov/geo/maps-data/maps/2010ref/st06_tract.html
		// find the city of interest, then download the "District to Map Sheet File"
//		File tractsForCityFile = new File(cityDir, "san_bernardino_DC10CT_C06071_CT2MS.txt");
		// or null to use region
		File tractsForCityFile = null;
		// portfolio file to filter
		File portInFile = new File(ealDir, "Porter-22-May-14-CA-ppty-90pct-Wills.txt");
//		File portInFile = new File(ealDir, "Porter-02-Jun-16-CA-ppty-90pct-Wald.txt");
//		File portOutFile = new File(cityDir, "san_bernardino_"+portInFile.getName());
//		File portOutFile = new File(cityDir, "coachella_valley_"+portInFile.getName());
		File portOutFile = new File("/tmp/coachella_valley_"+portInFile.getName());
		
		HashSet<String> tractNames = new HashSet<String>();
		
		if (tractsForCityFile != null) {
			for (String line : Files.readLines(tractsForCityFile, Charset.defaultCharset())) {
				line = line.trim();
				if (line.startsWith("TRACT")) {
					String[] split = line.split(";");
					String name = split[1];
					Preconditions.checkState(!tractNames.contains(name));
					tractNames.add(name);
				}
			}
		} else {
			// or do by region
//			LocationList border = new LocationList();
////			border.add(new Location(33.867970, -116.655201));
////			border.add(new Location(34.037557, -116.612410));
////			border.add(new Location(33.737092, -115.970720));
////			border.add(new Location(33.536610, -116.252075));
//			border.add(new Location(34.02706,-116.63163));
//			border.add(new Location(33.92798,-116.70888));
//			border.add(new Location(33.84569,-116.6094));
//			border.add(new Location(33.72848,-116.53602));
//			border.add(new Location(33.58517,-116.25467));
//			border.add(new Location(33.66443,-116.03133));
//			border.add(new Location(33.85317,-116.1884));
//			border.add(new Location(34.0028,-116.48967));
//			Region region = new Region(border, BorderType.MERCATOR_LINEAR);
			
			Region region = new Region(new Location(33.739320, -116.412738), 20d);
			
			Portfolio portfolio = Portfolio.createPortfolio(portInFile);
			for (Asset asset : portfolio.getAssetList()) {
				if (region.contains(asset.getLocation()))
					tractNames.add(asset.getParameterList().getParameter(String.class, "AssetName").getValue());
			}
		}
		
		System.out.println("Loaded "+tractNames.size()+" tracts");
		
		FileWriter outFW = new FileWriter(portOutFile);
		int matching = 0;
		for (String line : Files.readLines(portInFile, Charset.defaultCharset())) {
			line = line.trim();
			if (line.startsWith("AssetGroup")) {
				outFW.write(line+"\n");
				String[] split = line.split(",");
				Preconditions.checkState(split[2].equals("AssetName"));
			} else {
				String[] split = line.split(",");
				String name = split[2];
				if (tractNames.contains(name)) {
					matching++;
					outFW.write(line+"\n");
				}
			}
		}
		
		System.out.println("Wrote "+matching+" assets");
		
		outFW.close();
	}

}
