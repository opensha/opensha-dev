package scratch.kevin.ucerf3;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.util.FileUtils;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.observedEarthquake.parsers.UCERF3_CatalogParser;
import org.opensha.sha.earthquake.observedEarthquake.parsers.ngaWest.NGAWestParser;
import org.opensha.sha.faultSurface.CompoundSurface;
import org.opensha.sha.faultSurface.EvenlyGriddedSurface;
import org.opensha.sha.faultSurface.RuptureSurface;

import com.google.common.base.Preconditions;

public class CatalogWriter {
	
	protected static void writeFiniteSurfaceFile(ObsEqkRupture rup, File file, double distCutoff) throws IOException {
		writeFiniteSurfaceFile(rup.getRuptureSurface(), file, distCutoff);
	}
	
	protected static void writeFiniteSurfaceFile(RuptureSurface rupSurf, File file, double distCutoff) throws IOException {
		ArrayList<RuptureSurface> surfs = new ArrayList<RuptureSurface>();
		
		if (rupSurf instanceof CompoundSurface) 
			surfs.addAll(((CompoundSurface)rupSurf).getSurfaceList());
		else
			surfs.add(rupSurf);
		
		FileWriter fw = new FileWriter(file);
		
		fw.write("# FORMAT: <lon>\t<lat>\t<depth>\n");
		
		ArrayList<Location> prevLocs = new ArrayList<Location>();
		for (int i=0; i<surfs.size(); i++) {
			RuptureSurface surf = surfs.get(i);
			
			int skips = 0;
			ArrayList<Location> locs = new ArrayList<Location>();
			
			for (Location loc : surf.getEvenlyDiscritizedListOfLocsOnSurface()) {
				if (distCutoff > 0 && distCutoff < Double.MAX_VALUE) {
					// make sure we don't have any duplicates with the previous surface
					boolean skip = false;
					for (Location prevLoc : prevLocs) {
						if (LocationUtils.linearDistanceFast(loc, prevLoc) < distCutoff) {
							skip = true;
							break;
						}
					}
					if (skip) {
						skips++;
						break;
					}
				}
				locs.add(loc);
			}
			
			fw.write("# SURFACE "+i+" ("+skips+" pts skipped)\n");
			for (Location loc : locs) {
				fw.write(loc.getLongitude()+"\t"+loc.getLatitude()+"\t"+loc.getDepth()+"\n");
			}
			
			if (skips > 0)
				System.out.println(skips+" skipped!");
			
			prevLocs.addAll(locs);
		}
		
		fw.close();
	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File polTllDir = new File("src"+File.separator+"resources"+File.separator+"data"+File.separator+"ngaWest");
//		File excelFile = new File(polTllDir, "EQ.V8.xls");
		File excelFile = new File(polTllDir, "EQ.V8_modUCERF3.xls");
		File coulombDir = new File("/home/kevin/OpenSHA/UCERF3/rollins_stein_finite");
		
		File outputDir = new File("/home/kevin/OpenSHA/UCERF3/finiteCatalog");
			outputDir.mkdir();
		
		double gridSpacing = 1.0;
		double distCutoff = gridSpacing * 0.5;
		
		ObsEqkRupList ngaRups = NGAWestParser.loadNGAWestFiles(excelFile, polTllDir, gridSpacing);
		ObsEqkRupList finiteRups = new ObsEqkRupList();
		ArrayList<File> finiteFiles = new ArrayList<File>();
		
		for (ObsEqkRupture rup : ngaRups) {
			if (rup.getRuptureSurface() != null) {
				finiteRups.add(rup);
				File file = new File(outputDir, "ngaRupSurface_"+rup.getEventId()+".txt");
				writeFiniteSurfaceFile(rup, file, distCutoff);
				Preconditions.checkState(!finiteFiles.contains(file));
				finiteFiles.add(file);
			}
		}
		
		for (ObsEqkRupture rup : CoulombFileParser.loadObsEqkRups(coulombDir)) {
			finiteRups.add(rup);
			File file = new File(outputDir, rup.getEventId());
			writeFiniteSurfaceFile(rup, file, distCutoff);
			Preconditions.checkState(!finiteFiles.contains(file));
			finiteFiles.add(file);
		}
		
		File catalogFile = new File(outputDir.getParent(), "UCERF3CatalogVersion1.txt");
		
		ObsEqkRupList catalog = UCERF3_CatalogParser.loadCatalog(catalogFile);
		
		ArrayList<String> catalogFileLines = FileUtils.loadFile(catalogFile.getAbsolutePath());
		Preconditions.checkState(catalog.size() == catalogFileLines.size());
		
		for (int i=0; i<finiteRups.size(); i++) {
			ObsEqkRupture ngaRup = finiteRups.get(i);
			int ind = NGAWestParser.getBestMatchIndex(ngaRup, catalog);
			if (ind < 0) {
				System.out.println("No match for rupture: "+ngaRup.getEventId());
				continue;
			}
			String line = catalogFileLines.get(ind).trim();
			if (line.contains("txt")) {
				System.out.println("Conflict...rupture "+ngaRup.getEventId()+" also best match for:\n"+line);
				continue;
			}
			catalogFileLines.set(ind, line+" "+finiteFiles.get(i).getName());
		}
		
		File outputFile = new File(outputDir, "catalog.txt");
		
		FileWriter fw = new FileWriter(outputFile);
		
		fw.write("# This catalog is a modified version of Karen Feltzer's UCERF3 Catalog Version 1\n" +
			 	 "# (retrieved via e-mail from Ned Field, 12/13/2011) to include references to finite\n" +
				 "# fault ruptures provided by the NGA west project (distributed via e-mail by Paul\n" +
				 "# Spudich on 6/23/11).\n");
		
		for (String line : catalogFileLines)
			fw.write(line + "\n");
		fw.close();
	}

}
