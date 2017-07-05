package scratch.kevin.hazus;

import java.io.File;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.sha.calc.hazardMap.HazardDataSetLoader;
import org.opensha.sha.calc.hazus.parallel.HazusJobWriter;

import com.google.common.base.Preconditions;

public class HazusStitch {
	
	private static int copyExcept(LocationList badLocs, File from, File to) throws IOException {
		File[] files = from.listFiles();
		
		int skipped = 0;
		
		for (File file : files) {
			if (HazardDataSetLoader.shouldSkip(file))
				continue;
			
			// recursively parse dirs
			if (file.isDirectory()) {
				File newTo = new File(to, file.getName());
				if (!newTo.exists())
					newTo.mkdir();
				skipped += copyExcept(badLocs, file, newTo);
				continue;
			}
			
			Location loc = HazardDataSetLoader.decodeFileName(file.getName());
			Preconditions.checkNotNull(loc, "error decoding filename");
			
			if (badLocs != null && badLocs.contains(loc))
				skipped++;
			else {
				FileUtils.copyFileToDirectory(file, to);
			}
		}
		
		return skipped;
	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File newLocsFile = new File("/home/kevin/OpenSHA/hazus/grid_fix/new.csv");
		File oldLocsFile = new File("/home/kevin/OpenSHA/hazus/grid_fix/old.csv");
		
		File origCurves = new File("/home/kevin/OpenSHA/hazus/grid_fix/indep/curves_orig");
		File newCurves = new File("/home/kevin/OpenSHA/hazus/grid_fix/indep/curves_new");
		File mergeCurves = new File("/home/kevin/OpenSHA/hazus/grid_fix/indep/curves_merged");
		
		LocationList newLocs = HazusJobWriter.loadCSV(newLocsFile);
		LocationList oldLocs = HazusJobWriter.loadCSV(oldLocsFile);
		
		// first remove all of the bad locs
		if (mergeCurves.exists()) {
			System.out.println("NOT copying originals, since it already exists!");
		} else {
			mergeCurves.mkdir();
			int skipped = copyExcept(oldLocs, origCurves, mergeCurves);
			Preconditions.checkState(skipped > 0, "didn't skip anything!");
			int expectedSkip = oldLocs.size()*4;
			Preconditions.checkState(skipped == expectedSkip, "should have skipped "
					+expectedSkip+" ("+oldLocs.size()+" * 4), but actual: "+skipped);
		}
		
		// now copy new locs over
		copyExcept(null, newCurves, mergeCurves);
	}

}
