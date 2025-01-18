package scratch.kevin.prvi25.figures;

import static scratch.kevin.prvi25.figures.PRVI_Paths.FIGURES_DIR;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.PoliticalBoundariesData;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;

import com.google.common.base.Preconditions;

import scratch.kevin.latex.LaTeXUtils;

public class COV_StatsTexWriter {

	public static void main(String[] args) throws IOException {
		File outputDir = new File(FIGURES_DIR, "logic_tree_hazard");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		FileWriter fw = new FileWriter(new File(outputDir, "cov_stats.tex"));
		
		Region mapReg = PRVI25_RegionLoader.loadPRVI_MapExtents();
		XY_DataSet[] polBounds = PoliticalBoundariesData.loadDefaultOutlines(mapReg);
		Region[] polRegions = new Region[polBounds.length];
		
		for (int i=0; i<polBounds.length; i++) {
			LocationList list = new LocationList();
			for (Point2D pt : polBounds[i])
				list.add(new Location(pt.getY(), pt.getX()));
			polRegions[i] = new Region(list, BorderType.MERCATOR_LINEAR);
		}
		
		for (boolean gmm : new boolean[] {false,true}) {
			if (gmm)
				System.out.println("GMM");
			else
				System.out.println("ERF");
			String prefix = "Hazard";
			CSVFile<String> csv;
			double gridSpacing;
			if (gmm) {
				csv = CSVFile.readFile(new File(outputDir, "gmm_1.0s_TWO_IN_50.csv"), true);
				gridSpacing = 0.025;
				prefix = "GMM"+prefix;
			} else {
				csv = CSVFile.readFile(new File(outputDir, "1.0s_TWO_IN_50.csv"), true);
				gridSpacing = 0.1;
			}
			
			MinMaxAveTracker overallTrack = new MinMaxAveTracker();
			MinMaxAveTracker insideTrack = new MinMaxAveTracker();
			double halfSpacing = 0.5*gridSpacing;
			for (int row=1; row<csv.getNumRows(); row++) {
				double lat = csv.getDouble(row, 1);
				double lon = csv.getDouble(row, 2);
				double cov = csv.getDouble(row, 7);
				overallTrack.addValue(cov);
				
				Location center = new Location(lat, lon);
				if (mapReg.contains(center)) {
					Location[] testLocs = {
							center,
							new Location(lat-halfSpacing, lon-halfSpacing),
							new Location(lat-halfSpacing, lon+halfSpacing),
							new Location(lat+halfSpacing, lon-halfSpacing),
							new Location(lat+halfSpacing, lon+halfSpacing),
					};
					boolean inside = false;
					for (Location loc : testLocs) {
						for (Region reg : polRegions) {
							if (reg.contains(loc)) {
								inside = true;
								break;
							}
						}
						if (inside)
							break;
					}
					if (inside)
						insideTrack.addValue(cov);
				}
			}
			System.out.println("Overall stats: "+overallTrack);
			System.out.println("Inside stats: "+insideTrack);

			fw.write(LaTeXUtils.defineValueCommand(prefix+"Min",
					LaTeXUtils.numberExpFormatFixedDecimal(overallTrack.getMin(), 2))+"\n");
			fw.write(LaTeXUtils.defineValueCommand(prefix+"COVMax",
					LaTeXUtils.numberExpFormatFixedDecimal(overallTrack.getMax(), 2))+"\n");
			fw.write(LaTeXUtils.defineValueCommand(prefix+"COVAvg",
					LaTeXUtils.numberExpFormatFixedDecimal(overallTrack.getAverage(), 2))+"\n");
			fw.write(LaTeXUtils.defineValueCommand(prefix+"LandCOVMin",
					LaTeXUtils.numberExpFormatFixedDecimal(insideTrack.getMin(), 2))+"\n");
			fw.write(LaTeXUtils.defineValueCommand(prefix+"LandCOVMax",
					LaTeXUtils.numberExpFormatFixedDecimal(insideTrack.getMax(), 2))+"\n");
			fw.write(LaTeXUtils.defineValueCommand(prefix+"LandCOVAvg",
					LaTeXUtils.numberExpFormatFixedDecimal(insideTrack.getAverage(), 2))+"\n");
		}
		
		fw.close();
	}

}
