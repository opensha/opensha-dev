package scratch.kevin.prvi25.figures;

import static scratch.kevin.prvi25.figures.PRVI_Paths.FIGURES_DIR;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

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
		List<File> inputDirs = new ArrayList<>();
		List<String> names = new ArrayList<>();
		List<String> texPrefixes = new ArrayList<>();
		List<Double> gridSpacings = new ArrayList<>();
		
		File outputDir = new File(FIGURES_DIR, "logic_tree_hazard");
		
		inputDirs.add(outputDir);
		names.add("ERF");
		texPrefixes.add("");
		gridSpacings.add(0.1);
		
		inputDirs.add(new File(outputDir, "gmm_csvs"));
		names.add("GMM");
		texPrefixes.add("GMM");
		gridSpacings.add(0.025);
		
		inputDirs.add(new File(outputDir, "erf_plus_gmm_csvs"));
		names.add("ERF+GMM");
		texPrefixes.add("Comb");
		gridSpacings.add(0.025);
		
//		inputDirs.add(new File("/tmp/cov_test/new_model_test_full"));
//		names.add("New test");
//		texPrefixes.add("Test");
//		gridSpacings.add(0.1);
		
//		inputDirs.add(new File("/tmp/cov_test/mue_car_corr"));
//		names.add("New test correlated Mue&Car");
//		texPrefixes.add("Test");
//		gridSpacings.add(0.1);
		
		for (File dir : inputDirs)
			Preconditions.checkState(dir.exists());
		
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
		
		double[] periods = { 0d, 1d };
		
		for (double period : periods) {
			String csvName, texPrefix;
			if (period == 0d) {
				System.out.println("PGA");
				csvName = "pga_TWO_IN_50.csv";
				texPrefix = "HazardPGA";
			} else if (period == 1d) {
				System.out.println("1s SA");
				csvName = "1.0s_TWO_IN_50.csv";
				texPrefix = "HazardOneSec";
			} else {
				throw new IllegalStateException();
			}
			for (int i=0; i<inputDirs.size(); i++) {
				File inputDir = inputDirs.get(i);
				String myTexPrefix = texPrefixes.get(i)+texPrefix;
				System.out.println(names.get(i));
				CSVFile<String> csv = CSVFile.readFile(new File(inputDir, csvName), true);
				double gridSpacing = gridSpacings.get(i);
				
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
				System.out.println("Overall:\t"+(float)overallTrack.getAverage()+" ["+
						(float)overallTrack.getMin()+", "+(float)overallTrack.getMax()+"]");
				System.out.println("Inside stats:\t"+(float)insideTrack.getAverage()+" ["+
						(float)insideTrack.getMin()+", "+(float)insideTrack.getMax()+"]");

				fw.write(LaTeXUtils.defineValueCommand(myTexPrefix+"Min",
						LaTeXUtils.numberExpFormatFixedDecimal(overallTrack.getMin(), 2))+"\n");
				fw.write(LaTeXUtils.defineValueCommand(myTexPrefix+"COVMax",
						LaTeXUtils.numberExpFormatFixedDecimal(overallTrack.getMax(), 2))+"\n");
				fw.write(LaTeXUtils.defineValueCommand(myTexPrefix+"COVAvg",
						LaTeXUtils.numberExpFormatFixedDecimal(overallTrack.getAverage(), 2))+"\n");
				fw.write(LaTeXUtils.defineValueCommand(myTexPrefix+"LandCOVMin",
						LaTeXUtils.numberExpFormatFixedDecimal(insideTrack.getMin(), 2))+"\n");
				fw.write(LaTeXUtils.defineValueCommand(myTexPrefix+"LandCOVMax",
						LaTeXUtils.numberExpFormatFixedDecimal(insideTrack.getMax(), 2))+"\n");
				fw.write(LaTeXUtils.defineValueCommand(myTexPrefix+"LandCOVAvg",
						LaTeXUtils.numberExpFormatFixedDecimal(insideTrack.getAverage(), 2))+"\n");
			}
		}
		
		fw.close();
	}

}
