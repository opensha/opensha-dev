package scratch.kevin.pointSources;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.util.DataUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.FocalMech;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

public class TestFullGridSourceFileWriter {

	public static void main(String[] args) throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged_gridded.zip"));
		GridSourceProvider gridProv = sol.getGridSourceProvider();
		FaultGridAssociations associations = sol.getRupSet().requireModule(FaultGridAssociations.class);
		
		File outputDir = new File("/tmp/grid_source_tests");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		double minMag = 0d;
		outputDir = new File(outputDir, "updated");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
//		double minMag = 5d;
//		outputDir = new File(outputDir, "updated_m5");
//		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		CSVFile<String> gridCSV = new CSVFile<>(true);
		gridCSV.addLine("Grid Index", "Latitude", "Longitude");
		
		CSVFile<String> rupCSV = new CSVFile<>(false);
		rupCSV.addLine(
				"Grid Index",
				"Magnitude",
				"Annual Rate",
				"Rake",
				"Dip",
				"Dip Uncertainty",
				"Strike",
				"Strike Uncertainty",
				"Upper Depth (km)",
				"Lower Depth (km)",
				"Length (km)",
				"Tectonic Regime",
				"Associated Section Index 1",
				"Fraction Associated 1",
				"Associated Section Index N",
				"Fraction Associated N");
		
		WC1994_MagLengthRelationship WC94 = new WC1994_MagLengthRelationship();
		
		int locRoundScale = 3;
		int magRoundScale = 3;
		int mechRoundSigFigs = 3;
		int depthRoundSigFigs = 3;
		int lenRoundSigFigs = 3;
		int rateRoundSigFigs = 6;
		
		GriddedRegion gridReg = gridProv.getGriddedRegion();
		for (int i=0; i<gridReg.getNodeCount(); i++) {
			Location loc = gridReg.getLocation(i);
			
			gridCSV.addLine(i+"", getFixedPrecision(loc.lat, locRoundScale), getFixedPrecision(loc.lon, locRoundScale));
			
			double fractSS = gridProv.getFracStrikeSlip(i);
			double fractN = gridProv.getFracNormal(i);
			double fractR = gridProv.getFracReverse(i);
			
			IncrementalMagFreqDist mfd = gridProv.getMFD(i);
			if (mfd == null)
				continue;
			IncrementalMagFreqDist mfdAssoc = gridProv.getMFD_SubSeisOnFault(i);
			Map<Integer, Double> nodeFractAssociations = null;
			if (mfdAssoc != null) {
				nodeFractAssociations = new HashMap<>(associations.getScaledSectFracsOnNode(i));
				Preconditions.checkState(!nodeFractAssociations.isEmpty());
				// turn it into a fractional: scale to 1 if not already
				double sumFracts = 0d;
				for (double fract : nodeFractAssociations.values())
					sumFracts += fract;
				if ((float)sumFracts != 1f) {
					for (int sectIndex : new ArrayList<>(nodeFractAssociations.keySet()))
						nodeFractAssociations.put(sectIndex, nodeFractAssociations.get(sectIndex)/sumFracts);
				}
			}
			for (int m=0; m<mfd.size(); m++) {
				double mag = mfd.getX(m);
				if (mag < minMag)
					continue;
				double totRate = mfd.getY(m);
				if (totRate == 0d)
					continue;
				double associatedFract = 0d;
				if (mfdAssoc != null && mfdAssoc.size() > m) {
					Preconditions.checkState((float)mfdAssoc.getX(m) == (float)mag);
					double assocRate = mfdAssoc.getY(mag);
					associatedFract = assocRate/totRate;
					Preconditions.checkState((float)associatedFract <= 1f, "Bad associatedFract = %s / %s = %s",
							assocRate, totRate, associatedFract);
				}
				for (FocalMech mech : FocalMech.values()) {
					double mechRate;
					switch (mech) {
					case STRIKE_SLIP:
						mechRate = totRate*fractSS;
						break;
					case NORMAL:
						mechRate = totRate*fractN;
						break;
					case REVERSE:
						mechRate = totRate*fractR;
						break;

					default:
						throw new IllegalStateException();
					}
					if (mechRate == 0d)
						continue;
					
					double dipRad = Math.toRadians(mech.dip());
					
					double depth = (float)mag < 6.5f ? 5d : 1d;
					double length = WC94.getMedianLength(mag);
					double aspectWidth = length / 1.5;
					double ddWidth = (14.0 - depth) / Math.sin(dipRad);
					ddWidth = Math.min(aspectWidth, ddWidth);
					double lower = depth + ddWidth * Math.sin(dipRad);
					
					List<String> line = new ArrayList<>();
					line.add(i+"");
					line.add(getFixedPrecision(mag, magRoundScale));
					line.add(getSigFigs(mechRate, rateRoundSigFigs));
					line.add(getSigFigs(mech.rake(), mechRoundSigFigs));
					line.add(getSigFigs(mech.dip(), mechRoundSigFigs));
					line.add(getSigFigs(Double.NaN, mechRoundSigFigs));
					line.add(getSigFigs(Double.NaN, mechRoundSigFigs));
					line.add(getSigFigs(Double.NaN, mechRoundSigFigs));
					line.add(getSigFigs(depth, depthRoundSigFigs));
					line.add(getSigFigs(lower, depthRoundSigFigs));
					line.add(getSigFigs(length, lenRoundSigFigs));
					line.add(TectonicRegionType.ACTIVE_SHALLOW.name());
					if (associatedFract > 0) {
						for (int sectIndex : nodeFractAssociations.keySet()) {
							double fract = associatedFract*nodeFractAssociations.get(sectIndex);
							line.add(sectIndex+"");
							line.add(getSigFigs(fract, rateRoundSigFigs)+"");
						}
					}
					rupCSV.addLine(line);
				}
			}
		}
		
		gridCSV.writeToFile(new File(outputDir, "grid_source_locations.csv"));
		rupCSV.writeToFile(new File(outputDir, "grid_sources.csv"));
		Feature.write(gridReg.toFeature(), new File(outputDir, "grid_region.geojson"));
	}
	
	private static String getFixedPrecision(double val, int scale) {
		if (Double.isNaN(val))
			return "";
		if (!Double.isFinite(val))
			return val+"";
		if (val == Math.floor(val))
			return (int)val+"";
		return DataUtils.roundFixed(val, scale)+"";
	}
	
	private static String getSigFigs(double val, int sigFigs) {
		if (Double.isNaN(val))
			return "";
		if (!Double.isFinite(val))
			return val+"";
		if (val == Math.floor(val))
			return (int)val+"";
		return DataUtils.roundSigFigs(val, sigFigs)+"";
	}

}
