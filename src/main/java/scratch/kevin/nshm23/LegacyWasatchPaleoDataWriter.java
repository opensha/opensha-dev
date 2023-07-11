package scratch.kevin.nshm23;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.GeoJSONFaultReader;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.data.SegRateConstraint;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;

import scratch.alessandro.WasatchInversion;
import scratch.alessandro.logicTreeEnums.WasatchSlipRatesEnum;
import ucar.ma2.MAMath.MinMax;

public class LegacyWasatchPaleoDataWriter {

	public static void main(String[] args) throws IOException {
		WasatchInversion wasatch = new WasatchInversion(WasatchSlipRatesEnum.UNIFORM_MEAN);
		
		ArrayList<FaultSectionPrefData> sects = wasatch.getFaultSectionDataList();
		ArrayList<SegRateConstraint> paleo = wasatch.getSectionRateConstraints();
		
		CSVFile<String> csv = new CSVFile<>(true);
		
		/*
		 * * Paleo Site Name
* Latitude
* Longitude
* Paleo Rate (or recurrence interval)
* Lower 95% confidence bound
* Upper 95% confidence bound
* Standard Deviation (optional, otherwise I will assume it is (upper -
lower)/4)
		 */
		
		csv.addLine("Name", "Latitude", "Longitude", "Paleoseismic Rate", "Rate Std Dev", "Lower 95%", "Upper 95%");
		
		for (SegRateConstraint segRate : paleo) {
			int segIndex = segRate.getSegIndex();
			FaultSectionPrefData sect = sects.get(segIndex);
			
			System.out.println("Paleo for fault "+segRate.getFaultName());
			System.out.println("\tSection: "+sect.getSectionName());
			Location firstLoc = sect.getFaultTrace().first();
			Location lastLoc = sect.getFaultTrace().last();
			double lat = 0.5*(firstLoc.getLatitude()+lastLoc.getLatitude());
			double lon = 0.5*(firstLoc.getLongitude()+lastLoc.getLongitude());
			Location loc = new Location(lat, lon);
			System.out.println("\tLocation: "+loc);
			
			System.out.println("\tMean rate: "+(float)segRate.getMean()+" +/- "+(float)segRate.getStdDevOfMean());
			System.out.println("\t95% bounds: "+(float)segRate.getLower95Conf()+", "+(float)segRate.getUpper95Conf());
			System.out.println("\tCalc std. dev: "+(float)0.25d*(segRate.getUpper95Conf()-segRate.getLower95Conf()));
			
			List<String> line = new ArrayList<>();
			line.add(segRate.getFaultName());
			line.add(lat+"");
			line.add(lon+"");
			line.add(segRate.getMean()+"");
			line.add(segRate.getStdDevOfMean()+"");
			line.add(segRate.getLower95Conf()+"");
			line.add(segRate.getUpper95Conf()+"");
			csv.addLine(line);
		}
		
		csv.writeToFile(new File("/tmp/wasatch_legacy_data.csv"));
		
		// write wasatch sections
		List<GeoJSONFaultSection> geoSects = new ArrayList<>();
		for (FaultSectionPrefData sect : sects) {
			GeoJSONFaultSection geoSect = new GeoJSONFaultSection(sect);
			if (!Float.isFinite(geoSect.getDipDirection()))
				geoSect.setDipDirection((float)(geoSect.getFaultTrace().getAveStrike()+90d));
			geoSects.add(geoSect);
		}
		
		GeographicMapMaker mapMaker = new RupSetMapMaker(geoSects, RupSetMapMaker.buildBufferedRegion(geoSects));
		
		HashSet<FaultSection> highlights = new HashSet<>();
		for (SegRateConstraint segRate : paleo)
			highlights.add(geoSects.get(segRate.getSegIndex()));
		
		mapMaker.setSectHighlights(highlights, new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		mapMaker.setWriteGeoJSON(true);
		
		mapMaker.plot(new File("/tmp"), "wasatch_sects.geojson", "Legacy Wasatch Sections");
		
//		GeoJSONFaultReader.writeFaultSections(new File("/tmp/wasatch_sects.geojson"), geoSects);
	}

}
