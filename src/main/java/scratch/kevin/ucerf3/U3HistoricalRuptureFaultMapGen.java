package scratch.kevin.ucerf3;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.List;

import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.mapping.PoliticalBoundariesData;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.faultSurface.FaultSection;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.LastEventData;

public class U3HistoricalRuptureFaultMapGen {

	public static void main(String[] args) throws IOException {
		FaultModels fm = FaultModels.FM3_1;
		List<? extends FaultSection> subSects = fm.getDefaultDeformationModel().build(fm);
		LastEventData.populateSubSects(subSects, LastEventData.load());
		
		Region relm = new CaliforniaRegions.RELM_TESTING();
		Region ca = new Region(new Location(relm.getMinLat(), relm.getMinLon()),
				new Location(relm.getMaxLat(), relm.getMaxLon()));
		XY_DataSet[] caOutlines = PoliticalBoundariesData.loadCAOutlines();
		GeographicMapMaker mapMaker = new GeographicMapMaker(ca, caOutlines);
		
		mapMaker.setFaultSections(subSects);
		CPT yearCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().reverse().rescale(1800, 2023);
		yearCPT.setNanColor(Color.LIGHT_GRAY);
		
		double[] years = new double[subSects.size()];
		for (int s=0; s<years.length; s++) {
			FaultSection sect = subSects.get(s);
			long lastEvent = sect.getDateOfLastEvent();
			if (lastEvent > Long.MIN_VALUE) {
				Date date = new Date(lastEvent);
				GregorianCalendar cal = new GregorianCalendar();
				cal.setTime(date);
				years[s] = cal.get(GregorianCalendar.YEAR);
			} else {
				years[s] = Double.NaN;
			}
		}
		
		mapMaker.plotSectScalars(years, yearCPT, "Year of Last Event");
//		mapMaker.setSkipNaNs(true);
		mapMaker.setSectOutlineChar(null);
		
		mapMaker.plot(new File("/tmp"), "faults_date_last", " ");
	}

}
