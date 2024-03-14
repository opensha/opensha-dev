package scratch.kevin.nshm23.timeDependence;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.GregorianCalendar;
import java.util.List;

import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.faultSurface.FaultSection;

import scratch.UCERF3.utils.LastEventData;

public class U3_DOLE_MapGen {

	public static void main(String[] args) throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/FM3_1_branch_averaged.zip"));
		LastEventData.populateSubSects(sol.getRupSet().getFaultSectionDataList(), LastEventData.load());
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(sol.getRupSet().getFaultSectionDataList());
		mapMaker.setWriteGeoJSON(true);
		
		List<? extends FaultSection> subSects = sol.getRupSet().getFaultSectionDataList();
		double[] sectYears = new double[subSects.size()];
		for (int s=0; s<sectYears.length; s++) {
			long dole = subSects.get(s).getDateOfLastEvent();
			if (dole == Long.MIN_VALUE) {
				sectYears[s] = Double.NaN;
				continue;
			}
			GregorianCalendar cal = new GregorianCalendar();
			cal.setTimeInMillis(dole);
			int year = cal.get(GregorianCalendar.YEAR);
			sectYears[s] = year;
		}
		
		double minYear = 0;
		double maxYear = 2024;
		CPT yearCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(minYear, maxYear);
		yearCPT.setPreferredTickInterval(200);
		mapMaker.plotSectScalars(sectYears, yearCPT, "Year of Last Event");
		mapMaker.setSkipNaNs(true);
		mapMaker.setSectNaNChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(0, 0, 0, 80)));
		mapMaker.setSectOutlineChar(null);
		mapMaker.plot(new File("/tmp"), "u3_dole_year_map", " ");
	}

}
