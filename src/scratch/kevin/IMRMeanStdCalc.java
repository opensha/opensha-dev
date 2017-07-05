package scratch.kevin;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.siteData.OrderedSiteDataProviderList;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.data.siteData.SiteDataValueList;
import org.opensha.commons.data.siteData.impl.CVM4BasinDepth;
import org.opensha.commons.data.siteData.impl.CVMHBasinDepth;
import org.opensha.commons.data.siteData.impl.USGSBayAreaBasinDepth;
import org.opensha.commons.data.siteData.impl.WaldAllenGlobalVs30;
import org.opensha.commons.data.siteData.impl.WillsMap2006;
import org.opensha.commons.geo.Location;
import org.opensha.commons.param.Parameter;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.CB_2008_AttenRel;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGD_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceRupParameter;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.util.SiteTranslator;

import com.google.common.collect.Lists;

public class IMRMeanStdCalc {

	public static void main(String[] args) throws IOException {
		CSVFile<String> rupCSV = CSVFile.readFile(new File("/tmp/Sourcerepresentation_Mw79.csv"), false);
		double dip = 90;
		double rake = 180;
		double upper_depth = 0d;
		double lower_depth = 10d;
		double mag = 7.9;
		
		ScalarIMR imr = new CB_2008_AttenRel(null);
		imr.setParamDefaults();
		
		File siteFile = new File("/tmp/coords.txt");
		List<Location> locs = Lists.newArrayList();
		for (String line : FileUtils.readLines(siteFile)) {
			line = line.trim();
			if (line.isEmpty())
				continue;
			String[] split = line.split("\t");
			double lat = Double.parseDouble(split[0]);
			double lon = Double.parseDouble(split[1]);
			locs.add(new Location(lat, lon));
		}
		System.out.println("Loaded "+locs.size()+" sites");
		
		ArrayList<SiteData<?>> provList = new ArrayList<SiteData<?>>();
		provList.add(new WillsMap2006());
		// swap these two lines to switch between CVM4 and CVMH
//		provList.add(new CVM4BasinDepth(SiteData.TYPE_DEPTH_TO_2_5));
		provList.add(new CVMHBasinDepth(SiteData.TYPE_DEPTH_TO_2_5));
		ArrayList<Site> sites = new ArrayList<Site>();
		for (Location loc : locs) {
			Site site = new Site(loc);
			
			site.addParameter((Parameter)imr.getParameter(Vs30_Param.NAME).clone());
			site.addParameter((Parameter)imr.getParameter(DepthTo2pt5kmPerSecParam.NAME).clone());
			
			sites.add(site);
		}
		OrderedSiteDataProviderList provs = new OrderedSiteDataProviderList(provList);
		System.out.println("Fetching data!");
		ArrayList<SiteDataValueList<?>> vals  = provs.getAllAvailableData(sites);
		System.out.println("done.");
		
		SiteTranslator trans = new SiteTranslator();
		for (int i=0; i<sites.size(); i++) {
			ArrayList<SiteDataValue<?>> siteVals = new ArrayList<SiteDataValue<?>>();
			for (SiteDataValueList<?> valList : vals) {
				siteVals.add(valList.getValue(i));
			}
			Iterator<Parameter<?>> it = sites.get(i).getParametersIterator();
			while (it.hasNext()) {
				trans.setParameterValue(it.next(), siteVals);
			}
		}
		
		FaultSectionPrefData sect = new FaultSectionPrefData();
		FaultTrace trace = new FaultTrace("");
		for (int row=1; row<rupCSV.getNumRows(); row++) {
			List<String> line = rupCSV.getLine(row);
			double lat = Double.parseDouble(line.get(0));
			double lon = Double.parseDouble(line.get(1));
			trace.add(new Location(lat, lon));
		}
		sect.setFaultTrace(trace);
		sect.setAveDip(dip);
		sect.setAveRake(rake);
		sect.setAveLowerDepth(lower_depth);
		sect.setAveUpperDepth(upper_depth);
		RuptureSurface surf = sect.getStirlingGriddedSurface(1d);
		
		EqkRupture rup = new EqkRupture(mag, rake, surf, null);
		
		CSVFile<String> csv = new CSVFile<String>(true);
		csv.addLine("Lat", "Lon", "Vs30 (m/s)", "Z2.5 (km)", "DistRup (km)", "PGA Mean (g)", "PGA Std Dev",
				"PGV Mean (cm/s)", "PGV Std Dev", "PGD Mean (cm)", "PGD Std Dev");
		imr.setEqkRupture(rup);
		String[] imts = {PGA_Param.NAME, PGV_Param.NAME, PGD_Param.NAME};
		for (Site site : sites) {
			imr.setSite(site);
			List<String> line = Lists.newArrayList();
			
			line.add(site.getLocation().getLatitude()+"");
			line.add(site.getLocation().getLongitude()+"");
			line.add(site.getParameter(Vs30_Param.NAME).getValue()+"");
			line.add(site.getParameter(DepthTo2pt5kmPerSecParam.NAME).getValue()+"");
			line.add(surf.getDistanceRup(site.getLocation())+"");
			for (String imt : imts) {
				imr.setIntensityMeasure(imt);
				double mean = Math.exp(imr.getMean());
				double stdDev = Math.exp(imr.getStdDev());
				line.add(mean+"");
				line.add(stdDev+"");
			}
			csv.addLine(line);
		}
		csv.writeToFile(new File("/tmp/cb_2008_results.csv"));
	}

}
