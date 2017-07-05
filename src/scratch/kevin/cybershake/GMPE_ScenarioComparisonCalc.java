package scratch.kevin.cybershake;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.siteData.OrderedSiteDataProviderList;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.geo.Location;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.ParameterList;
import org.opensha.sha.cybershake.HazardCurveFetcher;
import org.opensha.sha.cybershake.db.CybershakeSite;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.db.MeanUCERF2_ToDB;
import org.opensha.sha.cybershake.plot.HazardCurvePlotter;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.SiteParams.DepthTo1pt0kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.util.SiteTranslator;

import com.google.common.collect.Lists;

public class GMPE_ScenarioComparisonCalc {

	public static void main(String[] args) throws IOException {
		DBAccess db = Cybershake_OpenSHA_DBApplication.getDB();
		
//		Double forceVs30 = null;
		Double forceVs30 = 760d;
		
		List<ScalarIMR> imrs = Lists.newArrayList();
		imrs.add(AttenRelRef.ASK_2014.instance(null));
		imrs.add(AttenRelRef.BSSA_2014.instance(null));
		imrs.add(AttenRelRef.CB_2014.instance(null));
		imrs.add(AttenRelRef.CY_2014.instance(null));
		
		AbstractERF erf = MeanUCERF2_ToDB.createUCERF2ERF();
//		int sourceID = 64;
//		int rupID = 3;
//		int sourceID = 184;
//		int rupID = 6;
//		int sourceID = 158;
//		int rupID = 6;
//		int sourceID = 242;
//		int rupID = 26;
		int sourceID = 265;
		int rupID = 94;
		ProbEqkRupture rup = erf.getRupture(sourceID, rupID);
		
		System.out.println("Scenario: M"+rup.getMag()+" on "+erf.getSource(sourceID).getName());
		
		double[] periods = { 2,3,5 };
		
		String outPrefix = "gmpe_amps_source"+sourceID+"_rup"+rupID;
		String vs30Name = "Wills 2006 Vs30 (m/s)";
		if (forceVs30 != null) {
			outPrefix += "_forceVs30"+forceVs30;
			vs30Name = "Hardcoded Vs30 (m/s)";
		}
		
		ParameterList siteParams = new ParameterList();
		for (ScalarIMR imr : imrs) {
			imr.setParamDefaults();
			imr.setIntensityMeasure(SA_Param.NAME);
			for (Parameter<?> param : imr.getSiteParams())
				if (!siteParams.containsParameter(param))
					siteParams.addParameter((Parameter<?>) param.clone());
		}
		
		int datasetID = 57;
		int velModelID = 5;
		int imTypeID = 21; // doesn't matter
		HazardCurveFetcher fetch = new HazardCurveFetcher(db, datasetID, imTypeID);
		
		List<CybershakeSite> sites = fetch.getCurveSites();
		System.out.println("Found "+sites.size()+" sites");
		
		Collections.sort(sites, new Comparator<CybershakeSite>() {

			@Override
			public int compare(CybershakeSite o1, CybershakeSite o2) {
				return new Integer(o1.id).compareTo(o2.id);
			}
		});
		
		List<String> header = Lists.newArrayList("Site Short Name", "Site ID",
				vs30Name, "CVM-S4.26 Z1.0 (km)", "CVM-S4.26 Z2.5 (km)");
		for (ScalarIMR imr : imrs)
			header.add(imr.getShortName());
		
		List<CSVFile<String>> csvs = Lists.newArrayList();
		for (int i=0; i<periods.length; i++) {
			CSVFile<String> csv = new CSVFile<String>(true);
			csv.addLine(header);
			csvs.add(csv);
		}
		
		OrderedSiteDataProviderList provs = HazardCurvePlotter.createProviders(velModelID, false);
		SiteTranslator trans = new SiteTranslator();
		
		for (CybershakeSite csSite : sites) {
			if (csSite.type_id == CybershakeSite.TYPE_TEST_SITE)
				continue;
			System.out.println("Calculating: "+csSite.short_name);
			
			Site site = new Site(csSite.createLocation());
			ArrayList<SiteDataValue<?>> datas = provs.getBestAvailableData(site.getLocation());
			for (Parameter<?> param : siteParams) {
				if (forceVs30 != null && param.getName().equals(Vs30_Param.NAME))
					trans.setParameterValue(param, new SiteDataValue<Double>(
							SiteData.TYPE_VS30, SiteData.TYPE_FLAG_INFERRED, forceVs30));
				else
					trans.setParameterValue(param, datas);
				site.addParameter(param);
			}
			
			double vs30 = site.getParameter(Double.class, Vs30_Param.NAME).getValue();
			double z1p0 = site.getParameter(Double.class, DepthTo1pt0kmPerSecParam.NAME).getValue()/1000d;
			double z2p5 = site.getParameter(Double.class, DepthTo2pt5kmPerSecParam.NAME).getValue();
			
			List<String> lineTemplate = Lists.newArrayList(csSite.short_name, csSite.id+"", vs30+"", z1p0+"", z2p5+"");
			List<List<String>> lines = Lists.newArrayList();
			for (int i=0; i<periods.length; i++)
				lines.add(Lists.newArrayList(lineTemplate));
			for (ScalarIMR imr : imrs) {
				imr.setSite(site);
				imr.setEqkRupture(rup);
				for (int i=0; i<periods.length; i++) {
					double period = periods[i];
					SA_Param.setPeriodInSA_Param(imr.getIntensityMeasure(), period);
					lines.get(i).add(Math.exp(imr.getMean())+"");
				}
			}
			for (int i=0; i<csvs.size(); i++)
				csvs.get(i).addLine(lines.get(i));
		}
		
		for (int i=0; i<periods.length; i++) {
			File outputFile = new File("/tmp", outPrefix+"_"+(int)periods[i]+"s.csv");
			csvs.get(i).writeToFile(outputFile);
		}
		
		db.destroy();
		
		// now write rupture surface
		RuptureSurface surf = rup.getRuptureSurface();
		String rupPrefix = "rup_surf_"+sourceID+"_"+rupID;
		FileWriter fw = new FileWriter(new File("/tmp", rupPrefix+"_full.txt"));
		for (Location loc : surf.getEvenlyDiscritizedListOfLocsOnSurface())
			fw.write(loc.getLatitude()+"\t"+loc.getLongitude()+"\t"+loc.getDepth()+"\n");
		fw.close();
		fw = new FileWriter(new File("/tmp", rupPrefix+"_upper.txt"));
		for (Location loc : surf.getUpperEdge())
			fw.write(loc.getLatitude()+"\t"+loc.getLongitude()+"\t"+loc.getDepth()+"\n");
		fw.close();
	}

}
