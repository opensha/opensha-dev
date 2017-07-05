package scratch.kevin.cybershake;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.impl.CVM4i26BasinDepth;
import org.opensha.commons.data.siteData.impl.CVM_Vs30;
import org.opensha.commons.data.siteData.impl.WillsMap2006;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.data.siteData.impl.CVM_Vs30.CVM;
import org.opensha.sha.cybershake.db.CachedPeakAmplitudesFromDB;
import org.opensha.sha.cybershake.db.CybershakeIM;
import org.opensha.sha.cybershake.db.CybershakeIM.CyberShakeComponent;
import org.opensha.sha.cybershake.db.CybershakeIM.IMType;
import org.opensha.sha.cybershake.db.CybershakeRun;
import org.opensha.sha.cybershake.db.CybershakeSite;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.db.ERF2DB;
import org.opensha.sha.cybershake.db.MeanUCERF2_ToDB;
import org.opensha.sha.cybershake.db.PeakAmplitudesFromDB;
import org.opensha.sha.cybershake.db.Runs2DB;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.faultSurface.RuptureSurface;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Range;
import com.google.common.primitives.Doubles;

public class ChristineMetadataGen {

	public static void main(String[] args) {
		DBAccess db = null;
		DBAccess dbForAmps = null;
		try {
			Range<Double> distRange = Range.open(0d, 70d);
			Range<Double> magRange = Range.open(6d, 6.8d);
			
			CyberShakeComponent comp = CyberShakeComponent.RotD50;
			// if we don't have curves for the former this is needed
			CyberShakeComponent compForID = CyberShakeComponent.GEOM_MEAN;
			double[] periods = { 0.5, 2, 5 };
			
			int hazardDatasetID = 61; // study 15.12
			String[] siteNames = { "LADT" };
			AbstractERF erf = MeanUCERF2_ToDB.createUCERF2ERF();
			
			List<String> siteDataNames = Lists.newArrayList();
			List<SiteData<Double>> siteDataProvs = Lists.newArrayList();
			siteDataProvs.add(new WillsMap2006());
			siteDataNames.add("Wills 2006 Vs30 (m/s)");
			siteDataProvs.add(new CVM_Vs30(CVM.CVMS4i26));
			siteDataNames.add("CVMS-4.26 Vs30 (m/s)");
			siteDataProvs.add(new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_1_0));
			siteDataNames.add("CVMS-4.26 Z1.0 (km)");
			siteDataProvs.add(new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_2_5));
			siteDataNames.add("CVMS-4.26 Z2.5 (km)");
			
			db = Cybershake_OpenSHA_DBApplication.getDB();
//			PeakAmplitudesFromDB amps2db = new PeakAmplitudesFromDB(db);
			dbForAmps = new DBAccess("moment.usc.edu", Cybershake_OpenSHA_DBApplication.DATABASE_NAME);
			PeakAmplitudesFromDB amps2db = new CachedPeakAmplitudesFromDB(
					dbForAmps, new File("/home/kevin/CyberShake/MCER/.amps_cache"), erf);
			ERF2DB erf2db = new ERF2DB(db);
			Runs2DB runs2db = new Runs2DB(db);
			List<CybershakeIM> imsForID = amps2db.getIMs(Doubles.asList(periods), IMType.SA, compForID);
			List<CybershakeIM> ims = amps2db.getIMs(Doubles.asList(periods), IMType.SA, comp);
			Preconditions.checkState(ims.size() == periods.length);
			for (int i=0; i<periods.length; i++)
				System.out.println("Mapped period "+periods[i]+" to: "+ims.get(i));
			
			RunIDFetcher runFetch = new RunIDFetcher(db, hazardDatasetID, imsForID.get(0).getID());
			
			CSVFile<String> csv = new CSVFile<String>(true);
			
			List<String> header = Lists.newArrayList(
					"Site Short Name",
					"Site ID",
					"Site Lat",
					"Site Lon",
					"Run ID");
			header.addAll(siteDataNames);
			header.addAll(Lists.newArrayList(
					"Source Name",
					"Source ID",
					"Rupture ID",
					"Magnitude",
					"Distance Rup (km)",
					"Distance X (km)",
					"Distance JB (km)",
					"Ztor (km)",
					"Rup Length (km)",
					"Rup Width (km)",
					"Dip",
					"Rake",
					"RV ID",
					"Hypo Lat",
					"Hypo Lon",
					"Hypo Depth (km)",
					"Epicenter Distance (km)"));
			for (CybershakeIM im : ims)
				header.add(PeakAmplitudesFromDB.getCleanedCS_Period(im.getVal())+" s SA, "
						+im.getComponent()+" ("+im.getUnits()+")");
			csv.addLine(header);
			
			for (String siteName : siteNames) {
				Integer runID = runFetch.getRunID(siteName);
				System.out.println(siteName+", run "+runID);
				Preconditions.checkNotNull(runID, "No run found for %s", siteName);
				CybershakeRun run = runs2db.getRun(runID);
				CybershakeSite site = runFetch.getSite(siteName);
				
				List<String> lineTemplate = Lists.newArrayList(
						site.short_name, site.id+"", site.lat+"", site.lon+"", runID+"");
				Location loc = site.createLocation();
				for (int i=0; i<siteDataProvs.size(); i++)
					lineTemplate.add(siteDataProvs.get(i).getValue(loc)+"");
				
				for (int sourceID=0; sourceID<erf.getNumSources(); sourceID++) {
					String sourceName = erf.getSource(sourceID).getName();
					System.out.println("Source "+sourceID+": "+sourceName);
					for (int rupID=0; rupID<erf.getNumRuptures(sourceID); rupID++) {
						ProbEqkRupture rup = erf.getRupture(sourceID, rupID);
						
						double mag = rup.getMag();
						if (!magRange.contains(mag))
							continue;
						
						RuptureSurface surf = rup.getRuptureSurface();
						
						double rRup = surf.getDistanceRup(loc);
						if (!distRange.contains(rRup))
							continue;
						
//						"Source Name",
//						"Source ID",
//						"Rupture ID",
//						"Magnitude",
//						"Distance Rup (km)",
//						"Distance X (km)",
//						"Distance JB (km)",
//						"Ztor (km)",
//						"Rup Length (km)",
//						"Rup Width (km)",
//						"Dip",
//						"Rake",
//						"RV ID",
//						"Hypo Lat",
//						"Hypo Lon",
//						"Hypo Depth (km)",
//						"Epicenter Distance (km)"
						List<String> rupLine = Lists.newArrayList(lineTemplate);
						rupLine.add(sourceName);
						rupLine.add(sourceID+"");
						rupLine.add(rupID+"");
						rupLine.add(mag+"");
						rupLine.add(rRup+"");
						rupLine.add(surf.getDistanceX(loc)+"");
						rupLine.add(surf.getDistanceJB(loc)+"");
						rupLine.add(surf.getAveRupTopDepth()+"");
						rupLine.add(surf.getAveLength()+"");
						rupLine.add(surf.getAveWidth()+"");
						rupLine.add(surf.getAveDip()+"");
						rupLine.add(rup.getAveRake()+"");
						List<List<Double>> imVals = Lists.newArrayList();
						for (CybershakeIM im : ims)
							imVals.add(amps2db.getIM_Values(runID, sourceID, rupID, im));
						Map<Integer, Location> hypos =
								erf2db.getHypocenters(run.getERFID(), sourceID, rupID, run.getRupVarScenID());
						for (int rvID=0; rvID<imVals.get(0).size(); rvID++) {
							Location hypo = hypos.get(rvID);
							Preconditions.checkNotNull(hypo);
							
							List<String> line = Lists.newArrayList(rupLine);
							line.add(rvID+"");
							line.add(hypo.getLatitude()+"");
							line.add(hypo.getLongitude()+"");
							line.add(hypo.getDepth()+"");
							line.add(LocationUtils.horzDistanceFast(loc, hypo)+"");
							
							for (int i=0; i<ims.size(); i++)
								line.add(imVals.get(i).get(rvID)+"");
							csv.addLine(line);
						}
					}
				}
			}
			
			System.out.println("Writing "+(csv.getNumRows()-1)+" records");
			csv.writeToFile(new File("/tmp/cybershake_data.csv"));
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} finally {
			if (db != null)
				db.destroy();
			if (dbForAmps != null && dbForAmps != db)
				dbForAmps.destroy();
		}
	}

}
