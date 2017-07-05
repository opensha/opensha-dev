package scratch.kevin.cybershake.ucerf3.safWallToWallTests;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.dom4j.DocumentException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.sha.cybershake.ERF_Rupture_File_Writer;
import org.opensha.sha.cybershake.db.CybershakeSite;
import org.opensha.sha.cybershake.db.CybershakeSiteInfo2DB;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.db.ERF2DB;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.faultSurface.RuptureSurface;

import com.google.common.collect.Lists;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.utils.FaultSystemIO;

public class FSSToCyberShakeDB {

	public static void main(String[] args) {
		DBAccess db = null;
		try {
			FaultSystemSolution fss = FaultSystemIO.loadSol(
					new File("/home/kevin/CyberShake/ucerf3/saf_wall_to_wall_tests/saf_subset_sol.zip"));
			
			String erfName = "UCERF3 SAF Downsampling Test ERF, 200m";
			String erfDescription = "Test UCERF3 ERF which contains only ruptures that are entirely on the SAF and at least"
					+ " partially within Southern California";
//			double gridSpacing = 0.2d;
			double gridSpacing = 1d;
			
			db = Cybershake_OpenSHA_DBApplication.getAuthenticatedDBAccess(true, true);
			if (db.isReadOnly())
				db.setIgnoreInserts(true);
			
			ERF2DB erf2db = new FSS_ERF2DB(fss, db, true, gridSpacing);
			File erfDir = new File("/home/kevin/CyberShake/rupSurfaces/37");
			erf2db.setFileBased(erfDir, true);
			System.out.println("ERF has "+erf2db.getERF_Instance().getNumSources()+" sources");
//			writeTestRup(erf2db.getERF_Instance());
			ERF_Rupture_File_Writer.writeRuptureFile(
					erf2db.getERF_Instance().getRupture(8182, 0), 8182, 0, new File("/tmp"), false);
//			erf2db.insertForecaseInDB(erfName, erfDescription);
//			int erfID = erf2db.getInserted_ERF_ID(erfName);
//			erf2db.insertSrcRupInDB(erfID, null, 10430, 0);
			
//			CybershakeSiteInfo2DB sites2db = new CybershakeSiteInfo2DB(db);
//			
//			ArrayList<CybershakeSite> sites = Lists.newArrayList();
//			sites.add(sites2db.getSiteFromDB(18)); // USC
//			
//			Cybershake_OpenSHA_DBApplication cs = new Cybershake_OpenSHA_DBApplication();
//			cs.setSiteInfoObject(sites2db);
//			boolean forceAdd = false; // can't remember what this is
//			cs.insertNewERFForSites(sites, erf2db, erfName, erfDescription, forceAdd);
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (db != null) {
				db.destroy();
				try {
					Thread.sleep(2000);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		System.exit(0);
	}
	
	private static void writeTestRup(ERF erf) throws IOException {
		RuptureSurface maxSurf = null;
		double maxMag = 0d;
		
		for (ProbEqkSource source : erf) {
			for (ProbEqkRupture rup : source) {
				double mag = rup.getMag();
				if (mag > maxMag) {
					maxMag = mag;
					maxSurf = rup.getRuptureSurface();
				}
			}
		}
		
		File outputFile = new File("/tmp/surf.txt");
		FileWriter fw = new FileWriter(outputFile);
		
		for (Location loc : maxSurf.getEvenlyDiscritizedListOfLocsOnSurface())
			fw.write(loc.getLatitude()+" "+loc.getLongitude()+" "+loc.getDepth()+"\n");
		
		fw.close();
	}

}
