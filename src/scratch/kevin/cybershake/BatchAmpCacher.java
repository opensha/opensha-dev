package scratch.kevin.cybershake;

import java.io.File;
import java.sql.SQLException;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.util.ClassUtils;
import org.opensha.sha.cybershake.HazardCurveFetcher;
import org.opensha.sha.cybershake.db.CachedPeakAmplitudesFromDB;
import org.opensha.sha.cybershake.db.CybershakeIM;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.db.HazardCurve2DB;
import org.opensha.sha.cybershake.db.MeanUCERF2_ToDB;
import org.opensha.sha.earthquake.ERF;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;
import com.google.common.primitives.Ints;

public class BatchAmpCacher {

	public static void main(String[] args) {
		File cacheDir = null;
		int datasetID = -1;
		int[] imTypeIDs = null;
		int killTimerMins = -1;
		if (args.length == 1 && args[0].equals("TEST_ECLIPSE")) {
			cacheDir = new File("/home/kevin/CyberShake/MCER/.amps_cache/");
			datasetID = 57;
//			imTypeIDs = new int[] { 151, 146, 142, 136 };
			imTypeIDs = new int[] { 151, 146, 144, 142, 138, 136 };
//			killTimerMins = 1;
		} else if (args.length == 3 || args.length == 4) {
			cacheDir = new File(args[0]);
			Preconditions.checkState(cacheDir.exists() || cacheDir.mkdir());
			datasetID = Integer.parseInt(args[1]);
			String[] imTypesStr = args[2].split(",");
			imTypeIDs = new int[imTypesStr.length];
			for (int i=0; i<imTypesStr.length; i++)
				imTypeIDs[i] = Integer.parseInt(imTypesStr[i]);
			if (args.length == 4)
				killTimerMins = Integer.parseInt(args[3]);
		} else {
			System.err.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(BatchAmpCacher.class)
					+" <cache-dir> <dataset-id> <im-type-ids> [<kill-timer-mins>]");
			System.exit(2);
		}
		System.out.println("Caching results for dataset "+datasetID+", IM Type IDs "
				+Joiner.on(",").join(Ints.asList(imTypeIDs))+" to "+cacheDir.getAbsolutePath());
		if (killTimerMins > 0) {
			System.out.println("Auto kill after "+killTimerMins+" minutes");
			new KillTimerThread(killTimerMins).start();
		}
		
		DBAccess db = null;
		int exitCode = 0;
		try {
			db = Cybershake_OpenSHA_DBApplication.getDB();
			
			ERF erf = MeanUCERF2_ToDB.createUCERF2ERF();
			
			CachedPeakAmplitudesFromDB.DD = true;
			CachedPeakAmplitudesFromDB amps2db = new CachedPeakAmplitudesFromDB(db, cacheDir, erf);
			HazardCurve2DB curves2db = new HazardCurve2DB(db);
			
			for (int imTypeID : imTypeIDs) {
				CybershakeIM im = curves2db.getIMFromID(imTypeID);
				
				System.out.println("Processing For IM: "+im);
				
				HazardCurveFetcher fetcher = new HazardCurveFetcher(db, datasetID, imTypeID);
				List<Integer> runIDs = fetcher.getRunIDs();
				
				int cached = 0;
				
				for (int i=0; i<runIDs.size(); i++) {
					int runID = runIDs.get(i);
					
					if (amps2db.isFileCached(runID, im))
						continue;
					
					System.out.println("Loading for run: "+runID+" ("+i+"/"+runIDs.size()+")");
					Stopwatch watch = Stopwatch.createStarted();
					amps2db.getAllIM_Values(runID, im);
					watch.stop();
					long secs = watch.elapsed(TimeUnit.SECONDS);
					String timeStr;
					if (secs > 60)
						timeStr = (float)((double)secs/60d)+" m";
					else
						timeStr = secs+" s";
					System.out.println("Done caching in "+timeStr);
					cached++;
				}
				System.out.println("Done caching amps for "+runIDs.size()+" runs ("+cached+" new)");
			}
		} catch (Throwable t) {
			t.printStackTrace();
			exitCode = 1;
		} finally {
			try {
				if (db != null)
					db.destroy();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		System.exit(exitCode);
	}
	
	private static class KillTimerThread extends Thread {
		
		private long millis;
		
		public KillTimerThread(long mins) {
			millis = mins * 60 * 1000;
		}
		
		@Override
		public void run() {
			try {
				Thread.sleep(millis);
			} catch (InterruptedException e) {}
			System.out.println("Killed after "+millis+" ms");
			System.out.flush();
			System.exit(1);
		}
	}

}
