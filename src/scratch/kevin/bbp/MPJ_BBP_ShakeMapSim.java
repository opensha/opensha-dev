package scratch.kevin.bbp;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FileUtils;
import org.opensha.commons.util.XMLUtils;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.kevin.bbp.BBP_Module.Method;
import scratch.kevin.bbp.BBP_Module.VelocityModel;

public class MPJ_BBP_ShakeMapSim extends MPJTaskCalculator {
	
	private VelocityModel vm;
	private Method method;
	
	private File srcFile;
	private File srfFile;
	private File mainOutputDir;
	
	private File resultsDir;
	private boolean doHF = true;
	
	private List<BBP_Site> sites;
	
	private File bbpDataDir = null;
	
	private ExecutorService exec;
	
	private boolean keepSRFs = false;
	
	private static final double DIST_BUFFER_DEFAULT = 200d;
	
	private static final int BUNDLE_SIZE = 1000;

	public MPJ_BBP_ShakeMapSim(CommandLine cmd) throws IOException {
		super(cmd);
		
		vm = VelocityModel.valueOf(cmd.getOptionValue("vm"));
		method = Method.valueOf(cmd.getOptionValue("method"));
		srcFile = new File(cmd.getOptionValue("src-file"));
		Preconditions.checkState(srcFile.exists());
		if (cmd.hasOption("srf-file")) {
			srfFile = new File(cmd.getOptionValue("srf-file"));
			Preconditions.checkState(srfFile.exists());
			if (rank == 0)
				debug("Using input SRF file: "+srfFile.getAbsolutePath());
		}
		mainOutputDir = new File(cmd.getOptionValue("output-dir"));
		resultsDir = new File(mainOutputDir, "results");
		if (rank == 0) {
			Preconditions.checkState((mainOutputDir.exists() && mainOutputDir.isDirectory()) || mainOutputDir.mkdir());
			Preconditions.checkState((resultsDir.exists() && resultsDir.isDirectory()) || resultsDir.mkdir());
		}
		doHF = !cmd.hasOption("no-hf");
		if (cmd.hasOption("bbp-data-dir")) {
			bbpDataDir = new File(cmd.getOptionValue("bbp-data-dir"));
			if (rank == 0)
				debug("BBP data dir: "+bbpDataDir.getAbsolutePath());
		}
		
		// define region
		double buffer = DIST_BUFFER_DEFAULT;
		if (cmd.hasOption("buffer-dist"))
			buffer = Double.parseDouble(cmd.getOptionValue("buffer-dist"));
		double spacing = Double.parseDouble(cmd.getOptionValue("spacing"));
		BBP_SourceFile src = BBP_SourceFile.readFile(srcFile);
		Location[] topPoints = Arrays.copyOf(src.getSurface().getRectangle(), 2);
		double[] azimuths = { 0, Math.PI/2d, Math.PI, 3d*Math.PI/2d };
		MinMaxAveTracker latTrack = new MinMaxAveTracker();
		MinMaxAveTracker lonTrack = new MinMaxAveTracker();
		for (Location loc : topPoints) {
			for (double az : azimuths) {
				Location oLoc = LocationUtils.location(loc, az, buffer);
				latTrack.addValue(oLoc.getLatitude());
				lonTrack.addValue(oLoc.getLongitude());
			}
		}
		double minLat = latTrack.getMin();
		double maxLat = latTrack.getMax();
		double minLon = lonTrack.getMin();
		double maxLon = lonTrack.getMax();
		// make them nicer
		minLat = Math.floor(minLat*10d)/10d;
		maxLat = Math.ceil(maxLat*10d)/10d;
		minLon = Math.floor(minLon*10d)/10d;
		maxLon = Math.ceil(maxLon*10d)/10d;
		GriddedRegion reg = new GriddedRegion(new Location(minLat, minLon),
				new Location(maxLat, maxLon), spacing, GriddedRegion.ANCHOR_0_0);
		
		if (rank == 0)
			debug("Building "+reg.getNodeCount()+" sites for ["+minLat+","+maxLat+"], ["+minLon+","+maxLon+"] with spacing="+spacing);
		sites = new ArrayList<>(reg.getNodeCount());
		for (int i=0; i<reg.getNodeCount(); i++)
			sites.add(new BBP_Site("s"+i, reg.getLocation(i), vm.getVs30(), 0.5, 100));
		if (rank == 0) {
			debug("Writing sites XML and STL files");
			File sitesXMLFile = new File(mainOutputDir, "sites.xml");
			XMLUtils.writeObjectToXMLAsRoot(reg, sitesXMLFile);
			File sitesSTLFile = new File(mainOutputDir, "sites.stl");
			BBP_Site.writeToFile(sitesSTLFile, sites);
		}
		
		exec = Executors.newFixedThreadPool(getNumThreads());
		
		if (rank == 0) {
			postBatchHook = new MasterZipHook();
			
			// initialize bundle dirs
			for (int i=0; i<=getNumTasks(); i+=BUNDLE_SIZE)
				getRunBundleDir(i, true);
		}
	}
	
	private File getRunBundleDir(int index, boolean create) {
		int siteBase = BUNDLE_SIZE*(int)(index / BUNDLE_SIZE);
		File parentDir = new File(resultsDir, "sites_"+siteBase+"_"+(siteBase+BUNDLE_SIZE));
		if (create)
			MPJ_BBP_Utils.waitOnDir(parentDir, 10, 2000);
		return parentDir;
	}
	
	private File getRunDir(int index, boolean create) {
		File runDir = new File(getRunBundleDir(index, create), "site_"+index);
		if (create)
			MPJ_BBP_Utils.waitOnDir(runDir, 10, 2000);
		return runDir;
	}
	
	private File getZipFile(int index) {
		File runDir = getRunDir(index, false);
		return new File(runDir.getParentFile(), runDir.getName()+".zip");
	}

	@Override
	protected int getNumTasks() {
		return sites.size();
	}
	
	private class MasterZipHook extends MPJ_BBP_Utils.MasterZipHook {

		public MasterZipHook() {
			super(new File(resultsDir.getParentFile(), resultsDir.getName()+".zip"),
					new File(resultsDir.getParentFile(), resultsDir.getName()+"_rotD.zip"));
		}

		@Override
		protected void debug(String message) {
			MPJ_BBP_ShakeMapSim.this.debug(message);
		}

		@Override
		protected void abortAndExit(int status) {
			MPJTaskCalculator.abortAndExit(status);
		}

		@Override
		protected File getSimZipFile(int index) {
			return getZipFile(index);
		}
		
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		ArrayList<Future<?>> futures = new ArrayList<>();
		for (int index : batch) {
			Task task = new Task(index);
			futures.add(exec.submit(task));
			Thread.sleep(100);
		}
		for (Future<?> f : futures)
			f.get();
	}
	
	private class Task implements Runnable {
		
		private int index;

		public Task(int index) {
			this.index = index;
		}

		@Override
		public void run() {
			try {
				File runDir = getRunDir(index, true);
				File zipFile = getZipFile(index);
				if (zipFile.exists()) {
					debug(index+" is already done, skipping");
					return;
				}
				// write site
				BBP_Site site = sites.get(index);
				File sitesFile = new File(runDir, site.getName()+".stl");
				List<BBP_Site> mySites = new ArrayList<>(1);
				mySites.add(site);
				BBP_Site.writeToFile(sitesFile, mySites);
				
				Preconditions.checkNotNull(srcFile, "No src file in %s", runDir.getAbsolutePath());
				BBP_Wrapper wrapper = new BBP_Wrapper(vm, method, srcFile, null, srfFile, sitesFile, runDir);
				wrapper.setDoHF(doHF);
				wrapper.setDataOnly(true);
				
				// run BBP
				debug("running BBP for "+index);
				wrapper.setDoHF(doHF);
				wrapper.setDataOnly(true);
				wrapper.setBBPDataDir(bbpDataDir);
				wrapper.run();
				
				if (!keepSRFs && srfFile == null) {
					for (File file : runDir.listFiles())
						if (file.getName().endsWith(".srf"))
							file.delete();
				}
				
				// zip it
				FileUtils.createZipFile(zipFile, runDir, true);
				FileUtils.deleteRecursive(runDir);
			} catch (IOException e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
		
	}

	@Override
	protected void doFinalAssembly() throws Exception {
		exec.shutdown();
		if (rank == 0) {
			debug("waiting for any post batch hook operations to finish");
			((MasterZipHook)postBatchHook).shutdown();
			debug("post batch hook done");
		}
	}
	
	public static Options createOptions() {
		Options ops = MPJ_BBP_Utils.addCommonOptions(MPJTaskCalculator.createOptions(), false, true, true, false);
		
		Option dist = new Option("dist", "buffer-dist", true,
				"Buffer distance around the surface for region. Default: "+DIST_BUFFER_DEFAULT+" km");
		dist.setRequired(false);
		ops.addOption(dist);
		
		Option discr = new Option("spc", "spacing", true, "Region discretization in degrees");
		discr.setRequired(true);
		ops.addOption(discr);
		
		return ops;
	}
	
	public static void main(String[] args) {
		args = MPJTaskCalculator.initMPJ(args);
		
		try {
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_BBP_ShakeMapSim.class);
			
			MPJ_BBP_ShakeMapSim driver = new MPJ_BBP_ShakeMapSim(cmd);
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
