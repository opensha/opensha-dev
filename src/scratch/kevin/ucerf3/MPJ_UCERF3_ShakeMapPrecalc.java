package scratch.kevin.ucerf3;

import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.nio.channels.FileChannel;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.dom4j.Document;
import org.dom4j.DocumentException;
import org.dom4j.Element;
import org.opensha.commons.data.Site;
import org.opensha.commons.geo.GeoTools;
import org.opensha.commons.hpc.mpj.taskDispatch.MPJTaskCalculator;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.XMLUtils;
import org.opensha.sha.cybershake.gui.util.AttenRelSaver;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.utils.FaultSystemIO;

/**
 * Calculates ShakeMaps (mean/std dev) for each grid node that are within a cutoff of each FaultSystemSolution rupture.
 * 
 * Calculation is organized by grid node, so each node is calculated and stored for each rupture at once. This is for efficiency
 * with distance calculations
 * 
 * @author kevin
 *
 */
public class MPJ_UCERF3_ShakeMapPrecalc extends MPJTaskCalculator {
	
	private double distCutoff;
	private FaultSystemSolution sol;
	private ScalarIMR[] gmpes;
	private List<Site> sites;
	private List<String> imts;
	private List<Double> saPeriods; // items null if corresponding IMT not SA
	private List<String> imtFileNames;
	
	private File outputDir;
	private File nodeOutputDir;
	private int numDigits;
	
	private FaultSystemSolutionERF erf;
	private ProbEqkSource[] sourcesForFSSRuptures;

	public MPJ_UCERF3_ShakeMapPrecalc(CommandLine cmd) throws IOException, DocumentException, InvocationTargetException {
		super(cmd);
		
		distCutoff = Double.parseDouble(cmd.getOptionValue("distance-cutoff"));
		
		File solFile = new File(cmd.getOptionValue("solution-file"));
		debug("Loading FSS from "+solFile.getAbsolutePath());
		sol = FaultSystemIO.loadSol(solFile);
		
		debug("Instantiating ERF");
		erf = new FaultSystemSolutionERF(sol);
		erf.updateForecast();
		
		// organize by FSS index
		sourcesForFSSRuptures = new ProbEqkSource[sol.getRupSet().getNumRuptures()];
		for (int sourceID=0; sourceID<erf.getNumFaultSystemSources(); sourceID++) {
			int fssIndex = erf.getFltSysRupIndexForSource(sourceID);
			sourcesForFSSRuptures[fssIndex] = erf.getSource(sourceID);
		}
		if (rank == 0)
			debug("Created "+erf.getNumFaultSystemSources()+"/"+sourcesForFSSRuptures.length+" sources");
		
		debug("Loading GMPEs");
		File gmpeFile = new File(cmd.getOptionValue("gmpe-file"));
		gmpes = new ScalarIMR[getNumThreads()];
		for (int i=0; i<gmpes.length; i++) {
			gmpes[i] = AttenRelSaver.LOAD_ATTEN_REL_FROM_FILE(gmpeFile.getAbsolutePath());
			gmpes[i].setParamDefaults();
		}
		
		debug("Loading sites");
		File sitesFile = new File(cmd.getOptionValue("sites-file"));
		Document siteDoc;
		try {
			siteDoc = XMLUtils.loadDocument(sitesFile);
		} catch (Exception e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		Element sitesRoot = siteDoc.getRootElement();
		ArrayList<Parameter<?>> paramsToAdd = new ArrayList<Parameter<?>>();
		for (Parameter<?> param : gmpes[0].getSiteParams())
			paramsToAdd.add(param);
		sites = Site.loadSitesFromXML(sitesRoot.element(Site.XML_METADATA_LIST_NAME), paramsToAdd);
		
		numDigits = ((getNumTasks()-1)+"").length();
		
		outputDir = new File(cmd.getOptionValue("output-dir"));
		Preconditions.checkState(rank != 0 || outputDir.exists() || outputDir.mkdir(),
				"Couldn't create main output dir: %s", outputDir.getAbsoluteFile());
		nodeOutputDir = new File(outputDir, "site_results");
		Preconditions.checkState(rank != 0 || nodeOutputDir.exists() || nodeOutputDir.mkdir(),
				"Couldn't create node output dir: %s", nodeOutputDir.getAbsoluteFile());
		
		String imtsOption = cmd.getOptionValue("imts");
		imts = Lists.newArrayList();
		imtFileNames = Lists.newArrayList();
		saPeriods = Lists.newArrayList();
		for (String imtCode : imtsOption.split(",")) {
			imtCode = imtCode.trim();
			if (imtCode.isEmpty())
				continue;
			imtCode = imtCode.toUpperCase();
			
			String imt, imtFileName;
			Double period;
			if (imtCode.equals(PGA_Param.NAME)) {
				imt = PGA_Param.NAME;
				imtFileName = "pga";
				period = null;
			} else if (imtCode.equals(PGV_Param.NAME)) {
				imt = PGV_Param.NAME;
				imtFileName = "pgv";
				period = null;
			} else {
				try {
					period = Double.parseDouble(imtCode);
				} catch (NumberFormatException e) {
					throw new IllegalStateException("Unknown IMT or couldn't parse SA period: "+imtCode+" (arg: "+imtsOption+")");
				}
				imt = SA_Param.NAME;
				imtFileName = "sa_"+period.floatValue()+"s";
			}
			
			imts.add(imt);
			imtFileNames.add(imtFileName);
			saPeriods.add(period);
		}
	}

	@Override
	protected int getNumTasks() {
		return sites.size();
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		for (int index : batch) {
			String prefix = index+"";
			while (prefix.length() < numDigits)
				prefix = "0"+prefix;
			prefix = "site_"+prefix;
			
			if (isAlreadyDone(index, prefix)) {
				debug(index+" is already done, skipping!");
				continue;
			}
			
			Site site = sites.get(index);
			
			List<SiteRuptureTask> siteRups = Lists.newArrayList();
			for (int fssIndex=0; fssIndex<sourcesForFSSRuptures.length; fssIndex++)
				siteRups.add(new SiteRuptureTask(fssIndex));
			
			ArrayDeque<SiteRuptureTask> deque = new ArrayDeque<SiteRuptureTask>(siteRups);
			
			if (gmpes.length == 1) {
				// single threaded
				new CalcRunnable(site, gmpes[0], deque).run();
			} else {
				// threaded
				List<Thread> threads = Lists.newArrayList();
				
				for (int i=0; i<gmpes.length; i++) {
					Thread thread = new Thread(new CalcRunnable(site, gmpes[i], deque));
					thread.start();
					threads.add(thread);
				}
				
				for (Thread thread : threads)
					thread.join();
			}
			
			// now write output
			int numRups = 0;
			for (SiteRuptureTask rup : siteRups)
				if (rup.means != null)
					numRups++;
			
			for (int i=0; i<imts.size(); i++) {
				String fileName = prefix+"_"+imtFileNames.get(i)+".bin";
				
				File file = new File(nodeOutputDir, fileName);
				
				DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file)));

				out.writeInt(index);
				out.writeDouble(site.getLocation().getLatitude());
				out.writeDouble(site.getLocation().getLongitude());
				out.writeInt(numRups);

				for (SiteRuptureTask rup : siteRups) {
					if (rup.means == null)
						continue;
					out.writeInt(rup.fssIndex);
					out.writeDouble(rup.means[i]);
					out.writeDouble(rup.stdDevs[i]);
				}

				out.close();
			}
		}
	}
	
	// site index (int), lat (double), lon (double), num rups (int)
	private static final long site_file_header_len = 4l + 8l + 8l + 4l;
	
	// fss index (int), mean (double), stdDev (double)
	private static final long file_len_per_rup = 4l + 8l + 8l;
	
	/**
	 * Checks if this site is already done by verifying that:
	 * * each IMT file exists and is non-empty
	 * * each IMT file has the same number of ruptures (and that number is within the allowable range)
	 * * each IMT file is the correct length given the number of ruptures
	 * * no I/O exceptions
	 * @param index
	 * @param sitePrefix
	 * @return
	 */
	private boolean isAlreadyDone(int index, String sitePrefix) {
		int numRups = -1;
		
		for (String imtName : imtFileNames) {
			File file = new File(nodeOutputDir, sitePrefix+"_"+imtName+".bin");
			long length = file.length();
			if (!file.exists() || length == 0)
				return false;
			
			// file exists
			if (length % 4 != 0)
				return false;
			
			DataInputStream in = null;
			try {
				in = new DataInputStream(new FileInputStream(file));
				
				int siteIndex = in.readInt();
				double lat = in.readDouble();
				double lon = in.readDouble();
				int myNumRups = in.readInt();
				in.close();
				
				if (siteIndex != index)
					return false;
				try {
					GeoTools.validateLat(lat);
					GeoTools.validateLat(lon);
				} catch (RuntimeException e) {
					return false;
				}
				
				if (myNumRups < 0 || myNumRups > sourcesForFSSRuptures.length)
					return false;
				
				if (numRups < 0)
					numRups = myNumRups;
				else if (myNumRups != numRups)
					return false;
				
				long calcLen = site_file_header_len + file_len_per_rup*myNumRups;
				if (calcLen != length)
					return false;
			} catch (IOException e) {
				try {
					if (in != null)
						in.close();
				} catch (IOException e1) {}
				return false;
			}
		}
		return true;
	}
	
	private class SiteRuptureTask {
		int fssIndex;
		double minDist;
		
		double[] means;
		double[] stdDevs;
		
		public SiteRuptureTask(int fssIndex) {
			this.fssIndex = fssIndex;
		}
		
		public void calc(Site site, ScalarIMR gmpe) {
			Preconditions.checkState(means == null);
			ProbEqkSource source = sourcesForFSSRuptures[fssIndex];
			if (source == null)
				return;
			Preconditions.checkState(source.getNumRuptures() == 1, "Must be a single rupture source");
			ProbEqkRupture rup = source.getRupture(0);
			
			minDist = source.getMinDistance(site);
			if (minDist > distCutoff)
				return;
			
			means = new double[imts.size()];
			stdDevs = new double[imts.size()];
			
			gmpe.setSite(site);
			gmpe.setEqkRupture(rup);
			
			for (int i=0; i<imts.size(); i++) {
				// set IMT
				String imt = imts.get(i);
				gmpe.setIntensityMeasure(imt);
				if (imt.equals(SA_Param.NAME))
					SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), saPeriods.get(i));
				
				means[i] = gmpe.getMean();
				stdDevs[i] = gmpe.getStdDev();
			}
		}
	}
	
	private static SiteRuptureTask popRuptureTask(Deque<SiteRuptureTask> rupTasks) {
		synchronized (rupTasks) {
			if (rupTasks.isEmpty())
				return null;
			return rupTasks.pop();
		}
	}
	
	private class CalcRunnable implements Runnable {
		
		private Site site;
		private ScalarIMR gmpe;
		private Deque<SiteRuptureTask> rupTasks;
		
		public CalcRunnable(Site site, ScalarIMR gmpe, Deque<SiteRuptureTask> rupTasks) {
			this.site = site;
			this.gmpe = gmpe;
			this.rupTasks = rupTasks;
		}

		@Override
		public void run() {
			while (true) {
				SiteRuptureTask task = popRuptureTask(rupTasks);
				if (task == null)
					break;
				
				task.calc(site, gmpe);
			}
		}
		
	}

	@Override
	protected void doFinalAssembly() throws Exception {
		if (rank != 0)
			return;
		
		// write combined output file
		
		for (int i=0; i<imts.size(); i++) {
			String imtName = imtFileNames.get(i);
			
			File outputFile = new File(outputDir, "results_"+imtName+".bin");
			
			// first write header with DataOutputStream
			DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
			out.writeInt(getNumTasks());
			out.close();
			
			// now write site results
			FileOutputStream outStream = new FileOutputStream(outputFile, true);
			FileChannel outChannel = outStream.getChannel();
			
			for (int index=0; index<getNumTasks(); index++) {
				if (index % 100 == 0)
					debug("Post processing site "+index);
				String prefix = index+"";
				while (prefix.length() < numDigits)
					prefix = "0"+prefix;
				prefix = "site_"+prefix;
				File siteFile = new File(nodeOutputDir, prefix+"_"+imtName+".bin");
				Preconditions.checkState(siteFile.exists(), "Site file doesn't exist! %s", siteFile.getAbsolutePath());
				
				FileInputStream inStream = new FileInputStream(siteFile);
				FileChannel inChannel = inStream.getChannel();
				
				outChannel.transferFrom(inChannel, outChannel.size(), inChannel.size());
				
				inChannel.close();
				inStream.close();
			}
			
			outChannel.close();
			outStream.close();
		}
	}
	
	public static Options createOptions() {
		Options ops = MPJTaskCalculator.createOptions();
		
		Option distCutoff = new Option("d", "distance-cutoff", true, "Distance cutoff in km");
		distCutoff.setRequired(true);
		ops.addOption(distCutoff);
		
		Option solFile = new Option("sol", "solution-file", true, "FaultSystemSolution file");
		solFile.setRequired(true);
		ops.addOption(solFile);
		
		Option gmpeFile = new Option("g", "gmpe-file", true, "GMPE XML file");
		gmpeFile.setRequired(true);
		ops.addOption(gmpeFile);
		
		Option siteFile = new Option("sites", "sites-file", true, "Sites XML file");
		siteFile.setRequired(true);
		ops.addOption(siteFile);
		
		Option imts = new Option("i", "imts", true, "IMTs. Comma separated values, options are 'PGA', 'PGV', and <sa-period>. "
				+ "Example: PGA,0.1,1.0 for PGA and 0.1s/1.0s Sa");
		imts.setRequired(true);
		ops.addOption(imts);
		
		Option outputDir = new Option("o", "output-dir", true, "Directory to store results");
		outputDir.setRequired(true);
		ops.addOption(outputDir);
		
		return ops;	
	}
	
	public static void main(String[] args) {
		args = MPJTaskCalculator.initMPJ(args);
		
		try {
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_UCERF3_ShakeMapPrecalc.class);
			
			MPJ_UCERF3_ShakeMapPrecalc calc = new MPJ_UCERF3_ShakeMapPrecalc(cmd);
			calc.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}
}
