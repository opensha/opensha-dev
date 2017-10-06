package scratch.kevin.bbp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.dom4j.Document;
import org.dom4j.Element;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.XMLUtils;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import scratch.kevin.bbp.BBP_Module.Method;
import scratch.kevin.bbp.BBP_Module.VelocityModel;

public class BBP_Wrapper implements Runnable {
	
	private static final boolean D = false;
	
	private static final String RESULTS_PREFIX = "You can find results in ";
	
	private VelocityModel vm;
	private Method method;
	
	private File srcFile;
	private Long seedOverride;
	private File srfFile;
	private File sitesFile;
	private File outputDir;
	
	private int maxRetries = 5;
	private boolean doHF = true;
	private boolean doRotD50 = true;
	private boolean doFAS = true;
	
	private String xmlFileName = "inputs.xml";
	private String scriptFileName = "run.sh";

	public BBP_Wrapper(VelocityModel vm, Method method, File srcFile, Long seedOverride, File srfFile, File sitesFile,
			File outputDir) {
		super();
		this.vm = vm;
		this.method = method;
		this.srcFile = srcFile;
		this.seedOverride = seedOverride;
		this.srfFile = srfFile;
		this.sitesFile = sitesFile;
		this.outputDir = outputDir;
	}
	
	public void run() {
		try {
			int retries = 0;
			boolean success = false;
			while (retries < maxRetries) {
				success = doRun();
				if (success)
					break;
				retries++;
				System.out.println("FAILED! Rety "+retries);
			}
			Preconditions.checkState(success, "BBP run failed with %s retries!", retries);
					
		} catch (IOException e) {
			ExceptionUtils.throwAsRuntimeException(e);
		}
	}
	
	static File writeSrcFileNewSeed(File srcFile, long seed, File outputDir) throws IOException {
		String fileName = srcFile.getName();
		if (fileName.endsWith(".src"))
			fileName = fileName.substring(0, fileName.indexOf(".src"));
		fileName += "_seed_"+seed+".src";
		File newSrcFile = new File(outputDir, fileName);
		FileWriter fw = new FileWriter(newSrcFile);
		
		boolean found = false;
		for (String line : Files.readLines(srcFile, Charset.defaultCharset())) {
			line = line.trim();
			if (line.startsWith("SEED")) {
				found = true;
				line = "SEED = "+seed;
			}
			fw.write(line+"\n");
		}
		
		fw.close();
		Preconditions.checkState(found, "SEED not found in input file");
		return newSrcFile;
	}

	void setScriptFileName(String scriptFileName) {
		this.scriptFileName = scriptFileName;
	}

	public void setMaxRetries(int maxRetries) {
		this.maxRetries = maxRetries;
	}

	public void setDoHF(boolean doHF) {
		this.doHF = doHF;
	}

	public void setDoRotD50(boolean doRotD50) {
		this.doRotD50 = doRotD50;
	}

	public void setDoFAS(boolean doFAS) {
		this.doFAS = doFAS;
	}

	public void setXMLFileName(String xmlFileName) {
		this.xmlFileName = xmlFileName;
	}

	private boolean doRun() throws IOException {
		Preconditions.checkState(srcFile.exists(), "Source file doesn't exist: %s", srcFile.getAbsolutePath());
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		if (seedOverride != null)
			srcFile = writeSrcFileNewSeed(srcFile, seedOverride, outputDir);
		
		long simID = System.currentTimeMillis();
		
		File bbpRunDir = doRun(outputDir, scriptFileName, xmlFileName, simID,
				srcFile, srfFile, sitesFile, vm, method, doHF, doRotD50, doFAS);
		
//		File inputsFile = new File(outputDir, inputsFileName);
//		
//		File srfFile = this.srfFile;
//		if (disableHF) {
//			String newSRFPrefix = inputsFileName;
//			if (newSRFPrefix.endsWith(".txt"))
//				newSRFPrefix = newSRFPrefix.substring(0, newSRFPrefix.indexOf(".txt"));
//			File newSRF = new File(outputDir, newSRFPrefix+"_noHF.srf");
//			if (srfFile != null) {
//				// easy, just remove the header from the SRF file
//				removeSRF_PlaneHeader(srfFile, newSRF);
//				srfFile = newSRF;
//			} else {
//				// have to run the platform to generate the SRF, then re-run
//				File slipImputFile = new File(outputDir, "genslip_"+inputsFileName);
//				writeBBP_InputsFile(slipImputFile, vmIndex, methodIndex, srcFile, srfFile, sitesFile);
//				File slipScriptfile = new File(outputDir, "genslip_"+scriptFileName);
//				writeBBP_Script(slipScriptfile, slipImputFile, System.currentTimeMillis(), "GenSlip");
//				
//				System.out.println("Running GenSlip only BBP with outputDir: "+outputDir.getAbsolutePath());
//				File bbpRunDir = runBBP_Script(slipScriptfile);
//				if (bbpRunDir == null)
//					return false;
//				String srcPrefix = srcFile.getName();
//				srcPrefix = srcPrefix.substring(0, srcPrefix.indexOf(".src"));
//				File genSRF = new File(bbpRunDir, srcPrefix+".srf");
//				Preconditions.checkState(genSRF.exists());
//				removeSRF_PlaneHeader(genSRF, newSRF);
//				srfFile = newSRF;
//			}
//		}
//		
//		writeBBP_InputsFile(inputsFile, vmIndex, methodIndex, srcFile, srfFile, sitesFile);
//		
//		File scriptFile = new File(outputDir, scriptFileName);
//		writeBBP_Script(scriptFile, inputsFile, System.currentTimeMillis());
//		
//		System.out.println("Running BBP with outputDir: "+outputDir.getAbsolutePath());
//		File bbpRunDir = runBBP_Script(scriptFile);
		if (bbpRunDir == null)
			return false;
		
		System.out.println("Detected BBP dir, moving data: "+bbpRunDir.getAbsolutePath());
		Preconditions.checkState(bbpRunDir.exists());
		for (File file : bbpRunDir.listFiles())
			move(file, new File(outputDir, file.getName()));
		bbpRunDir.delete();
		
		return true;
	}
	
	private static File doRun(File outputDir, String scriptFileName, String xmlFileName, long simID,
			File srcFile, File srfFile, File siteFile, VelocityModel vm, Method method,
			boolean doHF, boolean doRotD50, boolean doFAS) throws IOException {
		Preconditions.checkState(method == Method.GP, "Only GP supported currently");
		List<BBP_Module> modules = new ArrayList<>();
		if (srfFile == null)
			modules.add(BBP_Module.buildGenSlip(vm, srcFile));
		modules.add(BBP_Module.buildJBSim(vm, srcFile, srfFile, siteFile));
		if (doHF) {
			modules.add(BBP_Module.buildHFSims(vm, srcFile, srfFile, siteFile));
			modules.add(BBP_Module.buildMatch(vm, siteFile));
		} else {
			String preScriptName = "lf_only_"+scriptFileName;
			String preXMLFileName = "lf_only_"+xmlFileName;
			File dir = doRun(outputDir, preScriptName, preXMLFileName, simID, modules, siteFile);
			if (dir == null)
				return null;
			File tempDir = new File(new File(dir.getParentFile().getParentFile(), "tmpdata"), dir.getName());
			Preconditions.checkState(tempDir.exists(), "Temp dir doesn't exist: %s", tempDir.getAbsolutePath());
			if (D) System.out.println("Processing temp dir for velocity seismograms: "+tempDir.getAbsolutePath());
			for (File file : tempDir.listFiles()) {
				String name = file.getName();
				if (name.endsWith("-lf.bbp")) {
					if (D) System.out.println("Processing "+file.getAbsolutePath());
					// it's the velocity low frequency seismogram
					Preconditions.checkState(!name.contains("acc"));
					String prefix = name.substring(0, name.indexOf("-lf.bbp"));
					
					// first copy it
					File velDest = new File(dir, prefix+".vel.bbp");
					if (D) System.out.println("Copying "+file.getAbsolutePath()+" to "+velDest.getAbsolutePath());
					Files.copy(file, velDest);
					
					// now create acceleration file
					File accelDest = new File(dir, prefix+".acc.bbp");
					if (D) System.out.println("Creating acceleration file: "+accelDest.getAbsolutePath());
					velToAccel(file, accelDest);
				}
			}
			if (D) System.out.println("DONE processing velocity seismograms");
			modules.clear();
			scriptFileName = "pp_only_"+scriptFileName;
			xmlFileName = "pp_only_"+xmlFileName;
		}
		modules.add(BBP_Module.buildPlotMap(srcFile, srfFile, siteFile));
		modules.add(BBP_Module.buildPlotSeis(siteFile, srcFile));
		if (doRotD50)
			modules.add(BBP_Module.buildRotD50(siteFile));
		if (doFAS)
			modules.add(BBP_Module.buildFAS(siteFile));
		modules.add(BBP_Module.buildGenHTML(vm, method, siteFile, srcFile));
		
		return doRun(outputDir, scriptFileName, xmlFileName, simID, modules, siteFile);
	}
	
	private static File doRun(File outputDir, String scriptFileName, String xmlFileName, long simID,
			List<BBP_Module> modules, File siteFile) throws IOException {
		File scriptFile = new File(outputDir, scriptFileName);
		File xmlFile = new File(outputDir, xmlFileName);
		writeBBP_XML_File(xmlFile, siteFile, modules);
		writeBBP_Script(scriptFile, xmlFile, simID);
		
		System.out.println("Running BBP with simulation ID "+simID);
		return runBBP_Script(scriptFile);
	}
	
	private static void writeBBP_Script(File scriptFile, File xmlFile, long simID) throws IOException {
		FileWriter fw = new FileWriter(scriptFile);
		fw.write("#!/bin/bash"+"\n");
		fw.write("\n");
		fw.write("if [ -f ~/.bash_profile ]; then\n");
		fw.write("    . ~/.bash_profile\n");
		fw.write("fi\n");
		fw.write("\n");
		fw.write("if [ -z ${BBP_DATA_DIR+x} ];then\n");
		fw.write("    echo \"BBP_DATA_DIR is undefined\"\n");
		fw.write("    exit 2\n");
		fw.write("fi\n");
		fw.write("if [ ! -e $BBP_DATA_DIR ];then\n");
		fw.write("    echo \"BBP_DATA_DIR doesn't exist\"\n");
		fw.write("    exit 2\n");
		fw.write("fi\n");
		fw.write("\n");
		fw.write("export BBP_DATA_DIR=${BBP_DATA_DIR}/$HOST\n");
		fw.write("if [ ! -e $BBP_DATA_DIR ];then\n");
		fw.write("    mkdir $BBP_DATA_DIR\n");
		fw.write("fi\n");
		fw.write("\n");
		String command = "run_bbp.py --sim-id "+simID+" --xml-file="+xmlFile.getAbsolutePath();
//		if (endModule != null)
//			command += " --end "+endModule;
		fw.write(command+"\n");
		fw.close();
		scriptFile.setExecutable(true);
	}
	
	private static File runBBP_Script(File bbpScript) throws IOException {
		String[] command = { "/bin/bash", "-c", "bash "+bbpScript.getAbsolutePath() };
		
		Process p = Runtime.getRuntime().exec(command);
		int exit;
		try {
			exit = p.waitFor();
		} catch (InterruptedException e) {
			e.printStackTrace();
			return null;
		}
		List<String> stdout = streamToLines(p.getInputStream());
		List<String> stderr = streamToLines(p.getErrorStream());
		
		if (D) {
			System.out.println("======= STDOUT =======");
			for (String line : stdout)
				System.out.println(line);
			System.out.println("======================");
			System.out.println();
		}
		if (D || exit != 0) {
			System.out.println();
			System.out.println("======= STDERR =======");
			for (String line : stderr)
				System.out.println(line);
			System.out.println("======================");
		}
		
		System.out.println("Exit: "+exit);
		
		for (int i=stdout.size(); --i>=0;) {
			String line = stdout.get(i).trim();
			if (line.startsWith(RESULTS_PREFIX)) {
				File dir = new File(line.substring(RESULTS_PREFIX.length()).trim());
				return dir;
			}
		}
		return null;
	}
	
	private static void removeSRF_PlaneHeader(File inputSRF, File outputSRF) throws IOException {
		List<String> srf = Files.readLines(inputSRF, Charset.defaultCharset());
		List<String> modSRF = new ArrayList<>();
		boolean isPlane = false;
		for (String line : srf) {
			if (line.startsWith("PLANE"))
				isPlane = true;
			else if (line.startsWith("POINTS"))
				isPlane = false;
			if (!isPlane)
				modSRF.add(line);
		}
		System.out.println("Stripped out "+(srf.size() - modSRF.size())+" plane header lines");
		FileWriter fw = new FileWriter(outputSRF);
		
		for (String line : modSRF)
			fw.write(line+"\n");
		
		fw.close();
	}
	
	private static void move(File from, File to) throws IOException {
		if (from.isDirectory()) {
			Preconditions.checkState((to.exists() && to.isDirectory()) || to.mkdir());
			for (File file : from.listFiles())
				move(file, new File(to, file.getName()));
		} else {
			Files.move(from, to);
		}
	}
	
	private static List<String> streamToLines(InputStream stream) throws IOException {
		List<String> lines = new ArrayList<>();
		
		BufferedReader b = new BufferedReader(new InputStreamReader(stream));
		String line;
		while ((line = b.readLine()) != null)
			lines.add(line);
		
		return lines;
	}
	
	private static void writeBBP_XML_File(File xmlFile, File siteFile, List<BBP_Module> modules) throws IOException {
		Document doc = XMLUtils.createDocumentWithRoot("BBP_Run_Specification");
		Element root = doc.getRootElement();
		Element runEl = root.addElement("Scenario_Run");
		runEl.addAttribute("input_station_file", siteFile.getAbsolutePath());
		runEl.addAttribute("version", BBP_Module.VERSION);
		Element modsEl = root.addElement("BBP_Modules");
		for (BBP_Module module : modules)
			module.toXMLMetadata(modsEl);
		XMLUtils.writeDocumentToFile(xmlFile, doc);
	}
	
	private static final DecimalFormat expDF = new DecimalFormat("0.000000E00");
	
	private static void velToAccel(File bbpVelFile, File bbpAccelFile) throws IOException {
		List<Double> times = new ArrayList<>();
		List<double[]> velocities = new ArrayList<>();
		
		for (String line : Files.readLines(bbpVelFile, Charset.defaultCharset())) {
			if (line.startsWith("#") || line.isEmpty())
				continue;
			String[] split = line.split("\\s+");
			double t = Double.parseDouble(split[0]);
			double[] vals = new double[split.length-1];
			for (int i=0; i<vals.length; i++)
				vals[i] = Double.parseDouble(split[i+1]);
			if (!times.isEmpty())
				Preconditions.checkState(t > times.get(times.size()-1), "times not in order");
			times.add(t);
			if (!velocities.isEmpty())
				Preconditions.checkState(velocities.get(0).length == vals.length, "Lengths inconsistent");
			velocities.add(vals);
		}
		
		Preconditions.checkState(times.size() > 1, "Only one time?");
		
		FileWriter fw = new FileWriter(bbpAccelFile);
		fw.write("#    time(sec)      N-S(cm/s/s)      E-W(cm/s/s)      U-D(cm/s/s)\n");
		
		for (int i=0; i<velocities.size(); i++) {
			double t = times.get(i);
			double[] vels = velocities.get(i);
			
			double dt;
			double[] prevVels;
			if (i == 0) {
				dt = times.get(i+1) - t;
				prevVels = new double[vels.length];
			} else {
				dt = t - times.get(i-1);
				prevVels = velocities.get(i-1);
			}
			String line = bbpNumFormat(t);
			for (int j=0; j<vels.length; j++) {
				double acc = (vels[j] - prevVels[j])/dt;
				line += "\t"+bbpNumFormat(acc);
			}
			fw.write(line.replaceAll("E", "e")+"\n");
		}
		
		fw.close();
	}
	
	private static String bbpNumFormat(double val) {
		String str = expDF.format(val).replaceAll("E", "e");
		if (!str.contains("e-"))
			str = str.replaceAll("e", "e+");
		return str;
	}
	
//	public static void writeBBP_InputsFile(File inputsFile, int vmIndex, int methodIndex, File srcFile, File srfFile, File sitesFile)
//			throws IOException {
//		FileWriter fw = new FileWriter(inputsFile);
//		
//		fw.write("n\n"); // not validation run
//		fw.write(vmIndex+"\n"); // velocity model
//		fw.write(methodIndex+"\n"); // method
//		fw.write("2\n"); // full path to SRC
//		fw.write(srcFile.getAbsolutePath()+"\n");
//		if (srfFile == null) {
//			// run rupture generator
//			fw.write("y\n");
//		} else {
//			// use SRF
//			fw.write("n\n");
//			fw.write("2\n"); // full path to SRF
//			fw.write(srfFile.getAbsolutePath()+"\n");
//		}
//		fw.write("2\n"); // full path to sites
//		fw.write(sitesFile.getAbsolutePath()+"\n");
//		fw.write("n\n"); // no site response
//		fw.write("y\n"); // generate vel seismograms
//		fw.write("y\n"); // generate accel seismograms
//		
//		fw.close();
//	}

	public static void main(String[] args) throws IOException {
		File eventDir = new File("/data/kevin/simulators/catalogs/rundir2194_long/event_srfs");
		File srcFile = new File(eventDir, "event_136704.src");
//		File srfFile = new File(eventDir, "event_136704_0.05s_ADJ_VEL.srf");
		File srfFile = null;
		File sitesFile = new File("/home/kevin/bbp/bbp_data/run/stations_cs_sites.stl");
		File outputDir = new File("/tmp/bbp_test6");
		
//		runBBP(2, 1, srcFile, null, srfFile, sitesFile, outputDir);
//		removeSRF_PlaneHeader(new File("/tmp/test.srf"));
		
//		BBP_Wrapper wrapper = new BBP_Wrapper(2, 1, srcFile, null, null, sitesFile, outputDir);
		BBP_Wrapper wrapper = new BBP_Wrapper(VelocityModel.LA_BASIN, Method.GP, srcFile, null, srfFile, sitesFile, outputDir);
		wrapper.maxRetries = 1;
		wrapper.doHF = true;
		wrapper.run();
		
//		velToAccel(new File("/data/kevin/bbp/bbp_data/outdata/1507151299716/1507151299716.SBSM.vel.bbp"),
//				new File("/data/kevin/bbp/bbp_data/outdata/1507151299716/1507151299716.SBSM.acc.bbp.2"));
	}

}
