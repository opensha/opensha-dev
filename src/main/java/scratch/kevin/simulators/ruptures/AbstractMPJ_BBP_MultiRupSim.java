package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FileUtils;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.LogicalOrRupIden;
import org.opensha.sha.simulators.iden.RegionIden;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator.SRFInterpolationMode;
import org.opensha.sha.simulators.srf.SRF_PointData;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;
import com.google.common.io.Files;
import com.google.common.primitives.Ints;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import mpi.MPI;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.kevin.bbp.BBP_Module.Method;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;
import scratch.kevin.bbp.BBP_Wrapper;
import scratch.kevin.bbp.MPJ_BBP_Utils;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Loader;

public abstract class AbstractMPJ_BBP_MultiRupSim extends AbstractMPJ_BBP_Sim {
	
	private double dt;
	private SRFInterpolationMode interp;
	
	private double timeScalarFactor = 1d;
	private boolean velocityScale = false;
	
	protected final RSQSimCatalog catalog;

	public AbstractMPJ_BBP_MultiRupSim(CommandLine cmd) throws IOException {
		super(cmd);
		
		File catalogDir = new File(cmd.getOptionValue("catalog-dir"));
		if (!catalogDir.exists()) {
			String catalogDirName = catalogDir.getName();
			if (rank == 0)
				debug("Catalog dir doesn't exist, searching for: "+catalogDirName);
			catalogDir = RSQSimCatalog.locateCatalog(catalogDirName);
			Preconditions.checkNotNull(catalogDir, "Couldn't locate catalog: "+catalogDirName);
			if (rank == 0)
				debug("Located catalog: "+catalogDir.getAbsolutePath());
		}
		Preconditions.checkState(catalogDir.exists());
		
		dt = Double.parseDouble(cmd.getOptionValue("time-step"));
		interp = SRFInterpolationMode.valueOf(cmd.getOptionValue("srf-interp"));
		
		if (nodeScratch != null) {
			// copy catalog data over to node scratch
			catalogDir = copyCatalogDir(catalogDir, nodeScratch);
		} else if (sharedScratch != null) {
			// copy catalog data over to shared scratch
			if (rank == 0)
				copyCatalogDir(catalogDir, sharedScratch);
			resultsScratchDir = new File(sharedScratch, "results_tmp");
			if (rank == 0)
				MPJ_BBP_Utils.waitOnDir(resultsScratchDir, 10, 2000);
			MPI.COMM_WORLD.Barrier();
			catalogDir = new File(sharedScratch, catalogDir.getName());
		}
		
		// load the catalog
		if (cmd.hasOption("utm-zone") || cmd.hasOption("utm-band")) {
			Preconditions.checkState(cmd.hasOption("utm-zone") && cmd.hasOption("utm-band"),
					"Must supply both UTM zone and band if either is supplied.");
			int zone = Integer.parseInt(cmd.getOptionValue("utm-zone"));
			char band = cmd.getOptionValue("utm-band").trim().charAt(0);
			catalog = new RSQSimCatalog(catalogDir, catalogDir.getName(),
					null, null, null, null, null, zone, band);
		} else {
			catalog = new RSQSimCatalog(catalogDir, catalogDir.getName(),
					null, null, null, FaultModels.FM3_1, DeformationModels.GEOLOGIC); // TODO make option
		}

		if (cmd.hasOption("time-scalar"))
			timeScalarFactor = Double.parseDouble(cmd.getOptionValue("time-scalar"));
		if (cmd.hasOption("velocity-scale"))
			velocityScale = true;
	}
	
	protected List<SRF_PointData> getSRFPoints(int index) throws IOException {
		RSQSimEvent event = eventForIndex(index);
		RSQSimEventSlipTimeFunc func = catalog.getSlipTimeFunc(event);
		
		if (timeScalarFactor != 1d)
			func = func.getTimeScaledFunc(timeScalarFactor, velocityScale);
		
		// write SRF
		debug("bulding/writing SRF for "+index+", event "+event.getID());
		return RSQSimSRFGenerator.buildSRF(func, event.getAllElements(), dt, interp);
	}
	
	protected BBP_SourceFile getBBPSource(int index) {
		RSQSimEvent event = eventForIndex(index);
		debug("bulding/writing SRC for "+index);
		BBP_PlanarSurface bbpSurface;
		if (numRG > 0 && RSQSimBBP_Config.U3_SURFACES)
			bbpSurface = RSQSimBBP_Config.planarEquivalentU3Surface(catalog, event);
		else
			bbpSurface = RSQSimBBP_Config.estimateBBP_PlanarSurface(event);
		return RSQSimBBP_Config.buildBBP_Source(event, bbpSurface, 12345);
	}
	
	private File copyCatalogDir(File catDir, File scratchDir) throws IOException {
		File destDir = new File(scratchDir, catDir.getName());
		if (!destDir.exists()) {
			destDir.mkdir();
			MPJ_BBP_Utils.waitOnDir(destDir, 10, 2000);
		}
		List<File> filesToCopy = new ArrayList<>();
		for (File f : catDir.listFiles()) {
			String name = f.getName().toLowerCase();
			if (name.endsWith("list"))
				filesToCopy.add(f);
			else if (name.startsWith("trans.") && name.endsWith(".out"))
				filesToCopy.add(f);
			else if (name.startsWith("transv.") && name.endsWith(".out"))
				filesToCopy.add(f);
			else if (name.endsWith(".in"))
				filesToCopy.add(f);
			else if (name.endsWith(".flt"))
				filesToCopy.add(f);
		}
		if (filesToCopy.size() < 7) {
			String error = "Need at least 7 files: 4 list, trans, input, geom. Have "+filesToCopy.size()+":";
			for (File f : filesToCopy)
				error += "\n\t"+f.getAbsolutePath();
			throw new IllegalStateException(error);
		}
		
		debug("copying "+filesToCopy.size()+" catalog files to "+destDir.getAbsolutePath());
		for (File f : filesToCopy) {
			File destFile = new File(destDir, f.getName());
			if (destFile.exists() && destFile.length() == f.length())
				// skip copy, already exists
				continue;
			Files.copy(f, destFile);
		}
		
		return destDir;
	}
	
	protected abstract RSQSimEvent eventForIndex(int index);
	
	public static Options addCommonOptions(Options ops) {
		AbstractMPJ_BBP_Sim.addCommonOptions(ops);
		Option dt = new Option("dt", "time-step", true, "SRF time step");
		dt.setRequired(true);
		ops.addOption(dt);
		
		Option interp = new Option("interp", "srf-interp", true, "SRF interpolation mode");
		interp.setRequired(true);
		ops.addOption(interp);
		
		Option catalogDir = new Option("cdir", "catalog-dir", true, "RSQSim catalog dir");
		catalogDir.setRequired(true);
		ops.addOption(catalogDir);
		
		Option timeScalar = new Option("ts", "time-scalar", true,
				"Time scalar for slip/time functions. 1.5 will increase rupture propagation velocities by 50%");
		timeScalar.setRequired(false);
		ops.addOption(timeScalar);
		
		Option velScale = new Option("vs", "velocity-scale", false,
				"Works with time scalar. If supplied, velocities will also be scaled (each slip event is faster and shorter). "
				+ "Otherwise slip events are sooner, but the same speed/duration.");
		velScale.setRequired(false);
		ops.addOption(velScale);
		
		ops.addOption("uz", "utm-zone", true, "UTM zone (integer). Optional, but must also supply --utm-band if used.");
		ops.addOption("ub", "utm-band", true, "UTM band (character). Optional, but must also supply --utm-zone if used.");
		
		return ops;
	}

}
