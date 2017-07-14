package scratch.kevin.simulators;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Region;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;
import org.opensha.commons.util.XMLUtils;
import org.opensha.sha.calc.hazardMap.components.BinaryCurveArchiver;
import org.opensha.sha.calc.hazardMap.components.CalculationInputsXMLFile;
import org.opensha.sha.calc.hazardMap.components.CalculationSettings;
import org.opensha.sha.calc.hazardMap.components.CurveResultsArchiver;
import org.opensha.sha.calc.hazardMap.mpj.MPJHazardCurveDriver;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.AbstractRuptureIdentifier;
import org.opensha.sha.simulators.iden.MagRangeRuptureIdentifier;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.parsers.RSQSimFileReader;
import org.opensha.sha.simulators.utils.RSQSimUtils;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Files;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.utils.FaultSystemIO;

public class HazardMapCompareScriptGen {

	public static void main(String[] args) throws IOException, DocumentException {
		String runName = "2017_07_14-bruce2142";
		File catalogDir = new File("/home/kevin/Simulators/catalogs/bruce/rundir2142");
		File geomFile = new File(catalogDir, "zfault_Deepen.in");
		
//		String runName = "2017_07_14-bruce2194";
//		File catalogDir = new File("/home/kevin/Simulators/catalogs/bruce/rundir2194");
//		File geomFile = new File(catalogDir, "zfault_Deepen.in");
		
//		String runName = "2017_06_28-jacqui_slipWeakening_calibrated_2-vs-ucerf3";
//		File catalogDir = new File("/home/kevin/Simulators/catalogs/baseCatalog_slipWeakening_calibrated_2");
//		String runName = "2017_06_28-jacqui_shortTestCatalog-vs-ucerf3";
//		File catalogDir = new File("/home/kevin/Simulators/catalogs/shortTestCatalog");
//		File geomFile = new File(catalogDir, "UCERF3.D3.1.1km.tri.2.flt");
		
		boolean doUCERF3 = true;
		boolean u3SupraMinMag = true;
		double minMag = 6.5d;
		double minFractForInclusion = 0.1;
		if (u3SupraMinMag)
			runName += "-matchU3supra";
		else
			runName += "-m"+(float)minMag;
		if (minFractForInclusion > 0)
			runName += "-sectArea"+(float)minFractForInclusion;
		FaultModels fm = FaultModels.FM3_1;
		DeformationModels dm = DeformationModels.GEOLOGIC;
		
		// load/process U3 catalog
		FaultSystemSolution u3Sol = null;
		if (doUCERF3 || u3SupraMinMag) {
			System.out.println("Loading/filtering U3 solution");
			File localU3File = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
					+ "FM3_1_GEOL_MEAN_BRANCH_AVG_SOL.zip");
			u3Sol = FaultSystemIO.loadSol(localU3File);
			if (!u3SupraMinMag) {
				FaultSystemRupSet u3RupSet = u3Sol.getRupSet();
				double[] modRates = new double[u3RupSet.getNumRuptures()];
				for (int r=0; r<modRates.length; r++)
					if (u3RupSet.getMagForRup(r) >= minMag)
						modRates[r] = u3Sol.getRateForRup(r);
				u3Sol = new FaultSystemSolution(u3RupSet, modRates);
			}
		}
		
		// load RSQSim catalog
		System.out.println("Loading RSQSim catalog");
		List<SimulatorElement> geom = RSQSimFileReader.readGeometryFile(geomFile, 11, 'S');
		List<RuptureIdentifier> rupIdens = Lists.newArrayList();
		if (u3SupraMinMag) {
			int minSubSectIndex = Integer.MAX_VALUE;
			for (SimulatorElement elem : geom)
				if (elem.getSectionID() < minSubSectIndex)
					minSubSectIndex = elem.getSectionID();
			rupIdens.add(new U3SupraSeisRuptureIdentifier(u3Sol.getRupSet(), minSubSectIndex));
		} else {
			rupIdens.add(new MagRangeRuptureIdentifier(minMag, 11d));
		}
		List<RSQSimEvent> events = RSQSimFileReader.readEventsFile(catalogDir, geom, rupIdens);
		double length = events.get(events.size()-1).getTimeInYears() - events.get(0).getTimeInYears();
		System.out.println("Duration: "+(float)length+" years");
		System.out.println("Building RSQSim FSS");
		FaultSystemSolution rsqsimSol = RSQSimUtils.buildFaultSystemSolution(
				RSQSimUtils.getUCERF3SubSectsForComparison(fm, dm), geom, events, 0d, minFractForInclusion);
		
		File localMainDir = new File("/home/kevin/Simulators/hazard");
		File remoteMainDir = new File("/home/scec-02/kmilner/simulators/hazard");
		
		Region region = new CaliforniaRegions.RELM_TESTING();
		double spacing = 0.02;
		String imt = PGA_Param.NAME;
		double period = 0d;
		ScalarIMR imr = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null);
		double duration = 1d;
		
		File localDir = new File(localMainDir, runName);
		File remoteDir = new File(remoteMainDir, runName);
		System.out.println(localDir.getAbsolutePath()+" => "+remoteDir.getAbsolutePath());
		Preconditions.checkState(localDir.exists() || localDir.mkdir());
		
		int mins = 24*60;
		int nodes = 18;
		int ppn = 20;
		String queue = "scec";
		
		GriddedRegion gridded = new GriddedRegion(region, spacing, null);
		
		imr.setParamDefaults();
		imr.setIntensityMeasure(imt);
		if (period > 0)
			SA_Param.setPeriodInSA_Param(imr.getIntensityMeasure(), period);
		
		List<Site> sites = Lists.newArrayList();
		for (int i=0; i<gridded.getNodeCount(); i++) {
			Site site = new Site(gridded.locationForIndex(i));
			site.addParameterList(imr.getSiteParams());
			sites.add(site);
		}
		
		ArbitrarilyDiscretizedFunc xValues = new IMT_Info().getDefaultHazardCurve(imt);
		double maxSourceDistance = 200;
		
		Map<String, DiscretizedFunc> xValsMap = Maps.newHashMap();
		xValsMap.put("curves", xValues);
		CalculationSettings calcSettings = new CalculationSettings(xValues, maxSourceDistance);
		
		File javaBin = USC_HPCC_ScriptWriter.JAVA_BIN;
		File jarFile = new File(remoteDir, "OpenSHA_complete.jar");
		
		List<File> classpath = Lists.newArrayList();
		classpath.add(jarFile);
		
		MPJExpressShellScriptWriter mpj = new MPJExpressShellScriptWriter(javaBin, 60000, classpath,
				USC_HPCC_ScriptWriter.MPJ_HOME);
		
		List<Map<TectonicRegionType, ScalarIMR>> imrMaps = Lists.newArrayList();
		
		HashMap<TectonicRegionType, ScalarIMR> map = Maps.newHashMap();
		map.put(TectonicRegionType.ACTIVE_SHALLOW, imr);
		imrMaps.add(map);
		
		boolean[] bools;
		if (doUCERF3)
			bools = new boolean[] { true, false };
		else
			bools = new boolean[] { true };
		for (boolean rsqsim : bools) {
			String subDirName;
			if (rsqsim)
				subDirName = "rsqsim";
			else
				subDirName = "ucerf3";
			System.out.println("Processing "+subDirName);
			File localSubDir = new File(localDir, subDirName);
			Preconditions.checkState(localSubDir.exists() || localSubDir.mkdir());
			File remoteSubDir = new File(remoteDir, subDirName);
			
			File localSolFile;
			if (rsqsim) {
				// write solution
				localSolFile = new File(localSubDir, "rsqsim_solution.zip");
				FaultSystemIO.writeSol(rsqsimSol, localSolFile);
			} else {
				localSolFile = new File(localSubDir, "ucerf3_sol_filtered.zip");
				FaultSystemIO.writeSol(u3Sol, localSolFile);
			}
			File remoteSolFile = new File(remoteSubDir, localSolFile.getName());
			
			FaultSystemSolutionERF erf = new FaultSystemSolutionERF();
			erf.setParameter(FaultSystemSolutionERF.FILE_PARAM_NAME, remoteSolFile);
			erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
			erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
			erf.getTimeSpan().setDuration(duration);
			
			File curveDir = new File(remoteSubDir, "curves");
//			CurveResultsArchiver archiver = new AsciiFileCurveArchiver(
//				curveDir.getAbsolutePath()+File.separator, true, false);
			CurveResultsArchiver archiver = new BinaryCurveArchiver(curveDir, sites.size(), xValsMap);
			
			CalculationInputsXMLFile inputs = new CalculationInputsXMLFile(erf, imrMaps, sites, calcSettings, archiver);
			
			File localInputsFile = new File(localSubDir, "inputs.xml");
			File remoteInputsFile = new File(remoteSubDir, "inputs.xml");
			XMLUtils.writeObjectToXMLAsRoot(inputs, localInputsFile);
			
//			String cliArgs = "--max-dispatch 1000 --mult-erfs "+inputsFile.getAbsolutePath();
			String cliArgs = "--max-dispatch 1000 "+remoteInputsFile.getAbsolutePath();
			
			List<String> script = mpj.buildScript(MPJHazardCurveDriver.class.getName(), cliArgs);
			USC_HPCC_ScriptWriter writer = new USC_HPCC_ScriptWriter();
			
			script = writer.buildScript(script, mins, nodes, ppn, queue);
			
			File pbsFile = new File(localSubDir, subDirName+".pbs");
			JavaShellScriptWriter.writeScript(pbsFile, script);
		}
	}
	
	private static class U3SupraSeisRuptureIdentifier extends AbstractRuptureIdentifier {
		
		private double[] subSectMinMags;
		private int subSectOffset;
		
		public U3SupraSeisRuptureIdentifier(FaultSystemRupSet u3RupSet, int subSectOffset) {
			subSectMinMags = new double[u3RupSet.getNumSections()];
			this.subSectOffset = subSectOffset;
			
			for (int s=0; s<subSectMinMags.length; s++) {
				double minMag = Double.POSITIVE_INFINITY;
				for (int r : u3RupSet.getRupturesForSection(s))
					minMag = Math.min(minMag, u3RupSet.getMagForRup(r));
				subSectMinMags[s] = minMag;
			}
		}

		@Override
		public boolean isMatch(SimulatorEvent event) {
			double rupMag = event.getMagnitude();
			for (SimulatorElement elem : event.getAllElements()) {
				int subSectID = elem.getSectionID() - subSectOffset;
				if (rupMag >= subSectMinMags[subSectID])
					return true;
			}
			return false;
		}

		@Override
		public String getName() {
			return "U3 SupraSeismogenic Min Mag filter";
		}
		
	}

}
