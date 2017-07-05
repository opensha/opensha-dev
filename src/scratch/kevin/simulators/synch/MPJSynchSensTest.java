package scratch.kevin.simulators.synch;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;

import mpi.MPI;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.hpc.mpj.taskDispatch.MPJTaskCalculator;
import org.opensha.commons.util.ClassUtils;
import org.opensha.commons.util.threads.Task;
import org.opensha.commons.util.threads.ThreadedTaskComputer;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.ElementMagRangeDescription;
import org.opensha.sha.simulators.iden.LogicalAndRupIden;
import org.opensha.sha.simulators.iden.LogicalOrRupIden;
import org.opensha.sha.simulators.iden.RuptureIdentifier;

import scratch.kevin.markov.EmpiricalMarkovChain;
import scratch.kevin.simulators.MarkovChainBuilder;
import scratch.kevin.simulators.PeriodicityPlotter;
import scratch.kevin.simulators.SimAnalysisCatLoader;
import scratch.kevin.simulators.catBuild.RandomCatalogBuilder;
import scratch.kevin.simulators.dists.RandomDistType;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;

public class MPJSynchSensTest extends MPJTaskCalculator {
	
	private Map<Double, double[][]> gBarsMap;
	private Map<Double, List<double[][]>> gBarsTrialsMap;
	
	private int numTrials;
	private List<? extends SimulatorEvent> events;
	private List<RuptureIdentifier> rupIdens;
	
	private RandomDistType dist = RandomDistType.ACTUAL;
	private double distSpacing = 10d;
	
	private int nDims;
	
//	private static int[] lags = SynchParamCalculator.rangeInclusive(-30, 30);
	
	private File outputDir;
	
	private SensTestType sensType;
	
	private static final int maxTrialsPerTask = 20;
	
	public enum SensTestType {
		MAG("Min Mag", doubleRange(6.5, 7.5, 0.1)),
		BIN_SIZE("Bin Size (yr)", doubleRange(1d, 20d, 1)),
		ELEM_LOC("Elem. Loc.", doubleRange(0d, 19d, 1)),
		ELEM_VOL("Elem. VOL.", doubleRange(-50d, 50d, 5d)),
		START_SHIFT("Start Time Shift (yr)", doubleRange(0d, 10d, 0.5)),
		// decimal value is fraction of catalog, integer value is sequence number
		CAT_FRACT("Cat. Fract.", 1d, 0.5, 1.5, 0.25, 1.25, 2.25, 3.25,
				0.1, 1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1, 9.1);
		
		private double[] values;
		private String displayName;
		
		private SensTestType(String displayName, double... values) {
			// make sure there are no duplicates
			HashSet<Double> set = new HashSet<Double>();
			for (double value : values)
				set.add(value);
			Preconditions.checkState(set.size() == values.length);
			
			this.displayName = displayName;
			this.values = values;
		}
		
		public String getDisplayName() {
			return displayName;
		}
	}
	
	private static class ComputeTask {
		double sensValue;
		boolean doOrig = false;
		int numTrials;
		
		public ComputeTask(double sensValue,
				boolean doOrig, int numTrials) {
			this.sensValue = sensValue;
			this.doOrig = doOrig;
			this.numTrials = numTrials;
		}
	}
	
	private List<ComputeTask> tasks;
	
	/**
	 * Inclusive!
	 * @param min
	 * @param max
	 * @param delta
	 * @return
	 */
	private static double[] doubleRange(double min, double max, double delta) {
		List<Double> vals = Lists.newArrayList();
		for (double x=min; x<=max; x+=delta)
			vals.add(x);
		return Doubles.toArray(vals);
	}
	
	/**
	 * Inclusive!
	 * @param min
	 * @param max
	 * @return
	 */
	private static int[] intRange(int min, int max) {
		Preconditions.checkState(max > min);
		int[] vals = new int[max - min + 1];
		int index = 0;
		for (int i=min; i<=max; i++)
			vals[index++] = i;
		return vals;
	}
	
	private List<List<RuptureIdentifier>> elemNeighborIdens;
	
	private Map<Integer, SimulatorElement> elemIDMap;

	public MPJSynchSensTest(CommandLine cmd, File outputDir) throws IOException {
		super(cmd);
		
		this.outputDir = outputDir;
		
		Preconditions.checkArgument(cmd.hasOption("trials"));
		numTrials = Integer.parseInt(cmd.getOptionValue("trials"));
		
		sensType = SensTestType.valueOf(cmd.getOptionValue("type"));
		
		double minMag = 7;
		double maxMag = 10d;
		
		int[] include_elems = {
				ElementMagRangeDescription.SAF_CHOLAME_ELEMENT_ID,
				ElementMagRangeDescription.SAF_CARRIZO_ELEMENT_ID,
				ElementMagRangeDescription.GARLOCK_WEST_ELEMENT_ID,
				ElementMagRangeDescription.SAF_MOJAVE_ELEMENT_ID,
				ElementMagRangeDescription.SAF_COACHELLA_ELEMENT_ID,
				ElementMagRangeDescription.SAN_JACINTO__ELEMENT_ID
				};
		
		rupIdens = Lists.newArrayList();
		List<Color> colors = Lists.newArrayList();
		
		SimAnalysisCatLoader.loadElemMagIdens(include_elems, rupIdens, colors, minMag, maxMag);
		
		if (sensType == SensTestType.ELEM_LOC || sensType == SensTestType.ELEM_VOL) {
			elemIDMap = Maps.newHashMap();
			List<SimulatorElement> elems = SimAnalysisCatLoader.loadGeomOnly();
			for (SimulatorElement elem : elems)
				elemIDMap.put(elem.getID(), elem);
		}
		
		List<RuptureIdentifier> loadIdens = rupIdens;
		if (sensType == SensTestType.ELEM_LOC) {
			// need to create and load events for neighbors as well
			elemNeighborIdens = Lists.newArrayList();
			
			loadIdens = Lists.newArrayList();
			for (RuptureIdentifier rupIden : rupIdens) {
				List<RuptureIdentifier> neighbors = getNeighboringIdens(rupIden, 15d, 2d, 5, elemIDMap);
				loadIdens.addAll(neighbors);
				elemNeighborIdens.add(neighbors);
			}
		} else if (sensType == SensTestType.MAG) {
			loadIdens = Lists.newArrayList();
			for (RuptureIdentifier rupIden : rupIdens) {
				loadIdens.add(getAsDiffMinMag(rupIden, sensType.values[0]));
			}
		} else if (sensType == SensTestType.ELEM_VOL) {
			double max = Math.max(Math.abs(StatUtils.max(sensType.values)), Math.abs(StatUtils.min(sensType.values)));
			System.out.println("Max dist: "+max);
			// need to create and load events for neighbors as well
			elemNeighborIdens = Lists.newArrayList();

			loadIdens = Lists.newArrayList();
			for (RuptureIdentifier rupIden : rupIdens) {
				List<RuptureIdentifier> neighbors = getNeighboringIdens(rupIden, max, 2d, 5, elemIDMap);
				loadIdens.addAll(neighbors);
				elemNeighborIdens.add(neighbors);
			}
		}
		
		events = new SimAnalysisCatLoader(true, loadIdens, false).getEvents();
		
		nDims = rupIdens.size();
		
		tasks = Lists.newArrayList();
		
		gBarsMap = Maps.newHashMap();
		gBarsTrialsMap = Maps.newHashMap();
		
		for (double value : sensType.values) {
			
			if (numTrials == 0) {
				tasks.add(new ComputeTask(value, true, 0));
				continue;
			}
			
			int trialsAdded = 0;
			while (trialsAdded < numTrials) {
				int trialsToAdd = maxTrialsPerTask;
				if (trialsAdded + trialsToAdd > numTrials)
					trialsToAdd = numTrials - trialsAdded;
				
				tasks.add(new ComputeTask(value, trialsAdded == 0, trialsToAdd));
				
				trialsAdded += trialsToAdd;
			}
			Preconditions.checkState(trialsAdded == numTrials, "expected "+numTrials+", got "+trialsAdded);
			
			gBarsTrialsMap.put(value, new ArrayList<double[][]>());
		}
	}
	
	private List<RuptureIdentifier> getNeighboringIdens(RuptureIdentifier rupIden, double maxHorzDist, double maxVertDist,
			int num, Map<Integer, SimulatorElement> elemIDMap) {
		Preconditions.checkState(rupIden instanceof ElementMagRangeDescription);
		ElementMagRangeDescription elemIden = (ElementMagRangeDescription)rupIden;
		Preconditions.checkState(elemIden.getElementIDs().size() == 1);
		
		int origID = elemIden.getElementIDs().get(0);
		// find that elem
		SimulatorElement elem = elemIDMap.get(origID);
		Preconditions.checkNotNull(elem);
		
		int sectID = elem.getSectionID();
		
		List<SimulatorElement> matchingElems = Lists.newArrayList();
		
		for (SimulatorElement testElem : elemIDMap.values()) {
			if (testElem.getSectionID() == sectID) {
				// same section
				double vertDist = LocationUtils.vertDistance(elem.getCenterLocation(), testElem.getCenterLocation());
				if (vertDist > maxVertDist)
					continue;
				double horzDist = LocationUtils.horzDistanceFast(elem.getCenterLocation(), testElem.getCenterLocation());
				if (horzDist > maxHorzDist)
					continue;
				matchingElems.add(testElem);
			}
		}
		Preconditions.checkState(matchingElems.size() > 1);
		
		List<RuptureIdentifier> idens = Lists.newArrayList();
		
		if (num <= 0) {
			for (SimulatorElement matchingElem : matchingElems) {
				idens.add(new ElementMagRangeDescription(rupIden.getName(), matchingElem.getID(),
						elemIden.getMinMag(), elemIden.getMaxMag()));
			}
		} else {
			Collections.shuffle(matchingElems);
			
			int index = 0;
			for (int i=0; i<num; i++,index++) {
				if (index == matchingElems.size()) {
					index = 0;
					System.out.println("WARNING: not enough matching elems, resetting to the start");
				}
				idens.add(new ElementMagRangeDescription(rupIden.getName(), matchingElems.get(index).getID(),
						elemIden.getMinMag(), elemIden.getMaxMag()));
			}
		}
		
		return idens;
	}
	
	private ElementMagRangeDescription getAsDiffMinMag(RuptureIdentifier rupIden, double minMag) {
		Preconditions.checkState(rupIden instanceof ElementMagRangeDescription);
		ElementMagRangeDescription elemIden = (ElementMagRangeDescription)rupIden;
		return new ElementMagRangeDescription(elemIden.getName(), elemIden.getElementIDs(), minMag, elemIden.getMaxMag());
	}

	@Override
	protected int getNumTasks() {
		return tasks.size();
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		debug("Calculating batch of length "+batch.length+" ("+getMemoryDebug()+")");
		for (int item : batch) {
			ComputeTask task = tasks.get(item);
			
			int numThreads = getNumThreads();
			
			if (sensType == SensTestType.BIN_SIZE) {
				if (task.sensValue <= 2d)
					numThreads = 1;
				else if (task.sensValue <= 4d)
					numThreads = 2;
				else if (task.sensValue <= 6d)
					numThreads /= 4;
				else if (task.sensValue < 10d)
					numThreads /= 2;
				if (numThreads < 1)
					numThreads = 1;
			} else if (sensType == SensTestType.MAG) {
				if (task.sensValue < 7d)
					numThreads /= 2;
			}
			if (numThreads > getNumThreads())
				numThreads = getNumThreads();
			
			List<SynchCalc> tasks = Lists.newArrayList();
			
			if (task.doOrig) {
				double[][] gBars = new double[nDims][nDims];
				
				tasks.add(new SynchCalc(gBars, task.sensValue, false));
				
				gBarsMap.put(task.sensValue, gBars);
			}
			
			for (int trial=0; trial<task.numTrials; trial++) {
				double[][] gBars = new double[nDims][nDims];

				tasks.add(new SynchCalc(gBars, task.sensValue, true));

				gBarsTrialsMap.get(task.sensValue).add(gBars);
			}
			
			ThreadedTaskComputer comp = new ThreadedTaskComputer(tasks);
			debug("Computing item "+item+" (sensVal="+task.sensValue+") with "+tasks.size()+" tasks, "
					+numThreads+" threads ("+getMemoryDebug()+")");
			comp.computeThreaded(numThreads);
			debug("Done with item "+item+" ("+getMemoryDebug()+")");
		}
		debug("Done with batch  ("+getMemoryDebug()+")");
	}
	
	private EmpiricalMarkovChain buildChain(double sensVal, boolean random) {
		List<RuptureIdentifier> rupIdens = this.rupIdens;
		double distSpacing = this.distSpacing;
		List<SimulatorEvent> events = Lists.newArrayList(this.events);
		double startTimeShift = 0d;
		
		switch (sensType) {
		case BIN_SIZE:
			distSpacing = sensVal;
			break;
		case MAG:
			rupIdens = Lists.newArrayList();
			for (RuptureIdentifier rupIden : this.rupIdens)
				rupIdens.add(getAsDiffMinMag(rupIden, sensVal));
			break;
		case ELEM_LOC:
			Random r = new Random();
			rupIdens = Lists.newArrayList();
			for (int i=0; i<this.rupIdens.size(); i++) {
				List<RuptureIdentifier> neighbors = elemNeighborIdens.get(i);
				rupIdens.add(neighbors.get(r.nextInt(neighbors.size())));
			}
			break;
		case ELEM_VOL:
			boolean and = sensVal >= 0d;
			double maxDist = Math.abs(sensVal);
			rupIdens = Lists.newArrayList();
			for (int i=0; i<this.rupIdens.size(); i++) {
				RuptureIdentifier origIden = this.rupIdens.get(i);
				SimulatorElement elem = elemIDMap.get(((ElementMagRangeDescription)origIden).getElementIDs().get(0));
				List<RuptureIdentifier> neighbors = elemNeighborIdens.get(i);
				
				List<RuptureIdentifier> combined = Lists.newArrayList(origIden);
				for (RuptureIdentifier neighbor : neighbors) {
					SimulatorElement neighborElem = elemIDMap.get(((ElementMagRangeDescription)neighbor).getElementIDs().get(0));
					double hDist = LocationUtils.horzDistanceFast(elem.getCenterLocation(), neighborElem.getCenterLocation());
					if (hDist <= maxDist)
						combined.add(neighbor);
				}
				debug("Combined has "+combined.size()+" idens for "+origIden.getName());
				if (and)
					rupIdens.add(new LogicalAndRupIden(combined));
				else
					rupIdens.add(new LogicalOrRupIden(combined));
			}
			break;
		case CAT_FRACT:
			if (sensVal != 1d) {
				double origStart = events.get(0).getTime();
				double origEnd = events.get(events.size()-1).getTime();
				double origLen = origEnd - origStart;
				
				double fract = sensVal - Math.floor(sensVal);
				double position = Math.floor(sensVal);
				
				double newLen = fract*origLen;
				double newStart = origStart + newLen*position;
				double newEnd = newStart + newLen;
				
				events = Lists.newArrayList();
				for (SimulatorEvent event : this.events) {
					if (event.getTime() < newStart)
						continue;
					if (event.getTime() > newEnd)
						break;
					events.add(event);
				}
				Preconditions.checkState(events.size() > 1, "Too few events! newStart="+newStart+", newEnd="+newEnd
						+", sensVal="+sensVal+", fract="+fract+", pos="+position);
			}
			break;
		case START_SHIFT:
			startTimeShift = sensVal;
			break;

		default:
			throw new IllegalStateException("Unknown sensType: "+sensType);
		}
		
		if (random)
			events = RandomCatalogBuilder.getRandomResampledCatalog(events, rupIdens, dist, true, 1);
		
		return MarkovChainBuilder.build(distSpacing, events, rupIdens, startTimeShift);
	}
	
	private String getMemoryDebug() {
		System.gc();
		Runtime rt = Runtime.getRuntime();
		long totalMB = rt.totalMemory() / 1024 / 1024;
		long freeMB = rt.freeMemory() / 1024 / 1024;
		long usedMB = totalMB - freeMB;
		return "mem t/u/f: "+totalMB+"/"+usedMB+"/"+freeMB;
	}
	
	class SynchCalc implements Task {
		
		private double[][] gBars;
		private double sensValue;
		private boolean random;

		public SynchCalc(double[][] gBars, double sensValue, boolean random) {
			super();
			this.gBars = gBars;
			this.sensValue = sensValue;
			this.random = random;
		}

		@Override
		public void compute() {
			EmpiricalMarkovChain chain = buildChain(sensValue, random);
			
			for (int m=0; m<nDims; m++) {
				for (int n=m; n<nDims; n++) {
					double gBar = SynchParamCalculator.calcGBar(chain, m, n, 0);
					
					gBars[m][n] = gBar;
					gBars[n][m] = gBar;
				}
			}
			
			if (numTrials <= 100)
				debug("Done with task ("+getMemoryDebug()+")");
		}
		
	}

	@SuppressWarnings({ "rawtypes", "unchecked" })
	@Override
	protected void doFinalAssembly() throws Exception {
		Map[] sendbuf = new Map[] {gBarsMap};
		Map[] recvbuf1 = null; // gBars
		Map[] recvbuf2 = null; // gBarsTrials
		
		if (rank == 0) {
			recvbuf1 = new Map[size];
			recvbuf2 = new Map[size];
		}
		
		debug("Gather 1  ("+getMemoryDebug()+")");
		MPI.COMM_WORLD.Gather(sendbuf, 0, 1, MPI.OBJECT, recvbuf1, 0, 1, MPI.OBJECT, 0);
		if (numTrials > 0) {
			sendbuf = new Map[] {gBarsTrialsMap};
			debug("Gather 2  ("+getMemoryDebug()+")");
			MPI.COMM_WORLD.Gather(sendbuf, 0, 1, MPI.OBJECT, recvbuf2, 0, 1, MPI.OBJECT, 0);
		}
		
		if (rank == 0) {
			debug("Assembling/writing results ("+getMemoryDebug()+")");
			for (int i=1; i<size; i++) {
				Map<Double, double[][]> o1 = recvbuf1[i];
				gBarsMap.putAll(o1);
				
				Map<Double, List<double[][]>> o2 = recvbuf2[i];
				for (Double value : o2.keySet()) {
					List<double[][]> vals = o2.get(value);
					gBarsTrialsMap.get(value).addAll(vals);
				}
			}
			
			Preconditions.checkState(sensType.values.length == gBarsMap.size());
			for (Double value : gBarsTrialsMap.keySet())
				Preconditions.checkState(numTrials == gBarsTrialsMap.get(value).size(),
				"Trials of wrong length for sensVal="+value+". Expected "+numTrials
				+", got "+gBarsTrialsMap.get(value).size());
			
			File writeDir = new File(outputDir, SynchParamCalculator.getDirName());
			if (!writeDir.exists())
				writeDir.mkdir();
			
			File sensDir = new File(writeDir, "sens_tests");
			if (!sensDir.exists())
				sensDir.mkdir();
			
			File subSensDir = new File(sensDir, sensType.name());
			if (!subSensDir.exists())
				subSensDir.mkdir();
			
			for (int m=0; m<nDims; m++) {
				String name1 = rupIdens.get(m).getName();
				for (int n=m+1; n<nDims; n++) {
					CSVFile<String> csv = new CSVFile<String>(true);
					String name2 = rupIdens.get(n).getName();
					File csvFile = new File(subSensDir, PeriodicityPlotter.getFileSafeString(name1)
							+"_"+PeriodicityPlotter.getFileSafeString(name2)+"_"+numTrials+"trials.csv");
					
					List<String> line = Lists.newArrayList();
					line.add("Value");
					for (double value : sensType.values)
						line.add((float)value+"");
					csv.addLine(line);
					
					line = Lists.newArrayList();
					line.add("Gbar");
					for (double value : sensType.values)
						line.add(gBarsMap.get(value)[m][n]+"");
					csv.addLine(line);
					
					// now add trials
					for (int i=0; i<numTrials; i++) {
						line = Lists.newArrayList();
						line.add("Trial "+i);
						for (double value : sensType.values) {
							double[][] trialGBars = gBarsTrialsMap.get(value).get(i);
							line.add(trialGBars[m][n]+"");
						}
						csv.addLine(line);
					}
					
					csv.writeToFile(csvFile);
				}
			}
			debug("Done writing results ("+getMemoryDebug()+")");
		}
	}
	
	public static Options createOptions() {
		Options ops = MPJTaskCalculator.createOptions();
		
		Option trialsOption = new Option("trials", "num-trials", true,
				"Number of random trials");
		trialsOption.setRequired(true);
		ops.addOption(trialsOption);
		
		Option typeOption = new Option("type", "test-type", true, "Sensitivity test type");
		typeOption.setRequired(true);
		ops.addOption(typeOption);
		
		return ops;
	}
	
	public static void main(String[] args) {
		args = MPJTaskCalculator.initMPJ(args);
		
		try {
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJSynchSensTest.class);
			
			args = cmd.getArgs();
			
			if (args.length != 1) {
				System.err.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(MPJSynchSensTest.class)
						+" [options] <output-dir>");
				abortAndExit(2);
			}
			
			File outputDir = new File(args[0]);
			if (!outputDir.exists() && MPI.COMM_WORLD.Rank() == 0)
				Preconditions.checkState(outputDir.mkdir());
			
			MPJSynchSensTest calc = new MPJSynchSensTest(cmd, outputDir);
			calc.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}

