package scratch.kevin.simulators;

import java.awt.Color;
import java.awt.Dialog.ModalityType;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Random;

import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JOptionPane;

import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.gui.plot.GraphWidget;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.DataUtils;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.ElementMagRangeDescription;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.parsers.EQSIMv06FileReader;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;

import scratch.kevin.simulators.catBuild.RandomCatalogBuilder;
import scratch.kevin.simulators.dists.CompoundDistReturnPeriodProvider;
import scratch.kevin.simulators.dists.LogNormalDistReturnPeriodProvider;
import scratch.kevin.simulators.dists.RandomDistType;
import scratch.kevin.simulators.dists.RandomReturnPeriodProvider;

import com.google.common.collect.Lists;

public class RIDistPlot {

	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/Simulators");
		File geomFile = new File(dir, "ALLCAL2_1-7-11_Geometry.dat");
		System.out.println("Loading geometry...");
		General_EQSIM_Tools tools = new General_EQSIM_Tools(geomFile);
//		File eventFile = new File(dir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.barall");
		File eventFile = new File(dir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.long.barall");
		System.out.println("Loading events...");
		List<? extends SimulatorEvent> events = EQSIMv06FileReader.readEventsFile(eventFile, tools.getElementsList());
		
		File writeDir = new File(dir, "period_plots");
		
		List<RuptureIdentifier> rupIdens = Lists.newArrayList();
		
		rupIdens.add(new ElementMagRangeDescription("SAF Cholame 7+",
				ElementMagRangeDescription.SAF_CHOLAME_ELEMENT_ID, 7d, 10d));
		
		rupIdens.add(new ElementMagRangeDescription("SAF Carrizo 7+",
				ElementMagRangeDescription.SAF_CARRIZO_ELEMENT_ID, 7d, 10d));
		
		rupIdens.add(new ElementMagRangeDescription("Garlock 7+",
				ElementMagRangeDescription.GARLOCK_WEST_ELEMENT_ID, 7d, 10d));
		
		rupIdens.add(new ElementMagRangeDescription("SAF Mojave 7+",
				ElementMagRangeDescription.SAF_MOJAVE_ELEMENT_ID, 7d, 10d));
		
		rupIdens.add(new ElementMagRangeDescription("SAF Coachella 7+",
				ElementMagRangeDescription.SAF_COACHELLA_ELEMENT_ID, 7d, 10d));
		
		rupIdens.add(new ElementMagRangeDescription("San Jacinto 7+",
				ElementMagRangeDescription.SAN_JACINTO__ELEMENT_ID, 7d, 10d));
		
		List<List<SimulatorEvent>> resampledCatalogs = Lists.newArrayList();
		RandomDistType[] types = RandomDistType.values();
		for (RandomDistType type : types)
			resampledCatalogs.add(RandomCatalogBuilder.getRandomResampledCatalog(events, rupIdens, type, true));
		
		for (int i=0; i<rupIdens.size(); i++) {
			RuptureIdentifier rupIden = rupIdens.get(i);
			String name = rupIden.getName();
			
			File plotFile = new File(writeDir, "ri_dists_"+PeriodicityPlotter.getFileSafeString(name));
			
			System.out.println("Plotting: "+name);
			
			List<Color> colors = GraphWindow.generateDefaultColors();
			
			List<HistogramFunction> funcs = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			
			List<? extends SimulatorEvent> matches = rupIden.getMatches(events);
			funcs.add(getHist(matches, null));
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, colors.remove(0)));
			
			for (int j=0; j<resampledCatalogs.size(); j++) {
				RandomDistType type = types[j];
				matches = rupIden.getMatches(resampledCatalogs.get(j));
				funcs.add(getHist(matches, type));
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, colors.remove(0)));
			}
			
			GraphWindow gw = new GraphWindow(funcs, name+" RI Dists", chars);
			gw.saveAsPDF(plotFile.getAbsolutePath()+".pdf");
			gw.saveAsPNG(plotFile.getAbsolutePath()+".png");
			gw.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		}
		
		CSVFile<String> distParamsCSV = new CSVFile<String>(true);
		distParamsCSV.addLine("Name", "True Mean", "Median", "Variance", "Std Dev", "COV",
				"Trimmed Mean", "Median", "Variance", "Std Dev", "COV", "Scale: ln(mean)", "Shape");
		for (int i=0; i<rupIdens.size(); i++) {
			RuptureIdentifier rupIden = rupIdens.get(i);
			String name = rupIden.getName();
			
			List<String> line = Lists.newArrayList();
			
			List<? extends SimulatorEvent> matches = rupIden.getMatches(events);
			
			double[] rps = PeriodicityPlotter.getRPs(matches);
			double mean = StatUtils.mean(rps);
			double var = StatUtils.variance(rps, mean);
			double sd = Math.sqrt(var);
			double cov = sd/mean;
			line.add(name);
			line.add(mean+"");
			line.add(DataUtils.median(rps)+"");
			line.add(var+"");
			line.add(sd+"");
			line.add(cov+"");
			double[] trimmedRPs = LogNormalDistReturnPeriodProvider.getTrimmedRPs(rps, mean);
			mean = StatUtils.mean(trimmedRPs);
			var = StatUtils.variance(trimmedRPs, mean);
			sd = Math.sqrt(var);
			cov = sd/mean;
			line.add(mean+"");
			line.add(DataUtils.median(rps)+"");
			line.add(var+"");
			line.add(sd+"");
			line.add(cov+"");
			line.add(Math.log(mean)+"");
			line.add(cov+"");
			
			distParamsCSV.addLine(line);
		}
		
		distParamsCSV.writeToFile(new File(writeDir, "dist_params.csv"));
		
		events = null;
		System.gc();
	}
	
	private static HistogramFunction getHist(List<? extends SimulatorEvent> matches, RandomDistType type) {
		double[] rps = PeriodicityPlotter.getRPs(matches);
		HistogramFunction func = getHistFunc(rps);
		String name;
		if (type == null)
			name = "Actual Distribution";
		else
			name = "RANDOMIZED "+type.getName();
		name += " ("+matches.size()+" events)";
		func.setName(name);
		return func;
	}
	
	static HistogramFunction getHistFunc(double[] rps) {
		HistogramFunction func = new HistogramFunction(5d, 100, 10d);
//		HistogramFunction func = new HistogramFunction(5d, 200, 5d);
		double max = func.getMaxX()+5d;
		for (double rp : rps)
			if (rp < max)
				func.add(rp, 1d);
		return func;
	}
	
	static RandomReturnPeriodProvider invertForSJDist(HistogramFunction func,
			CompoundDistReturnPeriodProvider initial, int numIterations,
			double paramSweepFract, int numSweepsPerParam) {
		final boolean D = true;
		Random r = new Random();
		
		int numRPs = (int)func.calcSumOfY_Vals();
		
		int printMod = numIterations / 100;
		
		List<RandomReturnPeriodProvider> prevProvs = initial.getProvs();
		List<Double> prevWeights = initial.getWeights();
		double prevE = calcEnergy(func, numRPs, initial);
		
		if (numSweepsPerParam > 0 && paramSweepFract > 0) {
			ParamSweep sweep = new ParamSweep(func, prevProvs, prevWeights, paramSweepFract, numSweepsPerParam, numRPs, prevE);
			sweep.doSweep();
			if (sweep.bestE < prevE) {
				prevE = sweep.bestE;
				prevProvs = sweep.bestProvs;
				prevWeights = sweep.bestWeights;
			}
		}
		
		List<RandomReturnPeriodProvider> bestProvs = initial.getProvs();
		List<Double> bestWeights = initial.getWeights();
		double bestE = prevE;
		
		int numKept = 0;
		
		double pFact = 0.03;
		boolean multiP = true;
		
		for (int i=0; i<numIterations; i++) {
			// first select provider to work on
			int ind = r.nextInt(prevProvs.size());
			RandomReturnPeriodProvider prov = prevProvs.get(ind);
			double weight = prevWeights.get(ind);
			
			List<Double> curWeights;
			List<RandomReturnPeriodProvider> curProvs;
			
			// perturb
			boolean perturbWeight = prevProvs.size() > 1 && r.nextDouble() < 0.3;
			if (perturbWeight) {
				double newWeight = getPerturbed(r, weight, pFact);
				
				curWeights = Lists.newArrayList(prevWeights);
				curWeights.set(ind, newWeight);
				curProvs = prevProvs;
			} else {
				// perturb a distribution
				if (prov instanceof LogNormalDistReturnPeriodProvider) {
					LogNormalDistribution l = ((LogNormalDistReturnPeriodProvider)prov).getDist();
					double shape = l.getShape();
					double scale = l.getScale();
					
					boolean pShape, pScale;
					if (multiP) {
						pShape = true;
						pScale = true;
					} else {
						pShape = r.nextBoolean();
						pScale = !pShape;
					}
					
					if (pShape) {
						// perturb shape
						shape = getPerturbed(r, shape, pFact);
					}
					if (pScale) {
						// perturb scale in linear space (mean)
						double mean = Math.exp(scale);
						mean = getPerturbed(r, mean, pFact);
						scale = Math.log(mean);
					}
					curWeights = prevWeights;
					curProvs = Lists.newArrayList(prevProvs);
					curProvs.set(ind, new LogNormalDistReturnPeriodProvider(scale, shape));
				} else {
					throw new IllegalStateException("Not yet supported");
				}
			}
			
			// now evaluate
			// create synthetic dist
			CompoundDistReturnPeriodProvider newProv = new CompoundDistReturnPeriodProvider(curProvs, curWeights);
			double e = calcEnergy(func, numRPs, newProv);
			double T = 1/Math.log(i + 1); // classical SA cooling schedule (Geman and Geman, 1984) (slow but ensures convergence)
			double P;
			if (e < prevE) {
				P = 1; // Always keep new model if better
			} else {
				// Sometimes keep new model if worse (depends on T)
				P = Math.exp(((prevE - e)) / (double) T); 
			}
			// Use transition probability to determine (via random number draw) if solution is kept
			if (P > r.nextDouble()) {
				// we're keeping it
				prevE = e;
				prevProvs = curProvs;
				prevWeights = curWeights;
				
				// is it a new best?
				if (e < bestE) {
					bestE = e;
					bestProvs = curProvs;
					bestWeights = curWeights;
				}
				
				numKept++;
				
				if (D && i % printMod == 0){
					System.out.println("iter "+i+" kept, e="+e+" prevE="+prevE+" bestE="+bestE);
					System.out.println("P="+P+", T="+T);
					for (int n=0; n<curProvs.size(); n++) {
						RandomReturnPeriodProvider prov2 = curProvs.get(n);
						System.out.println("\t"+prov2+" (weight="+curWeights.get(n)+")");
					}
				}
			} else {
				if (D && i % printMod == 0){
					System.out.println("iter "+i+" dropped, e="+e+" prevE="+prevE+" bestE="+bestE);
					System.out.println("PREV");
					for (int n=0; n<curProvs.size(); n++) {
						RandomReturnPeriodProvider prov2 = prevProvs.get(n);
						System.out.println("\t"+prov2+" (weight="+prevWeights.get(n)+")");
					}
					System.out.println("CUR");
					for (int n=0; n<curProvs.size(); n++) {
						RandomReturnPeriodProvider prov2 = curProvs.get(n);
						System.out.println("\t"+prov2+" (weight="+curWeights.get(n)+")");
					}
				}
			}
		}
		
		if (D) {
			System.out.println("Inverted Dist ("+numKept+"/"+numIterations+" kept):");
			System.out.println("BEST");
			for (int i=0; i<bestProvs.size(); i++) {
				RandomReturnPeriodProvider prov = bestProvs.get(i);
				System.out.println("\t"+prov+" (weight="+bestWeights.get(i)+")");
			}
			System.out.println("PREV");
			for (int i=0; i<prevProvs.size(); i++) {
				RandomReturnPeriodProvider prov = prevProvs.get(i);
				System.out.println("\t"+prov+" (weight="+prevWeights.get(i)+")");
			}
		}
		
//		return new CompoundDistReturnPeriodProvider(bestProvs, bestWeights);
		// TODO replace
		return new CompoundDistReturnPeriodProvider(prevProvs, prevWeights);
	}
	
	private static final boolean graphic_debug = false;
	private static int graphic_cnt = 0;
	
	private static double calcEnergy(HistogramFunction func, int numRPs, CompoundDistReturnPeriodProvider newProv) {
		double[] rps = new double[numRPs];
		for (int n=0; n<rps.length; n++)
			rps[n] = newProv.getReturnPeriod();
		HistogramFunction synFunc = getHistFunc(rps);
		HistogramFunction eFunc = null;
		if (graphic_debug && graphic_cnt < 100)
			eFunc = new HistogramFunction(synFunc.getMinX(), synFunc.size(), synFunc.getDelta());
		double e = 0;
		for (int n=0; n<func.size(); n++) {
			double data = func.getY(n);
			double syn = synFunc.getY(n);
			double absDiff = Math.abs(data - syn);
			e += absDiff * absDiff;
//			e += absDiff;
			if (graphic_debug && eFunc != null)
				eFunc.set(n, absDiff*absDiff);
		}
		if (graphic_debug) {
			List<HistogramFunction> funcs = Lists.newArrayList(func, synFunc);
			List<PlotCurveCharacterstics> chars = Lists.newArrayList(
					new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK),
					new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
//			funcs.add(e)
			PlotSpec distChar = new PlotSpec(funcs, chars, "Dist (e="+e+")", "Years", "Number");
			funcs = Lists.newArrayList();
			funcs.add(eFunc);
			chars = Lists.newArrayList(
					new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 2f, Color.BLACK));
//			funcs.add(e)
			PlotSpec eChar = new PlotSpec(funcs, chars, "Energies", "Years", "Energy");
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			List<PlotSpec> specs = Lists.newArrayList(distChar, eChar);
			gp.setCombinedOnYAxis(false);
			gp.drawGraphPanel(specs, false, false, null, null);
			gp.getChartPanel().setSize(1000, 800);
			try {
				String numStr = graphic_cnt+"";
				if (graphic_cnt<10)
					numStr = "0"+graphic_cnt;
				gp.saveAsPNG("/tmp/energy_"+numStr+".png");
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			graphic_cnt++;
		}
		return e;
	}
	
	private static double getPerturbed(Random r, double orig, double maxPerturbFract) {
		double rand = r.nextDouble();
		// now in range -0.5, 0.5
		rand -= 0.5;
		return orig + orig*rand*maxPerturbFract;
	}
	
	private static class ParamSweep {
		
		private HistogramFunction func;
		private List<RandomReturnPeriodProvider> provs;
		private List<Double> weights;
		private double paramSweepFract;
		private int numSweepsPerParam;
		private int numRPs;
		private double startE;
		
		private double bestE;
		private List<RandomReturnPeriodProvider> bestProvs;
		private List<Double> bestWeights;
		
		private int sweepCount = 0;
		
		public ParamSweep(HistogramFunction func,
				List<RandomReturnPeriodProvider> provs, List<Double> weights,
				double paramSweepFract, int numSweepsPerParam, int numRPs,
				double startE) {
			super();
			this.func = func;
			this.provs = provs;
			this.weights = weights;
			this.paramSweepFract = paramSweepFract;
			this.numSweepsPerParam = numSweepsPerParam;
			this.numRPs = numRPs;
			this.startE = startE;
			this.bestE = startE;
		}
		
		public void doSweep() {
			sweepRecursive(provs, weights, 0);
		}
		
		private void sweepRecursive(List<RandomReturnPeriodProvider> curProvs, List<Double> curWeights, int sweepIndex) {
			if (sweepIndex == curProvs.size()) {
				if (sweepCount % 1000 == 0)
					System.out.println("Calculating sweep "+sweepCount);
				// actually calculate
				CompoundDistReturnPeriodProvider newProv = new CompoundDistReturnPeriodProvider(curProvs, curWeights);
				double e = calcEnergy(func, numRPs, newProv);
				if (e < bestE) {
					bestProvs = Lists.newArrayList(curProvs);
					bestWeights = Lists.newArrayList(curWeights);
					bestE = e;
					System.out.println("new best e at sweep "+sweepCount);
				}
				sweepCount++;
				return;
			}
			RandomReturnPeriodProvider curProv = curProvs.get(sweepIndex);
			double origWeight = curWeights.get(sweepIndex);
			if (curProv instanceof LogNormalDistReturnPeriodProvider) {
				LogNormalDistribution l = ((LogNormalDistReturnPeriodProvider)curProv).getDist();
				double origScale = l.getScale();
				double origShape = l.getShape();
				
				double origMean = Math.exp(origScale);
				
				double meanMin = origMean*(1d-paramSweepFract);
				double meanSweepAmount = (2d*paramSweepFract*origMean)/(double)numSweepsPerParam;
				
				System.out.println("Mean: orig="+origMean+", min="+meanMin+", delta="+meanSweepAmount);
				System.exit(0);
				
				double shapeMin = origShape*(1d-paramSweepFract);
				double shapeSweepAmount = (2d*paramSweepFract*origShape)/(double)numSweepsPerParam;
				
				double weightMin = origWeight*(1d-paramSweepFract);
				double weightSweepAmount = (2d*paramSweepFract*origWeight)/(double)numSweepsPerParam;
				
				double mean = meanMin;
				for (int i=0; i<numSweepsPerParam; i++) {
					double scale = Math.log(mean);
					double shape = shapeMin;
					for (int j=0; j<numSweepsPerParam; j++) {
						LogNormalDistReturnPeriodProvider lProv = new LogNormalDistReturnPeriodProvider(scale, shape);
						
						curProvs.set(sweepIndex, lProv);
						double weight = weightMin;
						for (int k=0; k<numSweepsPerParam; k++) {
							curWeights.set(sweepIndex, weight);
							
							sweepRecursive(curProvs, curWeights, sweepIndex+1);
							
							weight += weightSweepAmount;
						}
						
						shape += shapeSweepAmount;
					}
					
					mean += meanSweepAmount; 
				}
			} else {
				// skip
				sweepRecursive(curProvs, curWeights, sweepIndex+1);
			}
		}
	}

}
