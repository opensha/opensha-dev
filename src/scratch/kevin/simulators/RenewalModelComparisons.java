package scratch.kevin.simulators;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.calc.recurInterval.BPT_DistCalc;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.iden.SectionIDIden;
import org.opensha.sha.simulators.parsers.EQSIMv06FileReader;

import com.google.common.collect.Lists;

import scratch.kevin.markov.EmpiricalMarkovChain;
import scratch.kevin.markov.MarkovChain;
import scratch.kevin.simulators.SynchIdens.SynchFaults;

public class RenewalModelComparisons {
	
	public static void main(String[] args) throws IOException {
		double distSpacing = 10d;
		boolean sectIden = true;
		List<RuptureIdentifier> rupIdens;
		if (sectIden) {
			rupIdens = Lists.newArrayList();
			File geomFile = new File("/home/kevin/Simulators/ALLCAL2_1-7-11_Geometry.dat");
			List<SimulatorElement> elems = EQSIMv06FileReader.readGeometryFile(geomFile);
			rupIdens.add(SectionIDIden.getALLCAL2_SSAF_Mojave(elems));
			rupIdens.add(SectionIDIden.getALLCAL2_SSAF_Coachella(elems));
			for (RuptureIdentifier iden : rupIdens)
				((SectionIDIden)iden).setMomentFractForInclusion(0.05);
		} else {
			rupIdens = SynchIdens.getIndividualFaults(7d, 10d, SynchFaults.SAF_MOJAVE, SynchFaults.SAF_COACHELLA);
		}
		List<? extends SimulatorEvent> events = new SimAnalysisCatLoader(true, rupIdens, true).getEvents();
		// generate Markov Chain
		EmpiricalMarkovChain chain = MarkovChainBuilder.build(distSpacing, events, rupIdens);
		chain.getCollapsedChain(0);
		
		File outputDir = new File("/home/kevin/Simulators/recurrence_example");
		
		double maxYears = 350d;
		double cov = 0.3;
		
		String name1 = rupIdens.get(0).getName();
		String name2 = rupIdens.get(1).getName();
		
		double ri1 = calcMeanRI(rupIdens.get(0), events);
		double ri2 = calcMeanRI(rupIdens.get(1), events);
		
		EvenlyDiscretizedFunc bpt1 = calcBPT(ri1, maxYears, cov, distSpacing);
		EvenlyDiscretizedFunc bpt2 = calcBPT(ri2, maxYears, cov, distSpacing);
		
		EvenlyDiscretizedFunc actual1 = calcCondProb(0, chain, maxYears);
		EvenlyDiscretizedFunc actual2 = calcCondProb(1, chain, maxYears);
		
		write1D(name1, ri1, bpt1, actual1, outputDir);
		write1D(name2, ri2, bpt2, actual2, outputDir);
		
//		EvenlyDiscrXYZ_DataSet indepBothBPT = calcIndepProbBothBPT(bpt1, bpt2);
//		plot2D(indepBothBPT, "BPT Prob Both", name1, name2, outputDir, "bpt_prob_both");
//		EvenlyDiscrXYZ_DataSet indepEitherBPT = calcIndepProbEither(bpt1, bpt2);
//		plot2D(indepEitherBPT, "BPT Prob Either", name1, name2, outputDir, "bpt_prob_either");
//		
//		EvenlyDiscrXYZ_DataSet indepBothRSQSim = calcIndepProbBothBPT(actual1, actual2);
//		plot2D(indepBothRSQSim, "RSQSim Independent Prob Both", name1, name2, outputDir, "rsqsim_indep_prob_both");
//		EvenlyDiscrXYZ_DataSet indepEitherRsqsim = calcIndepProbEither(actual1, actual2);
//		plot2D(indepEitherRsqsim, "RSQSim Independent Prob Either", name1, name2, outputDir, "rsqsim_indep_prob_either");
//		
//		EvenlyDiscrXYZ_DataSet actualBoth = calcMarkovProbBoth(chain, bpt1.size()+2);
//		plot2D(actualBoth, "RSQSim Prob Both", name1, name2, outputDir, "rsqsim_actual_prob_both");
//		EvenlyDiscrXYZ_DataSet actualEither = calcMarkovProbEither(chain, bpt1.size()+2);
//		plot2D(actualEither, "RSQSim Prob Either", name1, name2, outputDir, "rsqsim_actual_prob_either");
		
		EvenlyDiscrXYZ_DataSet indepOne = calcIndepProbOne(actual1, 0);
		plot2D(indepOne, "RSQSim Prob "+name1, name1, name2, outputDir, "rsqsim_indep_prob_"+name1.replaceAll(" ", "_"));
		EvenlyDiscrXYZ_DataSet indepTwo = calcIndepProbOne(actual2, 1);
		plot2D(indepTwo, "RSQSim Prob "+name2, name1, name2, outputDir, "rsqsim_indep_prob_"+name2.replaceAll(" ", "_"));
		
		EvenlyDiscrXYZ_DataSet actualOne = calcMarkovProbOneOrBoth(chain, bpt1.size()+2, 0);
		plot2D(actualOne, "RSQSim Prob "+name1, name1, name2, outputDir, "rsqsim_actual_prob_"+name1.replaceAll(" ", "_"));
		EvenlyDiscrXYZ_DataSet actualTwo = calcMarkovProbOneOrBoth(chain, bpt1.size()+2, 1);
		plot2D(actualTwo, "RSQSim Prob "+name2, name1, name2, outputDir, "rsqsim_actual_prob_"+name2.replaceAll(" ", "_"));
	}
	
	private static double calcMeanRI(RuptureIdentifier iden, List<? extends SimulatorEvent> events) {
		double years = 0d;
		List<? extends SimulatorEvent> matches = iden.getMatches(events);
		for (int i=1; i<matches.size(); i++)
			years += matches.get(i).getTimeInYears() - matches.get(i-1).getTimeInYears();
		years /= (double)(matches.size()-1);
		
		System.out.println("Mean RI for "+iden.getName()+": "+years);
		
		return years;
	}
	
	private static EvenlyDiscretizedFunc calcCondProb(int index, EmpiricalMarkovChain chain, double max) {
		double distSpacing = chain.getDistSpacing();
		int numPoints = (int)(max/distSpacing+0.5);
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(distSpacing*0.5, numPoints, distSpacing);
		EmpiricalMarkovChain collapsed = chain.getCollapsedChain(index);
		for (int i=0; i<numPoints; i++) {
			int[] fromState = {i};
			int[] toState = {0};
			double rupProb = collapsed.getTransitionProb(fromState, toState);
			if (Double.isNaN(rupProb))
				rupProb = 0d;
			func.set(i, rupProb);
		}
		return func;
	}
	
	private static EvenlyDiscretizedFunc calcBPT(double mean, double max, double cov, double distSpacing) {
		BPT_DistCalc dist = new BPT_DistCalc();
		int numPoints = (int)(max/distSpacing+0.5);
		dist.setAll(mean, cov, distSpacing, numPoints, distSpacing);
		return dist.getCondProbFunc();
//		return dist.getPDF();
	}
	
	private static void write1D(String name, double meanRI, EvenlyDiscretizedFunc bpt, EvenlyDiscretizedFunc actual,
			File outputDir) throws IOException {
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		bpt.setName("BPT");
		funcs.add(bpt);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		actual.setName("RSQSim");
		funcs.add(actual);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "1-D Renewal Models", "Time (years)", "Conditional Probability");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		gp.setBackgroundColor(Color.WHITE);
		
//		gp.setUserBounds(new Range(0d, constraints.size()), null);
		gp.drawGraphPanel(spec);
		gp.getChartPanel().setSize(1000, 500);
		String prefix = name.replaceAll(" ", "_");
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
		gp.saveAsTXT(new File(outputDir, prefix+".txt").getAbsolutePath());
	}
	
	private static EvenlyDiscrXYZ_DataSet calcIndepProbBothBPT(EvenlyDiscretizedFunc bpt1, EvenlyDiscretizedFunc bpt2) {
		double gridSpacing = bpt1.getDelta();
		EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(bpt1.size(), bpt2.size(), 0.5*gridSpacing, 0.5*gridSpacing, gridSpacing);
		
		for (int x=0; x<bpt1.size(); x++)
			for (int y=0; y<bpt2.size(); y++)
				xyz.set(x, y, bpt1.getY(x)*bpt2.getY(y));
		
		return xyz;
	}
	
	private static EvenlyDiscrXYZ_DataSet calcIndepProbEither(EvenlyDiscretizedFunc bpt1, EvenlyDiscretizedFunc bpt2) {
		double gridSpacing = bpt1.getDelta();
		EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(bpt1.size(), bpt2.size(), 0.5*gridSpacing, 0.5*gridSpacing, gridSpacing);
		
		for (int x=0; x<bpt1.size(); x++)
			for (int y=0; y<bpt2.size(); y++)
				xyz.set(x, y, 1d - (1d-bpt1.getY(x))*(1d-bpt2.getY(y)));
		
		return xyz;
	}
	
	private static EvenlyDiscrXYZ_DataSet calcIndepProbOne(EvenlyDiscretizedFunc bpt, int index) {
		double gridSpacing = bpt.getDelta();
		EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(bpt.size(), bpt.size(), 0.5*gridSpacing, 0.5*gridSpacing, gridSpacing);
		
		for (int x=0; x<bpt.size(); x++) {
			for (int y=0; y<bpt.size(); y++) {
				if (index == 0)
					xyz.set(x, y, bpt.getY(x));
				else
					xyz.set(x, y, bpt.getY(y));
			}
		}
		
		return xyz;
	}
	
	private static EvenlyDiscrXYZ_DataSet calcMarkovProbEither(MarkovChain chain, int size) {
		double gridSpacing = chain.getDistSpacing();
		EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(size, size, 0.5*gridSpacing, 0.5*gridSpacing, gridSpacing);
		
		for (int x=0; x<size; x++) {
			for (int y=0; y<size; y++) {
				int[] fromState = {x,y};
				int[] toState = {x+1,y+1};
				double prob = 1d - chain.getTransitionProb(fromState, toState);
				xyz.set(x, y, prob);
			}
		}
		
		return xyz;
	}
	
	private static EvenlyDiscrXYZ_DataSet calcMarkovProbBoth(MarkovChain chain, int size) {
		double gridSpacing = chain.getDistSpacing();
		EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(size, size, 0.5*gridSpacing, 0.5*gridSpacing, gridSpacing);
		
		for (int x=0; x<size; x++) {
			for (int y=0; y<size; y++) {
				int[] fromState = {x,y};
				int[] toState = {0,0};
				double prob = chain.getTransitionProb(fromState, toState);
				xyz.set(x, y, prob);
			}
		}
		
		return xyz;
	}
	
	private static EvenlyDiscrXYZ_DataSet calcMarkovProbOneOrBoth(MarkovChain chain, int size, int index) {
		double gridSpacing = chain.getDistSpacing();
		EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(size, size, 0.5*gridSpacing, 0.5*gridSpacing, gridSpacing);
		
		for (int x=0; x<size; x++) {
			for (int y=0; y<size; y++) {
				int[] fromState = {x,y};
				int[] toState = {0,0};
				double both = chain.getTransitionProb(fromState, toState);
				if (index == 0)
					toState[1] = y+1;
				else
					toState[0] = x+1;
				double one = chain.getTransitionProb(fromState, toState);
				xyz.set(x, y, both+one);
			}
		}
		
		return xyz;
	}
	
	private static void plot2D(EvenlyDiscrXYZ_DataSet xyz, String title, String name1, String name2, File outputDir, String prefix)
			throws IOException {
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, 1d);
		cpt.setNanColor(Color.WHITE);
		
		XYZPlotSpec xyzSpec = new XYZPlotSpec(xyz, cpt, title,
				name1+" OI", name2+" OI", null);
		
		int width = 600;
		int height = 680;
		
		Range range = new Range(0d, xyz.getMaxY()+0.5*xyz.getGridSpacingY());
		XYZGraphPanel panel = new XYZGraphPanel();
		panel.drawPlot(xyzSpec, false, false, range, range);
		
		// write plot
		panel.getChartPanel().setSize(width, height);
		File out = new File(outputDir, prefix); 
		panel.saveAsPNG(out.getAbsolutePath()+".png");
		panel.saveAsPDF(out.getAbsolutePath()+".pdf");
	}

}
