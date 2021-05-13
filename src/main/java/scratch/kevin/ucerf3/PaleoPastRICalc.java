package scratch.kevin.ucerf3;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.paleoRateConstraints.PaleoRateConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoProbabilityModel;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoRateConstraintFetcher;

public class PaleoPastRICalc {

	public static void main(String[] args) throws IOException, DocumentException {
		FaultSystemSolution sol = FaultSystemIO.loadSol(
			 new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
			 		+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		FaultSystemRupSet rupSet = sol.getRupSet();
		List<PaleoRateConstraint> constraints =
				 UCERF3_PaleoRateConstraintFetcher.getConstraints(rupSet.getFaultSectionDataList());
		List<HashSet<Integer>> rupsForConstraints = Lists.newArrayList();
		List<Double> meanRIs = Lists.newArrayList();
		List<Double> meanPaleoVisibleRIs = Lists.newArrayList();
		List<Double> paleoConstraintRIs = Lists.newArrayList();
		CSVFile<String> paleoCSV = new CSVFile<String>(true);
		paleoCSV.addLine("Site Name", "Subsection Index", "Subsection Name", "Paleo Constraint MRI",
				"U3 Paleo Visible MRI", "U3 Supra-seis MRI", "Current Open Interval");
		long now = System.currentTimeMillis();
		UCERF3_PaleoProbabilityModel paleoProbModel = UCERF3_PaleoProbabilityModel.load();
		
		int numPastMean = 0;
		int numPastPaleoVisible = 0;
		int numPastPaleoConstraint = 0;
		for (PaleoRateConstraint constr : constraints) {
			int sect = constr.getSectionIndex();
			HashSet<Integer> rups = new HashSet<Integer>(rupSet.getRupturesForSection(sect));
			rupsForConstraints.add(rups);
			double ri = 1d/sol.calcParticRateForSect(sect, 0d, 10d);
			meanRIs.add(ri);
			double paleoVisibleRI = 1d/sol.calcTotPaleoVisibleRateForSect(sect, paleoProbModel);
			meanPaleoVisibleRIs.add(paleoVisibleRI);
			double paleoRI = 1d/constr.getMeanRate();
			paleoConstraintRIs.add(paleoRI);
			long dateLast = rupSet.getFaultSectionData(sect).getDateOfLastEvent();
			String oiStr = "(unknown)";
			if (dateLast > Long.MIN_VALUE) {
				long delta = now - dateLast;
				double oi = (double)delta/ProbabilityModelsCalc.MILLISEC_PER_YEAR;
				if (oi > ri)
					numPastMean++;
				if (oi > paleoVisibleRI)
					numPastPaleoVisible++;
				if (oi > paleoRI)
					numPastPaleoConstraint++;
				oiStr = (float)oi+"";
			}
			paleoCSV.addLine(constr.getPaleoSiteName(), constr.getSectionIndex()+"", constr.getFaultSectionName(),
					paleoRI+"", paleoVisibleRI+"", ri+"", oiStr);
		}
		
		System.out.println(numPastMean+" currently past U3 Supra-Seismogenic MRI");
		System.out.println(numPastPaleoVisible+" currently past U3 Paleo-Visible MRI");
		System.out.println(numPastPaleoConstraint+" currently past Paleo Constraint MRI");
		
		File runDir = new File("/home/kevin/OpenSHA/UCERF3/time_dep_catalogs/2017_04_10-long-catalogs/batch0");
		paleoCSV.writeToFile(new File(runDir, "paleo_site_mris.csv"));
		
		List<List<Double>> mrisList = Lists.newArrayList();
		List<String> mrisTitles = Lists.newArrayList();
		List<String> mrisPrefixes = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		List<Integer> mrisCurrent = Lists.newArrayList();
		
		mrisList.add(paleoConstraintRIs);
		mrisTitles.add("Paleo Constraint MRI");
		mrisPrefixes.add("paleo_data");
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		mrisCurrent.add(numPastPaleoConstraint);
		
		mrisList.add(meanRIs);
		mrisTitles.add("U3 Supra-Seismogenic MRI");
		mrisPrefixes.add("u3_supra_seis");
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		mrisCurrent.add(numPastMean);
		
		mrisList.add(meanPaleoVisibleRIs);
		mrisTitles.add("U3 Paleo-Visible MRI");
		mrisPrefixes.add("u3_paleo_visible");
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
		mrisCurrent.add(numPastPaleoVisible);
		
		List<XY_DataSet> funcs = Lists.newArrayList();
		List<XY_DataSet> cumulativeFuncs = Lists.newArrayList();
		
		List<List<Integer>> catalogIndexes = Lists.newArrayList();
		List<List<Double>> catalogYears = Lists.newArrayList();
		
		for (File subDir : runDir.listFiles()) {
			if (!subDir.isDirectory())
				continue;
			File eventsFile = new File(subDir, "sampledEventsData.txt");
			if (!eventsFile.exists())
				continue;
			List<Integer> fssIndexes = Lists.newArrayList();
			List<Double> years = Lists.newArrayList();
			loadCatalog(eventsFile, fssIndexes, years);
			double delta = years.get(years.size()-1) - years.get(0);
			System.out.println("Loaded "+delta+" years from "+subDir.getName());
			
			catalogIndexes.add(fssIndexes);
			catalogYears.add(years);
		}
		
		for (int m=0; m<mrisList.size(); m++) {
			List<Double> mris = mrisList.get(m);
			String mrisTitle = mrisTitles.get(m);
			
			EvenlyDiscretizedFunc meanFunc = null;
			for (int c=0; c<catalogIndexes.size(); c++) {
				List<Integer> fssIndexes = catalogIndexes.get(c);
				List<Double> years = catalogYears.get(c);
				
				EvenlyDiscretizedFunc func = calc(mris, rupsForConstraints, fssIndexes, years);
				if (meanFunc == null)
					meanFunc = func;
				else
					for (int i=0; i<meanFunc.size(); i++)
						meanFunc.add(i, func.getY(i));
			}
			
			Preconditions.checkState(meanFunc != null);
			double totYears = meanFunc.calcSumOfY_Vals();
			meanFunc.scale(1d/totYears);
			
			meanFunc.setName(mrisTitle);
			funcs.add(meanFunc);
			
			EvenlyDiscretizedFunc cumulativeFunc =
					new EvenlyDiscretizedFunc(meanFunc.getMinX(), meanFunc.getMaxX(), meanFunc.size());
			for (int i=0; i<meanFunc.size(); i++) {
				double sumAtAbove = 0d;
				for (int j=i; j<meanFunc.size(); j++)
					sumAtAbove += meanFunc.getY(j);
				cumulativeFunc.set(i, sumAtAbove);
			}
			
			cumulativeFunc.setName(mrisTitle);
			cumulativeFuncs.add(cumulativeFunc);
		}
		
		// now add crosses
		for (int i=0; i<mrisList.size(); i++) {
			int num = mrisCurrent.get(i);
			double histVal = funcs.get(i).getY(num);
			double cumVal = cumulativeFuncs.get(i).getY(num);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.BOLD_X, 5f, chars.get(i).getColor()));
			DefaultXY_DataSet xy = new DefaultXY_DataSet();
			xy.set((double)num, histVal);
			funcs.add(xy);
			xy = new DefaultXY_DataSet();
			xy.set((double)num, cumVal);
			cumulativeFuncs.add(xy);
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "U3-TD Open Intervals", "# Sites Beyond MRI", "Fraction of Catalog");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		gp.setBackgroundColor(Color.WHITE);
		
		gp.setUserBounds(new Range(0d, constraints.size()), null);
		gp.drawGraphPanel(spec);
		gp.getChartPanel().setSize(1000, 800);
		String prefix = "open_intervals";
		gp.saveAsPNG(new File(runDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(runDir, prefix+".pdf").getAbsolutePath());
		gp.saveAsTXT(new File(runDir, prefix+".txt").getAbsolutePath());
		
		List<String> header = Lists.newArrayList("Num Sites");
		for (String title : mrisTitles)
			header.add(title);
		CSVFile<String> csv = new CSVFile<String>(true);
		csv.addLine(header);
		for (int i=0; i<constraints.size(); i++) {
			List<String> line = Lists.newArrayList(i+"");
			for (int f=0; f<mrisList.size(); f++)
				line.add(funcs.get(f).getY(i)+"");
			csv.addLine(line);
		}
		csv.writeToFile(new File(runDir, prefix+".csv"));
		
		spec = new PlotSpec(cumulativeFuncs, chars, "U3-TD Open Intervals", "N Sites", "Prob At Least N Sites Above MRI");
		spec.setLegendVisible(true);
		
		gp.setUserBounds(new Range(0d, constraints.size()), new Range(0d, 1d));
		gp.drawGraphPanel(spec);
		gp.getChartPanel().setSize(1000, 800);
		prefix = "open_intervals_cumulative";
		gp.saveAsPNG(new File(runDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(runDir, prefix+".pdf").getAbsolutePath());
		gp.saveAsTXT(new File(runDir, prefix+".txt").getAbsolutePath());
		
		csv = new CSVFile<String>(true);
		csv.addLine(header);
		for (int i=0; i<constraints.size(); i++) {
			List<String> line = Lists.newArrayList(i+"");
			for (int f=0; f<mrisList.size(); f++)
				line.add(cumulativeFuncs.get(f).getY(i)+"");
			csv.addLine(line);
		}
		csv.writeToFile(new File(runDir, prefix+".csv"));
	}
	
	private static final double delta_years = 1;
	
	private static EvenlyDiscretizedFunc calc(List<Double> meanRIs, List<HashSet<Integer>> rupsForConstraints,
			List<Integer> fssIndexes, List<Double> years) {
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(0d, rupsForConstraints.size()+1, 1d);
		
		List<List<Double>> paleoRuptureYears = Lists.newArrayList();
		for (int i=0; i<rupsForConstraints.size(); i++)
			paleoRuptureYears.add(new ArrayList<Double>());
		
		for (int y=0; y<years.size(); y++) {
			int fssIndex = fssIndexes.get(y);
			double year = years.get(y);
			
			for (int i=0; i<rupsForConstraints.size(); i++) {
				if (rupsForConstraints.get(i).contains(fssIndex)) {
					// it's a hit on this section
					paleoRuptureYears.get(i).add(year);
				}
			}
		}
		
		double startYear = Double.NEGATIVE_INFINITY;
		double endYear = Double.NEGATIVE_INFINITY;
		for (int i=0; i<rupsForConstraints.size(); i++) {
			Preconditions.checkState(paleoRuptureYears.size() > 1);
			List<Double> sitePaleoYears = paleoRuptureYears.get(i);
			double siteStartYear = sitePaleoYears.get(0);
			startYear = Double.max(startYear, siteStartYear);
			endYear = Double.max(endYear, sitePaleoYears.get(sitePaleoYears.size()-1));
		}
		
		int[] curIndexes = new int[paleoRuptureYears.size()];
		
		for (double year=startYear; year<endYear; year+=delta_years) {
//			System.out.println("Year: "+year);
			int numAbove = 0;
			for (int i=0; i<paleoRuptureYears.size(); i++) {
				List<Double> sitePaleoYears = paleoRuptureYears.get(i);
				while (curIndexes[i] < sitePaleoYears.size()-1) {
					double nextYear = sitePaleoYears.get(curIndexes[i]+1);
					if (nextYear <= year)
						curIndexes[i]++;
					else
						break;
				}
				double sitePrevYear = sitePaleoYears.get(curIndexes[i]);
				double targetYear = sitePrevYear + meanRIs.get(i);
				if (year > targetYear)
					numAbove++;
			}
			
			func.add(numAbove, delta_years);
		}
		
		return func;
	}
	
	private static void loadCatalog(File file, List<Integer> fssIndexes, List<Double> years) throws IOException {
		Preconditions.checkState(fssIndexes.isEmpty() && years.isEmpty());
		
		for (String line : Files.readLines(file, Charset.defaultCharset())) {
			line = line.trim();
			if (line.startsWith("nth"))
				// header
				continue;
			
			// nthRupIndex	fssRupIndex	year	epoch	normRI	mag	nthCatalog	timeToNextInYrs	utilizedPaleoSite
			String[] split = line.split("\t");
			if (split.length < 3) // partial line
				break;
			int fssIndex = Integer.parseInt(split[1]);
			Preconditions.checkState(fssIndex >= 0);
			double year = Double.parseDouble(split[2]);
			Preconditions.checkState(years.isEmpty() || year >= years.get(years.size()-1));
			
			fssIndexes.add(fssIndex);
			years.add(year);
		}
	}

}
