package scratch.kevin.ucerf3.eal;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotElement;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.jfreechart.tornado.TornadoDiagram;
import org.opensha.commons.util.ClassUtils;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;
import org.opensha.sha.imr.AttenRelRef;

import scratch.UCERF3.analysis.BranchSensitivityHistogram;
import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.enumTreeBranches.MomentRateFixes;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.logicTree.LogicTreeBranchNode;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Table;
import com.google.common.collect.Table.Cell;

public class UCERF3_EALHistGen {

	public static void main(String[] args) throws IOException {
		Table<AttenRelRef, String, File> filesTable = HashBasedTable.create();
		
//		File origRunDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_01_15-ucerf3-eal-calc-NGA2s-2013");
//		File cbRerunDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_03_19-ucerf3-eal-calc-CB2014-recalc");
//		File askRerunDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_04_07-ucerf3-eal-calc-ASK2014-recalc");
//		File origRunDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_05_05-ucerf3-eal-calc-wald-vs30");
//		File origRunDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_05_16-ucerf3-99percent-wills");
//		File origRunDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_05_28-ucerf3-fatality-smaller");
		
		Map<String, File> runDirs = Maps.newHashMap();
		String runCategoryName = "Vs30";
		runDirs.put("Wills", new File("/home/kevin/OpenSHA/UCERF3/eal/2014_05_28-ucerf3-99percent-wills-smaller"));
		runDirs.put("Wald", new File("/home/kevin/OpenSHA/UCERF3/eal/2016_06_06-ucerf3-90percent-wald"));
		double inflationMultiplier = 10d/9d;
		
		String units = "$ Billions";
		double multiplier = 1d/1e6;
		double delta = 0.1;
		
//		String units = "Fatalities";
//		double multiplier = 1;
//		double delta = 1;
		
		multiplier *= inflationMultiplier;
		
		File plotDir = new File("/tmp/eal_plot");
		if (!plotDir.exists())
			plotDir.mkdir();
		
		Map<AttenRelRef, Double> imrWeightsMap = Maps.newHashMap();
		imrWeightsMap.put(AttenRelRef.CB_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.CY_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.ASK_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.BSSA_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.IDRISS_2014, 0.12);
		
//		imrWeightsMap.put(AttenRelRef.BSSA_2014, 1d);
		
		MagDependentAperiodicityOptions[] covs = { null, MagDependentAperiodicityOptions.HIGH_VALUES,
				MagDependentAperiodicityOptions.MID_VALUES, MagDependentAperiodicityOptions.LOW_VALUES };
		
		String gmpeCategoryName = "GMPE";
		String probModelCategoryName = "ProbModel";
		
		BranchSensitivityHistogram hist = new BranchSensitivityHistogram("EAL ("+units+")");
		
		for (String runName : runDirs.keySet()) {
			File runDir = runDirs.get(runName);
			
			Map<AttenRelRef, File> imrResultsDirMap = Maps.newHashMap();
			imrResultsDirMap.put(AttenRelRef.CB_2014, runDir);
			imrResultsDirMap.put(AttenRelRef.CY_2014, runDir);
			imrResultsDirMap.put(AttenRelRef.ASK_2014, runDir);
			imrResultsDirMap.put(AttenRelRef.BSSA_2014, runDir);
			imrResultsDirMap.put(AttenRelRef.IDRISS_2014, runDir);
			
			for (AttenRelRef ref : imrResultsDirMap.keySet()) {
				File resultDir = imrResultsDirMap.get(ref);
				for (MagDependentAperiodicityOptions cov : covs) {
					String covName;
					if (cov == null)
						covName = "POISSON";
					else
						covName = cov.name();
					filesTable.put(ref, covName, new File(resultDir, ref.name()+"_"+covName+"_eals.csv"));
				}
			}
			
			for (Cell<AttenRelRef, String, File> cell : filesTable.cellSet()) {
				AttenRelRef imr = cell.getRowKey();
				String imrName = imr.name();
				String probModelName = cell.getColumnKey();
				File csvFile = cell.getValue();
				CSVFile<String> csv = CSVFile.readFile(csvFile, true);
				
				double imrWeight = imrWeightsMap.get(imr);
				MagDependentAperiodicityOptions cov;
				if (probModelName.equals("POISSON"))
					cov = null;
				else
					cov = MagDependentAperiodicityOptions.valueOf(probModelName);
				double probModelWeight = FaultSystemSolutionERF.getWeightForCOV(cov);
				
				for (int row=1; row<csv.getNumRows(); row++) {
					List<String> line = csv.getLine(row);
					LogicTreeBranch branch = fromLine(line);
					double weight = Double.parseDouble(line.get(1));
					double eal = Double.parseDouble(line.get(2)); // total eal
					
					eal *= multiplier;
					
					// IMR & prob model weights
					weight *= imrWeight * probModelWeight;
					
					String[] extraVals;
					if (runDirs.size() == 1)
						extraVals = new String[] {
								gmpeCategoryName, imrName, probModelCategoryName, probModelName};
					else
						extraVals = new String[] {
								gmpeCategoryName, imrName, probModelCategoryName, probModelName, runCategoryName, runName};
					
					hist.addValues(branch, eal, weight, extraVals);
				}
			}
		}
		
		Map<String, PlotSpec> specs = hist.getStackedHistPlots(true, delta);
		
		List<File> histPDFs = Lists.newArrayList();
		List<String> names = Lists.newArrayList();
		for (Class<? extends LogicTreeBranchNode<?>> clazz : LogicTreeBranch.getLogicTreeNodeClasses()) {
			if (clazz.equals(InversionModels.class) || clazz.equals(MomentRateFixes.class))
				continue;
			names.add(ClassUtils.getClassNameWithoutPackage(LogicTreeBranch.getEnumEnclosingClass(clazz)));
		}
		names.add(probModelCategoryName);
		names.add(gmpeCategoryName);
		if (runDirs.size() > 1)
			names.add(runCategoryName);
		for (String name : names) {
			PlotSpec histSpec = specs.get(name);
			if (histSpec == null) {
				System.out.println("WARNING: no spec found for "+name);
				continue;
			}
			Preconditions.checkNotNull(histSpec, "No plot found for: "+name);
			
			List<? extends PlotElement> elems = histSpec.getPlotElems();
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			CommandLineInversionRunner.setFontSizes(gp);
			EvenlyDiscretizedFunc f1 = (EvenlyDiscretizedFunc) elems.get(0);
			double min = f1.getMinX();
			int num = f1.size();
			double plotMaxY = 0d;
			for (int i=0; i<elems.size(); i++) {
				if (elems.get(i) instanceof DiscretizedFunc) {
					double max = ((XY_DataSet)elems.get(i)).getMaxY();
					if (max > plotMaxY)
						plotMaxY = max;
				}
			}
			plotMaxY *= 1.3;
			if (plotMaxY > 1d)
				plotMaxY = 1d;
			// pad by a delta
			double plotMinX = min-0.5*delta - 0.5*delta;
			double plotMaxX = min+(num-0.5)*delta + 0.5*delta;
			gp.setUserBounds(plotMinX, plotMaxX, 0, plotMaxY);
			if (plotMinX >= plotMaxX) {
				System.out.println("Data bounds: "+f1.getMinX()+" "+f1.getMaxX()
						+" "+f1.getMinY()+" "+f1.getMaxY());
				System.out.println("Plot bounds: "+plotMinX+" "+plotMaxX+" 0 "+plotMaxY);
				System.out.println("Delta="+delta);
				System.out.println(f1);
			}
			
			gp.drawGraphPanel(histSpec);
			File file = new File(plotDir, name);
			gp.getChartPanel().setSize(500, 400);
			gp.saveAsPDF(file.getAbsolutePath()+".pdf");
			gp.saveAsPNG(file.getAbsolutePath()+".png");
			histPDFs.add(new File(file.getAbsolutePath()+".pdf"));
		}
		
		try {
			FaultSysSolutionERF_Calc.combineBranchSensHists(histPDFs, new File(plotDir, "eal_branch_sens_hists.pdf"));
		} catch (com.lowagie.text.DocumentException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		hist.getStaticsticsCSV().writeToFile(new File(plotDir, "eal_branch_sens_stats.csv"));
		
		TornadoDiagram tornadoShift = hist.getTornadoDiagram("EAL Sensitivity", true);
		tornadoShift.getHeadlessPlot(500, 600).saveAsPDF(
				new File(plotDir, "eal_branch_sens_tornado_meanshift.pdf").getAbsolutePath());
		tornadoShift.getHeadlessPlot(500, 600).saveAsPNG(
				new File(plotDir, "eal_branch_sens_tornado_meanshift.png").getAbsolutePath());
		TornadoDiagram tornado = hist.getTornadoDiagram("EAL Sensitivity", false);
		tornado.getHeadlessPlot(500, 600).saveAsPDF(
				new File(plotDir, "eal_branch_sens_tornado_mean.pdf").getAbsolutePath());
		tornado.getHeadlessPlot(500, 600).saveAsPNG(
				new File(plotDir, "eal_branch_sens_tornado_mean.png").getAbsolutePath());
	}
	
	private static LogicTreeBranch fromLine(List<String> line) {
		int col = 5;
		List<LogicTreeBranchNode<?>> nodes = Lists.newArrayList();
		
		for (Class<? extends LogicTreeBranchNode<?>> clazz : LogicTreeBranch.getLogicTreeNodeClasses()) {
			String strVal = line.get(col++);
			LogicTreeBranchNode<?> node = null;
			for (LogicTreeBranchNode<?> choice : clazz.getEnumConstants()) {
				if (choice.getShortName().equals(strVal)) {
					node = choice;
					break;
				}
			}
			Preconditions.checkNotNull(node);
			nodes.add(node);
		}
		Preconditions.checkState(col == line.size());
		
		return LogicTreeBranch.fromValues(nodes);
	}

}
