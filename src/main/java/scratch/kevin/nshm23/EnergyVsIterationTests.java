package scratch.kevin.nshm23;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.AnnealingProgress;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;

import com.google.common.base.Preconditions;

public class EnergyVsIterationTests {

	public static void main(String[] args) throws IOException {
		File invDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");

		File mainDir = new File(invDir, "2021_12_17-nshm23_draft_branches-FM3_1-CoulombRupSet");
//		File mainDir = new File(invDir, "2021_12_17-nshm23_draft_branches-max_dist-FM3_1-CoulombRupSet-TotNuclRate");
//		File mainDir = new File(invDir, "2021_12_17-nshm23_draft_branches-no_seg-FM3_1-CoulombRupSet");
//		File mainDir = new File(invDir, "2021_12_17-u3_branches-coulomb-FM3_1-5h");
//		File mainDir = new File(invDir, "2022_01_07-nshm23_draft_branches-no_seg-reweighted_even_fit-FM3_1-CoulombRupSet-SubB1-175_samples");
//		File mainDir = new File(invDir, "2022_01_10-nshm23_draft_branches-no_seg-reweighted_even_fit-conserve-FM3_1-CoulombRupSet-SubB1-105_samples");
//		File mainDir = new File(invDir, "");
		File resultsFile = new File(mainDir, "results.zip");
		
		File outputDir = new File(mainDir, "energy_vs_iters");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		SolutionLogicTree slt = SolutionLogicTree.load(resultsFile);
		LogicTree<?> tree = slt.getLogicTree();
		ZipFile	zip = new ZipFile(resultsFile);
		
		int numRuptures = -1;
		
		DefaultXY_DataSet fractChangeNextRound = new DefaultXY_DataSet();
		DefaultXY_DataSet fractChangeRateLater = new DefaultXY_DataSet();
		DefaultXY_DataSet fractChangeLater = new DefaultXY_DataSet();
		
		for (LogicTreeBranch<?> branch : tree) {
			if (numRuptures < 0)
				numRuptures = slt.forBranch(branch).getRupSet().getNumRuptures();
			
			String entryName = "solution_logic_tree/";
			for (int i=0; i<branch.size(); i++) {
				LogicTreeLevel<?> level = branch.getLevel(i);
				if (level.affects(AnnealingProgress.PROGRESS_FILE_NAME, true))
					entryName += branch.getValue(i).getFilePrefix()+"/";
			}
			entryName += AnnealingProgress.PROGRESS_FILE_NAME;
			System.out.println("Loading "+entryName);
			ZipEntry entry = zip.getEntry(entryName);
			Preconditions.checkNotNull(entry, "Entry not found: %s", entryName);
			
			CSVFile<String> csv = CSVFile.readStream(zip.getInputStream(entry), true);
			AnnealingProgress progress = new AnnealingProgress(csv);
			
			ArbitrarilyDiscretizedFunc energyIterFunc = new ArbitrarilyDiscretizedFunc();
			for (int i=0; i<progress.size(); i++)
				energyIterFunc.set((double)progress.getIterations(i), progress.getEnergies(i)[0]);
			
			double lastIter = energyIterFunc.getX(energyIterFunc.size()-1);
			double lastE = energyIterFunc.getY(energyIterFunc.size()-1);
			for (int i=0; i<energyIterFunc.size()-1; i++) {
				double iter1 = energyIterFunc.getX(i);
				double iter2 = energyIterFunc.getX(i+1);
				double e1 = energyIterFunc.getY(i);
				double e2 = energyIterFunc.getY(i+1);
				
				double iters = iter2-iter1;
				double eChange = Math.abs(e1 - e2);
				double eFract = eChange/e1;
				
				double round1 = iter1/(double)numRuptures;
				
				double deltaRounds = iters/(double)numRuptures;
				
				double fractPerRound = eFract/deltaRounds;
				
				fractChangeNextRound.set(round1, fractPerRound);
				
				// now until end
				iters = lastIter-iter1;
				eChange = Math.abs(e1 - lastE);
				eFract = eChange/e1;
				
				fractChangeLater.set(round1, eFract);
				
				deltaRounds = iters/(double)numRuptures;
				
				fractPerRound = eFract/deltaRounds;
				
				fractChangeRateLater.set(round1, fractPerRound);
			}
		}
		
		List<DefaultXY_DataSet> scatters = new ArrayList<>();
		List<String> titles = new ArrayList<>();
		List<String> yAxisLabels = new ArrayList<>();
		List<String> prefixes = new ArrayList<>();
		
		scatters.add(fractChangeNextRound);
		titles.add("Improvement Next Round");
		prefixes.add("improve_next_round");
		yAxisLabels.add("Fractional Change Round");
		
		scatters.add(fractChangeRateLater);
		titles.add("Eventual Improvement Rate");
		prefixes.add("improve_rate_to_end");
		yAxisLabels.add("Fractional Change Per Remaining Round");
		
		scatters.add(fractChangeLater);
		titles.add("Eventual Net Improvement");
		prefixes.add("improve_net_to_end");
		yAxisLabels.add("Fractional Change To Final");
		
		for (int i=0; i<scatters.size(); i++) {
			DefaultXY_DataSet scatter = scatters.get(i);
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			funcs.add(scatter);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLACK));
			
			PlotSpec spec = new PlotSpec(funcs, chars, titles.get(i),
					"Iteration Round (iters/numRups)", yAxisLabels.get(i));
			
			List<PlotSpec> specs = List.of(spec, spec);
			List<Range> yRanges = List.of(new Range(0d, 1.1*scatter.getMaxY()), new Range(1e-8, 1));
			List<Boolean> yLogs = List.of(false, true);
			List<Boolean> xLogs = List.of(false);
			List<Range> xRanges = List.of(new Range(0, 1.05*scatter.getMaxX()));
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(specs, xLogs, yLogs, xRanges, yRanges);
			
			String prefix = prefixes.get(i);
			
			PlotUtils.writePlots(outputDir, prefix, gp, 800, 1200, true, false, false);
		}
		
		zip.close();
	}

}
