package scratch.kevin.ucerf3;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;

import com.google.common.base.CharMatcher;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.inversion.InversionFaultSystemRupSetFactory;
import scratch.UCERF3.utils.paleoRateConstraints.PaleoRateConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoRateConstraintFetcher;

public class PaleoMagHistGen {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		FaultSystemRupSet rupSet = InversionFaultSystemRupSetFactory.forBranch(DeformationModels.GEOLOGIC_PLUS_ABM);

		ArrayList<PaleoRateConstraint> paleoConsts = UCERF3_PaleoRateConstraintFetcher.getConstraints(rupSet.getFaultSectionDataList());

		CharMatcher match = CharMatcher.inRange('a', 'z').or(CharMatcher.inRange('A', 'Z'))
				.or(CharMatcher.inRange('0', '9')).or(CharMatcher.is('-')).precomputed();

		File dir = new File("/tmp/paleo");
		if (!dir.exists())
			dir.mkdir();

		for (PaleoRateConstraint paleo : paleoConsts) {
			String paleoName = paleo.getPaleoSiteName();
			Preconditions.checkNotNull(paleoName, "Paleo name is null!!!!");
			String name = match.retainFrom(paleoName);
			
			System.out.println("Processing "+paleoName+" ("+name+")");

			File magCSVFile = new File(dir, name+"_mags.csv");
			File magTXTFile = new File(dir, name+"_mags.txt");

			HistogramFunction hist = new HistogramFunction(5.05, 40, 0.1);
			
			CSVFile<String> csv = new CSVFile<String>(true);
			csv.addLine(Lists.newArrayList("Rupture ID", "Mag"));

			int num = 0;
			for (int rupID : rupSet.getRupturesForSection(paleo.getSectionIndex())) {
				double mag = rupSet.getMagForRup(rupID);
				hist.add(mag, 1.0);
				num++;
				csv.addLine(rupID+"", mag+"");
			}
			
//			hist.normalizeBySumOfY_Vals();
			hist.setName("Mag Histogram for rups involving "+paleo.getPaleoSiteName()+" (Sect #"+paleo.getSectionIndex()+")");
			hist.setInfo("(based on "+num+" ruptures)");
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			
			ArrayList<DiscretizedFunc> funcs = new ArrayList<DiscretizedFunc>();
			funcs.add(hist);
			ArrayList<PlotCurveCharacterstics> chars =
					Lists.newArrayList(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
			
			gp.setBackground(Color.WHITE);
			gp.drawGraphPanel("Magnitude", "Number", funcs, chars, hist.getName());
			gp.getChartPanel().setSize(1000, 800);
			gp.saveAsPNG(new File(dir, name+"_hist.png").getAbsolutePath());
			
			csv.writeToFile(magCSVFile);
			csv.writeToTabSeparatedFile(magTXTFile, 1);
		}
	}

}
