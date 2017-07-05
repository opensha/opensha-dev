package scratch.kevin.ucerf3.inversion;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.util.ClassUtils;
import org.opensha.commons.gui.plot.GraphWindow;

import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.MatrixIO;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class CompoundSolAboveWaterlevelCalc {
	
	private static Map<String, boolean[]> getAbovesForDir(File dir) throws IOException {
		Map<String, boolean[]> abovesForPrefix = Maps.newHashMap();
		
		getAbovesForDir(dir, abovesForPrefix);
		
		System.out.println("Loaded "+abovesForPrefix.size()+" branches!");
		
		return abovesForPrefix;
	}
	
	private static void getAbovesForDir(File dir, Map<String, boolean[]> abovesForPrefix) throws IOException {
//		System.out.println("getAbovesForDir: "+dir.getName());
		for (File file : dir.listFiles()) {
			if (file.getName().startsWith("."))
				continue;
			if (file.isDirectory())
				getAbovesForDir(file, abovesForPrefix);
			
			String name = file.getName();
			if (!name.endsWith("_noMinRates.bin"))
				continue;
			
			String prefix = name.substring(0, name.indexOf("_noMinRates.bin"));
			if (prefix.contains("_run"))
				prefix = name.substring(0, name.indexOf("_run"));
			
			double[] rates = MatrixIO.doubleArrayFromFile(file);
			
			boolean[] aboves = abovesForPrefix.get(prefix);
			
			if (aboves == null) {
				aboves = new boolean[rates.length];
				abovesForPrefix.put(prefix, aboves);
			}
			
			for (int i=0; i<rates.length; i++)
				if (rates[i] > 0)
					aboves[i] = true;
		}
	}
	
	private static Map<String, boolean[]> getAbovesForZip(File zipFile) throws IOException {
		return getAbovesForZip(zipFile, -1);
	}
	
	private static Map<String, boolean[]> getAbovesForZip(File zipFile, int maxSols) throws IOException {
		ZipFile zip = new ZipFile(zipFile);
		Map<String, boolean[]> abovesForPrefix = Maps.newHashMap();
		
		Map<String, Integer> counts;
		if (maxSols > 0)
			counts = Maps.newHashMap();
		else
			counts = null;
		
		Enumeration<? extends ZipEntry> entries = zip.entries();
		while (entries.hasMoreElements()) {
			ZipEntry entry = entries.nextElement();
			
			String name = entry.getName();
			if (!name.endsWith("_noMinRates.bin"))
				continue;
			
			String prefix = name.substring(0, name.indexOf("_noMinRates.bin"));
			if (prefix.contains("_run"))
				prefix = name.substring(0, name.indexOf("_run"));
			
			if (counts != null) {
				Integer count = counts.get(prefix);
				if (count == null)
					count = 0;
				if (count == maxSols)
					continue;
				count++;
				counts.put(prefix, count);
			}
			
			double[] rates = MatrixIO.doubleArrayFromInputStream(zip.getInputStream(entry), entry.getSize());
			
			boolean[] aboves = abovesForPrefix.get(prefix);
			
			if (aboves == null) {
				aboves = new boolean[rates.length];
				abovesForPrefix.put(prefix, aboves);
			}
			
			for (int i=0; i<rates.length; i++)
				if (rates[i] > 0)
					aboves[i] = true;
		}
		
		System.out.println("Loaded "+abovesForPrefix.size()+" branches!");
		
		return abovesForPrefix;
	}
	
	private static int[] getAboveWaterlevelCounts(Map<String, boolean[]> abovesForPrefix) throws IOException {
		int[] aboves = null;
		for (boolean[] prefixAboves : abovesForPrefix.values()) {
			if (aboves == null)
				aboves = new int[prefixAboves.length];
			else
				Preconditions.checkState(aboves.length == prefixAboves.length);
			
			for (int i=0; i<aboves.length; i++)
				if (prefixAboves[i])
					aboves[i] = aboves[i] + 1;
		}
		
		return aboves;
	}
	
	private static void plotAboves(int[] aboves , File dir) throws IOException {
		List<int[]> abovesList = Lists.newArrayList(aboves);
		plotAboves(abovesList, dir);
	}
	
	private static void plotAboves(List<int[]> abovesList, File dir) throws IOException {
		System.out.println(abovesList.size()+" aboves");
		ArrayList<DiscretizedFunc> hists = Lists.newArrayList();
		for (int i=0; i<abovesList.size(); i++) {
			int[] aboves = abovesList.get(i);
			
			int max = 0;
			for (int above : aboves)
				if (above > max)
					max = above;
			
			if (i == abovesList.size()-1) {
				FileWriter fw = new FileWriter(new File(dir, "zeros_indexes.txt"));
				for (int rupIndex=0; rupIndex<aboves.length; rupIndex++) {
					if (aboves[rupIndex] == 0)
						fw.write(rupIndex+"\n");
				}
				fw.close();
			}
			
			HistogramFunction hist = new HistogramFunction(0d, max+1, 1d);
			
			for (int above : aboves)
				hist.add(above, 1d);
			
			hists.add(hist);
			
			EvenlyDiscretizedFunc cml = getCmlGreaterOrEqual(hist);
			
			ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
			funcs.add(hist);
//			funcs.add(cml);
			
			ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 5f, Color.BLACK));
//			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			CommandLineInversionRunner.setFontSizes(gp);
			
			String title = "Rups Above Waterlevel";
			
			String xAxisName = "# Solutions With Rup Above Waterlevel";
			String yAxisName = "# Ruptures";
			
			gp.setBackgroundColor(Color.WHITE);
			gp.drawGraphPanel(xAxisName, yAxisName, funcs, chars, title);
			
			String nameAdd;
			if (abovesList.size() > 1)
				nameAdd = "_"+i;
			else
				nameAdd = "";
			
			File file = new File(dir, "rups_above_waterlevel"+nameAdd);
			
			gp.getChartPanel().setSize(1000, 800);
			gp.saveAsPDF(file.getAbsolutePath()+".pdf");
			gp.saveAsPNG(file.getAbsolutePath()+".png");
			gp.saveAsTXT(file.getAbsolutePath()+".txt");
			file = new File(file.getAbsolutePath()+"_small");
			gp.getChartPanel().setSize(500, 400);
			gp.saveAsPDF(file.getAbsolutePath()+".pdf");
			gp.saveAsPNG(file.getAbsolutePath()+".png");
		}
		
		if (abovesList.size() > 1) {
			ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
			List<Color> colors = GraphWindow.generateDefaultColors();
			for (int i=0; i<abovesList.size(); i++)
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, colors.get(i)));
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			CommandLineInversionRunner.setFontSizes(gp);
			
			String title = "Rups Above Waterlevel";
			
			String xAxisName = "# Solutions With Rup Above Waterlevel";
			String yAxisName = "# Ruptures";
			
			double maxY = 0;
			for (DiscretizedFunc func : hists) {
				double y = func.getY(0);
				if (y > maxY)
					maxY = y;
			}
			
			gp.setUserBounds(0d, 10d, 0d, maxY+10d);
			gp.setBackgroundColor(Color.WHITE);
			gp.drawGraphPanel(xAxisName, yAxisName, hists, chars, title);
			
			File file = new File(dir, "rups_above_waterlevel_combined");
			
			gp.getChartPanel().setSize(1000, 800);
			gp.saveAsPDF(file.getAbsolutePath()+".pdf");
			gp.saveAsPNG(file.getAbsolutePath()+".png");
			gp.saveAsTXT(file.getAbsolutePath()+".txt");
			file = new File(file.getAbsolutePath()+"_small");
			gp.getChartPanel().setSize(500, 400);
			gp.saveAsPDF(file.getAbsolutePath()+".pdf");
			gp.saveAsPNG(file.getAbsolutePath()+".png");
		}
		
		// now number of above waterlevel as a function of the number of runs
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(1d, abovesList.size(), 1d);
		
		boolean[] aboves = null;
		for (int i=0; i<abovesList.size(); i++) {
			int[] runCounts = abovesList.get(i);
			if (aboves == null)
				aboves = new boolean[runCounts.length];
			for (int r=0; r<runCounts.length; r++)
				aboves[r] = aboves[r] || runCounts[r] > 0;
			
			int cnt = 0;
			for (boolean above : aboves)
				if (above)
					cnt++;
			
			func.set(i, cnt);
		}
		
		ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
		funcs.add(func);
//		GraphWindow gw = new GraphWindow(funcs, "Num Non Zeros");
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		CommandLineInversionRunner.setFontSizes(gp);
		gp.setBackgroundColor(Color.WHITE);
		gp.setUserBounds(0, func.getMaxX(), 0, aboves.length);
		gp.drawGraphPanel("# Runs", "# Ruptures", funcs, chars,
				"Ruptures Above Waterlevel");
		File file = new File("/tmp/compound_rups_above_waterlevel");
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		gp.saveAsTXT(file.getAbsolutePath()+".txt");
		file = new File(file.getAbsolutePath()+"_small");
		gp.getChartPanel().setSize(500, 400);
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		gp.saveAsPNG(file.getAbsolutePath()+".png");
	}
	
	private static EvenlyDiscretizedFunc getCmlGreaterOrEqual(EvenlyDiscretizedFunc func) {
		EvenlyDiscretizedFunc cml = new EvenlyDiscretizedFunc(func.getMinX(), func.size(), func.getDelta());
		
		double tot = 0d;
		for (int i=func.size(); --i>=0;) {
			tot += func.getY(i);
			cml.set(i, tot);
		}
		
		return cml;
	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		if (args.length != 1) {
			System.out.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(CompoundSolAboveWaterlevelCalc.class)+" <dir>");
			System.exit(2);
		}
		int num = 6;
		File dir = new File(args[0]);
		List<int[]> abovesList;
		if (dir.getName().endsWith(".zip")) {
			File zipFile = dir;
			dir = zipFile.getParentFile();
			Preconditions.checkArgument(zipFile.exists(), "Zip file "+zipFile.getAbsolutePath()+" doesn't exist!");
			abovesList = Lists.newArrayList();
			if (num > 0) {
				for (int i=0; i<num; i++) {
					Map<String, boolean[]> abovesMap = getAbovesForZip(zipFile, i+1);
					abovesList.add(getAboveWaterlevelCounts(abovesMap));
				}
			} else {
				Map<String, boolean[]> abovesMap = getAbovesForZip(zipFile, -1);
				abovesList.add(getAboveWaterlevelCounts(abovesMap));
			}
		} else {
			Preconditions.checkArgument(dir.exists(), "Dir "+dir.getAbsolutePath()+" doesn't exist!");
			Preconditions.checkArgument(dir.isDirectory(), dir.getAbsolutePath()+" isn't a directory!");
			Map<String, boolean[]> abovesMap = getAbovesForDir(dir);
			int[] aboves = getAboveWaterlevelCounts(abovesMap);
			abovesList = Lists.newArrayList(aboves);
		}
		
		plotAboves(abovesList, dir);
	}

}
