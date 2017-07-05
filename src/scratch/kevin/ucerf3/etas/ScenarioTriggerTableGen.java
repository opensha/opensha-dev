package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.util.ClassUtils;
import org.opensha.commons.util.FileNameComparator;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Table;

import scratch.UCERF3.erf.ETAS.ETAS_Simulator.TestScenario;
import scratch.UCERF3.erf.ETAS.ETAS_Params.U3ETAS_ProbabilityModelOptions;

public class ScenarioTriggerTableGen {

	public static void main(String[] args) throws IOException {
		if (args.length < 2 || args.length > 4) {
			System.err.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(ScenarioTriggerTableGen.class)+
					" <input-dir> <output-dir> [<start date: yyyy_MM_dd> <end date: yyyy_MM_dd>]");
			System.exit(2);
		}
		
		File inputDir = new File(args[0]);
		File outputDir = new File(args[1]);
		
		Preconditions.checkState(inputDir.exists() && inputDir.isDirectory(),
				"Input directory doesn't exist or isn't a directory: %s", inputDir.getAbsolutePath());
		Preconditions.checkState(outputDir.exists() && outputDir.isDirectory(),
				"Output directory doesn't exist or isn't a directory: %s", inputDir.getAbsolutePath());
		
		Date startDate = null;
		Date endDate = null;
		if (args.length > 2) {
			try {
				startDate = MPJ_ETAS_SimulatorScriptGen.df.parse(args[2]);
				if (args.length == 4)
					endDate = MPJ_ETAS_SimulatorScriptGen.df.parse(args[3]);
			} catch (ParseException e) {
				System.err.println("Couldn't parse date: "+e.getMessage());
				System.exit(1);
			}
		}
		
		File[] subDirs = inputDir.listFiles();
		Arrays.sort(subDirs, new FileNameComparator());
		
		Table<TestScenario, U3ETAS_ProbabilityModelOptions, EvenlyDiscretizedFunc> allMeanTable = HashBasedTable.create();
		Table<TestScenario, U3ETAS_ProbabilityModelOptions, EvenlyDiscretizedFunc> primaryTable = HashBasedTable.create();
		Table<TestScenario, U3ETAS_ProbabilityModelOptions, EvenlyDiscretizedFunc> fractWithAboveTable = HashBasedTable.create();
		
		for (File subDir : subDirs) {
			if (!subDir.isDirectory())
				continue;
			String name = subDir.getName();
			if (!name.startsWith("20") || !name.contains("-")) {
				System.out.println("(Not a job directory, skipping: "+name+")");
				continue;
			}
			if (startDate != null) {
				// check date range
				try {
					Date date = MPJ_ETAS_SimulatorScriptGen.df.parse(name.substring(0, name.indexOf("-")));
					if (date.getTime() < startDate.getTime() || (endDate != null && date.getTime() > endDate.getTime())) {
						System.out.println("(Outside of date range, skipping: "+name+")");
						continue;
					}
				} catch (ParseException e) {
					System.out.println("(Couldn't parse date from dir, skipping: "+name+")");
					continue;
				}
			}
			// now parse scenario
			TestScenario scenario = null;
			String nameLower = name.toLowerCase();
			for (TestScenario candidate : TestScenario.values()) {
				if (nameLower.contains(candidate.name().toLowerCase())) {
					scenario = candidate;
					break;
				}
			}
			if (scenario == null) {
				System.out.println("(No scenario detected, skipping: "+name+")");
				continue;
			}
			U3ETAS_ProbabilityModelOptions probModel = null;
			for (U3ETAS_ProbabilityModelOptions candidate : U3ETAS_ProbabilityModelOptions.values()) {
				if (nameLower.contains(candidate.name().toLowerCase())) {
					probModel = candidate;
					break;
				}
			}
			if (probModel == null) {
				System.out.println("(No prob model detected, skipping: "+name+")");
				continue;
			}
			System.out.println("Detected "+scenario.name()+" "+probModel.name());
			
			File tableFile = new File(new File(subDir, "plots"), "consolidated_aftershocks_mag_num_cumulative.csv");
			if (!tableFile.exists()) {
				System.out.println("\tTable not present, skipping");
				continue;
			}
			
			CSVFile<String> csv = CSVFile.readFile(tableFile, true);
			
			if (allMeanTable.contains(scenario, probModel))
				System.out.println("Already exists, overwriting with newer");
			
			allMeanTable.put(scenario, probModel, loadColumn(csv, "Mean", 1));
			primaryTable.put(scenario, probModel, loadColumn(csv, "Primary", 11));
			fractWithAboveTable.put(scenario, probModel, loadColumn(csv, "Fract With ≥ Mag", 10));
		}
		
		System.out.println("Loaded "+allMeanTable.size()+" scenarios");
		
		List<TestScenario> scenarios = Lists.newArrayList(allMeanTable.rowKeySet());
		Collections.sort(scenarios, new EnumNameComparator());
		List<U3ETAS_ProbabilityModelOptions> probModels = Lists.newArrayList(allMeanTable.columnKeySet());
		Collections.sort(probModels, new EnumNameComparator());
		
		CSVFile<String> allMeanCSV = new CSVFile<String>(true);
		CSVFile<String> primaryCSV = new CSVFile<String>(true);
		CSVFile<String> fractWithAboveCSV = new CSVFile<String>(true);
		
		// N≥5, N≥5.5, N≥6.3, N≥6.7, N≥7.0, N≥7.4, N≥7.8, N≥8.2
		double[] mags = { 5, 5.5, 6.3, 6.7, 7, 7.4, 7.8, 8.2 };
		List<String> numHeader = Lists.newArrayList("Scenario", "Prob Model");
		List<String> fractHeader = Lists.newArrayList("Scenario", "Prob Model");
		for (double mag : mags) {
			numHeader.add("N≥"+(float)mag);
			fractHeader.add("Fract≥"+(float)mag);
		}
		allMeanCSV.addLine(numHeader);
		primaryCSV.addLine(numHeader);
		fractWithAboveCSV.addLine(fractHeader);
		
		for (TestScenario scenario : scenarios) {
			for (U3ETAS_ProbabilityModelOptions probModel : probModels) {
				EvenlyDiscretizedFunc allFunc = allMeanTable.get(scenario, probModel);
				if (allFunc == null)
					continue;
				allMeanCSV.addLine(getCSVLine(scenario, probModel, mags, allFunc));
				primaryCSV.addLine(getCSVLine(scenario, probModel, mags, primaryTable.get(scenario, probModel)));
				fractWithAboveCSV.addLine(getCSVLine(scenario, probModel, mags, fractWithAboveTable.get(scenario, probModel)));
			}
		}
		
		allMeanCSV.writeToFile(new File(outputDir, "scenarios_consolidated_table_all_children.csv"));
		primaryCSV.writeToFile(new File(outputDir, "scenarios_consolidated_table_primary.csv"));
		fractWithAboveCSV.writeToFile(new File(outputDir, "scenarios_consolidated_table_fract_with_above.csv"));
	}
	
	private static class EnumNameComparator implements Comparator<Enum<?>> {

		@Override
		public int compare(Enum<?> o1, Enum<?> o2) {
			return o1.name().compareTo(o2.name());
		}
		
	}
	
	private static EvenlyDiscretizedFunc loadColumn(CSVFile<String> csv, String expectedName, int col) {
		String actualName = csv.get(0, col);
		Preconditions.checkState(actualName.equals(expectedName),
				"Bad column index. Expected %s, found %s", expectedName, actualName);
		
		double minMag = Double.parseDouble(csv.get(1, 0));
		double maxMag = Double.parseDouble(csv.get(csv.getNumRows()-1, 0));
		int num = csv.getNumRows()-1;
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(minMag, maxMag, num);
		
		for (int i=0; i<num; i++) {
			int row = i+1;
			double val = Double.parseDouble(csv.get(row, col));
			func.set(i, val);
		}
		
		return func;
	}
	
	private static List<String> getCSVLine(TestScenario scenario, U3ETAS_ProbabilityModelOptions probModel,
			double[] mags, EvenlyDiscretizedFunc func) {
		List<String> line = Lists.newArrayList(scenario.toString(), probModel.toString());
		
		for (double mag : mags) {
			int index = func.getClosestXIndex(mag);
			Preconditions.checkState((float)mag == (float)func.getX(index), "Mag mismatch: %s != %s", mag, func.getX(index));
			line.add((float)func.getY(index)+"");
		}
		
		return line;
	}

}
